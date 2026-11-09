# DToL Centromere Atlas — data specification (v1)

This describes **everything needed to regenerate the atlas data automatically**:

1. how to supply a new species (one proteome + one manifest row),
2. how the species tree is rebuilt with FastSpeciesTree,
3. the exact `tree_data.json` format the web app reads.

Files in this folder:

| File | What it is |
|---|---|
| `DATA_SPEC.md` | this document |
| `tree_data.schema.json` | JSON Schema (draft 2020-12) for `tree_data.json`. Validate with any validator, e.g. `check-jsonschema --schemafile tree_data.schema.json tree_data.json` |
| `species_manifest.template.tsv` | the per-species input table (one row per species) with 3 real example rows |
| `tsv/` | **the fill-in template:** `species.tsv`, `cenpa_sequences.tsv`, `sat_families.tsv` (pre-filled with the 325 species), `fields.tsv` (column dictionary), `tsv_to_json.py` (converter). Start here. |
| `frozen/` | the frozen paper tree (`paper_325sp_v1.nwk`) and `FROZEN.md` (provenance, md5s, how to deploy a new tree) |
| `validate_tree_data.py` | standard-library-only checker: schema-level checks **plus** newick ↔ species-key consistency. `python3 validate_tree_data.py ../tree_data.json` |

---

## 0. The one rule: a single stable species ID

Everything joins on **`species_id` = ToLID + assembly version**, e.g. `fTauBub2.1`, `drHerSpho1.hap1.1`, `pEimTen1.1`.

- The same string is used as the proteome filename stem, the tree tip label, the manifest key and (lower-cased) the JSON key.
- Allowed characters: `A–Z a–z 0–9 . _ -`. **No spaces, `:`, `(`, `)`, `,`, `;`, `'`** — they break Newick.
- Never reuse or rename an ID. If a tip is found to be mislabelled (this happened: `ddHotPalu2` was actually `drHerSpho1`, *Heracleum sphondylium*), fix it **in the manifest**, which is the only source of truth, then regenerate.

Why this matters: in the current 325-species build, the tree tips went through three different naming schemes
(`ucPycProv1_helixer.pep` → BUSCO-style lowercase code → `ucPycProv1.1.fa`), with hand-written exception lists
(`rosCan_S27_v1`, `daTanVulg1.hap2.1` vs `hap1.1`). A single ID removes those exceptions.

---

## Pipeline (since 2026-11-09): TSV-first

```
Piotr's tables (November2026_finalsubmission_files/)            data_spec/tsv/*.tsv  (source of truth)
  sat_bp_gc_count_data.csv, cen_families_Nov2026.csv  ──refresh_from_authoritative.py──▶  sat_families.tsv + species.tsv (sat_pct, HOR)
                                                                                                    │ tsv_to_json.py (+ frozen tree)
                                                                                                    ▼
                                                                         tree_data.json ── validate_tree_data.py (must pass)
```

One command: `python3 build_atlas_data.py` (in `web_prototype/`). It runs all three steps and restores the previous `tree_data.json` if validation fails. Edit the TSVs, never the JSON. A species Piotr marks `satellite == FALSE` gets **no** satellite data: no families, and `sat_pct`/`rep_seq` are cleared.

## 1. Input per species

### 1a. Manifest row (`species_manifest.tsv`, tab-separated, UTF-8)

| Column | Required | Example | Notes |
|---|---|---|---|
| `species_id` | yes | `fTauBub2.1` | primary key, see §0 |
| `scientific_name` | yes | `Taurulus bubalis` | binomial, no authority |
| `ncbi_taxid` | yes | `61643` | used for TimeTree/NCBI lookups and calibration |
| `assembly_accession` | yes | `GCA_910589615.1` | INSDC accession, versioned |
| `taxa1` | yes | `Actinopterygii` | display subgroup (see §3 for the current vocabulary) |
| `broad` | yes | `Vertebrates` | one of `Vertebrates`, `Invertebrates`, `Viridiplantae`, `Fungi`, `Protists` |
| `proteome_uri` | yes | `s3://<bucket>/proteomes/fTauBub2.1.faa` | where the proteome lives |
| `proteome_md5` | yes | `9e107d9d…` | md5 of the proteome file; lets a pipeline detect changes |
| `annotation_method` | yes | `helixer` | the 325 set is all Helixer — keep new species on the same method if possible |
| `include_in_tree` | yes | `yes` | `no` = keep in database but leave out of the tree (e.g. failed QC) |
| `notes` | no | | free text |

### 1b. Proteome file (input to FastSpeciesTree)

- **One protein FASTA per species**, named exactly **`<species_id>.faa`**.
- Amino acids only; **one protein per gene** (longest isoform); unique header IDs; strip a trailing `*` stop character.
- Not gzipped when FastSpeciesTree runs (decompress while staging).

FastSpeciesTree details that affect how files are laid out (checked against `FastSpeciesTree.py` v1.0):

- `-f <dir>` treats **every non-hidden file in the directory** as a proteome (it only skips `.dmnd`). The directory must contain **only** the `.faa` files: no manifest, no README, no `.gz`.
- It writes DIAMOND databases (`<file>.dmnd`) **into the input directory**, so the directory must be writable. Stage a fresh local copy from S3 for each run.
- **The tip label is the input filename**, extension included (`fTauBub2.1.faa`). Strip `.faa` when post-processing.

---

## 2. Rebuilding the tree (current recipe)

```
# stage: download every manifest row with include_in_tree == yes into proteomes/
python3 FastSpeciesTree.py -f proteomes/ -o fst_out/ -t 32 -s sensitive
# -> fst_out/Results/*.treefile   (IQ-TREE, sensitive mode)
```

What the 325-species paper tree used:

- FastSpeciesTree v1.0 (Holmes 2026, bioRxiv 10.64898/2026.01.20.700630) in `-s sensitive` mode (IQ-TREE).
- 194 single-copy ODB12 orthologues, DIAMOND pseudo-alignment, concatenated.
- Partitioned ModelFinder Plus over {LG, JTT, Q.BIRD, Q.MAMMAL, Q.INSECT, Q.PLANT, Q.YEAST} × {I, G, I+G}.
- 1000 UFBoot + 1000 SH-aLRT.

Post-processing, in order:

1. **Strip the extension** from tip labels (`.faa`, or legacy `_helixer.pep` / `.fa` / `.fasta`).
2. **Root.** The paper tree is **midpoint-rooted** (`phytools::midpoint.root`, then `ape::multi2di`); there is no explicit outgroup. On the 325-species tree the midpoint happens to fall between the two *Eimeria* (Alveolata) and everything else. Adding species can move a midpoint root, so **a new tree should state its rooting explicitly**, e.g. a named outgroup set, and the choice should be agreed with Jacob. The placement of Alveolata/Excavata is a known soft spot, accepted for the paper.
3. **Time-calibrate.** The atlas uses `ape::chronos(model = "correlated", lambda = 0.1)` with the 62-constraint calibration set (`01_species_tree/over_calib.tsv`; root = Eukaryota 1085–1671 Ma). Generator: `01_species_tree/calibrate_chronos_correlated_325sp.R <constraints.tsv> <out.nwk>`. Its λ is hard-coded to 1, and the paper tree used λ = 0.1, so expose λ as an argument before automating this. Calibration points are defined as the MRCA of named tip pairs, so they survive adding species. A **new clade with no calibration point is only interpolated**, so add calibration points when a new lineage arrives. The branch lengths in `newick` are in **millions of years**.
4. **Write `tree_data.json`** (§3) and run `validate_tree_data.py`.

Important caveats for automating this:

- **Adding species is a full rebuild, not an incremental insert.** FastSpeciesTree re-selects genes and re-infers everything. Its pseudo-alignment step is cheap enough that this is the intended mode, but **topology can shift**, so node-level results are tied to one tree version.
- **The paper tree is frozen** (`frozen/paper_325sp_v1.nwk`, see `frozen/FROZEN.md`) until a new tree is deliberately deployed, in some months. All paper analyses (ASR, Pagel's λ, calibration QC) refer to it. Until then, generators must **copy `newick` and `meta` from the frozen tree unchanged** and only regenerate `species` records. `build_atlas_data.py` refuses to run if the frozen tree's md5 changes, and `validate_tree_data.py` checks it too.
- A new tree gets a new `tree_version` and a new file in `frozen/`; the paper file is never edited.

---

## 3. `tree_data.json` (v1, what the web app reads today)

```jsonc
{
  "meta":    { "tree_version": "paper_325sp_v1", "frozen": true,
               "newick_md5": "4b403fa34ad997413d975c9b93d1517b",   // md5 of the newick string
               "tree": "FastSpeciesTree … midpoint-rooted … chronos λ=0.1 …", "schema": "…" },
  "newick":  "((pEimTen1.1.fa:737.30,pxEimPrae1.1.fa:737.30):737.30, ... );",
  "species": {
    "ftaubub2.1": { "name": "Taurulus bubalis", "taxa1": "Actinopterygii", ... },
    ...
  }
}
```

`meta` is ignored by the app. It records which frozen tree the file was built on. While a tree is frozen, `md5(newick)` must equal `meta.newick_md5`.

**Join rule (exact JS used by the app):**

```js
sp(tipName) = species[ tipName.replace(/\.fa$/i, '').toLowerCase() ]
```

So the species key is the **tip label, lower-cased, with any trailing `.fa` removed**. Every tip must have a key, and every key must be a tip (the build drops extra keys). The current tips end in `.fa` for historical reasons. New builds can drop the extension entirely; the same rule still works.

`null` means "not measured / not applicable". A field can also be **absent**: 2 species (`roscan_s27_v1`, `datanvulg1.hap2.1`) currently lack the legacy annotation fields. Treat absent and `null` the same way. Generators should **always emit every key, using `null` when there is no value**.

### 3a. Identity and taxonomy

| Field | Type | Unit / values | Used by app | Source (current build) |
|---|---|---|---|---|
| `name` | string | scientific name | label, dossier title, search | species table (legacy); manifest `scientific_name` going forward |
| `taxa1` | string | `Actinopterygii`, `Aves`, `Mammalia`, `Reptilia`, `Coleoptera`, `Diptera`, `Hymenoptera`, `Lepidoptera`, `Dicots`, `Monocots`, `Bryophyta`, `Algae`, `Fungi`, `Alveolata`, `Discoba`, … (36 values today) | dossier subtitle | species table |
| `broad` | string enum | `Vertebrates` `Invertebrates` `Viridiplantae` `Fungi` `Protists` | clade colouring, view filters (`?view=metazoa` etc.) | species table |
| `branch_col` | hex colour | must match `broad`: Vertebrates `#1565C0`, Invertebrates `#EF6C00`, Viridiplantae `#2E7D32`, Fungi `#6A1B9A`, Protists `#C62828` | branch and label colour | derived from `broad` |

### 3b. Centromere architecture and CENP-A

| Field | Type | Unit / values | Used by app | Source |
|---|---|---|---|---|
| `arch` | string enum | `Satellite`, `Transposon`, `Satellite/transposon`, `Holocentric`, `Monocentric sequence unknown`, `Unknown` | tip marker, badge | master `DTOL_327_master_March.xlsx` → `classification` |
| `arch_col` | hex colour | Satellite `#ff006e`, Transposon `#3a86ff`, Satellite/transposon `#8338ec`, Holocentric `#fb5609`, Monocentric sequence unknown / Unknown `#808080` | tip marker, badge | derived from `arch` |
| `n_cenpa` | int ≥ 0 | number of CENP-A/CENH3 paralogues | ring + badge | CENP-A copy-number table (unique loci, ≥40 % reference identity); equals `len(cenpa_sequences)` for every species today |
| `cenpa_sequences` | array of `{id: string, aa: string, length: int}` | protein sequence(s); `length == len(aa)` | dossier | CENP-A rescue pipeline (`17_cenpa_phylogeny_with_rescue`) |
| `satellite1` | string \| null | **free-text curator note** about the centromere call (e.g. `"Copia-Bryco"`, `"Hi-C identifies clearly!"`, `"fTauBub.185 [186_1]"`) | dossier | master table notes |

### 3c. Satellites and HORs

| Field | Type | Unit | Used by app | Source |
|---|---|---|---|---|
| `sat_pct` | number \| null | % = 100 × Σ family `bp` / chromosome length; null if no centromeric satellite | stats counter, outer bar, dossier | computed from Piotr's `sat_bp_gc_count_data.csv` |
| `cen_mb` | number \| null | Mb, centromere size | dossier | centromere-length table |
| `n_sat_fam` | int \| null | number of centromeric satellite families | dossier | **derived** = `len(sat_families)` |
| `sat_families` | array of `{fam: string, n: int, bp: int, width: int, gc: number, trash_fam: string}` | **all** centromeric families (`satellite == TRUE`), sorted by `n` descending (first = dominant). `fam` = `trash_fam` = TRASH name; `n` copies; `width` mean monomer bp; `bp` = n × width; `gc` % | dossier table (★ = dominant) | Piotr's `sat_bp_gc_count_data.csv` (Nov 2026) |
| `rep_seq` | string \| null | DNA, 600 bp: a representative satellite stretch, **not a monomer**; null if no centromeric satellite | dossier | legacy TRASH output (not in Piotr's tables) |
| `sat_monomer` | int \| null | bp, dominant monomer length | ring | **derived** = dominant family `width` |
| `sat_totbp` | int \| null | bp, total array size of dominant family | ring | **derived** = dominant family `bp` |
| `sat_gc` | number \| null | %, GC of dominant family (0–100) | dossier | **derived** = dominant family `gc` |
| `hor_regimentation` | number \| null | 0–100 | ring, dossier | Piotr's `cen_families_Nov2026.csv` → `regimentation_score`, dominant family (highest `count_total`); **null for Holocentric** |
| `hor_score` | number \| null | % | ring, dossier | Piotr's `cen_families_Nov2026.csv` → `HOR_score`; **null for Holocentric** |

"Dominant family" for species with several satellite families means **the family with the highest copy number**. This follows the paper's convention (Ian's); it is not a weighted mean.

**The species-level satellite numbers are derived, never entered.** `n_sat_fam`, `sat_monomer`, `sat_totbp` and `sat_gc` are computed from `sat_families`: sort by `n` descending (ties broken by `bp`), take the first family. The converter and `build_atlas_data.py` both compute them, and `validate_tree_data.py` raises an error if a file disagrees. Before 2026-11-06 the JSON kept only the top 3 families **by total bp** while the summaries used copy number, so the two could disagree (e.g. *iyLasVill1.1*). The HOR values (`hor_regimentation`, `hor_score`) are not derived this way: they come from the HOR table's dominant family, joined through `trash_fam`.

### 3d. Genome and annotation

| Field | Type | Unit | Used by app | Source |
|---|---|---|---|---|
| `genome_mb` | number | Mb (1 decimal) | ring, dossier | master `genome.bp` / 1e6 |
| `chrs` | int | chromosome count | ring, dossier | master `chrs` |
| `genome_gc` | number | % (2 decimals) | ring, dossier | master `genome.gc` |
| `te_mb`, `te_pct` | number \| null | Mb, % of genome | outer bar, dossier | EDTA annotation |
| `gene_mb`, `gene_pct` | number \| null | Mb, % of genome | outer bar, dossier | Helixer gene annotation |
| `trash_mb`, `trash_pct` | number \| null | Mb, % of genome, all TRASH satellites | *not shown* | TRASH |
| `stack_pct` | number \| null | % (≈ te + gene + trash; can exceed 100 because of overlaps) | *not shown* | derived |

---

## 4. Known quirks in v1 (fix in v2, don't copy them)

1. **Display colours are stored in the data** (`branch_col`, `arch_col`). They are pure functions of `broad`/`arch` and belong in the app.
2. **The key is a lower-cased tip label** instead of an explicit `species_id` field. ToLIDs are case-sensitive in meaning (`fTauBub2.1`), so the original casing is only recoverable from the newick.
3. **No accession, taxid or provenance per species.** The dossier cannot link to ENA/NCBI or a genome note yet.
4. `satellite1` mixes curator notes with family IDs.
5. **One monolithic file.** At 325 species it is ~420 KB, mostly CENP-A sequences and satellite tables. At thousands of species it should be split.

### Proposed v2 layout (for the AWS-hosted version, not implemented yet)

```
tree/<tree_version>.json        { "tree_version", "built", "fst_version", "calibration", "newick" }   # tips = species_id
species_index.json              { "<species_id>": { name, taxa1, broad, arch, n_cenpa, genome_mb, chrs,
                                                     genome_gc, sat_monomer, sat_totbp, sat_pct,
                                                     hor_regimentation, hor_score, te_pct, gene_pct } }
species/<species_id>.json       full dossier: everything above + cenpa_sequences, sat_families, rep_seq,
                                assembly_accession, ncbi_taxid, notes, provenance{source, version, date}
```

The app would load the tree and the small index up front, then fetch `species/<id>.json` when a tip is clicked. Each file can be served straight from S3/CloudFront, with CORS enabled for the atlas origin. Until the app is switched over, the generator should keep emitting the v1 `tree_data.json` described in §3.
