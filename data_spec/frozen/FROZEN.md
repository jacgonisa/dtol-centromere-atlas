# Frozen trees

The atlas serves exactly one species tree at a time. A tree is **frozen** once it is deployed:
its newick never changes. A new tree is a new file with a new `tree_version`. It is never an edit.

## `paper_325sp_v1` — the published 325-species tree (frozen 2026-11-06, current)

| | |
|---|---|
| File | `paper_325sp_v1.nwk` |
| md5 of the newick string (whitespace stripped) | `4b403fa34ad997413d975c9b93d1517b`. This equals `tree_data.json → meta.newick_md5` |
| md5 of the file on disk (with trailing newline) | `a5e967bf6408b5f1aa0005e37792d873` |
| Tips | 325, labelled `<species_id>.fa` |
| Topology | FastSpeciesTree v1.0, `-s sensitive` (IQ-TREE): 194 ODB12 single-copy orthologues, partitioned MFP, 1000 UFBoot + 1000 SH-aLRT. Raw tree: `2026_trees/fast_species_tree_iqtree_only325/Results/fast_species_tree_iqtree_only325.treefile` |
| Rooting | **Midpoint** (`phytools::midpoint.root`, then `multi2di`). No explicit outgroup. On this tree the midpoint falls between the two *Eimeria* (Alveolata) and all other species. |
| Dating | `ape::chronos(model = "correlated", lambda = 0.1)`, 62 calibration constraints (`01_species_tree/over_calib.tsv`, md5 `fc215dcc6cfa3316026e8727a36c710a`, root = Eukaryota 1085–1671 Ma). Generator: `01_species_tree/calibrate_chronos_correlated_325sp.R` (its default is λ = 1; the paper tree used λ = 0.1). Corroborated by MCMCTree. Branch lengths are in millions of years. |
| Upstream file | `01_species_tree/outputs/full_325sp_calibrated_correlatedlambda01.nwk`. Identical topology and branch lengths (checked 2026-11-06). |
| Only change from upstream | tip `ddHotPalu2.hap1.1.fa` → `drHerSpho1.hap1.1.fa`. That assembly is *Heracleum sphondylium* (GCA_965233475.1), mislabelled as *Hottonia palustris* in the master table. |

The paper analyses (ASR, Pagel's λ, half-lives, calibration QC) all refer to this tree.

## Deploying a new tree (planned, months away)

1. Build it following `../DATA_SPEC.md` §2. Choose and document the rooting **explicitly**: with more species, a midpoint root can move.
2. Save it as `frozen/<new_version>.nwk` and add a section above with md5s and the method.
3. Point `FROZEN` in `build_atlas_data.py` at the new file and md5. Leave `paper_325sp_v1.nwk` in place.
4. Rebuild, run `validate_tree_data.py`, and deploy.
