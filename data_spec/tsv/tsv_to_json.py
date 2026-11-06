#!/usr/bin/env python3
"""TSV tables -> tree_data.json for the DToL Centromere Atlas.  Standard library only.

    python3 tsv_to_json.py <dir with species.tsv, cenpa_sequences.tsv, sat_families.tsv> \
                           <frozen tree .nwk> <out tree_data.json>

Column types come from fields.tsv (same folder as this script). Empty cell = null.
Derived here (do NOT put them in the TSVs): branch_col (from broad), arch_col (from arch),
cenpa length (from aa), the JSON key (species_id lower-cased), the meta block, and the
satellite summaries from sat_families.tsv: families are sorted by copy number (n, then bp);
dominant = first; n_sat_fam = number of families; sat_monomer/sat_totbp/sat_gc = dominant's
width/bp/gc.
"""
import csv, hashlib, json, re, sys
from pathlib import Path

BROAD_COL = {"Vertebrates": "#1565C0", "Invertebrates": "#EF6C00", "Viridiplantae": "#2E7D32",
             "Fungi": "#6A1B9A", "Protists": "#C62828"}
ARCH_COL = {"Satellite": "#ff006e", "Transposon": "#3a86ff", "Satellite/transposon": "#8338ec",
            "Holocentric": "#fb5609", "Monocentric sequence unknown": "#808080", "Unknown": "#808080"}
FROZEN_VERSION = {"4b403fa34ad997413d975c9b93d1517b": "paper_325sp_v1"}
TREE_DESC = ("FastSpeciesTree v1.0 (IQ-TREE sensitive), midpoint-rooted, ape::chronos correlated "
             "lambda=0.1, 62 calibrations (over_calib.tsv); branch lengths in million years")

HERE = Path(__file__).resolve().parent
FIELDS = {}
for r in csv.DictReader(open(HERE / "fields.tsv"), delimiter="\t"):
    FIELDS.setdefault(r["table"], {})[r["column"]] = r
CAST = {"int": int, "float": float, "string": str, "enum": str}


def read(path, table):
    spec, rows = FIELDS[table], []
    with open(path, newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        missing = set(spec) - set(rd.fieldnames or [])
        if missing: sys.exit(f"{path}: missing columns {sorted(missing)}")
        for i, r in enumerate(rd, start=2):
            out = {}
            for col, s in spec.items():
                v = r[col]
                if v == "":
                    if s["required"] == "yes": sys.exit(f"{path} line {i}: required '{col}' is empty")
                    out[col] = None
                    continue
                try: out[col] = CAST[s["type"]](v.replace("\\n", "\n").replace("\\t", "\t"))
                except ValueError: sys.exit(f"{path} line {i}: '{col}'={v!r} is not {s['type']}")
            rows.append(out)
    return rows


def main(src, tree, out):
    src = Path(src)
    newick = Path(tree).read_text().strip()
    md5 = hashlib.md5(newick.encode()).hexdigest()
    tips = [re.sub(r"\.fa$", "", t, flags=re.I) for t in re.findall(r"[(,]\s*([^():,;\s]+)\s*(?=[:,)])", newick)]

    species = {}
    for r in read(src / "species.tsv", "species"):
        sid = r.pop("species_id"); k = sid.lower()
        if k in species: sys.exit(f"duplicate species_id {sid}")
        if r["broad"] not in BROAD_COL: sys.exit(f"{sid}: broad '{r['broad']}' not in {list(BROAD_COL)}")
        if r["arch"] not in ARCH_COL: sys.exit(f"{sid}: arch '{r['arch']}' not in {list(ARCH_COL)}")
        r["branch_col"], r["arch_col"] = BROAD_COL[r["broad"]], ARCH_COL[r["arch"]]
        r["cenpa_sequences"], r["sat_families"] = [], []
        species[k] = r
    for r in read(src / "cenpa_sequences.tsv", "cenpa_sequences"):
        k = r.pop("species_id").lower()
        if k not in species: sys.exit(f"cenpa_sequences: unknown species_id {k}")
        r["length"] = len(r["aa"]); species[k]["cenpa_sequences"].append(r)
    for r in read(src / "sat_families.tsv", "sat_families"):
        k = r.pop("species_id").lower()
        if k not in species: sys.exit(f"sat_families: unknown species_id {k}")
        species[k]["sat_families"].append(r)

    for r in species.values():                       # satellite summaries derive from the families
        f = sorted(r["sat_families"], key=lambda x: (-x["n"], -x["bp"]))
        r["sat_families"] = f
        top = f[0] if f else {}
        r["n_sat_fam"] = len(f) or None
        r["sat_monomer"], r["sat_totbp"], r["sat_gc"] = top.get("width"), top.get("bp"), top.get("gc")

    keys = [t.lower() for t in tips]
    if set(keys) != set(species):
        sys.exit(f"tree/species mismatch: no row for {sorted(set(keys) - set(species))[:10]}; "
                 f"not in tree {sorted(set(species) - set(keys))[:10]}")
    D = {"meta": {"tree_version": FROZEN_VERSION.get(md5, "unversioned"), "frozen": md5 in FROZEN_VERSION,
                  "newick_md5": md5, "tree": TREE_DESC, "schema": "data_spec/tree_data.schema.json (v1)"},
         "newick": newick, "species": {k: species[k] for k in keys}}
    json.dump(D, open(out, "w"))
    print(f"{len(keys)} species, {sum(len(s['cenpa_sequences']) for s in species.values())} CENP-A seqs, "
          f"{sum(len(s['sat_families']) for s in species.values())} sat families, tree {D['meta']['tree_version']} -> {out}")


if __name__ == "__main__":
    if len(sys.argv) != 4: sys.exit(__doc__)
    main(*sys.argv[1:])
