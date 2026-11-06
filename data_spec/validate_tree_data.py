#!/usr/bin/env python3
"""Validate a DToL Centromere Atlas tree_data.json (v1).  Standard library only.

    python3 validate_tree_data.py path/to/tree_data.json [--manifest species_manifest.tsv]

Checks
  * top-level shape: {"newick": str, "species": {key: record}}
  * every newick tip maps to a species key with the app's join rule
    (tip.lower() with a trailing ".fa" removed) and every key is a tip
  * required fields, types, enums, ranges and colour consistency (see DATA_SPEC.md §3)
  * n_cenpa == len(cenpa_sequences), cenpa length == len(aa), holocentric => no HOR values
  * frozen tree: md5(newick) == meta.newick_md5 (paper tree = 4b403fa34ad997413d975c9b93d1517b)
  * optional: every tip is a species_id in the manifest (include_in_tree == yes)
Exit code 0 = valid (warnings allowed), 1 = errors found.
"""
import csv, hashlib, json, re, sys

BROAD_COL = {"Vertebrates": "#1565C0", "Invertebrates": "#EF6C00", "Viridiplantae": "#2E7D32",
             "Fungi": "#6A1B9A", "Protists": "#C62828"}
ARCH_COL = {"Satellite": "#ff006e", "Transposon": "#3a86ff", "Satellite/transposon": "#8338ec",
            "Holocentric": "#fb5609", "Monocentric sequence unknown": "#808080", "Unknown": "#808080"}
REQUIRED = {"name": str, "taxa1": str, "broad": str, "branch_col": str, "arch": str,
            "arch_col": str, "genome_mb": (int, float), "chrs": int, "genome_gc": (int, float)}
NUM = (int, float)
OPTIONAL = {  # field: (allowed types, min, max)   None = unbounded
    "n_cenpa": (int, 0, None), "sat_pct": (NUM, 0, 100), "cen_mb": (NUM, 0, None),
    "n_sat_fam": (int, 0, None), "sat_monomer": (int, 1, None), "sat_totbp": (int, 0, None),
    "sat_gc": (NUM, 0, 100), "hor_regimentation": (NUM, 0, 100), "hor_score": (NUM, 0, None),
    "te_mb": (NUM, 0, None), "te_pct": (NUM, 0, None), "gene_mb": (NUM, 0, None),
    "gene_pct": (NUM, 0, None), "trash_mb": (NUM, 0, None), "trash_pct": (NUM, 0, None),
    "stack_pct": (NUM, 0, None), "satellite1": (str, None, None), "rep_seq": (str, None, None)}


def tip_names(newick):
    # leaf = label directly after '(' or ',' ; strip quotes and branch length
    return [m.strip().strip("'\"") for m in re.findall(r"[(,]\s*([^():,;\s]+)\s*(?=[:,)])", newick)]


def key_of(tip):
    return re.sub(r"\.fa$", "", tip, flags=re.I).lower()


def main(argv):
    if len(argv) < 2:
        print(__doc__); return 1
    D = json.load(open(argv[1]))
    errs, warns = [], []
    E, W = errs.append, warns.append

    if not isinstance(D.get("newick"), str) or not isinstance(D.get("species"), dict):
        print("ERROR: need {'newick': str, 'species': object}"); return 1
    S = D["species"]
    meta = D.get("meta")
    if not isinstance(meta, dict): W("no 'meta' block: cannot verify which frozen tree this was built on")
    else:
        md5 = hashlib.md5(D["newick"].strip().encode()).hexdigest()
        if meta.get("frozen") and md5 != meta.get("newick_md5"):
            E(f"newick md5 {md5} != meta.newick_md5 {meta.get('newick_md5')} (tree '{meta.get('tree_version')}' is frozen)")
    tips = tip_names(D["newick"])
    keys = [key_of(t) for t in tips]
    dup = {k for k in keys if keys.count(k) > 1}
    if dup: E(f"duplicate tips after key normalisation: {sorted(dup)[:10]}")
    for t, k in zip(tips, keys):
        if k not in S: E(f"tip '{t}' -> key '{k}' has no species record")
    for k in set(S) - set(keys): E(f"species key '{k}' is not a tip in newick")
    for t in tips:
        if re.search(r"[\s:(),;']", t): E(f"tip '{t}' contains newick-unsafe characters")

    for k, r in S.items():
        if not re.fullmatch(r"[a-z0-9._-]+", k): E(f"{k}: key must be lower-case [a-z0-9._-]")
        for f, t in REQUIRED.items():
            if f not in r or r[f] is None: E(f"{k}: missing required '{f}'")
            elif not isinstance(r[f], t) or isinstance(r[f], bool): E(f"{k}: '{f}' wrong type {type(r[f]).__name__}")
        if r.get("broad") not in BROAD_COL: E(f"{k}: broad '{r.get('broad')}' not in {list(BROAD_COL)}")
        elif r.get("branch_col") != BROAD_COL[r["broad"]]: W(f"{k}: branch_col {r.get('branch_col')} != {BROAD_COL[r['broad']]} for {r['broad']}")
        if r.get("arch") not in ARCH_COL: E(f"{k}: arch '{r.get('arch')}' not in {list(ARCH_COL)}")
        elif r.get("arch_col") != ARCH_COL[r["arch"]]: W(f"{k}: arch_col {r.get('arch_col')} != {ARCH_COL[r['arch']]}")
        if isinstance(r.get("genome_gc"), NUM) and not 0 <= r["genome_gc"] <= 100: E(f"{k}: genome_gc out of 0-100")

        for f, (t, lo, hi) in OPTIONAL.items():
            if f not in r: W(f"{k}: '{f}' absent (emit null instead)"); continue
            v = r[f]
            if v is None: continue
            if isinstance(v, bool) or not isinstance(v, t): E(f"{k}: '{f}'={v!r} wrong type"); continue
            if lo is not None and v < lo: E(f"{k}: '{f}'={v} < {lo}")
            if hi is not None and v > hi: E(f"{k}: '{f}'={v} > {hi}")

        seqs = r.get("cenpa_sequences", [])
        if not isinstance(seqs, list): E(f"{k}: cenpa_sequences must be a list")
        else:
            for s in seqs:
                if not {"id", "aa", "length"} <= set(s): E(f"{k}: cenpa entry missing id/aa/length"); continue
                if s["length"] != len(s["aa"]): E(f"{k}: {s['id']} length {s['length']} != len(aa) {len(s['aa'])}")
                if not re.fullmatch(r"[A-Za-z*]+", s["aa"]): E(f"{k}: {s['id']} aa has non-letter characters")
            if isinstance(r.get("n_cenpa"), int) and r["n_cenpa"] != len(seqs):
                W(f"{k}: n_cenpa={r['n_cenpa']} but {len(seqs)} sequences")
        for fam in r.get("sat_families", []) or []:
            if not {"fam", "n", "bp", "width", "gc"} <= set(fam): E(f"{k}: sat_families entry missing fam/n/bp/width/gc")
        if r.get("arch") == "Holocentric" and (r.get("hor_regimentation") is not None or r.get("hor_score") is not None):
            E(f"{k}: Holocentric species must have null HOR values")

    if "--manifest" in argv:
        path = argv[argv.index("--manifest") + 1]
        rows = list(csv.DictReader(open(path), delimiter="\t"))
        ids = {row["species_id"].lower() for row in rows if row.get("include_in_tree", "yes").lower() == "yes"}
        for k in keys:
            if k not in ids: E(f"tip key '{k}' not in manifest (include_in_tree=yes)")
        for i in ids - set(keys): W(f"manifest species '{i}' not in tree")

    for w in warns[:40]: print("WARN ", w)
    if len(warns) > 40: print(f"WARN  … {len(warns) - 40} more warnings")
    for e in errs[:60]: print("ERROR", e)
    if len(errs) > 60: print(f"ERROR … {len(errs) - 60} more errors")
    print(f"\n{len(tips)} tips, {len(S)} species records, {len(errs)} errors, {len(warns)} warnings")
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
