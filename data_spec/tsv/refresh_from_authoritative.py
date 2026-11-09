#!/usr/bin/env python3
"""refresh_from_authoritative.py — refresh the model TSVs from the authoritative source tables.

Rebuilds the satellite + HOR data in species.tsv / sat_families.tsv from Piotr's Nov-2026
authoritative tables (the single source of truth for centromeric-satellite & HOR metrics):
  - sat_bp_gc_count_data.csv : per-family satellite count / monomer length / GC  (satellite==TRUE only)
  - cen_families_Nov2026.csv : HOR score + regimentation (dominant family = highest count_total)

Updates, per species (matched by ToLID prefix):
  sat_families.tsv  : one row per centromeric satellite family (n, bp, width, gc, trash_fam)
                      -> the converter derives sat_monomer / sat_totbp / sat_gc / n_sat_fam
  species.tsv       : sat_pct  = 100 * sum(family bp) / chromosome length  (satellite segment of the
                                 stacked-genome ring); species Piotr marks satellite==FALSE get
                                 sat_pct and rep_seq cleared (no centromeric satellite)
                      hor_regimentation, hor_score  (dominant family; blank for Holocentric)

Everything else (name, taxonomy, genome size/GC, architecture, CENP-A, TE%/gene% from Charla's
EDTA+Helixer annotation) is left untouched — those are not part of Piotr's tables.

Then run:  python3 tsv_to_json.py . ../frozen/paper_325sp_v1.nwk ../../tree_data.json
"""
import csv, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOV  = Path("/home/jg2070/Desktop/dtol_review_August/November2026_finalsubmission_files")
SAT  = NOV/"sat_bp_gc_count_data.csv"
CF   = NOV/"cen_families_Nov2026.csv"
pfx  = lambda s: re.sub(r'[0-9].*$', '', str(s).lower())
ALIAS = {"drroscani": "roscan_s"}        # Rosa canina: new data drRosCani -> atlas id rosCan_S27_v1 (pfx roscan_s)
remap = lambda s: ALIAS.get(pfx(s), pfx(s))

# ── satellite families + chromosome length, grouped by species prefix ──
fams, chrlen, nosat = {}, {}, set()
for r in csv.DictReader(open(SAT)):
    p = remap(r["fasta_name"])
    if r["satellite"] != "TRUE": nosat.add(p); continue
    n = int(round(float(r["satellite.count"]))); w = float(r["satellite.mean.length"])
    fams.setdefault(p, []).append({"fam": r["trash.name"], "n": n,
        "bp": int(round(n*w)), "width": int(round(w)), "gc": round(float(r["satellite.GC"]), 1),
        "trash_fam": r["trash.name"]})
    chrlen[p] = float(r["Chromosome.length.bp."])

# ── HOR: dominant family (highest count_total) per species prefix ──
hor = {}
for r in csv.DictReader(open(CF)):
    if not r["fasta"] or r["fasta"].lower().startswith("gca_"): continue   # skip human benchmark
    p = remap(r["fasta"])
    try: c = float(r["count_total"])
    except ValueError: c = 0.0
    if p not in hor or c > hor[p][0]:
        hor[p] = (c, r["HOR_score"], r["regimentation_score"], r["is_holocentric"])
numf = lambda v: "" if v in ("", "NA", None) else round(float(v), 2)

# ── patch species.tsv ──
sp = list(csv.DictReader(open(HERE/"species.tsv"), delimiter="\t"))
cols = list(sp[0].keys())
n_sat = n_hor = n_clear = 0
for row in sp:
    p = pfx(row["species_id"])
    if p in fams:                                       # sat_pct from new data (ChrOnly denom)
        tot = sum(f["bp"] for f in fams[p])
        row["sat_pct"] = round(100*tot/chrlen[p], 2) if chrlen.get(p) else row["sat_pct"]
        n_sat += 1
    elif p in nosat:                                    # Piotr: no centromeric satellite -> clear
        if row["sat_pct"] or row["rep_seq"]: n_clear += 1
        row["sat_pct"] = row["rep_seq"] = ""            # legacy values (would draw a ghost satellite segment)
    row["hor_regimentation"] = row["hor_score"] = ""    # all HOR originates from cen_families -> reset
    if p in hor and row["arch"] != "Holocentric":
        _, hs, rs, _ = hor[p]
        row["hor_regimentation"], row["hor_score"] = numf(rs), numf(hs); n_hor += 1

with open(HERE/"species.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(sp)

# ── rebuild sat_families.tsv ──
sid_by_pfx = {pfx(r["species_id"]): r["species_id"] for r in sp}
sf_cols = ["species_id","fam","n","bp","width","gc","trash_fam"]
with open(HERE/"sat_families.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=sf_cols, delimiter="\t"); w.writeheader()
    nrows = 0
    for p, fl in fams.items():
        sid = sid_by_pfx.get(p)
        if not sid: continue
        for fam in sorted(fl, key=lambda x: (-x["n"], -x["bp"])):
            w.writerow({"species_id": sid, **{k: fam[k] for k in ("fam","n","bp","width","gc","trash_fam")}})
            nrows += 1

print(f"species.tsv: sat_pct set for {n_sat} spp, HOR set for {n_hor} spp (non-holo, dominant family)")
print(f"species.tsv: cleared legacy sat_pct/rep_seq for {n_clear} spp Piotr marks satellite==FALSE")
print(f"sat_families.tsv: {nrows} family rows across {len(fams)} species")
