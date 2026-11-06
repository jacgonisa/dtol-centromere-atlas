# Atlas data as TSV: fill these, then convert to JSON

Three tab-separated tables, **pre-filled with the current 325 species**. Use them as the template: edit rows, or append a row per new species. `fields.tsv` lists every column with its type, unit, allowed values and whether it is required.

| File | One row per | Key |
|---|---|---|
| `species.tsv` | species (22 columns) | `species_id` |
| `cenpa_sequences.tsv` | CENP-A/CENH3 protein (0..n per species) | `species_id` + `id` |
| `sat_families.tsv` | satellite family: list **all** of them (0..n per species) | `species_id` + `fam` |
| `fields.tsv` | column of the three tables above | `table` + `column` |

Rules:
- `species_id` = ToLID + assembly version, **exact case** (`fTauBub2.1`), and it must match a tip of the tree. While the paper tree is frozen, the set of species is exactly the 325 tips. New species arrive with the next tree (see `../frozen/FROZEN.md`).
- **Empty cell = no data (null).** Required columns cannot be empty.
- Write tabs and newlines inside text as `\t` and `\n`.
- **Don't add colour columns or a sequence-length column.** `branch_col`, `arch_col` and the CENP-A `length` are computed by the converter.
- **The satellite summaries are not in `species.tsv`.** `n_sat_fam`, `sat_monomer`, `sat_totbp` and `sat_gc` are computed from `sat_families.tsv`: dominant family = most copies (`n`). Just list every centromeric family, in any order.
- If editing in Excel, save as *UTF-8 tab-delimited text* and check that no ID or number got auto-reformatted.

Convert and check:

```
python3 tsv_to_json.py . ../frozen/paper_325sp_v1.nwk tree_data.json
python3 ../validate_tree_data.py tree_data.json
```

Converting these tables unchanged reproduces the live `tree_data.json` field-for-field, which was checked on 2026-11-06.
