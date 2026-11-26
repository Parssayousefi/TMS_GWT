# TMS evidence dashboard

Open `index.html` directly in a modern browser. All interactive data, styling, scripts and the exact control-workspace figure are embedded, so the main dashboard works offline. Keep `tables/` and the JSON receipts beside it for downloads. No framework, package installation, remote fonts or external service is required.

Rebuild from the pinned aggregate sources with:

```powershell
C:\Miniforge3\envs\gwt\python.exe setup/build_tms_dashboard.py
```

The builder only writes `dashboard/`. It does not recompute analyses, alter input tables, or change exclusions. `provenance.json` lists each input table, original source hash, export hash, row count, transformations and verification-receipt paths. `build_receipt.json` records checks at build time. Browser verification is separate in `browser_qa.json`.

All raw and normalized corrected core families, acute SUDMEX comparisons, all-post descriptions, clinical associations, motion/intensity/influence checks, seed stability, method/HRF/network/temporal sensitivities, the nine-test haemodynamic follow-up and every rank/omission sensitivity, corrected Aachen target flow and validated SUDMEX target flow are included. Families overlap and rows must not be counted as independent hypotheses. All original summary p-value columns are accessible in result details and exports.

Only SUDMEX rows are imported from the older mixed target-flow run. Historical Aachen results are invalidated. No patient-level clinical records or raw imaging are bundled; leave-one-out aggregate fits use anonymous omission labels. The combined control image is copied byte-for-byte.

The main view shows every reported row by default. Optional adjusted-p filtering uses the displayed, explicitly labeled correction family and does not apply a new dashboard-wide correction. Missing numeric values remain `null` and display as “—”. Method interactions, HRF associations and descriptive control overlap are distinguished from treatment comparisons.

The corrected target-flow release passed independent numerical checks and agreement with saved corrected production-MATLAB summaries. A fresh 24-case direct MATLAB atom check was unavailable because of license-server error 96; this limit remains explicit.

`template.html` and `DESIGN.md` are editable presentation sources. `index.html`, `data.json`, summary CSV copies and provenance are generated from them and the allowlisted source tables. Git publication is handled separately; this builder does not stage, commit, or push files.
