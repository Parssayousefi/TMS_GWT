# TMS evidence dashboard

The overview and tables labelled **QC refit** or **QC PID** use the FD-enforced
samples: Aachen13/SUDMEX33 under the configured rules,12/30 under the stricter
sensitivity. Current group models use no extra continuous FD covariate; SUDMEX
adjusts baseline outcome. Raw scores, the original normalization reference and
fresh QC-clean references are separately labelled. Earlier tables are explicitly
labelled **Before full FD QC**, and the descriptive control-workspace image
retains its historical sample.

Open `index.html` directly in a modern browser. All interactive data, styling, scripts and the exact control-workspace figure are embedded, so the main dashboard works offline. Keep `tables/` and the JSON receipts beside it for downloads. No framework, package installation, remote fonts or external service is required.

Rebuild from the pinned aggregate sources with:

```powershell
C:\Miniforge3\envs\gwt\python.exe setup/build_tms_dashboard.py
```

The builder writes `dashboard/` and a compact release receipt under
`PROVENANCE/fd_reanalysis_20261126/`. It does not recompute statistical tests,
alter input tables or choose exclusions. `provenance.json` lists each source
hash, export hash, row count, transformations and verification receipts.
`build_receipt.json` records build checks; browser verification is separate in
`browser_qa.json`. The generated comparison figure uses the same saved estimates
and marginal intervals.

All raw and normalized corrected core families, acute SUDMEX comparisons, all-post descriptions, clinical associations, motion/intensity/influence checks, seed stability, method/HRF/network/temporal sensitivities, the nine-test haemodynamic follow-up and every rank/omission sensitivity, corrected Aachen target flow and validated SUDMEX target flow are included. Families overlap and rows must not be counted as independent hypotheses. All original summary p-value columns are accessible in result details and exports.

Only SUDMEX rows are imported from the older mixed target-flow run. Historical Aachen results are invalidated. No patient-level clinical records or raw imaging are bundled; leave-one-out aggregate fits use anonymous omission labels. The combined control image is copied byte-for-byte.

The main view shows every reported row by default. Optional adjusted-p filtering uses the displayed, explicitly labeled correction family and does not apply a new dashboard-wide correction. Missing numeric values remain `null` and display as “—”. Method interactions, HRF associations and descriptive control overlap are distinguished from treatment comparisons.

The corrected target-flow release passed independent numerical checks and agreement with saved corrected production-MATLAB summaries. A fresh 24-case direct MATLAB atom check was unavailable because of license-server error 96; this limit remains explicit.

`template.html` and `DESIGN.md` are editable presentation sources. `index.html`, `data.json`, summary CSV copies and provenance are generated from them and the allowlisted source tables. Git publication is handled separately; this builder does not stage, commit, or push files.

The historical source inventory includes HRF associations, spatial contrasts,
baseline/arm sensitivities and omission fits. These are retained records, not the
selected publication narrative. The current overview displays the raw-versus-
normalized outside-workspace synergy comparison. Neither haemodynamic
associations nor reinterpretation of nonsignificant changes as stability is used
as the shared treatment story.

The complete six-mode PID follow-up is included: seven fixed pair scopes, both representations, both normalization seeds plus their combined reference, external regional change templates, all omission fits, and calibration diagnostics. Mode-level sums are formed before conditional normalization. No shared treatment effect was established; these are exploratory follow-ups, not a new preregistration.
