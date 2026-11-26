"""Refresh current-result pointers without replacing historical data or notebooks."""
from pathlib import Path
import json
import shutil
import nbformat
ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'output/corrected_tms_analysis/20261126T112054_284547Z'
MAN=ROOT/'manuscript/aachen_brief_report/corrected_20261126'
readme=ROOT/'README.md';old=readme.read_text(encoding='utf-8')
header='''# TMS_GWT

Information decomposition and workspace analyses of two open resting-state fMRI
datasets following prefrontal TMS. The manuscript cites the earlier psychedelic
paper; its private participant data and effect maps are not included.

## Current results 26 November 2026

The corrected Aachen analysis is complete: 15 participants, 60 scans, all five
information measures and 20 participation/role outcomes, four spatial scopes,
fresh conditional-null normalization, and four GWT-related expectations.
The earlier significant global Φ finding is superseded: corrected change
+0.001349 nats, 95% CI −0.001525 to +0.004224, maximum-t p5=.7441.
No GWT expectation survived its four-outcome family. Redundancy participation
increased nominally (p=.0295), but Holm4 p=.1178; the amount of whole-brain
persistent redundancy was essentially unchanged (p=.9543).

Exploratory haemodynamic associations are present. Narrower estimated HRFs
were associated with greater Φ changes across the whole brain (Holm9 p=.00849)
and outside the workspace (p=.000477). These are same-signal, model-dependent
associations and do not establish a vascular mechanism, a uniform TMS effect,
or validation of the superseded global increase.

- [Interactive analysis dashboard](dashboard/index.html) covers full results,
  correction families, clinical outcomes, targeted flow and diagnostics.
- [Current manuscript](manuscript/aachen_brief_report/corrected_20261126/manuscript.md)
  and [supplement](manuscript/aachen_brief_report/corrected_20261126/supplement.md).
- [Manuscript DOCX](manuscript/aachen_brief_report/corrected_20261126/Aachen_Brief_Research_Report_Corrected.docx)
  and [supplement DOCX](manuscript/aachen_brief_report/corrected_20261126/Supplementary_Material_Corrected.docx).
- [Supplementary tables and scalar outputs](manuscript/aachen_brief_report/corrected_20261126/TMS_Corrected_Supplementary_Data.zip).
- [Corrected analysis report](docs/TMS_CORRECTED_ANALYSES_20261126.md).

The requested control-workspace image is Supplementary Figure S1, with a
section citing Luppi and Yousefi. Dice overlaps are .776/.741 for matched SUDMEX
sham scans and .743 for the corrected Aachen paired-baseline mean. These are
descriptive partial overlaps, not evidence of identical masks or an intervention
effect. The fixed 235-region HCP reference was not changed.

## Input routing and history

Aachen anatomy was already skull-stripped; forced repeat extraction had clipped
brain tissue. The corrected branch uses `--skull-strip-t1w skip` and all other
production processing definitions are unchanged. Valid current Aachen inputs:
`reruns/aachen_skullstrip_skip_20261125/pipeline/results/`.
Original Aachen outputs under `pipeline/results/` are historical and must not
be mixed with the corrected branch. SUDMEX inputs remain in `pipeline/results/`.
The separate SUDMEX anatomical audit found no comparable extraction failure;
inferior signal coverage remains a limitation in both datasets.

The historical proposal and 25 November exploratory reports are retained for
traceability and do not describe the current corrected Aachen findings. Current
numerical results are in `output/corrected_tms_analysis/20261126T112054_284547Z`;
the current local notebook is `output/jupyter-notebook/TMS_Corrected_Analyses.ipynb`.

## Repository contents and reproduction

Git tracks analysis code, configuration, small group-level references,
manuscript documents, summary tables and the dashboard. Raw imaging, individual
time series, large intermediate matrices, licenses, credentials, local logs,
third-party software installations and private psychedelic results are excluded.
The supplementary ZIP includes deidentified public-dataset scalar outputs.
Open the dashboard locally in a browser; it needs no server or external service.

## Sharing with collaborators

The private repository is at https://github.com/Parssayousefi/TMS_GWT.
The owner can invite colleagues through **Settings → Collaborators → Add people**.
Each colleague must accept the invitation before opening the private repository.
See [GitHub's invitation instructions](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/repository-access-and-collaboration/inviting-collaborators-to-a-personal-repository).

For a results-only review, send `dashboard/index.html`; the interactive dashboard
and control figure are embedded and work offline. To retain the individual source
CSV downloads, send the whole `dashboard/` folder. GitHub displays HTML source;
download the repository ZIP, extract it and open `dashboard/index.html` in a browser.

The executed notebook is in [notebooks/TMS_Corrected_Analyses.ipynb](notebooks/TMS_Corrected_Analyses.ipynb).
Its stored results can be reviewed without running it. Rerunning requires the
documented local analysis inputs and updating the explicit project-root path.

## Reproducing the analyses

With the documented local inputs and the `gwt` Python environment available:

```powershell
python setup/analyse_corrected_tms.py
python setup/normalize_corrected_aachen.py --analysis <new-output-folder>
python setup/extend_corrected_tms_checks.py --analysis <new-output-folder>
python setup/run_corrected_aachen_sensitivities.py --analysis <new-output-folder>
```

Publication-build scripts point to the released analysis folder for reproducible
reporting. Update their explicit run constants only when intentionally releasing
a new analysis. They must never fall back to historical Aachen inputs.
The pipeline requires MATLAB/CONN/rsHRF and original data for production stages;
the statistical summaries alone cannot reconstruct participant imaging.
Fresh selected-pair MATLAB verification for the targeted flow extension was
unavailable because of a license-server error. Independent Python identities
and comparisons with the already computed MATLAB matrices passed.

'''
readme.write_text(header+old[old.index('## Rules'):],encoding='utf-8')
report='''# Corrected TMS analysis release

The previous global Aachen Φ increase does not survive anatomical revalidation.
The current whole-brain estimate is +0.001349 nats (95% CI −0.001525 to +0.004224;
maximum-t p5=.7441). The normalized estimate is +0.682 percentile points
(−0.543 to +1.908; p5=.6027). No outcome survives the complete forty-endpoint
family, and no four-outcome GWT expectation survives in either representation.

Redundancy participation increases nominally (+.004661; p=.02946; Holm4=.11782),
whereas mean persistent redundancy does not (+.00004679 nats; p=.9543).
These measure network distribution and information amount, respectively.

The haemodynamic follow-up retained nine associations and all participant
omission checks. Whole-brain Φ changes covary negatively with fitted HRF width
(−.004593 nats per SD, Holm9=.008486); outside-width and outside-peak associations
also survive (.000477 and .046324). All width omission slopes retain direction,
but corrected support persists in only 12/15 whole-brain and 14/15 outside
omissions. HRF and Φ share the same data and deconvolution model; no independent
physiological or causal mediation inference follows.

Target-flow primary families are null in the eight-person mask subset. Three
left-minus-right baseline-plus-motion models survive their Holm10 sensitivity
(.04762), and left persistent synergy decreases at 120% in an intensity-specific
sensitivity (Holm20=.03311). These findings do not replace the primary test.

SUDMEX's acute FPN contrast falls from +3.14 unadjusted to +1.09 after baseline
and motion adjustment (p=.498). No four GWT expectation or clinical association
survives the relevant family. Its all-post sensitivity retains raw synergy and
normalized workspace source-synergy changes, but later open-label visits prevent
interpreting them as randomized treatment effects.

Control maps have broad spatial agreement with the fixed HCP reference:
SUDMEX matched sham Dice .776/.741; Aachen paired-baseline Dice .743.
The supplied comparison figure is Supplementary Figure S1 with Luppi/Yousefi
citations. No reference mask was redefined from treatment results.

## Released artifacts

- Core: `output/corrected_tms_analysis/20261126T112054_284547Z`
- Target flow: `output/corrected_aachen_target_flow/20261126T113156_656159Z`
- Independent review: `output/corrected_tms_analysis_review`
- Manuscript: `manuscript/aachen_brief_report/corrected_20261126`
- Dashboard: `dashboard/index.html`
- Executed notebook: `output/jupyter-notebook/TMS_Corrected_Analyses.ipynb`

Scientific tables replace the historical Aachen results. Remaining manuscript
administration: Zain's affiliation, author contributions, funding, conflicts and
a public archival code/data identifier. Private GitHub storage does not itself
satisfy a public data-availability statement.
'''
(ROOT/'docs/TMS_CORRECTED_ANALYSES_20261126.md').write_text(report,encoding='utf-8')
(MAN.parent/'CURRENT.txt').write_text('corrected_20261126\n',encoding='utf-8')
(MAN.parent/'ANATOMICAL_REVIEW_PENDING.txt').write_text('HISTORICAL DRAFT SUPERSEDED\nAnatomical correction and statistical reanalysis are complete. The root-folder draft retains invalidated historical Aachen findings for traceability. Use corrected_20261126/ for the current manuscript and supplement. Do not submit historical files.\n',encoding='utf-8')
for name in ['TMS_Global_Workspace_Publication.ipynb','TMS_Fixed_Workspace_Reanalysis.ipynb','TMS_Exploratory_Findings.ipynb','TMS_Targeted_Information_Flow.ipynb']:
    path=ROOT/'output/jupyter-notebook'/name
    if not path.exists():continue
    nb=nbformat.read(path,as_version=4)
    notice='Historical results — Aachen treatment effects superseded\n\nUse TMS_Corrected_Analyses.ipynb and dashboard/index.html for the corrected 26 November results. SUDMEX results remain valid. Historical Aachen cells below are preserved for traceability and must not be used as current findings.'
    if nb.cells and nb.cells[0].metadata.get('corrected_release_notice'):nb.cells[0].source=notice
    else:nb.cells.insert(0,nbformat.v4.new_markdown_cell(notice,metadata={'corrected_release_notice':True}))
    nbformat.write(nb,path)
status=json.loads((RUN/'status.json').read_text());status.update(target_flow='complete; fresh MATLAB sentinel check unavailable; independent and production-output checks passed',haemodynamic_followup='complete',notebook='executed, zero errors')
status.setdefault('manuscript','updated; render review pending')
status.setdefault('dashboard','in progress')
(RUN/'status.json').write_text(json.dumps(status,indent=2)+'\n',encoding='utf-8')
print('Current pointers and historical notices refreshed')
