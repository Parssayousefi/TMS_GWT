"""Save the bounded spatial follow-up interpretation and extend the review notebook."""
from pathlib import Path
import json
import nbformat

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/aachen_spatial_followup'/(ROOT/'output/aachen_spatial_followup/LATEST.txt').read_text().strip()
assert json.loads((OUT/'followup_verification.json').read_text())['status']=='passed'
report='''# Aachen spatial follow-up and short-report interpretation

The outside-workspace Φ increase is well supported within the original deconvolved-signal analysis, but it is not a demonstrated preferential outside-workspace effect or a preprocessing-independent neural effect.

## What strengthened the finding

- Raw outside-workspace Φ increased in 13/15 participants; normalized Φ increased in 12/15.
- Omitting any participant retained a positive raw mean (0.002366–0.003095 nats), with all nominal p≤0.00505.
- Removing all edges touching each canonical network in turn retained the increase (eight-test maxT p≤0.03791 raw, ≤0.03418 normalized). No single network is necessary for the aggregate result.
- Baseline-plus-motion adjustment retained outside-workspace Φ increases: raw +0.002893 nats and normalized +1.191 percentile points. Holm correction over the six targeted adjustment models gives p=0.03562 and p=0.01970. These are diagnostic families; the earlier five-information-measure and 20/40-endpoint results remain visible.
- Using the first 294 time points in every deconvolved scan retained the outside-workspace estimate (+0.002848 nats). The one-volume scan-length differences are not a plausible explanation of the original magnitude.

## What limits the claim

- The direct outside-minus-inside contrast remains uncertain after correction. Even baseline-plus-motion models give six-model Holm p=0.09102 raw and p=0.07492 normalized. We cannot infer spatial selectivity from one region set being significant and another not.
- The outside-workspace estimate is larger at 90% rMT (+0.003722 nats) than 120% (+0.001896), but neither a corrected intensity contrast nor a first-versus-second-visit contrast establishes a difference. Individual visit changes correlate negatively (raw r=−0.461; normalized r=−0.544). Averaging is not independent replication.
- No individual network-pair block survives correction across the 35 estimable blocks. The complement includes visual, sensorimotor and attention regions, but also some DMN and FPN regions; it is not synonymous with a sensory system.
- Most importantly, the increase is not reproduced in denoised BOLD before HRF deconvolution. Whole-brain Φ changes by −0.000284 nats (nominal p=0.685) and outside-workspace Φ by +0.000200 (p=0.833). The direct deconvolution-by-time difference is corrected-significant for whole-brain Φ (maxT p=0.03607), but not outside-workspace Φ (p=0.07629). Therefore do not claim all spatial effects differ significantly between processing variants.

## What the information decomposition says

For pairs outside the workspace, the estimated group-mean changes are:

| Quantity | Change in nats |
|---|---:|
| Joint past-to-future information | −0.008918 |
| Sum of each region's own past-to-future information | −0.011727 |
| Φ = joint minus summed individual information | +0.002809 |

Both ingredients decrease in the group means; the individual-information sum decreases more. Thus the Φ increase is a relative shift in temporal dependence, not an increase in total joint information. None of the component-wise tests survives its nine-test family correction, so this identity should not be presented as proof of a biological mechanism. Normalized Φ remains positive, but normalization does not validate the haemodynamic model.

Neither the nine mean HRF-parameter tests nor six motion-adjusted HRF–Φ associations survives correction. This does not rule out regional, nonlinear or estimation effects of deconvolution. The remaining alternatives include neural dynamics revealed by deconvolution, haemodynamic/model sensitivity, and nonspecific pre/post state changes. A sham condition is absent.

## Suggested short report paragraph

Following left DLPFC iTBS, mean pairwise whole-minus-sum information increased in haemodynamically deconvolved resting-state signals. The clearest evidence concerned connections outside an independently defined global workspace, whereas workspace integration and role organization showed no corrected changes. The increase reflected a shift in the balance between joint and within-region temporal predictability, and persisted across participant, network, scan-length and motion checks. However, spatial selectivity was not established, visit-level responses were inconsistent, and the increase was not reproduced in denoised BOLD before deconvolution. The findings therefore identify an exploratory, model-dependent change in temporal information following prefrontal stimulation.

## Relation to TMS biology and earlier work

DLPFC stimulation can affect connected regions beyond its target; the original acquisition study reported intensity-dependent frontostriatal connectivity changes. This makes a distributed response plausible, but does not establish a DLPFC-specific causal mechanism in a pre/post dataset. A later entropy reanalysis of these same scans reported striatal effects; it is related evidence from the same data, not independent replication. Avoid claims of increased consciousness, enhanced cognitive control, increased synergy, therapeutic benefit or a validated psychedelic-like state.

Sources: [Alkhasli et al. 2019](https://www.frontiersin.org/journals/human-neuroscience/articles/10.3389/fnhum.2019.00190/full); [Liu et al. 2025](https://doi.org/10.1016/j.neurot.2025.e00556).

## Reproducibility

`setup/deepen_aachen_spatial_effect.py` performs the bounded follow-up and retains all tests. `setup/verify_aachen_spatial_followup.py` independently checks the covariance calculations, weighted network partition, direct method interaction sign flips and immutable hashes. Independent log-determinant calculations agreed within 1.3×10⁻¹⁴ nats. No production preprocessing or frozen HCP/psychedelic analysis was rerun.
'''
(ROOT/'docs/AACHEN_SPATIAL_FOLLOWUP_20261125.md').write_text(report,encoding='utf-8')
(OUT/'README.md').write_text(report,encoding='utf-8')
path=ROOT/'output/jupyter-notebook/TMS_Fixed_Workspace_Reanalysis.ipynb'
nb=nbformat.read(path,as_version=4)
nb.cells=[c for c in nb.cells if not c.metadata.get('spatial_followup')]
md=nbformat.v4.new_markdown_cell;code=nbformat.v4.new_code_cell
cells=[md('## Deeper spatial follow-up\n\nThe outside-workspace effect is distributed across participants and remains after each network omission, but does not establish preferential modulation outside the workspace. The key additional limitation is signal representation: the increase is not reproduced before haemodynamic deconvolution.'),
 code("FOLLOW = ROOT / 'output/aachen_spatial_followup' / '"+OUT.name+"'\nassert json.loads((FOLLOW / 'followup_verification.json').read_text())['status'] == 'passed'\ndisplay(pd.read_csv(FOLLOW / 'estimator_sensitivity_results.csv')[['metric','estimate','ci_low','ci_high','p','p_maxT_family']].round(6))\ndisplay(pd.read_csv(FOLLOW / 'method_interactions.csv')[['metric','estimate','p','p_maxT_family']].round(6))"),
 md('### Temporal-information balance\n\nOutside-workspace joint information decreased by 0.008918 nats, while summed within-region self-information decreased by 0.011727. Their difference is the +0.002809 Φ change. This describes group means; the separate components do not survive their nine-test correction.'),
 code("display(pd.read_csv(FOLLOW / 'temporal_component_results.csv').query(\"branch == 'deconvolved'\")[[ 'metric','estimate','ci_low','ci_high','p_maxT_family']].round(6))\ndisplay(pd.read_csv(FOLLOW / 'nonparametric_checks.csv').round(6))"),
 md('### Scope of interpretation\n\nThe short report should describe an exploratory, model-dependent change in temporal information following iTBS. The data do not establish increased synergy, workspace integration, consciousness, or clinical efficacy. The complete follow-up is in `docs/AACHEN_SPATIAL_FOLLOWUP_20261125.md`.')]
for cell in cells:cell.metadata['spatial_followup']=True
nb.cells.extend(cells);nbformat.validate(nb);nbformat.write(nb,path)
print('Saved follow-up report and notebook section')
