"""Reframe the Aachen draft after the fixed-workspace audit; archive its precursor."""
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'manuscript/aachen_brief_report'
path=FOLDER/'manuscript_source.md'
archive=FOLDER/'archive/manuscript_before_workspace_audit.md'
archive.parent.mkdir(exist_ok=True)
if not archive.exists():archive.write_bytes(path.read_bytes())
text=archive.read_text(encoding='utf-8')


def replace_section(start,end,content):
    global text
    first=text.index(start);last=text.index(end,first)
    text=text[:first]+content.strip()+'\n\n'+text[last:]


text=text.replace('# Global information dynamics following prefrontal intermittent theta burst stimulation',
                  '# Whole brain and workspace information dynamics following prefrontal theta burst stimulation',1)
replace_section('## Abstract','## Introduction', '''## Abstract

Prefrontal intermittent theta burst stimulation (iTBS) can alter distributed brain activity, but its effects on temporal information remain incompletely characterized. We reanalysed resting-state fMRI from 15 healthy adults with pre/post scans at two active stimulation intensities. Denoised, haemodynamically deconvolved signals from 454 regions underwent Gaussian integrated information decomposition at a 2-second lag. Five information measures were averaged across all region pairs and within an independently defined 235-region workspace. The two matched visit changes were averaged within participant. Whole-brain whole-minus-sum information, Φ, increased by 0.00179 nats (95% confidence interval 0.000515–0.003070; 8.5% of the baseline mean). Exact maximum-t correction gave p=0.0239 across five measures and p=0.0429 across both spatial scopes. Within the workspace, neither Φ nor redundancy-corrected ΦR survived five-measure correction (p=0.284 and 0.109). Conditional-null-normalized whole-brain Φ increased, but did not survive joint correction across both scopes (p=0.0782). A spatial follow-up identified increases in raw and normalized Φ among pairs outside the workspace, surviving correction across 20 information/scope combinations; the direct outside-versus-inside contrast did not survive correction. Whole-brain Φ remained positive under every leave-one-participant-out analysis, while motion adjustment weakened corrected support. These exploratory findings distinguish a global pairwise information change from evidence of increased workspace integration. The small sample, absence of sham and dependence on inferential scope limit causal and consciousness-related interpretations.''')
text=text.replace('We therefore examined whether mean pairwise Φ changed after left dorsolateral prefrontal cortex (DLPFC) iTBS in an open, repeated-visit fMRI dataset. The report emerged from an exploratory analysis of global and workspace measures and focuses on the visit-averaged Φ finding. We retained the complete endpoint inventory and examined intensity-specific changes, participant influence, motion adjustment and sensitivity to total mutual information. Because both visits involved active stimulation, the analysis estimates pre/post change following iTBS without isolating a stimulation-specific causal effect.',
'''We examined changes in mean pairwise Φ after left dorsolateral prefrontal cortex (DLPFC) iTBS in an open, repeated-visit fMRI dataset, asking whether a whole-brain finding extended to an independently defined workspace. This report emerged from exploratory analysis and a subsequent audit of spatial aggregation. We retained all five information measures and the original workspace-role outcomes, and examined intensity, participant influence, motion and total-information normalization. Both visits involved active stimulation, so the analysis estimates pre/post change without isolating a stimulation-specific causal effect.''')
text=text.replace('Together, these constituted five global information measures.','Together, these constituted one family of five information measures at each spatial scope.')
text=text.replace('For each scan, we averaged each measure across the 102,831 unique unordered pairs, excluding self-pairs. Every scan had complete pair coverage. “Global” thus denotes an equally weighted mean of pairwise quantities; it is neither a 454-region whole-system Φ estimate nor a sum of synergy alone. Secondary workspace participation and role measures used a fixed, independently derived 235-region HCP reference, without re-estimating the workspace from these intervention data. Supplementary Methods describe these outcomes, and Supplementary Table S1 retains all 25 neural endpoints.',
'''For each scan, whole-brain means used all 102,831 unique unordered region pairs, excluding self-pairs. Within-workspace means used the 27,495 pairs with both endpoints in the unchanged 235-region HCP mask used by Yousefi and colleagues [@yousefi2026]. Every scan had complete coverage. A spatial follow-up partitioned the graph into within-workspace, outside-workspace (23,871 pairs), and workspace-to-outside (51,465 pairs) connections. Global Φ is thus a mean of pairwise quantities, not a 454-region whole-system estimate or a sum of synergy alone.

Secondary participation coefficients (PCs) used interactions with all 454 regions, followed by averaging and rank-based gateway/broadcaster classification within the fixed workspace, as in the earlier GWT analysis [@yousefi2026]. The mask was not re-estimated from TMS data. The original 20 PC/role endpoints were retained (Supplementary Table S1). The within-workspace information averages extend those outcomes. They are also distinct from Luppi's edgewise network-based statistical analysis, which assessed the overlap of altered connections with the workspace [@luppi2024].''')
text=text.replace('Normalized changes are reported in percentile points, not nats.','The same pairwise reference was retained across spatial scopes. Normalized changes are reported in percentile points, not nats.')
replace_section('### Participant-level estimation and inference','## Results', '''### Participant-level estimation and inference

We calculated post-minus-pre change at each intensity, then averaged the two changes within participant. This gives equal weight to visits and retains n=15 independent people. Signals and covariance matrices were not pooled across visits. With complete pairs, the resulting change equals the average post-scan metric minus the average baseline metric.

Analyses were exploratory and not preregistered. The focus on Φ, correction framework and spatial follow-up developed during exploration. We used a two-sided single-step maximum-|t| test across all 2¹⁵=32,768 participant sign patterns. Each participant's sign was shared across endpoints, preserving their dependence. The five information measures formed one family at each scope. We additionally evaluated the ten whole-brain/workspace outcomes jointly, and all 30 endpoints when including the original 20 PC/role measures. For the four-scope spatial follow-up, we corrected across 20 information/scope combinations and provided a 40-endpoint sensitivity including PC/role measures. These corrections cover the specified families, not every historical exploratory analysis. Joint sign symmetry of null errors is assumed; this is not a sham-randomization test [@winkler2014].

We report Holm-adjusted parametric p values as an alternative [@holm1979]. Raw and normalized analyses are dependent method sensitivities, not independent replications. All 95% confidence intervals are marginal Student-t intervals, without simultaneous or selection adjustment. Paired Hedges' gz uses the SD of participant changes and a small-sample correction.

Sensitivity analyses included visit-specific estimates, their paired difference, participant influence, signed-rank tests and bootstrap intervals. Motion adjustment regressed participant-averaged change on an intercept and change in mean framewise displacement (FD), using HC3 uncertainty. The intercept estimates change at zero FD change. We report Holm correction over the five measures at each scope and, for spatial diagnostics, over all 20 scope/information combinations. We directly tested outside-minus-within-workspace changes across the five measures, rather than inferring a spatial difference from separate significance tests.''')

# The result and discussion text is supplied by the companion source fragment.
fragment=(FOLDER/'workspace_results_discussion.md').read_text(encoding='utf-8')
replace_section('## Results','## Data availability',fragment)
text=text.replace('Whole-brain Φ remained positive under every leave-one-participant-out analysis, while motion adjustment weakened corrected support.',
                  'Whole-brain Φ remained positive under every leave-one-participant-out analysis, but the increase was not reproduced in denoised BOLD without deconvolution.')
text=text.replace('## Results\n', 'Auxiliary checks decomposed Φ into joint time-delayed information minus summed self-information, examined network contributions, and compared the original deconvolved signals with an equal-length 294-sample variant and denoised BOLD before deconvolution. All retained the Gaussian estimator and 2-second lag. The direct method-by-time contrasts were corrected over eight scope/measure combinations (Φ and ΦR at four scopes); further complete families are described in Supplementary Methods.\n\n## Results\n',1)
replace_section('## Figures and tables','## References', '''## Figures and tables

{{FIGURE1}}

Figure 1. Participant-level Φ before and after iTBS. Whole-brain means (A, C) and within-workspace means (B, D) are shown in raw nats (A, B) and conditional-null percentiles (C, D). Each line links one participant's average across the two baseline scans to their average across the two post scans; red diamonds show group means. There are 15 independent participants. Conditional percentiles are not nats and do not preserve additive information identities.

{{FIGURE2}}

Figure 2. Spatial scope of pairwise information changes. Group means and marginal 95% confidence intervals are shown for raw Φ (A), raw ΦR (B), normalized Φ (C), and normalized ΦR (D). Whole-brain means include all three disjoint blocks below them: within-workspace, outside-workspace and workspace-to-outside pairs. Each point uses the same 15 participant-averaged changes. Confidence intervals do not incorporate multiplicity correction; corrected tests appear in the text and supplementary tables. The direct outside-minus-within contrast does not survive correction.

Table 1. Raw information changes across the whole brain and within the fixed HCP workspace.

{{TABLE1}}

Estimates and marginal 95% confidence intervals are in nats. Maximum-t p5 corrects the five information measures within each scope; p10 corrects all ten rows jointly. Each row uses the same 15 participant-averaged changes. Persistent synergy and total source-side synergy are distinct quantities. Holm corrections, all 30 endpoints, and further spatial analyses are retained in the supplementary data.

Table 2. Conditional-null-normalized information changes across the two scopes.

{{TABLE2}}

Estimates and marginal 95% confidence intervals are in percentile points. Maximum-t p5 and p10 have the same family definitions as Table 1. Normalization is a model-dependent sensitivity of the same scans, and normalized measures are not additive information components.''')
path.write_text(text,encoding='utf-8')
abstract=text.split('## Abstract\n\n')[1].split('\n## Introduction')[0]
body=text.split('## Introduction\n\n')[1].split('\n## Data availability')[0]
print('Abstract words:',len(abstract.split()),'Main body words:',len(body.split()))
assert len(abstract.split())<=250 and len(body.split())<4000
