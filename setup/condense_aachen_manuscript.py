"""Condense the canonical manuscript, preserving complete diagnostics in supplement."""
from pathlib import Path
from datetime import datetime,timezone
import re
import shutil

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'manuscript/aachen_brief_report'

ABSTRACT='''Prefrontal intermittent theta burst stimulation (iTBS) can alter distributed brain activity, but its effects on temporal information are less clear. We applied integrated information decomposition to resting-state fMRI from 15 healthy adults scanned before and after two active iTBS intensities. Five measures were calculated from 454 haemodynamically deconvolved regional signals, across the whole brain and within an independent 235-region workspace. Matched visit changes were averaged within participant. Whole-brain mean pairwise whole-minus-sum information, Φ, increased by 0.00179 nats (95% confidence interval 0.000515–0.003070; 8.5% of baseline), surviving maximum-t correction across five information measures (p=0.0239) and both spatial scopes (p=0.0429). No within-workspace measure survived correction. Conditional-null normalization supported the whole-brain increase within the five-measure family but not across both scopes. A spatial follow-up identified increased raw and normalized Φ among non-workspace pairs, although the direct outside-versus-inside contrast did not survive correction. Participant-influence checks retained positive estimates. These exploratory results distinguish a global pairwise information change from increased workspace integration. The small sample, absence of sham and sensitivity to analytical assumptions limit causal interpretation.'''

BODY='''## Introduction

Prefrontal intermittent theta burst stimulation (iTBS) can change relationships between distributed brain regions [@alkhasli2019]. Functional connectivity describes statistical dependence, whereas information decomposition distinguishes information retained by individual regions from information available through their joint dynamics. This distinction may clarify what changes after stimulation.

Whole-minus-sum information, Φ, compares joint past–future information with the summed past–future information of individual components [@barrett2011]. Integrated information decomposition (ΦID) further separates redundant, unique and synergistic contributions [@mediano2025]. Φ is therefore not equivalent to synergy alone. Previous fMRI studies identified a synergistic workspace across conscious states [@luppi2024] and examined workspace organization under psychedelics [@yousefi2026]. They motivate the present analysis without specifying an expected direction for iTBS.

We asked whether prefrontal iTBS was followed by a change in mean pairwise Φ, and whether that change extended to an independently defined global workspace. We retained five information measures, examined conditional-null normalization and tested robustness to participant influence, motion and spatial aggregation. The dataset contains two active stimulation visits and no sham condition; it can characterize pre/post changes but cannot isolate their cause.

## Methods

### Dataset and processing

We reanalysed OpenNeuro ds001832, version 1.0.1, collected by Alkhasli and colleagues [@alkhasli2019]. Participants received 600-pulse left dorsolateral prefrontal cortex iTBS at 90% and 120% of resting motor threshold on separate, counterbalanced visits. Ten-minute resting-state scans were acquired before stimulation and approximately seven minutes afterward; no pulses were delivered during fMRI. One of 16 participants was excluded because the deposited slice-timing metadata did not match acquisition dimensions. The final sample comprised 15 healthy adults, with both complete pre/post pairs and 60 retained scans.

We used the established GWT processing framework [@yousefi2026], including fMRIPrep, CONN denoising, parcellation into 400 Schaefer cortical and 54 Tian subcortical regions, and blind haemodynamic deconvolution with rsHRF [@wu2013]. Gaussian minimum-mutual-information ΦID followed the implementation used in prior workspace analyses [@luppi2024; @mediano2025], at one repetition time (2 seconds). Acquisition, denoising, deconvolution settings and numerical checks are provided in Supplementary Methods.

### Measures and spatial summaries

For each regional pair X,Y, we calculated:

ΦXY = I[(Xt,Yt);(Xt+τ,Yt+τ)] − I(Xt;Xt+τ) − I(Yt;Yt+τ).

Natural logarithms give units of nats. The first term measures joint temporal predictability; the remaining terms measure each region's own temporal predictability. Φ can be negative because the separate terms can count redundant information twice. Redundancy-corrected ΦR adds persistent redundancy. We also retained persistent synergy, persistent redundancy and total source-side synergy, forming one family of five information measures.

Whole-brain means included 102,831 unordered pairs, excluding self-pairs. Workspace means included the 27,495 pairs with both endpoints in the unchanged 235-region HCP mask used previously [@yousefi2026]. A spatial follow-up separated within-workspace, outside-workspace and cross-boundary pairs. These are means of pairwise quantities, not whole-system integrated information. The original 20 participation-coefficient and gateway/broadcaster outcomes were also retained (Supplementary Table S1). Scalar workspace averages differ from Luppi's edgewise network-based analysis [@luppi2024].

### Normalization and inference

A NuMIT-inspired sensitivity compared each observed measure with a model null matched to joint time-delayed mutual information [@liardi2025; @blackburne2026]. This ΦID adaptation used a shared bivariate autoregressive reference learned from both baseline visits, with two independently generated null sets. Conditional percentiles were averaged across references and pairs. Changes are reported in percentile points; normalized components are not additive. Supplementary Methods describe generation and calibration.

We calculated post-minus-pre changes separately at each intensity and averaged the two within participant, retaining 15 independent observations. Signals and covariance matrices were not pooled across visits. Exploratory inference used two-sided single-step maximum-|t| correction over all 32,768 shared participant sign patterns [@winkler2014]. The five information measures formed a family at each scope; joint ten-test correction covered whole-brain and workspace outcomes. The four-scope follow-up used 20 tests. Wider sensitivities included the 20 participation/role endpoints, producing 30- and 40-test families. These tests assume jointly sign-symmetric null errors and are not sham-randomization tests.

For spatial specificity, we directly tested outside-minus-inside changes, with the five information measures forming one family. We report uncorrected and corrected p values for Φ. Raw and normalized analyses are dependent sensitivities. Holm corrections [@holm1979], marginal 95% confidence intervals, visit-specific estimates, participant omission, equal record lengths, motion-adjusted HC3 models and signal-representation checks appear in the supplement. Reported families do not encompass every previous exploratory analytical choice.

## Results

### Whole-brain change and workspace comparison

Whole-brain Φ increased from 0.02100 to 0.02279 nats. The mean change was 0.001792 nats (95% CI 0.000515–0.003070), corresponding to 8.5% of baseline and paired Hedges' gz=0.735. Eleven of 15 participants showed an increase (Figure 1). Maximum-t p values were 0.02393 across five measures and 0.04291 across both scopes (Table 1). Whole-brain ΦR increased by 0.002220 nats, but did not survive five-measure correction (p=0.06580); persistent synergy had a negative point estimate.

Within-workspace Φ increased by 0.001299 nats (95% CI −0.000293–0.002892; five-measure p=0.28406), and ΦR by 0.001981 (0.000173–0.003790; corrected p=0.10883). No workspace information or role outcome survived its correction. Normalized whole-brain Φ increased by 0.712 percentile points (0.163–1.262), with maximum-t p=0.04828 across five measures and 0.07825 across both scopes. No normalized workspace measure survived correction (Table 2).

### Spatial follow-up and robustness

Among pairs outside the workspace, Φ increased by 0.002809 nats (95% CI 0.001198–0.004420), or 1.151 conditional percentile points (0.513–1.790). Both survived correction across 20 scope/measure combinations (raw p=0.00946; normalized p=0.01270) and the 40-endpoint sensitivity (p=0.03851 and 0.03479; Figure 2).

The direct outside-minus-inside Φ change was 0.001510 nats. Its uncorrected sign-flip p value was 0.03790, increasing to 0.16522 within the five-measure spatial-contrast family. The normalized contrast gave p=0.03723 uncorrected and 0.16809 corrected. Thus, separate significance outside the workspace does not establish a larger effect there.

Every leave-one-participant-out whole-brain estimate remained positive (0.001392–0.002030 nats). Both intensity-specific point estimates were positive, but their intervals included zero and individual visit changes correlated negatively (r=−0.346). Motion adjustment barely changed the whole-brain estimate (+0.001787 nats), although its five-measure Holm p value was 0.16647. The wider 30-endpoint maximum-t p value was 0.15936. Equal-length signals retained the positive outside-workspace estimate. The whole-brain Φ change was sensitive to signal representation (Supplementary Results; Table S9).

Outside-workspace Φ increased in 13 of 15 participants. Participant and network omission checks retained positive aggregate estimates, but no individual network-pair block survived correction. Component estimates indicated that summed self-information decreased more than joint temporal information; their difference yields the positive Φ change. Neither component survived its nine-test correction. Full influence, intensity, motion and component results are retained in Supplementary Tables S2–S10.

## Discussion

Mean pairwise Φ increased after prefrontal iTBS when matched changes from two active visits were averaged. This finding did not extend to corrected evidence for greater integration or altered roles within the independent workspace. Outside-workspace pairs showed a clearer aggregate change, but the direct comparison did not establish spatial specificity. The supported interpretation is a change in temporal information balance following stimulation.

An increase in Φ need not mean that the brain contains more joint information or more persistent synergy. Φ subtracts summed regional self-information from joint temporal information, so it can increase when the former decreases more than the latter [@mediano2025]. Our point estimates exhibited that pattern, but component-wise uncertainty prevents a specific physiological explanation. The current result also differs from changes in workspace participation and gateway/broadcaster organization under psychedelics [@yousefi2026], and from Luppi's edgewise ΦR findings [@luppi2024]. Shared methods facilitate comparison without making these outcomes interchangeable.

The original study reported frontostriatal connectivity changes [@alkhasli2019], consistent with distributed consequences of prefrontal stimulation. Our aggregate Φ finding does not identify the stimulated parcels as its origin or establish directed propagation from DLPFC. It is also not a measurement of consciousness or subjective experience. The healthy sample and absence of subjective or therapeutic endpoints preclude conclusions about antidepressant benefit.

The estimate was not attributable to a single participant, and mean-motion adjustment changed its magnitude little. Its support nevertheless depended on correction scope and model uncertainty. Mean framewise displacement cannot capture all motion-related contamination [@power2012]. Without sham, repeated scanning, arousal, expectancy and physiological changes remain plausible contributors. Negative correlation between visit changes further limits claims about reproducible individual responses; a nonsignificant intensity difference does not establish equivalence.

Conditional normalization provides a model-based check at matched total information [@liardi2025], rather than a universal bias correction or independent replication. The Gaussian estimator, baseline-derived reference, single atlas, short records and estimated haemodynamic model constrain interpretation. Deconvolution is part of the intended processing chain and does not turn fMRI into a direct neuronal measurement. Detailed processing sensitivities are reported in the supplement.

These exploratory data support a distinction between increased global pairwise Φ and increased integration within a fixed workspace. A sham-controlled replication with reliable individual target localization is needed to determine whether this information change is attributable to stimulation and where it arises.'''

def section(text,start,end):return text.split(start,1)[1].split(end,1)[0].strip()

def main():
    marker=FOLDER/'archive/manuscript_before_condensation_20261125.md'
    if not marker.exists():
        shutil.copyfile(FOLDER/'manuscript_source.md',marker)
        shutil.copyfile(FOLDER/'supplement_source.md',FOLDER/'archive/supplement_before_condensation_20261125.md')
    old=marker.read_text(encoding='utf-8')
    maintext=old.split('## Abstract',1)[0]+'## Abstract\n\n'+ABSTRACT+'\n\n'+BODY+'\n\n## Data availability'+old.split('## Data availability',1)[1]
    supplement=(FOLDER/'archive/supplement_before_condensation_20261125.md').read_text(encoding='utf-8')
    processing=section(old,'### Dataset and stimulation design','### Information measures')
    processing='### Acquisition and processing details\n\n'+processing.replace('### Image processing and signal extraction','### Signal processing')
    supplement=supplement.replace('## Additional methods','## Additional methods\n\n'+processing,1)
    extra='''## Supplementary results

### Signal representation and estimated haemodynamics

Deconvolution was the intended preprocessing for the primary analysis [@yousefi2026; @wu2013]. A sensitivity on denoised BOLD is a comparison of signal representations, not a requirement that the same result appear before preprocessing. Deconvolution changes temporal dependencies and targets a different signal representation; agreement is therefore not guaranteed. Conversely, differences alone do not establish which representation better reflects neuronal dynamics.

'''+section(old,'### Temporal decomposition and signal representation','### Participant influence, intensity and motion')+'''

Neither the nine estimated-HRF parameter contrasts nor six motion-adjusted HRF–Φ associations survived correction (Table S10). These checks do not independently validate a neuronal interpretation or exclude nonlinear haemodynamic contributions.

### Participant influence, visits, motion and correction scope

'''+section(old,'### Participant influence, intensity and motion','## Discussion')+'''

### Direct spatial contrast

The correction family contains five outside-minus-inside pre/post contrasts: Φ, ΦR, persistent synergy, persistent redundancy and total source-side synergy. Raw Φ had an uncorrected sign-flip p=0.03790 and five-measure maximum-t p=0.16522; normalized Φ had p=0.03723 and 0.16809, respectively. A single isolated contrast would not inherently require five-test correction, but all five measures were examined here. The raw Φ paired-t p was 0.04999. These are different inferential procedures, not interchangeable labels for the same p value.

'''
    supplement=supplement.replace('## Numerical audit',extra+'## Numerical audit',1)
    (FOLDER/'manuscript_source.md').write_text(maintext,encoding='utf-8')
    (FOLDER/'supplement_source.md').write_text(supplement,encoding='utf-8')
    words=lambda t:len(re.findall(r"\b[\w]+(?:[’'-][\w]+)*\b",t))
    print('Abstract:',words(ABSTRACT),'Body:',words(BODY))
    assert 'denoised BOLD' not in maintext and 'absent before deconvolution' not in maintext

if __name__=='__main__':main()
