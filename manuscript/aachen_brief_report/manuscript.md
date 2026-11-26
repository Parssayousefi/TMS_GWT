**INTERNAL DRAFT: Aachen anatomical reprocessing is under review; results require revalidation before submission.**

# Whole brain and workspace information dynamics following prefrontal theta burst stimulation

Zain Salamor†, Parsa Yousefi¹†*, Hector Taylor², Abigail Calder¹, Gregor Hasler¹

¹ Molecular Psychiatry Laboratory, Faculty of Science and Medicine, University of Fribourg, Villars-sur-Glâne, Switzerland

² Department of Adult Psychiatry and Psychotherapy, University Hospital of Psychiatry Zurich, University of Zurich, Zurich, Switzerland

[TODO: Add Zain Salamor’s institutional affiliation and corresponding superscript.]

† These authors contributed equally to this work and share first authorship.

*Correspondence: Parsa Yousefi, parsa.yousefi@unifr.ch

Article type: Brief Research Report

Running title: Information dynamics after prefrontal iTBS

Keywords: intermittent theta burst stimulation; resting state fMRI; integrated information decomposition; prefrontal cortex; neuromodulation; information theory

## Abstract

Prefrontal intermittent theta burst stimulation (iTBS) can alter distributed brain activity, but its effects on temporal information are less clear. We applied integrated information decomposition to resting-state fMRI from 15 healthy adults scanned before and after two active iTBS intensities. Five measures were calculated from 454 haemodynamically deconvolved regional signals, across the whole brain and within an independent 235-region workspace. Matched visit changes were averaged within participant. Whole-brain mean pairwise whole-minus-sum information, Φ, increased by 0.00179 nats (95% confidence interval 0.000515–0.003070; 8.5% of baseline), surviving maximum-t correction across five information measures (p=0.0239) and both spatial scopes (p=0.0429). No within-workspace measure survived correction. Conditional-null normalization supported the whole-brain increase within the five-measure family but not across both scopes. A spatial follow-up identified increased raw and normalized Φ among non-workspace pairs, although the direct outside-versus-inside contrast did not survive correction. Participant-influence checks retained positive estimates. These exploratory results distinguish a global pairwise information change from increased workspace integration. The small sample, absence of sham and sensitivity to analytical assumptions limit causal interpretation.

## Introduction

Prefrontal intermittent theta burst stimulation (iTBS) can change relationships between distributed brain regions (1). Functional connectivity describes statistical dependence, whereas information decomposition distinguishes information retained by individual regions from information available through their joint dynamics. This distinction may clarify what changes after stimulation.

Whole-minus-sum information, Φ, compares joint past–future information with the summed past–future information of individual components (2). Integrated information decomposition (ΦID) further separates redundant, unique and synergistic contributions (3). Φ is therefore not equivalent to synergy alone. Previous fMRI studies identified a synergistic workspace across conscious states (4) and examined workspace organization under psychedelics (5). They motivate the present analysis without specifying an expected direction for iTBS.

We asked whether prefrontal iTBS was followed by a change in mean pairwise Φ, and whether that change extended to an independently defined global workspace. We retained five information measures, examined conditional-null normalization and tested robustness to participant influence, motion and spatial aggregation. The dataset contains two active stimulation visits and no sham condition; it can characterize pre/post changes but cannot isolate their cause.

## Methods

### Dataset and processing

We reanalysed OpenNeuro ds001832, version 1.0.1, collected by Alkhasli and colleagues (1). Participants received 600-pulse left dorsolateral prefrontal cortex iTBS at 90% and 120% of resting motor threshold on separate, counterbalanced visits. Ten-minute resting-state scans were acquired before stimulation and approximately seven minutes afterward; no pulses were delivered during fMRI. One of 16 participants was excluded because the deposited slice-timing metadata did not match acquisition dimensions. The final sample comprised 15 healthy adults, with both complete pre/post pairs and 60 retained scans.

We used the established GWT processing framework (5), including fMRIPrep, CONN denoising, parcellation into 400 Schaefer cortical and 54 Tian subcortical regions, and blind haemodynamic deconvolution with rsHRF (6). Gaussian minimum-mutual-information ΦID followed the implementation used in prior workspace analyses (4, 3), at one repetition time (2 seconds). Acquisition, denoising, deconvolution settings and numerical checks are provided in Supplementary Methods.

### Measures and spatial summaries

For each regional pair X,Y, we calculated:

ΦXY = I[(Xt,Yt);(Xt+τ,Yt+τ)] − I(Xt;Xt+τ) − I(Yt;Yt+τ).

Natural logarithms give units of nats. The first term measures joint temporal predictability; the remaining terms measure each region's own temporal predictability. Φ can be negative because the separate terms can count redundant information twice. Redundancy-corrected ΦR adds persistent redundancy. We also retained persistent synergy, persistent redundancy and total source-side synergy, forming one family of five information measures.

Whole-brain means included 102,831 unordered pairs, excluding self-pairs. Workspace means included the 27,495 pairs with both endpoints in the unchanged 235-region HCP mask used previously (5). A spatial follow-up separated within-workspace, outside-workspace and cross-boundary pairs. These are means of pairwise quantities, not whole-system integrated information. The original 20 participation-coefficient and gateway/broadcaster outcomes were also retained (Supplementary Table S1). Scalar workspace averages differ from Luppi's edgewise network-based analysis (4).

### Normalization and inference

A NuMIT-inspired sensitivity compared each observed measure with a model null matched to joint time-delayed mutual information (7, 8). This ΦID adaptation used a shared bivariate autoregressive reference learned from both baseline visits, with two independently generated null sets. Conditional percentiles were averaged across references and pairs. Changes are reported in percentile points; normalized components are not additive. Supplementary Methods describe generation and calibration.

We calculated post-minus-pre changes separately at each intensity and averaged the two within participant, retaining 15 independent observations. Signals and covariance matrices were not pooled across visits. Exploratory inference used two-sided single-step maximum-|t| correction over all 32,768 shared participant sign patterns (9). The five information measures formed a family at each scope; joint ten-test correction covered whole-brain and workspace outcomes. The four-scope follow-up used 20 tests. Wider sensitivities included the 20 participation/role endpoints, producing 30- and 40-test families. These tests assume jointly sign-symmetric null errors and are not sham-randomization tests.

For spatial specificity, we directly tested outside-minus-inside changes, with the five information measures forming one family. We report uncorrected and corrected p values for Φ. Raw and normalized analyses are dependent sensitivities. Holm corrections (10), marginal 95% confidence intervals, visit-specific estimates, participant omission, equal record lengths, motion-adjusted HC3 models and signal-representation checks appear in the supplement. Reported families do not encompass every previous exploratory analytical choice.

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

An increase in Φ need not mean that the brain contains more joint information or more persistent synergy. Φ subtracts summed regional self-information from joint temporal information, so it can increase when the former decreases more than the latter (3). Our point estimates exhibited that pattern, but component-wise uncertainty prevents a specific physiological explanation. The current result also differs from changes in workspace participation and gateway/broadcaster organization under psychedelics (5), and from Luppi's edgewise ΦR findings (4). Shared methods facilitate comparison without making these outcomes interchangeable.

The original study reported frontostriatal connectivity changes (1), consistent with distributed consequences of prefrontal stimulation. Our aggregate Φ finding does not identify the stimulated parcels as its origin or establish directed propagation from DLPFC. It is also not a measurement of consciousness or subjective experience. The healthy sample and absence of subjective or therapeutic endpoints preclude conclusions about antidepressant benefit.

The estimate was not attributable to a single participant, and mean-motion adjustment changed its magnitude little. Its support nevertheless depended on correction scope and model uncertainty. Mean framewise displacement cannot capture all motion-related contamination (11). Without sham, repeated scanning, arousal, expectancy and physiological changes remain plausible contributors. Negative correlation between visit changes further limits claims about reproducible individual responses; a nonsignificant intensity difference does not establish equivalence.

Conditional normalization provides a model-based check at matched total information (7), rather than a universal bias correction or independent replication. The Gaussian estimator, baseline-derived reference, single atlas, short records and estimated haemodynamic model constrain interpretation. Deconvolution is part of the intended processing chain and does not turn fMRI into a direct neuronal measurement. Detailed processing sensitivities are reported in the supplement.

These exploratory data support a distinction between increased global pairwise Φ and increased integration within a fixed workspace. A sham-controlled replication with reliable individual target localization is needed to determine whether this information change is attributable to stimulation and where it arises.

## Data availability

The original imaging dataset is publicly available under CC0 at OpenNeuro, version 1.0.1, doi:10.18112/openneuro.ds001832.v1.0.1. Participant-level scalar outputs, complete endpoint tables and numerical sensitivity results accompany this report as Supplementary Data. No private psychedelic participant data or effect maps were included in the present analysis. [TODO: Insert the persistent public repository URL or DOI for the analysis code and submission materials before submission.]

## Ethics statement

The original acquisition study received approval from the ethics committee of the Faculty of Medicine, RWTH Aachen University, Germany, protocol 003/15, and all participants provided written informed consent, as reported by Alkhasli and colleagues (1). The present study is a secondary analysis of publicly released, deidentified data and involved no new recruitment or participant contact.

## Author contributions

Zain Salamor and Parsa Yousefi share first authorship. Gregor Hasler served as principal investigator and supervisor. [TODO: Confirm and complete the CRediT contributions for all authors; obtain each author’s approval of the manuscript and author order.]

## Funding

[TODO: Insert the funding statement for this secondary analysis, including grant numbers where applicable.]

## Conflict of interest

[TODO: Obtain and insert declarations from all five authors.]

## Acknowledgments

We thank Isabel Alkhasli, Katrin Sakreida, Felix M. Mottaghy and Ferdinand Binkofski for making the imaging dataset available, and the developers of the analysis tools used in this work.

## Generative AI statement

OpenAI Codex assisted with analysis-code preparation, reference checking and manuscript drafting. Numerical claims were checked against the saved analysis outputs. The authors are responsible for the final content.

## Figures and tables

![Figure 1](figures/Figure_1_participant_changes.png)

Figure 1. Participant-level Φ before and after iTBS. Whole-brain means (A, C) and within-workspace means (B, D) are shown in raw nats (A, B) and conditional-null percentiles (C, D). Each line links one participant's average across the two baseline scans to their average across the two post scans; red diamonds show group means. There are 15 independent participants. Conditional percentiles are not nats and do not preserve additive information identities.

![Figure 2](figures/Figure_2_spatial_scope.png)

Figure 2. Spatial scope of pairwise information changes. Group means and marginal 95% confidence intervals are shown for raw Φ (A), raw ΦR (B), normalized Φ (C), and normalized ΦR (D). Whole-brain means include all three disjoint blocks below them: within-workspace, outside-workspace and workspace-to-outside pairs. Each point uses the same 15 participant-averaged changes. Confidence intervals do not incorporate multiplicity correction; corrected tests appear in the text and supplementary tables. The direct outside-minus-within contrast does not survive correction.

Table 1. Raw information changes across the whole brain and within the fixed HCP workspace.

| scope | measure | estimate | ci_low | ci_high | p_maxT_family | p_maxT_whole_and_workspace10 |
| --- | --- | --- | --- | --- | --- | --- |
| Whole brain | Persistent synergy | -0.0076121321286168 | -0.0154097659070605 | 0.0001855016498268 | 0.157470703125 | 0.2337646484375 |
| Whole brain | Persistent redundancy | 0.0004282253231195 | -0.0004702054814039 | 0.0013266561276431 | 0.67926025390625 | 0.80560302734375 |
| Whole brain | Total source synergy | -0.0020708582657842 | -0.0052049733426909 | 0.0010632568111224 | 0.43829345703125 | 0.5694580078125 |
| Whole brain | Φ | 0.0017921716024821 | 0.0005148332551216 | 0.0030695099498425 | 0.02392578125 | 0.04290771484375 |
| Whole brain | ΦR | 0.0022203969256016 | 0.0003560113155407 | 0.0040847825356626 | 0.0657958984375 | 0.108154296875 |
| Within workspace | Persistent synergy | -0.0037330345682997 | -0.0125339461245755 | 0.005067876987976 | 0.73626708984375 | 0.86529541015625 |
| Within workspace | Persistent redundancy | 0.0006820133409015 | 6.969037475389173e-05 | 0.0012943363070491 | 0.0997314453125 | 0.14581298828125 |
| Within workspace | Total source synergy | -0.000136913642021 | -0.0038642211997759 | 0.0035903939157339 | 0.999755859375 | 1.0 |
| Within workspace | Φ | 0.0012992932715628 | -0.0002931174816175 | 0.0028917040247432 | 0.2840576171875 | 0.38763427734375 |
| Within workspace | ΦR | 0.0019813066124643 | 0.0001726109161415 | 0.0037900023087872 | 0.10882568359375 | 0.15740966796875 |

Estimates and marginal 95% confidence intervals are in nats. Maximum-t p5 corrects the five information measures within each scope; p10 corrects all ten rows jointly. Each row uses the same 15 participant-averaged changes. Persistent synergy and total source-side synergy are distinct quantities. Holm corrections, all 30 endpoints, and further spatial analyses are retained in the supplementary data.

Table 2. Conditional-null-normalized information changes across the two scopes.

| scope | measure | estimate | ci_low | ci_high | p_maxT_family | p_maxT_whole_and_workspace10 |
| --- | --- | --- | --- | --- | --- | --- |
| Whole brain | Persistent synergy | -0.48462559575609004 | -1.0888520744071 | 0.11960088289491 | 0.31280517578125 | 0.40814208984375 |
| Whole brain | Persistent redundancy | 0.5380342258395701 | -0.33593131931477 | 1.41199977099393 | 0.51690673828125 | 0.6312255859375 |
| Whole brain | Total source synergy | 0.21545565339088002 | -0.88606966630473 | 1.31698097308651 | 0.9681396484375 | 0.99407958984375 |
| Whole brain | Φ | 0.71230876840225 | 0.16296142750858 | 1.26165610929591 | 0.04827880859375 | 0.0782470703125 |
| Whole brain | ΦR | 0.79803423698238 | 0.12741696849453 | 1.46865150547022 | 0.076904296875 | 0.118408203125 |
| Within workspace | Persistent synergy | -0.42789283656309 | -1.01792007000916 | 0.16213439688296 | 0.3919677734375 | 0.49444580078125 |
| Within workspace | Persistent redundancy | 0.70387368547887 | 0.09582749360691 | 1.31191987735083 | 0.09674072265625 | 0.134521484375 |
| Within workspace | Total source synergy | 0.22646265867030999 | -1.07217449894249 | 1.52509981628313 | 0.9747314453125 | 0.9964599609375 |
| Within workspace | Φ | 0.44664291246960003 | -0.27948387459616 | 1.17276969953537 | 0.51947021484375 | 0.632080078125 |
| Within workspace | ΦR | 0.65828536514508 | -0.02402874976881 | 1.34059948005898 | 0.188232421875 | 0.25518798828125 |

Estimates and marginal 95% confidence intervals are in percentile points. Maximum-t p5 and p10 have the same family definitions as Table 1. Normalization is a model-dependent sensitivity of the same scans, and normalized measures are not additive information components.

## References

1. Alkhasli I, Sakreida K, Mottaghy FM, Binkofski F. Modulation of Fronto-Striatal Functional Connectivity Using Transcranial Magnetic Stimulation. Frontiers in Human Neuroscience. 2019;13:190. doi:10.3389/fnhum.2019.00190.

2. Barrett AB, Seth AK. Practical Measures of Integrated Information for Time-Series Data. PLoS Computational Biology. 2011;7:e1001052. doi:10.1371/journal.pcbi.1001052.

3. Mediano PAM, Rosas FE, Luppi AI, Carhart-Harris RL, Bor D, Seth AK, et al. Toward a unified taxonomy of information dynamics via Integrated Information Decomposition. Proceedings of the National Academy of Sciences. 2025;122:e2423297122. doi:10.1073/pnas.2423297122.

4. Luppi AI, Mediano PA, Rosas FE, Allanson J, Pickard J, Carhart-Harris RL, et al. A synergistic workspace for human consciousness revealed by Integrated Information Decomposition. eLife. 2024;12:RP88173. doi:10.7554/elife.88173.

5. Yousefi P, Shinozuka K, Vizeli P, Holze F, Becker AM, Muller F, et al. Psychedelics Restructure the Global Workspace of Consciousness. bioRxiv [Preprint]. 2026:version 1. doi:10.64898/2026.09.28.755172v1.

6. Wu GR, Liao W, Stramaglia S, Ding JR, Chen H, Marinazzo D. A blind deconvolution approach to recover effective connectivity brain networks from resting state fMRI data. Medical Image Analysis. 2013;17:365-374. doi:10.1016/j.media.2013.01.003.

7. Liardi A, Rosas FE, Carhart-Harris RL, Blackburne G, Bor D, Mediano PAM. Null models for comparing information decomposition across complex systems. PLOS Computational Biology. 2025;21:e1013629. doi:10.1371/journal.pcbi.1013629.

8. Blackburne G, Liardi A, Mallaroni P, McAlpine R, Carhart-Harris R, Timmermann C, et al. Psychedelics restructure information in the human brain. bioRxiv [Preprint]. 2026:version 2. doi:10.64898/2026.10.01.755915v2.

9. Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE. Permutation inference for the general linear model. NeuroImage. 2014;92:381-397. doi:10.1016/j.neuroimage.2014.01.060.

10. Holm S. A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics. 1979;6:65–70.

11. Power JD, Barnes KA, Snyder AZ, Schlaggar BL, Petersen SE. Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage. 2012;59:2142-2154. doi:10.1016/j.neuroimage.2011.10.018.
