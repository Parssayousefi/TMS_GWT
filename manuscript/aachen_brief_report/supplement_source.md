# Supplementary material

Whole brain and workspace information dynamics following prefrontal theta burst stimulation

Zain Salamor and Parsa Yousefi contributed equally. Correspondence: Parsa Yousefi, parsa.yousefi@unifr.ch.

## Additional methods

### Acquisition and processing details

We used OpenNeuro dataset ds001832, version 1.0.1 (doi:10.18112/openneuro.ds001832.v1.0.1), originally collected by Alkhasli and colleagues [@alkhasli2019]. The release contains 16 healthy adults. Each underwent left DLPFC iTBS on two separate visits, at 90% and 120% of resting motor threshold (rMT), in counterbalanced order. Stimulation comprised 600 pulses over 3 minutes 20 seconds, using a figure-of-eight coil and MRI-guided neuronavigation to a target corresponding to Talairach coordinates −45, 45, 35. The protocol consisted of 50-Hz triplets repeated at 5 Hz in intermittent trains [@huang2005]. Each visit included a 10-minute resting scan before stimulation and another beginning approximately 7 minutes afterward. Participants fixated with eyes open. Stimulation occurred outside the scanner; no pulses were delivered during fMRI acquisition. There was no sham visit.

We excluded one participant because all four runs contained 36 slices but the deposited timing metadata specified 34 slice times, preventing uniform slice-timing correction. The analysed sample comprised 15 participants, 7 women and 8 men, aged 19–44 years (mean 27.87, SD 7.12), with 60 scans. Seven received 90% rMT first and eight received 120% first. All analysed participants contributed both complete pre/post pairs.

### Signal processing

Functional images were acquired at 3 T with TR = 2 seconds, TE = 28 milliseconds, flip angle 77° and voxel dimensions 3 × 3 × 3.75 mm. Most runs contained 300 volumes; two contained 299. We used fMRIPrep 23.2.1 [@esteban2019] for motion correction, slice-timing correction, anatomical registration and normalization to MNI152NLin2009cAsym at 2-mm resolution. Susceptibility-distortion correction and FreeSurfer surface reconstruction were not performed. The first five volumes were removed.

Run-level quality control used translation and rotation thresholds of 3 mm and 3.15°, with a continuous-segment salvage rule requiring at least 8 minutes. All 60 analysed runs were retained. One run flagged relative to the motion-correction reference passed the segment-relative rule without losing additional volumes. Denoising used CONN 25b [@whitfield2012], including five white-matter and five cerebrospinal-fluid aCompCor components [@behzadi2007], six motion parameters and their first derivatives, linear detrending, and regression followed by 0.008–0.09-Hz filtering. Global-signal regression, volume scrubbing and additional spatial smoothing were not applied.

We spatially averaged the cleaned signal within 400 Schaefer cortical parcels [@schaefer2018] and 54 Tian subcortical parcels [@tian2020]. Regional haemodynamic responses were estimated using the rsHRF implementation of blind deconvolution [@wu2013]: a 1-SD event threshold, third-order gamma basis, 24-second response window, threefold temporal upsampling, and a 4–8-second onset search. Iterative Wiener deconvolution used 50 iterations. Outputs contained 295 time points for 58 scans and 294 for two scans. Numerical audits confirmed finite outputs, 454 nonconstant regional signals per scan, and deconvolved signals that differed from both their inputs and simple standardized inputs. These checks establish processing integrity, not direct validation of recovered neuronal activity.

### Fixed mask and information aggregation

The fixed HCP workspace contains 235 of 454 regions. For a symmetric pairwise matrix M with undefined diagonal, the within-workspace mean is the mean of Mij over i<j with both i and j in that mask. It contains 27,495 unique pairs. The complement contains 219 regions and 23,871 unique within-complement pairs. The cross-boundary block contains 235×219=51,465 pairs. These three disjoint blocks reconstruct the whole-brain mean when weighted by their pair counts; simply averaging their three means would give the wrong weighting. All 60 Aachen scans had finite estimates on every off-diagonal pair.

Workspace membership was inherited from the independent HCP reference used in the earlier GWT pipeline [@yousefi2026]. No HCP participant calculation or psychedelic participant analysis was rerun. TMS processing outputs were read without modification. The scalar within-workspace information mean is an additional analysis, distinct from both the earlier GWT participation/role outcomes and Luppi's edgewise network-based statistics [@luppi2024].

### Secondary participation and role endpoints

For each node, the weighted participation coefficient is PCi=1−Σm(sim/si)², where sim is the sum of edge weights to module m and si is their sum across modules. The seven cortical networks and a single subcortical module define eight modules. Self-edges are excluded. PCs use all 454 destinations; restricting destinations to the workspace would alter the original definition. Synergy PCs use persistent synergy and redundancy PCs use persistent redundancy.

The synergy and redundancy PCs are ranked separately among the 235 workspace nodes for each scan. A positive rank difference defines a gateway and a negative difference a broadcaster; ties are neither. The original 20 secondary endpoints comprise two mean workspace PCs, two total role counts, eight network gateway counts and eight network broadcaster counts. These outcomes retain the historical families: two PCs, two total counts, FPN gateways, the seven other network gateway counts, and eight broadcaster counts. The expanded 30- and 40-endpoint sensitivities avoid relying solely on these separate families. No regional role is inferred from the significance of a group-average matrix.

### Conditional null reference

Following the motivation of NuMIT [@liardi2025], we compared each observed pairwise measure against a model reference at matched joint time-delayed mutual information. This ΦID application is a sensitivity extension; it does not make conditional percentiles additive information atoms [@blackburne2026].

The Aachen coefficient pool used 1,024 randomly selected edges from each of the 30 baseline scans, producing 30,720 bivariate first-order autoregressive coefficient matrices. The sampling seed was 2026112503. The pool was shared across participants and conditions. Random positive-definite innovation covariances and rescaling of stable dynamics generated stationary model covariances at target mutual-information values. Reference generation retained 32,768 valid draws at each 0.005-nat grid point for each of two independent seeds, 2026112511 and 2026112512. Nonfinite draws and draws missing the target information by more than 10⁻⁶ nats were rejected consistently across measures and replenished.

For each empirical value, midrank empirical cumulative probabilities were evaluated in the adjacent reference distributions and linearly interpolated in target mutual information. The two seed estimates were averaged before edge aggregation. The same reference pool and distributions were retained for whole-brain, workspace, outside-workspace and cross-boundary summaries. No scope-specific prior was fitted. Normalized means were multiplied by 100 for reporting as conditional percentiles, and changes as percentile points.

Validation included comparison with the original MATLAB estimates, independent log-determinant calculations, interpolation checks, simulations at finite record lengths and reference-seed comparisons. These establish numerical consistency rather than unbiased estimation or biological validity. Finite-record calibration is retained in Table S7, without claiming that every mean null percentile equals 0.5. Reference-prior uncertainty was not propagated; inference is conditional on the shared baseline-derived reference.

### Sensitivity analyses and correction

All primary changes use the participant as the independent unit. Both visit changes must be available; the analysis does not silently average one visit when the other is missing. The exact maximum-|t| procedure applies each of the 32,768 sign vectors jointly across endpoints and compares observed statistics with the permutation maximum [@winkler2014]. Shared signs preserve dependence between Φ, ΦR and other measures. The reported marginal intervals are not simultaneous intervals.

Holm correction [@holm1979] is provided alongside the maximum-t results. Motion models use the uncentered FD change so that the intercept represents zero motion change; centering FD would instead force the intercept to the observed unadjusted mean. HC3 models also examine intensity differences with baseline differences, FD differences and counterbalanced visit order. Results from different spatial scopes and model variants are retained, not chosen according to which crosses a threshold.

Participant-bootstrap intervals used 10,000 resamples. The historical whole-brain bootstrap used seed 2026112515; the workspace audit used seed 2026112526. Leave-one-out analyses remove one person, including both visits, at a time. Signed-rank tests and bootstrap intervals are dependent sensitivities of the same observations and do not provide new confirmatory evidence. All such intervals and unadjusted p values are labelled accordingly in their tables.

### Bounded spatial and signal-representation follow-up

After examining the fixed-workspace results, we evaluated two targets: the outside-workspace Φ change and the direct outside-minus-inside change. Visit/intensity consistency used 12 jointly corrected contrasts: two targets each at 90%, 120%, their difference, first visit, second visit and their difference. All use 15 independent participants. Baseline and motion sensitivity comprised six models per raw/normalized branch: two targets with baseline, motion, or both. Baseline values were centered at their sample mean; FD change remained uncentered. Participant influence, signed-rank, sign-test and 10,000-resample bootstrap checks were reported as dependent diagnostics.

The outside-workspace graph was partitioned by the eight fixed network labels. Of 36 possible unordered network-pair blocks, 35 contained edges; the subcortical complement contains only one node, so its within-network block is empty. Joint maximum-t correction covered the 35 estimable blocks. A separate eight-test omission sensitivity removed all edges touching one network at a time. Block means reconstruct the outside-workspace mean only when weighted by their edge counts. A positive result after each omission indicates that no single network is necessary for the aggregate finding; it does not establish significance of any individual network block.

We decomposed mean Φ as mean joint time-delayed mutual information minus mean summed within-region self-information. For each scope, these terms were computed on exactly the same pairs. Nine component tests covered three quantities (joint information, self-information sum, persistent redundancy) at three scopes (within, outside, outside-minus-inside). The group-mean identity is exact, whereas the separate component hypotheses remain subject to correction.

Signal-representation sensitivities used the original deconvolved scans, their first 294 samples, and the denoised BOLD scans before deconvolution. The equal-length check retains each scan's previously estimated HRF; it does not rerun deconvolution on truncated inputs. The same Gaussian covariance-based Φ formula and 2-second lag were used throughout. The 12-test family comprises Φ at four scopes (whole brain, within, outside, their difference) across the three representations. A separate direct method-by-time comparison tests deconvolved-minus-denoised changes in Φ and ΦR at the four scopes, with eight-test correction. No normalized null was applied to the alternative signal representations.

Estimated HRF height, time to peak and full width at half maximum were averaged within and outside the mask; their differences formed a nine-test parameter family. Six motion-adjusted association models related the three HRF-parameter changes to outside-workspace Φ or the outside-minus-inside Φ contrast, with Holm correction. These checks cannot exclude regional or nonlinear haemodynamic contributions. Independent log-determinant calculations agreed with the denoised-BOLD temporal information calculation to within 1.3×10⁻¹⁴ nats.

## Supplementary results

### Signal representation and estimated haemodynamics

Deconvolution was the intended preprocessing for the primary analysis [@yousefi2026; @wu2013]. A sensitivity on denoised BOLD is a comparison of signal representations, not a requirement that the same result appear before preprocessing. Deconvolution changes temporal dependencies and targets a different signal representation; agreement is therefore not guaranteed. Conversely, differences alone do not establish which representation better reflects neuronal dynamics.

The outside-workspace Φ increase did not correspond to an increase in total joint temporal information. At the group-mean level, joint time-delayed mutual information decreased by 0.008918 nats, while the summed within-region self-information decreased by 0.011727 nats. Their difference is the observed +0.002809-nat Φ change. The individual component changes did not survive the nine-test mechanistic correction, so this identity describes the estimated balance rather than establishing a separately significant physiological mechanism.

Restricting every deconvolved scan to 294 samples preserved the outside-workspace Φ estimate (+0.002848 nats). Applying the same Gaussian Φ formula to denoised BOLD before deconvolution gave an outside-workspace estimate of +0.000200 nats (95% CI −0.001794–0.002194; uncorrected p=0.8330), and a whole-brain estimate of −0.000284 nats (p=0.6850). The direct method-by-time difference survived correction for whole-brain Φ (eight-test maximum-t p=0.03607), but not outside-workspace Φ (p=0.07629). The positive finding is therefore sensitive to signal representation and is not reproduced directly in denoised BOLD (Supplementary Table S9).

Neither the nine estimated-HRF parameter contrasts nor six motion-adjusted HRF–Φ associations survived correction (Table S10). These checks do not independently validate a neuronal interpretation or exclude nonlinear haemodynamic contributions.

### Participant influence, visits, motion and correction scope

All leave-one-participant-out whole-brain Φ estimates remained positive, ranging from 0.001392 to 0.002030 nats, with uncorrected p<0.021 in each omission (Supplementary Figure S1). The original 10,000-resample bootstrap interval was 0.000762–0.002978 nats and the Wilcoxon p value was 0.00537. In contrast, the workspace Φ bootstrap interval included zero (−0.000123–0.002717), consistent with its uncertain mean effect. These are dependent checks of the same sample.

Whole-brain Φ changes were positive at both intensities: 0.002043 nats at 90% rMT (95% CI −0.000254–0.004341) and 0.001541 at 120% (−0.000624–0.003706). Their paired difference was −0.000503 nats (−0.004164–0.003159). The individual visit changes correlated negatively (r=−0.346), so improved precision after averaging should not be interpreted as independent replication. Neither workspace Φ intensity-specific interval excluded zero; no corrected intensity contrast established a dose difference (Supplementary Table S3).

Mean FD increased from 0.1364 to 0.1415 mm, a difference of 0.0052 mm (95% CI −0.0183–0.0286). Motion-adjusted whole-brain Φ changed by 0.001787 nats, close to the unadjusted estimate, but its five-measure Holm p value was 0.16647. Motion-adjusted workspace Φ and ΦR intervals included zero (nominal p=0.18277 and 0.11110). Outside-workspace Φ retained five-measure Holm support after motion adjustment (raw p=0.03983; normalized p=0.02262), but not correction across all 20 scope/information models (p=0.15930 and 0.09048; Supplementary Table S2).

The outside-workspace increase was positive in 13 of 15 participants in raw units and 12 of 15 after normalization. Every leave-one-participant-out estimate remained positive (raw range 0.002366–0.003095 nats; uncorrected p≤0.00505). Removing each canonical network in turn also retained a positive corrected outside-workspace mean. However, none of the 35 estimable individual network-pair blocks survived correction, precluding a specific network localization. Baseline-plus-motion adjustment retained the outside-workspace estimate, but the adjusted outside-minus-inside contrast still did not survive the six-model sensitivity correction. These findings strengthen the within-sample aggregate result without establishing spatial specificity (Supplementary Table S9).

The inferential scope also mattered without motion adjustment. Treating the ten whole-brain/workspace information outcomes and 20 PC/role outcomes as one 30-endpoint family yielded whole-brain Φ maximum-t p=0.15936 raw and p=0.25073 normalized. Thus the global result was not retained by every correction or sensitivity analysis.

### Direct spatial contrast

The correction family contains five outside-minus-inside pre/post contrasts: Φ, ΦR, persistent synergy, persistent redundancy and total source-side synergy. Raw Φ had an uncorrected sign-flip p=0.03790 and five-measure maximum-t p=0.16522; normalized Φ had p=0.03723 and 0.16809, respectively. A single isolated contrast would not inherently require five-test correction, but all five measures were examined here. The raw Φ paired-t p was 0.04999. These are different inferential procedures, not interchangeable labels for the same p value.

## Numerical audit

An independent verifier reconstructed within-workspace means using direct indexing of both matrix axes, rather than the edge-mask routine used for extraction. Agreement for Φ and ΦR was within 1.1×10⁻¹⁶ nats. It also checked that averaged paired differences equal mean post minus mean pre, independently computed every workspace sign-flip statistic using SciPy, verified graph-partition identities, and checked immutable input hashes. All historical 25 scalar outcomes were reproduced across the completed TMS scans; the manuscript data package retains only the Aachen rows.

The expanded analyses use the same 15 people. All leave-one-out workspace Φ point estimates remain positive, but their uncorrected p values extend to 0.201. Positivity under omission alone therefore does not establish a reliable workspace effect. Figure S1 shows uncertainty for each omission rather than displaying only point estimates.

## Supplementary tables and participant outputs

{{SUPPLEMENT_TABLE}}

Raw estimates are nats. Normalized CSV estimates use conditional-CDF units on a 0–1 scale unless a table explicitly says percentile points; multiply these estimates and intervals by 100 for manuscript units. Integer counts refer to regions, not participants. P values are two-sided. The columns identify the correction family explicitly. Missing p values denote checks without a hypothesis test or unestimable models and must not be treated as zero.

## Supplementary figure

{{SUPPLEMENT_FIGURE}}

Figure S1. Leave-one-participant-out raw Φ changes across the whole brain and within the fixed workspace. Each point is the mean of the remaining 14 participants' visit-averaged changes; bars are marginal 95% Student-t intervals. The omitted participant is labelled by their original dataset identifier. Positive estimates are not equivalent to intervals excluding zero, and these dependent checks are not independent replications.

## References

References use the numbering of the main manuscript.

{{REFERENCES}}
