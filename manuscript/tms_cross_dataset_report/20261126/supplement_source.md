# Supplementary material

Haemodynamic covariation of information integration across two prefrontal stimulation datasets

## Acquisition and processing

### Aachen

We used OpenNeuro dataset ds001832, version 1.0.1 (doi:10.18112/openneuro.ds001832.v1.0.1), originally collected by Alkhasli and colleagues [@alkhasli2019]. The release contains 16 healthy adults. Each underwent left DLPFC iTBS on two separate visits, at 90% and 120% of resting motor threshold (rMT), in counterbalanced order. Stimulation comprised 600 pulses over 3 minutes 20 seconds, using a figure-of-eight coil and MRI-guided neuronavigation to a target corresponding to Talairach coordinates −45, 45, 35. The protocol consisted of 50-Hz triplets repeated at 5 Hz in intermittent trains [@huang2005]. Each visit included a 10-minute resting scan before stimulation and another beginning approximately 7 minutes afterward. Participants fixated with eyes open. Stimulation occurred outside the scanner; no pulses were delivered during fMRI acquisition. There was no sham visit.

We excluded one participant because all four runs contained 36 slices but the deposited timing metadata specified 34 slice times, preventing uniform slice-timing correction. The analysed sample comprised 15 participants, 7 women and 8 men, aged 19–44 years (mean 27.87, SD 7.12), with 60 scans. Seven received 90% rMT first and eight received 120% first. All analysed participants contributed both complete pre/post pairs.

### Signal processing

Functional images were acquired at 3 T with TR = 2 seconds, TE = 28 milliseconds, flip angle 77° and voxel dimensions 3 × 3 × 3.75 mm. Most runs contained 300 volumes; two contained 299. We used fMRIPrep 23.2.1 [@esteban2019] for motion correction, slice-timing correction, anatomical registration and normalization to MNI152NLin2009cAsym at 2-mm resolution. Because the deposited anatomical images were already skull-stripped, fMRIPrep used --skull-strip-t1w skip. Susceptibility-distortion correction and FreeSurfer surface reconstruction were not performed. The first five volumes were removed.

Run-level quality control used translation and rotation thresholds of 3 mm and 3.15°, with a continuous-segment salvage rule requiring at least 8 minutes. All 60 analysed runs were retained. One run flagged relative to the motion-correction reference passed the segment-relative rule without losing additional volumes. Denoising used CONN 25b [@whitfield2012], including five white-matter and five cerebrospinal-fluid aCompCor components [@behzadi2007], six motion parameters and their first derivatives, linear detrending, and regression followed by 0.008–0.09-Hz filtering. Global-signal regression, volume scrubbing and additional spatial smoothing were not applied.

We spatially averaged the cleaned signal within 400 Schaefer cortical parcels [@schaefer2018] and 54 Tian subcortical parcels [@tian2020]. Regional haemodynamic responses were estimated using the rsHRF implementation of blind deconvolution [@wu2013]: a 1-SD event threshold, third-order gamma basis, 24-second response window, threefold temporal upsampling, and a 4–8-second onset search. Iterative Wiener deconvolution used 50 iterations. Outputs contained 295 time points for 58 scans and 294 for two scans. Numerical audits confirmed finite outputs, 454 nonconstant regional signals per scan, and deconvolved signals that differed from both their inputs and simple standardized inputs. These checks establish processing integrity, not direct validation of recovered neuronal activity.

### SUDMEX

SUDMEX is a double-blind randomized trial of 5-Hz left-DLPFC rTMS in cocaine use disorder [@garza2021]. Our acute MRI comparison retained 40 complete pairs, 22 active and 18 sham, at baseline and two weeks. Models compared follow-up outcomes by arm with baseline adjustment and then baseline plus change in mean framewise displacement, using HC3 uncertainty. Sham maps in Figure S1 use the same 18 complete sham pairs. SUDMEX anatomy was audited separately; no Aachen-like repeat-extraction failure was found, although inferior coverage limitations remained.

The present cross-dataset analysis uses only the randomized baseline-to-two-week comparison. For both datasets, saved HRF parameters were read from the production MATLAB outputs, with 454 regions, three finite parameters and TR=2 seconds verified in every analysed scan. The 140 relevant scans comprise 60 Aachen scans and 80 SUDMEX scans. No processing was rerun for this follow-up.

## Spatial estimands and models

The independent HCP mask contains 235 regions and its complement 219. Information summaries average unique off-diagonal pairs; regional HRF parameters average regions, so the pair and region weights differ. The outside-minus-inside contrast is computed within each participant. Cross-boundary edges contribute to the whole-brain summary but belong to neither within-compartment mean. The direct contrast is not inferred by comparing significance in separate compartments.

Aachen's two complete within-visit contrasts are averaged before inference. Five direct spatial information contrasts use 32,768 exact joint sign flips and Holm correction. SUDMEX models the follow-up spatial difference as a function of treatment arm, its baseline value and motion change, using HC3 uncertainty. Unadjusted, baseline-only and arm-specific descriptive contrasts remain available. Every five-measure family is retained in both raw and normalized representations.

## Complete haemodynamic family

Fifteen associations are retained for each dataset, representation and model: height, peak latency and width at whole-brain, within-workspace, outside-workspace and outside-minus-inside scopes (twelve), plus those three whole-brain HRF predictors for the outside-minus-inside Φ outcome. The last three models test spatial slope differences using a common predictor. They differ from models that relate the HRF spatial difference to the Φ spatial difference. Holm15 controls each complete family, including nonsignificant tests.

Primary association models include motion change and, for SUDMEX, randomized arm. HRF change is standardized within dataset. Baseline sensitivities additionally include baseline Φ and the corresponding baseline HRF parameter. These are association models and do not adjust a treatment effect for an established mediator. Arm-by-HRF interaction models include both main effects and their product. Arm-specific slopes retain the full-cohort SD scale.

Same-direction cross-dataset conjunctions use the maximum of the two cohort p values for each endpoint, setting the conjunction p to one when slope directions disagree. Holm correction then covers all fifteen endpoints within each representation. This intersection-union procedure asks whether both cohorts support a same-direction association; it is not a meta-analysis or a test of equal effect sizes. The follow-up was motivated by already observed Aachen findings and remains exploratory.

## Robustness and numerical verification

The complete raw fifteen-model family was refitted after each participant omission: 225 Aachen and 600 SUDMEX fits. Correction was repeated within every omission. All width slopes at whole-brain, within-workspace and outside-workspace scopes retained their negative direction. Corrected whole-brain support remained in 12/15 Aachen and 38/40 SUDMEX omissions; outside-workspace support remained in 14/15 and 26/40. These dependent checks do not provide additional independent replications.

Rank sensitivities use twelve-decimal rounding to resolve numerical ties far below HRF measurement resolution. Whole-brain raw width correlations survive Holm15 in both datasets; outside-workspace raw rank support survives in Aachen but not SUDMEX. Normalized outside-workspace rank associations survive in both. Baseline adjustment preserves both raw outside-workspace associations but not the normalized Aachen association. The full tables retain these qualifications.

All 140 saved HRF parameter arrays were checked. The earlier 60 Aachen HRF summaries, ten spatial contrasts and nine haemodynamic associations were reproduced. Each newly fitted HC3 result was checked against a separate pseudoinverse-based sandwich-covariance calculation. The largest absolute discrepancy was below 2.4×10⁻¹². Input hashes were verified after analysis. No original matrices, masks or processing files were changed.

## Scope of interpretation

No direct spatial treatment contrast survived correction in either dataset. The shared HRF-width association does not establish preferential modulation outside the workspace, stimulation specificity, a physiological mechanism or a causal neuronal pathway. Nonsignificant interactions are not equivalence tests. The fitted response and information estimates share input signals and a deconvolution procedure. Existing denoised-versus-deconvolved Aachen checks are available in the full repository; they are method sensitivities, not independent physiological validation.

The repository retains the full earlier PID, participation, gateway/broadcaster, clinical, target-flow and all-post analyses. Later SUDMEX open-label visits are excluded from the present randomized comparison. Historical Aachen outputs preceding anatomical correction are clearly marked as superseded.

## Supplementary figure

{{SUPPLEMENT_FIGURE1}}

Figure S1. Control workspaces compared with the frozen HCP reference. The supplied image is unchanged. SUDMEX sham maps before and after two weeks use the same eighteen participants; Aachen maps average pretreatment visits within fifteen participants. Shared regions and Dice overlap document broad correspondence with previous work [@luppi2024; @yousefi2026], with incomplete overlap. These descriptive maps do not redefine the mask used in the statistical analyses.

## Supplementary data

The accompanying ZIP retains all seventy spatial result rows, 270 haemodynamic model/rank rows, thirty cross-dataset conjunction rows, 825 participant-omission fits and the influence summaries, together with the analysis specification and numerical verification. Rows from overlapping models and representations are not independent tests. Original coefficient and p-value precision is retained. No missing value is replaced by zero.

## References

{{REFERENCES}}
