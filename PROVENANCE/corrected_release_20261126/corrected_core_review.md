# Independent corrected-TMS analysis review

Completed 26 November 2026. Reviewed run: `output/corrected_tms_analysis/20261126T112054_284547Z`.

**No material numerical or manuscript-claim discrepancy was found in the bounded review.** The corrected Aachen route, current motion values, within-person pairing, matrix summaries, principal inference and normalized summaries reproduced independently. Both representations contain no Aachen endpoint surviving the 40-endpoint maximum-t correction and no expectation surviving the four-outcome correction. This does not establish equivalence or absence of an effect.

## Source-grounded expectation family

The initial inventory was written at 11:22:38 UTC before this reviewer inspected corrected effect outputs. Its exact text and SHA256 are preserved in `gwt_expectation_inventory_initial.md` and `expectation_inventory_history.json`. This was an analysis-review specification, not a preregistration; historical SUDMEX and anatomically superseded Aachen results already existed.

The four retained directions were increases in FPN gateways, DAN gateways, mean redundancy participation and mean synergy participation. The last expectation has weaker consistency in the cited bioRxiv v1 full text. The inventory documents precise source sections, rejects private psychedelic maps and superseded proposal plans, and distinguishes the paper's network/participation findings from an invented directional prediction for global Φ or ΦR.

`setup/test_corrected_gwt_expectations.py` generates `gwt_expectations_raw.csv` and `gwt_expectations_normalized.csv`. It uses two-sided Holm4 for each dataset, Aachen exact maxT4 over all 32,768 participant sign patterns, and SUDMEX HC3 follow-up models adjusted for arm, baseline outcome and FD change. The latter is a conditional sensitivity contrast. All endpoints remain in the family regardless of their direction or significance.

| Aachen expectation | Raw change | Raw Holm4 | Raw maxT4 | Normalized change | Normalized Holm4 | Normalized maxT4 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| FPN gateways | −0.1333 | 1.0000 | 0.9998 | −0.3667 | 1.0000 | 0.9824 |
| DAN gateways | −0.4000 | 0.3114 | 0.3888 | −0.4667 | 0.1527 | 0.1621 |
| Redundancy participation | +0.004661 | 0.1178 | 0.1141 | +0.001194 | 0.1527 | 0.1924 |
| Synergy participation | +0.000069 | 1.0000 | 0.9721 | −0.000005 | 1.0000 | 1.0000 |

All four SUDMEX Holm4 values are 1 in both branches. Normalization changes the graph weights used to derive participation and roles; participation coefficients and region counts retain their own units. Neither is reported in percentile points.

## Independent numerical checks

| Check | Verified extent |
| --- | --- |
| Raw scan summaries | 193 scans: 60 corrected Aachen and 133 original SUDMEX; all 965 matrices; 46 summaries per scan, including every scope mean, global mean, PC mean, network role count and ties |
| Direct matrix agreement | Maximum absolute summary discrepancy 3.33×10⁻¹⁶; symmetry, finite off-diagonal support and ΦR = Φ + persistent redundancy checked |
| Fixed edge sets | 102,831 whole-brain; 27,495 within workspace; 23,871 outside; 51,465 cross-boundary pairs |
| Routing and motion | All 360 Aachen matrix/summary provenance paths use the corrected rerun; all 60 FD values match current corrected QC |
| Source integrity | 1,166 raw input hashes and 7 statistics-source hashes match; counts are per receipt and include overlapping files |
| Pairing | Both branches: 30 Aachen visit contrasts averaged within the same 15 ordered participants; 40 acute SUDMEX pairs, 22 active and 18 sham; explicit t0/t1 and pre/post intensity identities |
| Aachen inference | Independent signed-observation enumeration of all 32,768 patterns, recomputing SD for each pattern; all nine endpoint subfamilies, joint whole/workspace 10, and 20/30/40 families agree; marginal t intervals and Holm values agree |
| SUDMEX acute inference | Per branch: 30 Welch contrasts, 30 baseline HC3 and 30 baseline-plus-motion HC3 models independently reproduced; all 90 results also match the historical fixed-workspace scope outputs |
| Clinical inference | Per branch: all 56 HC3 associations independently reproduced, including Holm56 and the nested 35-model family; unchanged from historical scope outputs |
| Original publication outputs | Original 25 motion-adjusted models and 35 clinical associations reproduced separately, preserving their original definitions |

The review uses independent PC aggregation and pandas ranks, independent explicit HC3 sandwich calculations, and direct signed-observation enumeration rather than calling the production summary or testing routines. Small numerical helper fixtures compare against SciPy/statsmodels. Saved CSV round-tripping can amplify submachine-scale changes in raw-scale clinical designs; fitted-summary comparisons use an absolute tolerance of 2×10⁻⁸. Matrix comparisons remain at much tighter tolerance.

Clinical records join one-to-one on dataset, participant and contrast, with arm agreement and source-QC masking. Complete-case sample sizes are 40 for craving VAS and CCQ-now, 38 for CCQ-general, and 39 for each BIS summary. The original 35-association table standardized neural change after removing clinical missingness; the current 56-association table standardizes across the 40 neural participants first. Both were reproduced under their respective definitions. The original whole-brain ΦR association is not interchangeable with the nested current within-workspace ΦR endpoint.

## Normalization and reuse

All 60 corrected Aachen normalized scans were independently reconstructed from the saved two-seed conditional CDF references and 300 observed matrices. Maximum summary discrepancy was 2.22×10⁻¹⁶. Independent covariance/log-determinant calculations checked every TDMI cache edge, with maximum discrepancy 1.70×10⁻¹³.

The baseline coefficient pool was independently rebuilt from exactly 30 corrected baseline scans and the specified 1,024 sampled edges per scan: 30,720 coefficient matrices, maximum discrepancy 4.44×10⁻¹⁶. It differs from the superseded Aachen pool. Both existing references contain 369 grid points and 32,768 sorted finite draws per point. No historical Aachen null was reused. Reference draws were read for verification, not regenerated.

All 129 reused normalized SUDMEX scan summaries match their historical source. The original wrapper checked 480 source hashes covering the 80 acute scans. This review additionally checked all 941 relevant SUDMEX hashes in the historical scope provenance, including later inputs and references; they match. The root agent subsequently added the same 941-file guard to the wrapper and saved `normalization/extended_sudmex_reuse_guard.json`, without rerunning numerical analyses. Removing only the new guard function and its invocation restores the byte-identical original generating script. This closes the future reuse-guard gap as well as validating the current reuse.

Maximum saved interpolation error is 0.00076075 CDF units. The finite-record calibration's maximum mean-percentile departure from 0.5 is 0.06313; these are model diagnostics, not universal bias correction. Reference-seed variability is nonzero, and the shared fitted prior's uncertainty is not propagated into participant confidence intervals. The normalized branch remains a dependent model sensitivity. Its percentiles are not additive ΦID atoms and cannot be assigned the raw ΦR identity.

## Representation, HRF and manuscript review

Verified all 120 raw/normalized outside-network profile rows against fixed ROI memberships and edge-count weighting, all 180 signal-summary rows and their 90 visit contrasts, all eight direct method-by-time tests with independent maxT enumeration, all 60 saved HRF parameter files, and all six HC3 HRF associations. The associated 250 follow-up input hashes, 60 HRF hashes and 5 routing-source hashes match.

Within-workspace Φ has a deconvolved-minus-denoised change interaction at maxT8 = 0.00250244. In the original six-association family, outside-workspace Φ change is associated with HRF-width change at Holm6 = 0.000317854 and HRF-peak change at Holm6 = 0.0330886. These are method/association results from shared signals. They do not establish a stimulation-induced increase, independent replication, causal mediation, or a vascular artifact. The current manuscript maintains these distinctions.

A subsequent user-requested follow-up added three whole-brain HRF associations after the initial six results were known. Its entire nine-model family and all 135 leave-one-person-out fits reproduced from the 60 HRF MAT files and independent core participant changes, with maximum numerical discrepancy 9.25×10⁻¹⁴. The weighted 235/454 within-workspace plus 219/454 outside mean agrees with the direct mean of all 454 regions. Whole-brain Φ–width slope is −0.00459281 nats per SD (95% CI −0.00692845 to −0.00225717; Holm9 = 0.00848555). Outside width and peak survive the expanded family at 0.000476781 and 0.0463241. Width slopes remain negative in all 15 omissions, with Holm9 support in 12/15 whole-brain and 14/15 outside fits. The direct outside-minus-inside width association does not survive Holm9.

Ancillary Spearman estimates initially exposed machine-precision near-ties. The declared amendment rounds both variables to 12 decimals solely before rank correlation. Independent direct-region reconstruction then reproduces all nine rank results and Holm9 values; HC3 and leave-one-out fits are unchanged. Whole-brain and outside width correlations are approximately −0.6953 and −0.8536. `global_hrf_review.json` binds the current follow-up source and outputs. These dependent secondary associations cannot explain away the uncertain global intervention mean or establish a neuronal/vascular mechanism; the revised manuscript explicitly discusses shared-signal and algorithmic coupling.

The bounded audit of `manuscript_source.md` and `supplement_source.md` found no material numerical or inferential mismatch. Raw/normalized estimates, marginal intervals, correction labels, pair counts, four-expectation outcomes and the distinction between redundancy participation and persistent redundancy match the current tables. Figure S1 is an exact copy of the previously reviewed control-workspace comparison, SHA256 `2445849762ff8931802f355e2f0737d2c6212f88396b5ee3a8d47f6f68687245`. Its descriptive broad-correspondence wording does not assert identical maps or independent duplicate HCP comparators.

The target release's 21 listed file hashes match. Its eight-person primary analysis has minimum maxT5 = 0.1796875. The supplement correctly discloses three left-minus-right baseline-plus-motion HC3 results at Holm10 = 0.0476172 and decreased left-target persistent synergy at 120% rMT, −0.0411986 nats, Holm20 = 0.0331098. These remain model/intensity-specific sensitivities with n=8 and five residual degrees of freedom in the motion models. The unavailable fresh MATLAB check is explicitly disclosed; existing production-MATLAB and independent Python verification are not described as a fresh MATLAB rerun.

The longitudinal SUDMEX statements match the saved 41-person, 40-endpoint, 200,000-sign-vector output and are clearly separated from the acute randomized contrast. This bounded manuscript audit checked those claims and the sign-sampling implementation, without repeating the full longitudinal Monte Carlo analysis. It did not review Word layout or re-audit original acquisition references.

## Inferential limits and provenance

Exact sign-flip p values depend on independent participants and the joint sign-symmetry assumptions stated in the manuscript; they are not randomized-sham tests for Aachen. HC3 is a finite-sample t approximation, not a guarantee against misspecification. Maximum leverage reached 0.467 in raw and 0.666 in normalized acute neural models, and approximately 0.517 in clinical models. FD adjustment conditions on a post-baseline quantity and does not establish causal mediation control. Marginal intervals are not simultaneous intervals. Historical BH `q_family` values must not be presented as current 30/40-endpoint FWER or Holm4.

Anatomical recovery, complete finite matrices and control-workspace overlap do not prove full regional signal quality. Inferior temporal/orbitofrontal coverage limitations remain relevant, as documented in the earlier independent anatomy and control-workspace reviews. They are not resolved by significant method or HRF associations.

`verify_corrected_core.py` contains the read-only reproduction checks. The per-stage JSON receipts retain the script hash when each phase completed; later additions only extended the review with other phases. `corrected_core_review.json` binds the final script, review files and reviewed outputs. Only review artifacts, the separate four-expectation script and its two authorized CSVs were written. Production data/code, the core analysis wrappers, manuscript sources, HCP and psychedelic computations were not modified by this reviewer.
