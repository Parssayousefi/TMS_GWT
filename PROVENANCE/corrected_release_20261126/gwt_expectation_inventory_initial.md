# GWT expectation inventory and SUDMEX analysis audit

This inventory was defined on 26 November 2026 before this reviewer opened corrected Aachen treatment-effect results. It is a dated analysis-review specification, not a retrospective preregistration. Existing SUDMEX results and the original, anatomically superseded Aachen analyses already exist. No private psychedelic arrays or maps are needed or authorized here.

The primary evidence is the local full text [yousefi2026.xml](D:/projects/TMS_GWT/manuscript/aachen_brief_report/sources/yousefi2026.xml), *Psychedelics Restructure the Global Workspace of Consciousness*, DOI 10.64898/2026.09.28.755172. Its metadata identify the supplied version as bioRxiv v1, posted 5 October 2026. “Reported” below refers to that source; no peer-review or TMS-registration status is inferred.

## Finite expected-effect family, defined before corrected effects

Use these four metrics, retain two-sided tests, and apply Holm correction across the four expectations separately for each dataset's principal estimand. For SUDMEX that estimand is the active-minus-sham difference in t1-minus-t0 changes; for Aachen it is the participant-level mean of the two post-minus-pre changes. The latter is an uncontrolled pre/post association, not a sham-adjusted treatment effect. FPN gateways have the strongest and most direct rationale. All four outcomes remain in the family regardless of their later direction or significance. An effect opposite to the expected direction is not support, even if its two-sided test is significant.

| Exact existing metric | Expected direction | Source and strength of expectation | Limits |
| --- | --- | --- | --- |
| `gateways_FPN` | Increase | XML Results `s2a`, paragraph 2: increases in all four classic-psychedelic conditions; mean +4.37 regions, individual changes +4.03 to +4.60. Paragraph 4: every FPN increase survives motion adjustment. Discussion `s3`, paragraph 4 interprets the shift from FPN broadcasters toward gateways. | A cross-intervention expectation, not a published TMS finding. Do not import the psychedelic magnitude as a TMS effect-size threshold. |
| `gateways_DAN` | Increase | XML Results `s2a`, paragraph 2: increases in all four classic conditions. Abstract and Conclusion `s4` emphasize FPN and DAN. | DAN belonged to the paper's exploratory six-network family, and MDMA also increased it. It is not a uniquely psychedelic marker. |
| `mean_redundancy_pc` | Increase | XML Results `s2a`, paragraph 3: increases in all four classic conditions and in MDMA and amphetamine. Methods `s5f`: mean across the fixed 235 workspace regions. | Generic redistribution across modules, not drug-class specificity or an increase in raw redundancy magnitude. |
| `mean_synergy_pc` | Increase, with weaker evidence | XML Results `s2a`, paragraph 3: increases for LSD-LPM, psilocybin and mescaline, but not LSD-LAM or the non-classic conditions. | Three of four classic conditions support the direction; never describe it as a universal replicated psychedelic result. The proposal's private robustness-based rejection does not override this public source finding. |

The existing core's family-based BH values are retained as originally defined and reported separately. They do not equal a four-expectation Holm family: the core separates FPN, other network gateways and PC outcomes. A later review overlay can calculate Holm from the four unchanged two-sided p-values without altering source metrics or analyses. Adjusted and unadjusted models are sensitivity analyses of the same endpoints, not additional opportunities to select a favorable result.

## Findings that must not become invented expectations

Broadcasters are workspace regions with redundancy-PC rank above synergy-PC rank; gateways use the opposite sign, with ranks calculated among the fixed workspace regions. Equal ranks belong to neither class (XML Methods `s5d` and `s5e`). Therefore `broadcasters_FPN` and `broadcasters_DAN` are supporting complementary role summaries: a gateway increase generally accompanies a broadcaster decrease, but ties prevent an exact one-for-one identity. They are not independent replications of gateway evidence. The source does not establish a uniform direction for `broadcasters_total`, `gateways_total`, or every network's broadcaster count.

The paper's original primary network family was FPN and DMN, not FPN and DAN (Methods `s5e`). DMN increased significantly in only two of four classic conditions; one of those did not retain significance after motion adjustment (`s2a`, paragraphs 2 and 4). DMN, VIS, SMN, VAN, LIM and SUB outcomes remain in the full exploratory inventory, rather than being promoted to the finite consistent-direction family. VIS/SMN decreases are shared with MDMA/amphetamine and are condition-dependent among classic psychedelics.

`global_phi` and `global_phi_r` are distinct existing exploratory endpoints. The code averages each matrix over 102,831 unique off-diagonal region pairs. The production identity is Phi_R = Phi + persistent redundancy. A change in Phi therefore neither implies the same result for Phi_R nor establishes a change in the gateway/PC signature. These pairwise means are not a 454-region whole-system integrated-information estimate. XML `s5c` specifies persistent synergy and redundancy atoms; it does not report a directional global-Phi or Phi_R drug result. No directional expectation for either is supplied by this paper.

The source's regional drug-profile correlations and phenomenology comparisons (`s2b`, `s2c`, `s5g`, `s5h`) do not supply participant-level TMS–craving or TMS–impulsivity predictions. Published across-drug 5D-ASC profile correspondence is not an individual clinical-change association. Such SUDMEX associations remain exploratory and corrected as a complete family.

## Superseded proposal components

[PROPOSAL.md](D:/projects/TMS_GWT/docs/PROPOSAL.md), “Rationale and hypotheses,” “Data” and “Analysis plan,” combines the public rationale with private Basel tier-1 results and plans that were not implemented as written. Discard private robustness/effect-size claims, private psychedelic effect maps, the H2 TMS-versus-psychedelic map correlation, MDMA/amphetamine map comparisons, and additions of new psychedelic cohorts. Do not treat the proposal's statement about OSF preregistration or masked labels as evidence that either occurred. No registration URL or completed registration record is established by these two local sources.

Retain the public FPN rationale, fixed independent HCP workspace, unchanged metric definitions and two-sided tests. The stated mixed model with age/sex is a proposal, not the existing implementation: current SUDMEX analysis uses Welch tests of changes and HC3 regression of change on arm, baseline outcome and change in FD. The private ±2-region equivalence bound, Bayes-factor plan, five-workspace checks, GSR, altered lags and other unexecuted robustness proposals must not be claimed as completed. Open-label later visits cannot enter the randomized acute t0/t1 estimand.

## Existing implementation map

All metric/function references below are in [tms_publication.py](D:/projects/TMS_GWT/pipeline/scripts/analysis/tms_publication.py).

| Responsibility | Existing functions / exact fields | Definition |
| --- | --- | --- |
| Fixed independent reference | `load_reference`, lines 53–78 | 235 HCP workspace regions; checks ROI order and the eight-module mapping. Does not estimate a TMS-specific workspace for intervention inference. |
| Cross-module spread | `participation_coefficient`, line 81; `metric_values`, line 108 | PC = 1 minus the sum of squared fractions of each region's strength assigned to the eight modules. Uses all 454 connection targets; means use the 235 workspace regions. |
| Gateway/broadcaster roles | `classify_roles`, line 95 | Average-tie ranks of synergy/redudancy PC within valid workspace; positive difference gateway, negative broadcaster, zero neither. |
| Four expectations | `METRIC_INFO`, lines 24–33; `metric_values`, lines 124–131 | `gateways_FPN`, `gateways_DAN`, `mean_synergy_pc`, `mean_redundancy_pc`. |
| Role companions | Same definitions | `broadcasters_FPN`, `broadcasters_DAN`, `gateways_total`, `broadcasters_total`, `ties_total`; eight network-specific counts per role. |
| Pairwise information | `MATRIX_NAMES`, line 21; `metric_values`, lines 111–114; `extract_scan_metrics`, line 205 | `global_synergy`, `global_redundancy`, `global_synergy_total`, `global_phi`, `global_phi_r`; checks Phi_R identity. |
| Pair construction | `paired_changes`, line 221; `pool_aachen`, line 251 | Explicit t0/t1 or paired intensity-session labels; same valid ROI support required; two Aachen changes averaged within each participant, with no single-visit fallback. |
| Main tests and existing multiplicity | `welch_summary`, line 282; `paired_summary`, line 307; `neural_tests`, line 331; `correct_families`, line 322 | Two-sided Welch SUDMEX contrast; two-sided participant-level Aachen test; BH within the original metric families. |
| Baseline/motion sensitivity | `fit_hc3`, line 344; `motion_sensitivity`, line 365 | Change ~ active + baseline outcome + delta-FD, HC3 covariance with t reference; complete cases and rank checks. |
| Clinical associations | `clinical_analyses`, line 377 | Clinical change ~ standardized neural change + arm + clinical baseline + neural baseline + delta-FD; HC3/t inference, complete-case QC and BH across the supplied clinical–neural family. |
| Public-paper comparison | `comparison_to_paper`, line 420 | Existing qualitative mapping covers FPN, DAN and redundancy PC only. It does not already implement the present four-expectation family, nor justify a global-Phi direction. |

## SUDMEX reproducibility audit

Pending bounded checks of the existing acute output tables and clinical joins. The expectation family above is fixed before reading those tables in this review; corrected Aachen treatment-effect outputs have not been inspected.
