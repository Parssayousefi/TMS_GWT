# Fixed HCP workspace audit and TMS reanalysis

The fixed HCP workspace was already used correctly for participation coefficients, gateways and broadcasters. The five global information averages used all 454 regions. This rerun adds means restricted to pairs whose **both endpoints** are in the unchanged 235-region HCP workspace (27,495 unique pairs), and compares them with the 102,831-pair whole-brain means.

The GWT publication notebook and Yousefi paper report participation and role outcomes, not a matching mean-Φ outcome. Luppi tested pairwise ΦR using network-based statistics and subsequently examined workspace overlap. This scalar workspace analysis is therefore an extension of the shared estimator, not an exact repetition of Luppi's inferential procedure.

## What was rerun

All 193 completed TMS scans were re-audited from their existing pairwise matrices. All 25 historical scalar outcomes were reproduced. Conditional normalization was reevaluated for the 189 scans contributing to either the acute comparisons or all-follow-up analysis, retaining the original shared reference distributions and both seeds. Preprocessing and ΦID estimation were unchanged. No HCP or psychedelic participant analysis was rerun.

Aachen retains 15 independent participants: each contributes the mean of their two matched post-minus-pre changes. SUDMEX retains the randomized two-week comparison (40 participants) separately from its all-post longitudinal sensitivity (41 participants). The latter averages each person's available follow-ups, including open-label treatment, and does not estimate a randomized active-versus-sham effect.

## Families and correction

The five information measures are **one family**, not five families. Exact joint sign flips and the maximum absolute t statistic preserve their dependence. Scope-specific five-test results, combined whole-brain/workspace ten-test results, and a 30-endpoint sensitivity including the original 20 role/PC endpoints are all retained. Four-scope diagnostics use a joint 20-test family, with all 40 neural endpoints also provided as a sensitivity. These nested analyses remain exploratory; they do not make the complete historical search confirmatory. Raw and normalized analyses are reported together as dependent method sensitivities.

## Aachen results

Raw changes are nats. Normalized entries below are conditional-CDF units; multiply estimates and intervals by 100 for percentile points. Confidence intervals are marginal.

| branch | metric | estimate | ci_low | ci_high | p | p_maxT_family | p_maxT_whole_and_workspace10 | p_maxT_all |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| raw | whole_brain__synergy | -0.00761213 | -0.0154098 | 0.000185502 | 0.0549634 | 0.157471 | 0.233765 | 0.611206 |
| raw | whole_brain__redundancy | 0.000428225 | -0.000470205 | 0.00132666 | 0.323987 | 0.67926 | 0.805603 | 0.997437 |
| raw | whole_brain__synergy_total | -0.00207086 | -0.00520497 | 0.00106326 | 0.178304 | 0.438293 | 0.569458 | 0.952454 |
| raw | whole_brain__phi | 0.00179217 | 0.000514833 | 0.00306951 | 0.00937804 | 0.0239258 | 0.0429077 | 0.159363 |
| raw | whole_brain__phi_r | 0.0022204 | 0.000356011 | 0.00408478 | 0.0229235 | 0.0657959 | 0.108154 | 0.320984 |
| raw | within_workspace__synergy | -0.00373303 | -0.0125339 | 0.00506788 | 0.378355 | 0.736267 | 0.865295 | 0.99939 |
| raw | within_workspace__redundancy | 0.000682013 | 6.96904e-05 | 0.00129434 | 0.0315324 | 0.0997314 | 0.145813 | 0.415833 |
| raw | within_workspace__synergy_total | -0.000136914 | -0.00386422 | 0.00359039 | 0.93832 | 0.999756 | 1 | 1 |
| raw | within_workspace__phi | 0.00129929 | -0.000293117 | 0.0028917 | 0.101994 | 0.284058 | 0.387634 | 0.827209 |
| raw | within_workspace__phi_r | 0.00198131 | 0.000172611 | 0.00379 | 0.0339976 | 0.108826 | 0.15741 | 0.440491 |
| normalized | whole_brain__synergy | -0.00484626 | -0.0108885 | 0.00119601 | 0.107404 | 0.312805 | 0.408142 | 0.847778 |
| normalized | whole_brain__redundancy | 0.00538034 | -0.00335931 | 0.01412 | 0.207888 | 0.516907 | 0.631226 | 0.980713 |
| normalized | whole_brain__synergy_total | 0.00215456 | -0.0088607 | 0.0131698 | 0.681206 | 0.96814 | 0.99408 | 1 |
| normalized | whole_brain__phi | 0.00712309 | 0.00162961 | 0.0126166 | 0.0147196 | 0.0482788 | 0.0782471 | 0.250732 |
| normalized | whole_brain__phi_r | 0.00798034 | 0.00127417 | 0.0146865 | 0.0230147 | 0.0769043 | 0.118408 | 0.348999 |
| normalized | within_workspace__synergy | -0.00427893 | -0.0101792 | 0.00162134 | 0.142157 | 0.391968 | 0.494446 | 0.915955 |
| normalized | within_workspace__redundancy | 0.00703874 | 0.000958275 | 0.0131192 | 0.0263264 | 0.0967407 | 0.134521 | 0.375488 |
| normalized | within_workspace__synergy_total | 0.00226463 | -0.0107217 | 0.015251 | 0.713994 | 0.974731 | 0.99646 | 1 |
| normalized | within_workspace__phi | 0.00446643 | -0.00279484 | 0.0117277 | 0.208251 | 0.51947 | 0.63208 | 0.980835 |
| normalized | within_workspace__phi_r | 0.00658285 | -0.000240287 | 0.013406 | 0.0575062 | 0.188232 | 0.255188 | 0.615723 |

No within-workspace information endpoint survives its five-measure family correction. Workspace Φ rises by 0.001299 nats (95% CI −0.000293 to 0.002892; maxT p=0.2841). Workspace ΦR rises by 0.001981 nats (0.000173 to 0.003790; uncorrected p=0.0340, maxT p=0.1088). Normalization does not strengthen these findings: the respective five-measure maxT p values are 0.5195 and 0.1882. Motion-adjusted workspace Φ and ΦR intervals include zero; the nominal p values are 0.1828 and 0.1111.

The original whole-brain Φ finding is reproduced. Its raw maxT p is 0.0239 across five information measures and 0.0429 across the ten whole-brain/workspace outcomes. The normalized equivalents are 0.0483 and 0.0782. Neither survives the expanded 30-endpoint correction. Participation and role outcomes are unchanged by this audit; the mask had already been applied to them.

## Spatial diagnosis

The 219 non-workspace regions define 23,871 outside–outside pairs; 51,465 pairs connect workspace to outside regions. The raw Φ change outside the workspace is 0.002809 nats (95% CI 0.001198 to 0.004420), and normalized Φ rises by 1.151 percentile points (0.513 to 1.790). These survive correction across all 20 information/scope combinations (raw p=0.00946; normalized p=0.01270) and the 40-endpoint sensitivity (0.03851 and 0.03479). Within the five outside-workspace information measures, motion-adjusted Holm p values are 0.03983 and 0.02262; across all 20 scope/information measures, they are 0.15930 and 0.09048.

However, the **direct difference in change outside versus within the workspace** does not survive correction across the five information measures (raw Φ maxT p=0.1652; normalized Φ p=0.1681). A significant effect outside and a nonsignificant effect inside are not evidence of a significant difference between locations. The defensible statement is that the clearest corrected evidence concerns pairs outside the fixed workspace; a selective or preferential outside-workspace effect has not been established.

## SUDMEX and clinical outcomes

For the randomized comparison, all five within-workspace information measures have Holm-adjusted p=1.00 in both the baseline-only and baseline-plus-motion models, for raw and normalized data. Raw baseline-adjusted workspace Φ is −0.000125 nats (95% CI −0.004117 to 0.003866), and ΦR is +0.000007 (−0.003930 to 0.003945). There is no evidence here for a randomized workspace-integration increase.

The all-post sensitivity shows raw increases in workspace persistent synergy and total source synergy (five-measure maxT p=0.00116 and 0.00110), but normalized total source synergy decreases (p=0.00813); normalized persistent synergy is not supported. This is a metric- and normalization-dependent longitudinal pattern, with mixed treatment stages and unequal follow-up availability. It should not replace the randomized comparison or be called a replicated Aachen effect. No clinical–neural association survives its declared family correction. Aachen provides no matched craving or impulsivity endpoint.

## Implication for the manuscript

Retain a clearly exploratory paper about information dynamics after iTBS. Do not describe the Aachen global Φ increase as increased workspace synergy, a replicated Luppi/Yousefi result, or improved consciousness. The new fixed-mask analysis belongs prominently in the main results. Report the outside-workspace decomposition as a diagnostic extension, including the nonsignificant direct spatial contrast. The study remains small, lacks sham in Aachen, and is not evidence of clinical efficacy.

## Verification and reproducibility

The independent verifier reconstructed workspace means by indexing both matrix axes, reproduced participant averaging and all 32,768 workspace sign-flip statistics with SciPy, and checked graph-partition identities and input hashes. The direct-indexing discrepancy was below 1.1×10⁻¹⁶ nats. All frozen inputs remained unchanged.

- Run directory: `output/tms_workspace_scope/20261125T181725_063043Z`.
- Analysis: `setup/analyse_tms_workspace_scope.py`.
- Verification: `setup/verify_tms_workspace_scope.py`.
- This report, expanded diagnostics and figures: `setup/build_tms_workspace_report.py`.
- Executable review notebook: `output/jupyter-notebook/TMS_Fixed_Workspace_Reanalysis.ipynb`.
- Source methods: [Yousefi et al.](https://www.biorxiv.org/content/10.64898/2026.09.28.755172v1.full), [Luppi et al.](https://andrealuppi.github.io/publication/luppi2023elife/luppi2023elife.pdf).
