# Corrected Aachen control-workspace QC

Both pretreatment scans are included with equal weight per participant. Separate visit maps provide a stability check. These are prestimulation controls; Aachen has no sham group. All source matrices are from the isolated corrected preprocessing branch.

| Condition | Participants | Workspace regions | Shared with HCP | Dice | Gradient rho |
| --- | ---: | ---: | ---: | ---: | ---: |
| baseline_mean | 15 | 239 | 176 | 0.743 | 0.596 |
| pre90 | 15 | 234 | 180 | 0.768 | 0.607 |
| pre120 | 15 | 229 | 158 | 0.681 | 0.512 |

Between-baseline Dice: 0.782; gradient Spearman rho: 0.716.

The fixed HCP mask and gradient are byte-identical to the GWT project reference. The independently derived Aachen workspaces are QC outputs only; they do not replace the fixed reference for treatment-effect inference. No HCP or psychedelic calculations were rerun.

Participant bootstraps, split-half outputs and uniform leave-one-participant-out maps quantify sampling stability and individual influence conditional on the pipeline. The full cohort remains the main map; no participant omission was selected to improve overlap. Spatial resemblance is descriptive and cannot itself certify registration or demonstrate a treatment effect.

Independent review identified uneven inferior frontal/temporal functional coverage. Reference contrast in subjects 10 and 14 was investigated: saved dataflow and sampled actual volumes distinguish reference enhancement from the unchanged 4D intensity input. Both anatomical releases passed with this evidence; spatial bias and temporal-quality limits remain documented in the cohort audit. Nonzero reference support alone does not establish usable temporal signal.

Outputs: `output/control_workspace_qc/aachen_corrected_20261125T234156_788094Z`.
