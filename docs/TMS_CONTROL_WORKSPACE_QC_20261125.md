# Control workspace quality check

Workspaces were independently derived using the unchanged GWT core: arithmetic means of persistent synergy/redundancy matrices, node-strength ranks, and synergy rank greater than redundancy rank. The HCP mask was not used to construct or tune these maps.

| Dataset | Condition | Participants | Workspace regions | Shared with HCP | Dice | Gradient rho |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| SUDMEX | All pretreatment scans | 41 | 227 | 179 | 0.775 | 0.676 |
| SUDMEX | Matched sham: pretreatment | 18 | 234 | 182 | 0.776 | 0.652 |
| SUDMEX | Matched sham: two weeks | 18 | 232 | 173 | 0.741 | 0.615 |
| Aachen | Two pretreatment visits averaged | 15 | 239 | 176 | 0.743 | 0.596 |
| Aachen | Before 90% session | 15 | 234 | 180 | 0.768 | 0.607 |
| Aachen | Before 120% session | 15 | 229 | 158 | 0.681 | 0.512 |

Dice measures mask overlap (1 means identical); rho is the Spearman correlation of the unthresholded synergy-minus-redundancy rank gradients.

The same 18 sham participants contribute at both SUDMEX visits: between-visit Dice 0.803, gradient rho 0.790. Later open-label treatment is excluded.

Aachen has no sham arm. Its two prestimulation visits are paired within the same 15 participants; between-visit Dice 0.782, gradient rho 0.716. The combined map averages matrices within participant before averaging participants, with no time-series concatenation.

For the combined Aachen map, the central 95% participant-bootstrap range of Dice against HCP is 0.670–0.766; omitting each participant in turn gives 0.720–0.748. These are stability diagnostics with the reference held fixed.

Our 235-region HCP mask and its gradient are byte-identical to the original GWT project files. They are one frozen reference, not two independent comparators. HCP and psychedelic calculations were not rerun.

Each condition includes 2,000 participant bootstraps and 1,000 split halves. Aachen also includes every leave-one-participant-out map. These describe sampling stability conditional on the pipeline; no spatial significance test, tuned overlap threshold, or selected participant exclusion was used.

Aachen uses the corrected branch only: all 15 participants and 60 runs passed anatomical review and output validation through deconvolution and PhiID. Native anatomical support retention improved from 44.2–52.9% to 99.63–100%. Regional functional coverage remains uneven, particularly in inferior frontal and temporal parcels. The SUDMEX audit found no comparable anatomical over-extraction, but also documented inferior functional-coverage limitations. Workspace resemblance does not certify the temporal signal quality of every parcel.

The previous Aachen treatment-effect estimates remain on hold pending statistical reanalysis of these corrected inputs. No treatment-effect inference is made in this control QC.

![Control workspaces compared with the frozen HCP mask](D:/projects/TMS_GWT/output/control_workspace_qc/combined_20261125T234309_541220Z/control_workspace_comparison.png)

SUDMEX source run: `output/control_workspace_qc/20261125T201530_636948Z`.
Corrected Aachen source run: `output/control_workspace_qc/aachen_corrected_20261125T234156_788094Z`.
Full anatomical review: `output/aachen_alignment_review/corrected_cohort/COHORT_FINDINGS.md`.
