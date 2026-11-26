# Publication narrative assessment

> This assessment predates full FD enforcement. The updated evidence and current
> model choices are in [the QC reanalysis report](TMS_FD_REANALYSIS_20261126.md).
> The newer raw Aachen outside-workspace synergy finding is a candidate
> within-cohort result; the shared two-dataset treatment story remains unsupported.

Reviewed 26 November 2026 against the corrected Aachen release and the current
SUDMEX aggregate results. This is a working editorial assessment, not a new
analysis or a submission manuscript. The intended format remains a short
exploratory research report.

**The existing two cohorts do not establish a shared positive PID treatment
effect.** The next useful paper question is whether information-processing
responses depend on the stimulation site and protocol. The newly downloaded
ds001927 can address the site comparison directly; its outcomes are not yet
known. Aachen and SUDMEX should retain their current estimates as separate tests
of generalizability, rather than being presented as positive replications.

Three decisions follow from the evidence:

1. Keep the scientific question about information processing after stimulation.
   Neither haemodynamic associations nor a claim of preserved workspace
   structure is needed for this question.
2. Keep the experimental contrasts intact: matched pre/post changes in Aachen,
   and the two-week active-versus-sham contrast in SUDMEX. Pooling later SUDMEX
   visits would change the causal question.
3. Give ds001927 a narrow, frozen analysis before its outcomes are inspected.
   Further changes to correction families cannot create a replicated effect
   where the second cohort has little unadjusted evidence.

## What the existing results support

Aachen contributes 15 independent participants, each with two matched
intensity-specific pre/post comparisons. SUDMEX contributes 40 acute-phase
completers, 22 active and 18 sham. The populations, protocols, timescales and
experimental controls differ; a pooled effect would not represent one common
intervention. In particular, a significant effect in one cohort and a
nonsignificant effect in the other would not itself establish a difference
between protocols.

| Candidate claim | Current evidence | Assessment |
| --- | --- | --- |
| TMS increases global Φ | Corrected Aachen change +0.001349 nats, 95% CI −0.001525 to +0.004224; maxT across five information measures p=.7441 | The previous positive finding is superseded. The interval leaves both decreases and increases possible. |
| TMS reduces DAN gateways | Normalized Aachen change −0.467 parcels, marginal CI −0.904 to −0.029; Holm4 p=.1527, exact sign-flip p=.0566. SUDMEX adjusted active–sham estimate +0.400, CI −0.599 to +1.399, p=.422 | No replicated reduction. The plotted Aachen interval is nominal, and the sham-controlled SUDMEX estimate points in the other direction. |
| TMS broadens redundant network participation | Raw Aachen coefficient change +.004661, CI +.000538 to +.008785; Holm4 p=.1178. SUDMEX adjusted estimate −.000295, CI −.010557 to +.009966, p=.954 | An exploratory Aachen signal. Participation describes distribution across networks; it is distinct from the amount of redundant information. |
| TMS preferentially changes information outside the workspace | Aachen raw outside-minus-inside source-synergy change −.004707 nats, exact p=.0323, maxT5 p=.1203. SUDMEX adjusted counterpart −.002580 nats, CI −.014122 to +.008963, p=.653 | Neither cohort establishes the spatial interaction after the stated family correction. Separate inside/outside significance cannot substitute for the interaction. |
| TMS changes the six processing modes consistently | Across 42 outcomes, the smallest SUDMEX baseline-and-motion-adjusted p is .3468 raw and .3504 normalized; every within-scope Holm6 value is 1 | The lack of a shared effect is not explained only by a large correction family. |
| A common regional change is hidden by global averaging | Six regional change templates learned in one cohort and tested in the other have conjunction p=1. Descriptive full-atlas correlations range from −.036 to −.006 | No common pattern was established. These small negative correlations also do not establish opposing responses. |
| Targeted information flow supplies the shared mechanism | SUDMEX's ten primary target tests all have maxT10 p≥.4057. Corrected Aachen's primary target families also do not pass in the eight-person localization subset | No replicated primary target effect. Pairwise predictive direction is not proof of a neuronal causal connection. |
| Clinical improvement links the shared neural effect to experience | All 56 clinical–neural tests per representation have Holm56 p=1; the smallest unadjusted p is .0299 | No supported clinical association in this analysis. The clinical outcomes are available in SUDMEX only, so even a future association would not replicate across both cohorts. |

Some secondary tests pass their stated family correction and should remain
visible. In the corrected eight-person Aachen target subset, three
left-minus-right baseline-plus-motion models pass Holm10 (p=.0476), and the
120%-intensity left-target persistent-synergy change is −.04120 nats
(Holm20 p=.0331). These depend on the adjustment or intensity branch and do not
replace the primary average-across-intensities test. They also do not have a
corresponding supported SUDMEX effect.

The SUDMEX all-post sensitivity also contains corrected longitudinal changes.
Within-workspace persistent synergy rises by .01616 nats and source synergy
by .007346 nats (maxT5 p=.00116 and .00110), whereas normalized source synergy
falls by .01664 on the 0–1 conditional-percentile scale (maxT5 p=.00813).
Raw amounts and conditional percentiles answer different questions. These
41-person pooled changes combine randomized and later open-label observations,
and therefore cannot be presented as an active-versus-sham treatment effect or
as support for a common increase in synergy.

The conclusions above concern failure to establish specific effects. They are
not equivalence tests, do not establish stability, and do not show that TMS has
no neural effects. The illustrative intervals are marginal 95% intervals;
family-wise p values are reported separately.

## The story we could write now

A defensible current title is **Information processing after prefrontal TMS:
a cross-dataset test**. Its arc is:

- **Question:** Does prefrontal stimulation produce an information-processing
  signature that generalizes across an acute healthy-volunteer experiment and
  a clinical randomized trial?
- **Test:** Apply the same estimator, fixed workspace reference and matched
  participant contrasts, including the complete bounded families and their
  normalization checks.
- **Finding:** Nominal changes in Aachen did not yield a shared treatment
  signature in SUDMEX, either in scalar measures or regional change patterns.
- **Contribution:** Report the estimated magnitudes and uncertainty, identify
  which claims the controlled comparison does not support, and provide a
  reproducible basis for an independent, better-controlled test.

Possible short pitch:

> We tested whether prefrontal TMS produces a reproducible signature of brain
> information processing across two open resting-state fMRI datasets. Using
> corrected anatomical processing, a fixed workspace reference and
> participant-level contrasts, we examined information measures, network roles
> and regional change patterns. Exploratory Aachen signals did not generalize
> to the randomized active-versus-sham SUDMEX comparison. Conditional
> normalization did not resolve this discrepancy. The results leave a shared
> response unestablished and motivate testing stimulation-site dependence
> directly, while retaining the uncertainty around each cohort's effect.

This would be an exploratory report with inconclusive treatment findings,
not the positive common mechanism the project currently seeks. Its suitability
for a particular journal requires an editorial judgment; no acceptance or
positive biological narrative is implied. A positive-mechanism manuscript
should wait for relevant results from the new controlled analysis.

## A distinct remaining test in the existing data

One narrowly motivated follow-up remains distinguishable from the analyses
already completed: **outgoing and incoming transfer entropy between the
localized left DLPFC target and a fixed striatal destination set**.

The completed target analysis averages all remote destinations and tests
outgoing-minus-incoming flow. Equal changes in both directions can cancel in
that difference. The six-mode FPN–striatum analysis averages the whole FPN,
which is a different source region. A fixed target–striatum test would address
a more specific circuit question, with both directional outcomes retained.
The original Aachen study directly investigated frontostriatal connectivity
and reported intensity-dependent changes. A later SUDMEX analysis studied
striatal functional topography and its treatment changes; it motivates the
circuit, but reuses the same participants and is not independent replication
of our results. [Alkhasli et al. (2019)](https://doi.org/10.3389/fnhum.2019.00190);
[Zhang et al. (2025)](https://doi.org/10.1038/s41386-025-02098-z).

If pursued, freeze the existing target definitions, all 20 atlas-labelled
striatal parcels, one native lag and the matched participant contrasts before
extracting the new outcomes. Test TEout and TEin together with Holm2 per cohort.
For a shared-effect claim, use the larger of the two cohort p values for each
direction, set inconsistent effect directions to p=1, and apply Holm2 to
those two conjunction tests. Report the ΦID components descriptively so that
a TE change is not automatically labelled
synergy. Further inferential tests of components, intensities or alternative
destinations require their own declared wider family and full reporting.
This is an additional exploratory family, not project-wide error control.

This is a lower-priority follow-up: deposited Aachen target masks with corrected
registration are available for only eight usable participants, SUDMEX's common target is
approximate, and the existing broader circuit and target tests are already
unconvincing. It should not become an unrestricted search over targets or
networks. No new target-circuit test was run for this assessment.

## The stronger next paper question

**Does the stimulated cortical site alter the balance of information storage
and integration?** This question can be tested within ds001927's 23 participants
under the same 1-Hz protocol at prefrontal, occipital and temporoparietal control
sites. CTR is active stimulation at a comparison site, not sham. The acquisition
paper reported site-dependent connectivity responses, making site dependence
a grounded hypothesis for information decomposition; it does not establish
that PID measures will show the same result.
[Castrillon et al. (2020)](https://doi.org/10.1126/sciadv.aay2739).

Keep the existing proposal small: a within-participant FRO-minus-CTR difference
in pre/post change for mean pairwise ΦID storage and redundancy-corrected
integration ΦR, with Holm2. Fix storage as `rtr + xtx + yty + sts`; it is not
interchangeable with a sum of single-region active information storage.
The distinctions among information measures follow the
[ΦID taxonomy](https://doi.org/10.1073/pnas.2423297122).

OCC provides a planned secondary site comparison with its full family reported.
If a result is confined to OCC, the paper must describe that result and its
secondary status; it cannot be rewritten as a preregistered frontal effect.
Anatomical quality, coverage, motion, usable length, preprocessing and physical
lag choices should be fixed before inspecting PID contrasts. Native one-sample
lags differ between ds001927 (1.25 s) and the existing datasets (2 s), so raw
magnitudes cannot silently be treated as measuring the same timescale.

The three datasets could then answer one coherent question about the
generality and site dependence of information-processing changes. Each must
retain its own estimate and design. A positive new result would stand on the
new comparison; it would not retroactively make Aachen and SUDMEX positive
replications. A site-specific or protocol-specific explanation is a testable
possibility, not the present conclusion.

## Evidence trail

Numbers were checked directly against current dashboard exports, which route
to the corrected Aachen release:

- `s00`, `s03`, `s04`, `s08`, `s09`: core and GWT-related contrasts.
- `s18`, `s19`: all-post workspace changes and clinical associations.
- `s34`–`s36`, `s38`: corrected Aachen and SUDMEX target-flow results.
- `s47`: direct outside-minus-inside interactions.
- `s52`–`s57`, `s59`: raw/normalized processing modes and independent-cohort
  pattern tests.

All are under `dashboard/tables/`. The read-only audit receipt and input hashes
are in `output/narrative_review_20261126/aggregate_review_receipt.json`.
No participant exclusions, production outputs, frozen reference data or
manuscript files were changed. No new hypothesis tests were run. This
assessment preserves the distinction between observed results, published
background evidence and proposed analyses.
