"""Explicit aggregate-only FD reanalysis sources for the existing dashboard."""
from pathlib import Path
import html
import json
import base64
import hashlib
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "output/tms_fd_reanalysis_20261126"


def attach(specs, add):
    calibration = json.loads((RUN / "recalibration/verification.json").read_text())
    assert calibration["passed"] and calibration["state"] == "complete"
    core_release = json.loads((RUN / "core_qc_only/verification.json").read_text())
    assert core_release["group_motion_covariates"] is False
    for policy in ["registered", "strict"]:
        assert (RUN / "core_qc_only" / policy / "normalized_recalibrated_aachen_all40.csv").exists()
    old_count = len(specs)
    for spec in specs:
        spec["label"] = "Before full FD QC · " + spec["label"]
        spec["group"] = "Before FD QC · " + spec["group"]
        spec["note"] = "Historical sample before the configured FD exclusions. See current QC refits for intervention inference. " + spec["note"]
    for policy, label, aa, active, sham in [("registered", "Configured FD rules", 13, 18, 15),
                                           ("strict", "Stricter FD sensitivity", 12, 17, 13)]:
        core = RUN / "core_qc_only" / policy
        note = f"{label}: Aachen n={aa}, complete two-visit means; acute SUDMEX n={active+sham} ({active} active/{sham} sham). All-post SUDMEX uses n={36 if policy == 'registered' else 33}. "
        for branch in ["raw", "normalized_frozen_reference", "normalized_recalibrated"]:
            scale_note = ("Original normalization reference retained as sensitivity. " if branch == "normalized_frozen_reference" else
                          "Normalization references exclude QC-failing baseline records. " if branch == "normalized_recalibrated" else "Raw information units. ")
            for family in ["all40", "main30", "information20"]:
                add(core, f"{branch}_aachen_{family}.csv", f"{label} · Aachen {family} · {branch}",
                    "QC refit · information and roles", "Within-person change", branch=branch,
                    note=note + scale_note + "Five information measures per spatial scope; original role families and wider corrections retained. No sham in Aachen.")
            add(core, f"{branch}_sudmex_acute_all40.csv", f"{label} · SUDMEX acute · {branch}",
                "QC refit · randomized comparison", "Randomized acute comparison", "SUDMEX", branch,
                note + scale_note + "Baseline-outcome-adjusted HC3 and unadjusted contrasts. No extra group-level FD covariate. Corrections apply within each original family/model.")
            add(core, f"gwt_expectations_{branch}.csv", f"{label} · GWT4 · {branch}",
                "QC refit · GWT expectations", "Four-outcome expectation family", "Mixed", branch,
                note + scale_note + "Four expectations tested two-sided; Holm4 and Aachen exact maxT4.")
        for filename, group, role, dataset, extra in [
            ("aachen_spatial_interactions.csv", "spatial differences", "Spatial interaction", "Aachen", "Direct outside-minus-inside changes, five measures per representation."),
            ("sudmex_spatial_interactions.csv", "spatial differences", "Spatial interaction", "SUDMEX", "Active-sham effect on outside-minus-inside contrast; baseline-outcome-adjusted HC3; Holm5."),
            ("aachen_intensity_results.csv", "intensity", "Intensity sensitivity", "Aachen", "Original five-measure by three-contrast families plus a separate order-and-baseline-adjusted five-measure interaction family; no extra FD covariate."),
            ("aachen_robustness.csv", "robustness", "Sensitivity", "Aachen", "Bootstrap and Wilcoxon; bootstrap intervals are marginal. No extra FD covariate."),
            ("aachen_information_leave_one_out.csv", "influence", "Influence diagnostic", "Aachen", "All20 information endpoints refit after each omission; corrections repeated within omission, not independent replications."),
            ("clinical_neural_associations.csv", "clinical associations", "Association", "SUDMEX", "All56 clinical-neural associations per representation. Flagged questionnaire values remain excluded."),
            ("sudmex_all_post_all40_results.csv", "all-post longitudinal", "Descriptive longitudinal", "SUDMEX", "Passing baseline and every passing follow-up averaged within person. Later open-label visits and unequal follow-up prevent a randomized treatment interpretation."),
            ("sudmex_all_post_arm_results.csv", "all-post longitudinal", "Descriptive longitudinal", "SUDMEX", "Original arms after later open-label treatment; not a continuing randomized contrast."),
            ("cross_dataset_conjunction40.csv", "shared effects", "Cross-dataset conjunction", "Mixed", "Maximum of both cohort p-values; same-direction flag additionally required. Unestimable members remain missing and count in correction burden."),
        ]:
            add(core, filename, f"{label} · {filename.removesuffix('.csv')}", "QC refit · " + group, role, dataset, "mixed", note + extra)
        modes = RUN / "modes_qc_only" / policy
        for sub, branch in [("raw", "raw"), ("normalized_fixed_reference", "normalized_frozen_reference"),
                            ("normalized_recalibrated", "normalized_recalibrated")]:
            for name, role in [("mode_treatment_results.csv", "Exploratory treatment comparison"),
                               ("cross_dataset_conjunction.csv", "Cross-dataset conjunction")]:
                add(modes / sub, name, f"{label} · six PID modes · {branch} · {role}", "QC PID modes", role, "Mixed", branch,
                    note + "Six modes in seven fixed pair groups; Holm6/maxT6 and Holm42/maxT42. Reference choice is explicit in the representation label.")
        for name, role in [("external_template_tests.csv", "Exploratory treatment comparison"),
                           ("cross_dataset_conjunction.csv", "Cross-dataset conjunction"), ("leave_one_out.csv", "Influence diagnostic")]:
            add(modes / "regional_templates", name, f"{label} · regional templates · {role}", "QC PID regional change", role, "Mixed", "raw",
                note + "Both external-cohort templates were refitted after exclusions. Full454/cortical400 branches; Holm6/Holm12. Participant-level tests.")
        targets = RUN / "targets_qc_only" / policy
        for name, role in [("inference.csv", "Target-flow comparison"),
                           ("aachen_intensity_descriptive.csv", "Sensitivity"), ("aachen_intensity_exact_sensitivity.csv", "Exact sensitivity"),
                           ("leave_one_out.csv", "Influence diagnostic")]:
            add(targets, name, f"{label} · targeted flow · {name.removesuffix('.csv')}", "QC target flow", role, "Aachen" if name.startswith("aachen_") else "Mixed", "raw",
                note + "Target-eligible subsets are smaller: see row N. Five flow measures/two scopes, original lag and equal-length variants. SUDMEX adjusts baseline outcome only.")
        sens = RUN / "sensitivities_qc_only" / policy
        for name in ["estimator_sensitivity_results.csv", "method_interactions.csv",
                     "network_results.csv", "nonparametric_checks.csv", "participant_influence.csv", "temporal_component_results.csv", "visit_consistency_results.csv"]:
            add(sens, name, f"{label} · method check · {name.removesuffix('.csv')}", "QC method diagnostics", "Method diagnostic", "Aachen", "mixed",
                note + "Cached corrected signals/subsets; original method/network families refitted. Processing interactions are not treatment effects. Cached normalized checks retain original calibration.")
    specs[:] = specs[old_count:] + specs[:old_count]


def headline():
    p = RUN / "core_qc_only/registered"
    a = pd.read_csv(p / "raw_aachen_all40.csv").set_index("metric")
    strict = pd.read_csv(RUN / "core_qc_only/strict/raw_aachen_all40.csv").set_index("metric")
    row = a.loc["outside_workspace__synergy_total"]
    phi = a.loc["whole_brain__phi"]
    return f'''<p class="intro">FD-enforced reanalysis: Aachen 13 participants and SUDMEX 33 (18 active, 15 sham). The stricter sensitivity retains 12 and 30. Some Aachen within-person changes survive correction; a shared intervention effect across the two datasets remains unsupported.</p>
<div class="notice"><p><strong>QC correction:</strong> displacement checks had run, but configured mean/peak FD exclusions had not. Every retained scan window has now been audited. Thresholds are not a preregistration. Raw, frozen-reference normalized and freshly recalibrated normalized results are all available below.</p></div>
<div class="two"><article class="surface"><h2>Strongest Aachen signal</h2><div class="bigNumber">{row.estimate:+.6f}</div><p>Total synergy between pairs outside the fixed workspace decreases. Exact maxT5 p={row.p_maxT_family:.4f}; full40 maxT p={row.p_maxT_all:.4f}.</p><p class="muted">The stricter sample still passes its five-measure family (p={strict.loc['outside_workspace__synergy_total','p_maxT_family']:.4f}), but not full40. Every participant omission retains the decrease and passes maxT5 in both QC samples. Normalization and the randomized SUDMEX comparison do not establish the same reduction.</p></article>
<article class="surface"><h2>Global Φ and replication</h2><p>Aachen global Φ change is {phi.estimate:+.6f} nats; 95% CI [{phi.ci_low:.6f}, {phi.ci_high:.6f}], maxT5 p={phi.p_maxT_family:.3f}.</p><p>No GWT expectation survives its four-test family. No corrected acute active–sham effect is established in SUDMEX. Neither modes nor independently refitted regional templates establish a common response.</p></article></div>
<article class="surface" style="margin-top:28px"><h2>Interpretation for the paper</h2><p>A candidate Aachen pre/post synergy reduction warrants qualification. It does not yet support a replicated TMS mechanism: Aachen has no sham, and the normalized analysis estimates a different conditional quantity. SUDMEX's pooled all-post changes remain longitudinal findings, not evidence of active treatment versus sham.</p><p>Current models use motion QC and the existing denoising, without an additional group-level FD covariate. SUDMEX adjusts for baseline outcome. <a href="../docs/TMS_FD_REANALYSIS_20261126.md">Read the QC reanalysis report</a>. Tables labelled “Before full FD QC” are historical snapshots.</p></article>
<figure class="surface" style="margin-top:28px"><img alt="Raw total synergy decreases in both QC samples, while normalized intervals include zero" src="{comparison_figure()}" style="width:100%;height:auto"><figcaption class="small muted">Error bars are marginal95% confidence intervals. Printed p-values use exact maxT across the five information measures within the outside-workspace scope. Normalized scores use QC-clean baseline reference pools.</figcaption></figure>
<button class="plain primary" data-open="results">Explore the full QC result ledger</button>'''


def comparison_figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")
    for ax, branch, title in zip(axes, ["raw", "normalized_recalibrated"],
                                 ["Raw total synergy (nats)", "Normalized total synergy (conditional percentile)"]):
        for i, (policy, color) in enumerate([("registered", "#17617b"), ("strict", "#8b589d")]):
            path = RUN / "core_qc_only" / policy / f"{branch}_aachen_all40.csv"
            d = pd.read_csv(path).set_index("metric").loc["outside_workspace__synergy_total"]
            ax.errorbar(d.estimate, i, xerr=[[d.estimate-d.ci_low], [d.ci_high-d.estimate]],
                        fmt="o", color=color, capsize=5, lw=2, markersize=7)
            ax.text(.03, .87-.45*i, f"maxT5 p={d.p_maxT_family:.4f}", transform=ax.transAxes, color=color, fontsize=10)
        ax.axvline(0, color="#909090", linestyle="--", lw=1)
        ax.set_yticks([0, 1], ["Configured QC", "Stricter QC"])
        ax.set_ylim(1.5, -.55)
        ax.set_title(title, fontsize=12, pad=13)
        ax.set_xlabel("Post minus pre")
        ax.spines[["right", "top"]].set_visible(False)
    fig.suptitle("Aachen: total synergy outside the fixed workspace", fontsize=15)
    path = ROOT / "dashboard/qc_synergy_comparison.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()


def publish_receipt():
    out = ROOT / "PROVENANCE/fd_reanalysis_20261126"
    out.mkdir(parents=True, exist_ok=True)
    audit = ROOT / "output/existing_tms_fd_audit_20261126/audit_receipt.json"
    record = json.loads(audit.read_text())
    source_paths = [audit, RUN / "core_qc_only/verification.json", RUN / "modes_qc_only/verification.json",
                    RUN / "targets_qc_only/current_model_receipt.json", RUN / "sensitivities_qc_only/current_model_receipt.json",
                    RUN / "recalibration/verification.json"]
    summary = {"scope": "Current QC-enforced TMS statistical reanalysis", "configured_thresholds": record["registered_thresholds"],
               "stricter_sensitivity": record["strict_sensitivity"], "counts": record["counts"],
               "group_motion_covariates": False, "sudmex_model": "followup ~ active + baseline outcome; HC3/t",
               "aachen_model": "mean of two matched intensity changes per participant; exact sign flips",
               "normalization": "QC-clean references for main normalized scores; frozen references explicitly labelled",
               "not_rerun": "Hemodynamic associations and frozen HCP/psychedelic cohorts",
               "receipt_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
               "report": "docs/TMS_FD_REANALYSIS_20261126.md", "dashboard": "dashboard/index.html"}
    (out / "analysis_release.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
