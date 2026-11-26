"""Build a static, self-contained dashboard from explicitly allowed aggregate tables.

No estimates or p-values are recomputed. Original summary columns are retained.
Historical Aachen rows in the older mixed target-flow run are excluded explicitly.
"""
from __future__ import annotations

import base64
import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "dashboard"
CORE = ROOT / "output/corrected_tms_analysis/20261126T112054_284547Z"
AA = ROOT / "output/corrected_aachen_target_flow/20261126T113156_656159Z"
SUD = ROOT / "output/targeted_information_flow/20261125T193623_492441Z"
QC = ROOT / "output/control_workspace_qc/combined_20261125T234309_541220Z"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_utf8(path, text):
    """Preserve LF output on the workspace's Python version and on Windows."""
    path.write_bytes(text.encode("utf-8"))


def value(s):
    if s is None or s.strip().lower() in ("", "nan", "none", "null", "na"):
        return None
    try:
        x = float(s)
        if not math.isfinite(x):
            return None
        return int(x) if x.is_integer() else x
    except ValueError:
        return s


def main():
    OUT.mkdir(exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    specs = []

    def add(root, name, label, group, role, dataset="Aachen", branch="raw", note="", diagnostic=False, sud_only=False):
        specs.append(dict(root=root, name=name, label=label, group=group, role=role,
                          dataset=dataset, branch=branch, note=note, diagnostic=diagnostic,
                          sud_only=sud_only))

    primary = "Aachen: participant mean of the two post-minus-pre contrasts; n=15, no sham. Exact joint sign flips. Families overlap across the 20-, 30-, and 40-endpoint tables."
    for branch in ("raw", "normalized"):
        for k, desc in ((40, "all spatial and role endpoints"), (30, "main endpoints"), (20, "four-scope information endpoints")):
            add(CORE, f"{branch}_aachen_{'all40' if k == 40 else 'main30' if k == 30 else 'information20'}.csv",
                f"Aachen · {k} {desc} · {branch}", "Global and spatial", "Planned comparison", branch=branch, note=primary)
        add(CORE, f"{branch}_sudmex_acute.csv", f"SUDMEX · acute active–sham · {branch}", "Acute active–sham", "Randomized acute comparison", "SUDMEX", branch,
            "Two-week active versus sham; 22 active and 18 sham. Unadjusted change difference and HC3 baseline / baseline-plus-motion models are all retained. Correction within model.")
        add(CORE, f"gwt_expectations_{branch}.csv", f"GWT four expectations · {branch}", "GWT expectations", "Planned comparison", "Mixed", branch,
            "Four directional expectations tested two-sided. Aachen uses paired changes; SUDMEX uses baseline-plus-motion HC3. Holm4 and Aachen maxT4 are both retained.")
    for name, label, group, role, note in [
        ("aachen_all_spatial_motion.csv", "Aachen · motion across four scopes", "Motion and intensity", "Sensitivity", "HC3 adjustment; Holm5 within spatial scope and Holm20 across four scopes, separately by branch."),
        ("aachen_intensity_results.csv", "Aachen · intensity contrasts", "Motion and intensity", "Sensitivity", "Holm15 for five metrics × three unadjusted intensity contrasts per scope and branch. The adjusted intensity contrast has its own Holm5 family. A within-intensity change is not an intensity interaction."),
        ("aachen_robustness.csv", "Aachen · motion, bootstrap and Wilcoxon", "Motion and intensity", "Sensitivity", "Motion correction uses Holm5 per scope and Holm10 across whole brain plus workspace. Bootstrap intervals and Wilcoxon checks are descriptive."),
        ("aachen_spatial_interactions.csv", "Aachen · outside minus within workspace", "Global and spatial", "Sensitivity", "Five spatial difference-in-change tests per branch; exact maxT5 and Holm5."),
        ("aachen_leave_one_out.csv", "Aachen · leave-one-out primary fits", "Influence diagnostics", "Influence diagnostic", "Aggregate refits, not independent replications. Omitted participant labels are anonymized; nominal p-values are not multiplicity corrected."),
        ("normalization_seed_tests.csv", "Aachen · normalization seed stability", "Normalization stability", "Sensitivity", "All 40 endpoints at each saved reference seed and the combined reference; seed variants are not independent replications."),
    ]:
        add(CORE, name, label, group, role, branch="normalized" if "seed" in name else "mixed", note=note)
    for name, label, note in [
        ("sudmex_all_post_all40_results.csv", "SUDMEX · all available post scans · 40 endpoints", "Within-person post average minus baseline; n=41. MaxT across 40 and per family. 200,000 Monte Carlo sign vectors including identity. Later open-label treatment and unequal follow-ups prevent a randomized treatment interpretation."),
        ("sudmex_all_post_arm_results.csv", "SUDMEX · all-post original-arm contrasts", "All available post scans averaged per person. Original assignment comparisons after later open-label treatment are descriptive, not continuing randomized exposure."),
        ("sudmex_all_post_workspace_results.csv", "SUDMEX · all-post workspace five", "Five within-workspace measures. Same descriptive all-post design; later open-label treatment and unequal follow-ups."),
    ]:
        add(CORE, name, label, "All-post longitudinal", "Descriptive longitudinal", "SUDMEX", "mixed", note)
    add(CORE, "clinical_neural_associations.csv", "SUDMEX · clinical–neural associations", "Clinical associations", "Association", "SUDMEX", "mixed",
        "All 56 clinical–neural associations per branch; Holm56 plus the original workspace35 subset. Aggregate model summaries only. Association does not identify a treatment mechanism.")
    for name, label, group, role, note in [
        ("baseline_motion_sensitivity.csv", "Aachen · baseline and motion follow-up", "Motion and intensity", "Sensitivity", "Complete six-test family per branch; HC3 models."),
        ("estimator_sensitivity_results.csv", "Aachen · estimator and equal-length variants", "Method and temporal checks", "Sensitivity", "All 12 estimator variants retained, one joint maxT12 family. Equal 294-volume windows and no-deconvolution variants reuse the same participants."),
        ("hrf_parameter_tests.csv", "Aachen · HRF parameter changes", "HRF associations", "Method diagnostic", "Nine HRF-parameter tests; a physiological/model diagnostic, not independent TMS evidence."),
        ("hrf_phi_associations.csv", "Aachen · original six HRF–Φ associations", "HRF associations", "Association", "Original six associations, Holm6, retained for full accounting. The expanded current nine-association family is in Haemodynamic follow-up with Holm9. These associations do not establish a TMS treatment effect or a causal mechanism."),
        ("method_interactions.csv", "Aachen · deconvolution method × time", "Method and temporal checks", "Method diagnostic", "Eight paired method-by-time interactions; exact maxT8. A method interaction means the estimated change differs between processing variants; it is not a TMS treatment effect."),
        ("network_results.csv", "Aachen · fixed-network contributions", "Network checks", "Sensitivity", "Fixed-network block contributions and all eight network omissions; joint maxT within each saved family. No network selected from p-values."),
        ("nonparametric_checks.csv", "Aachen · nonparametric checks", "Influence diagnostics", "Influence diagnostic", "Wilcoxon, sign tests, bootstrap intervals and visit correlation are diagnostic. No adjusted p-value supplied."),
        ("participant_influence.csv", "Aachen · spatial leave-one-out fits", "Influence diagnostics", "Influence diagnostic", "All participant omissions retained as aggregate refits with anonymized omission labels. Nominal p-values only."),
        ("temporal_component_results.csv", "Aachen · TDMI and self-information components", "Method and temporal checks", "Method diagnostic", "TDMI minus self-information identity and estimator variants; joint maxT9 within each processing branch. Same participants, not replications."),
        ("visit_consistency_results.csv", "Aachen · visit and intensity consistency", "Motion and intensity", "Sensitivity", "All 12 visit/intensity questions per branch; exact joint maxT12."),
    ]:
        add(CORE / "sensitivities", name, label, group, role, branch="mixed", note=note)
    for name, label, role, note in [
        ("association_results.csv", "Haemodynamics · all nine Φ–HRF associations", "Association", "Three new whole-brain associations join the previous six in Holm9. HC3 slope per SD of HRF change, adjusted for motion. HRF and deconvolved Φ are estimated from the same signals; coupling is not independent vascular causality or rescue of the historical global increase."),
        ("rank_sensitivity.csv", "Haemodynamics · all nine rank sensitivities", "Sensitivity", "Unadjusted Spearman sensitivity, all nine pairs retained with Holm9. The displayed estimate is rho, not an HC3 slope."),
        ("leave_one_out.csv", "Haemodynamics · all omission fits", "Influence diagnostic", "All nine associations refit after each of 15 participant omissions, with Holm9 repeated within each omission. Omission labels anonymized. Stability does not identify causality."),
    ]:
        add(CORE / "haemodynamic_followup", name, label, "Haemodynamic follow-up", role, note=note)
    add(CORE / "haemodynamic_followup", "influence_summary.csv", "Haemodynamics · influence summary", "Haemodynamic diagnostics", "QC", note="Aggregate ranges across all 15 omissions for all nine associations; no selection of omissions.", diagnostic=True)
    add(AA, "inference.csv", "Aachen · corrected target flow", "Target flow", "Target comparison", note="Eight participants, 16 deposited sides. Left and left-minus-right, five measures, three lag/length variants. maxT5 per scope, Holm10 per variant. No sham; ROI identity does not independently verify coil/peak-field location.")
    add(AA, "aachen_baseline_motion_sensitivity.csv", "Aachen target · baseline + motion", "Target flow", "Sensitivity", note="HC3 baseline and motion model; n=8, residual df=5. Holm10. High leverage (maximum .896) and model sensitivity limit interpretation.")
    add(AA, "aachen_intensity_descriptive.csv", "Aachen target · intensity-specific changes", "Target flow", "Descriptive longitudinal", note="All 20 intensity-specific changes; Holm20. A change within the 120% session is not a direct 120%-versus-90% contrast.")
    add(AA, "leave_one_out.csv", "Aachen target · leave-one-out", "Influence diagnostics", "Influence diagnostic", note="All 80 aggregate refits; anonymized omission labels; nominal p-values only.")
    add(SUD, "inference.csv", "SUDMEX · validated target flow", "Target flow", "Randomized acute comparison", "SUDMEX", "raw",
        "Only SUDMEX rows from the validated mixed run. Common baseline target n=40 (22/18); individual target subset n=19 (12/7). Both models and lag/length checks retained. Freedman–Lane permutation p-values; maxT5 by scope and maxT10 across both scopes.", sud_only=True)
    add(SUD, "leave_one_out.csv", "SUDMEX target · leave-one-out", "Influence diagnostics", "Influence diagnostic", "SUDMEX", "raw",
        "Only SUDMEX aggregate refits; omitted participant labels anonymized. Nominal p-values, not independent replications.", sud_only=True)
    for name, label in [
        ("monte_carlo_metric_stability.csv", "Normalization · Monte Carlo stability"),
        ("interpolation_validation.csv", "Normalization · interpolation validation"),
        ("synthetic_calibration_summary.csv", "Normalization · synthetic calibration"),
    ]:
        add(CORE / "normalization", name, label, "Calibration", "Calibration diagnostic", "Mixed", "normalized", "Diagnostic calibration only. KS p-values are not treatment-effect tests.", diagnostic=True)
    for name, label in [("network_node_counts.csv", "Fixed-network node counts"), ("network_pair_counts.csv", "Fixed-network pair counts")]:
        add(CORE / "sensitivities", name, label, "Network structure", "QC", note="Fixed descriptive network accounting; no inferential treatment test.", diagnostic=True)
    for name, label in [("motion_model_diagnostics.csv", "Aachen target · model leverage"), ("leave_one_out_summary.csv", "Aachen target · omission summary")]:
        add(AA, name, label, "Target diagnostics", "QC", note="Descriptive model diagnostics; no new exclusion or test selection.", diagnostic=True)

    sources, rows, diagnostics, original_hashes = [], [], [], {}
    for index, spec in enumerate(specs):
        p = spec["root"] / spec["name"]
        original_hashes[str(p)] = sha(p)
        with p.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            original = list(reader)
        selected = [r for r in original if not spec["sud_only"] or r.get("dataset") == "sudmex_tms"]
        assert selected, p
        # Leave-one-out model summaries may be public; participant identifiers are not.
        id_fields = [f for f in fields if f in ("omitted", "omitted_subject")]
        omitted_values = sorted({r[f] for r in selected for f in id_fields})
        aliases = {v: f"Omission {i+1}" for i, v in enumerate(omitted_values)}
        for row in selected:
            for f in id_fields:
                row[f] = aliases[row[f]]
        assert not set(fields) & {"subject", "participant", "participant_id", "subject_id", "prefix", "scan", "session"}, p
        ident = f"s{index:02}"
        export = OUT / "tables" / f"{ident}_{p.name}"
        with export.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(selected)
        source = {k: v for k, v in spec.items() if k != "root"}
        source.update(id=ident, path=str(p), sha256=original_hashes[str(p)], rows=len(selected),
                      source_rows=len(original), columns=fields, download=f"tables/{export.name}",
                      export_sha256=sha(export), transformations=[
                          *(["Keep dataset == sudmex_tms; discard all historical Aachen rows"] if spec["sud_only"] else []),
                          *(["Replace omitted participant identifiers with stable within-table omission numbers"] if id_fields else [])])
        sources.append(source)
        parsed = [{k: value(v) for k, v in r.items()} for r in selected]
        if spec["diagnostic"]:
            diagnostics.append(dict(source=ident, rows=parsed))
            continue
        for ix, raw in enumerate(parsed):
            dataset = {"aachen_itbs": "Aachen", "sudmex_tms": "SUDMEX"}.get(raw.get("dataset"), spec["dataset"])
            branch = raw.get("branch") or spec["branch"]
            if branch == "mixed":
                branch = "raw"
            scale = "normalized" if branch == "normalized" else "raw"
            endpoint = raw.get("metric") or raw.get("neural_metric") or raw.get("outcome") or (f"{raw['scope']}__phi / HRF {raw['parameter']}" if raw.get("parameter") else "Not supplied")
            detail = " · ".join(str(raw[k]) for k in ("scope", "comparison", "contrast", "clinical_metric", "hrf_parameter", "model", "check", "variant", "definition", "reference_seed", "omitted_subject", "omitted") if raw.get(k) is not None)
            rows.append(dict(id=f"{ident}_{ix}", source=ident, dataset=dataset, scale=scale,
                             branch=branch, endpoint=endpoint, detail=detail, role=spec["role"],
                             group=spec["group"], original=raw))

    # Remove only obsolete numbered exports previously generated in this owned folder.
    expected_exports = {(OUT / s["download"]).resolve() for s in sources}
    table_root = (OUT / "tables").resolve()
    for old in table_root.glob("s[0-9][0-9]_*.csv"):
        resolved = old.resolve()
        assert resolved.parent == table_root
        if resolved not in expected_exports:
            resolved.unlink()

    def lookup(name, metric, **where):
        return next(r["original"] for r in rows if sources[int(r["source"][1:])]["name"] == name
                    and r["endpoint"] == metric and all(r["original"].get(k) == v for k,v in where.items()))
    phi = lookup("raw_aachen_all40.csv", "whole_brain__phi")
    assert abs(phi["estimate"] - 0.0013491622181517832) < 1e-14
    assert all(r["dataset"] == "SUDMEX" for r in rows if sources[int(r["source"][1:])]["sud_only"])
    assert all(sha(Path(p)) == h for p, h in original_hashes.items())

    receipt_paths = [CORE / "analysis_specification.json", CORE / "statistics_source_sha256.json", CORE / "raw_input_sha256.json",
                     CORE / "extended_checks_receipt.json", CORE / "normalization/verification.json",
                     CORE / "sensitivities/analysis_scope.json", CORE / "sensitivities/followup_verification.json",
                     AA / "target_flow_release_receipt.json", AA / "independent_verification.json",
                     SUD / "independent_verification.json", SUD / "matlab_verification.json",
                     ROOT / "output/aachen_alignment_review/corrected_cohort/cohort_release_receipt.json",
                     CORE / "haemodynamic_followup/analysis_specification.json", CORE / "haemodynamic_followup/verification.json",
                     QC / "README.md", QC / "qa.json", QC / "provenance.json"]
    receipts = [dict(path=str(p), sha256=sha(p)) for p in receipt_paths]
    image_path = QC / "control_workspace_comparison.png"
    image_bytes = image_path.read_bytes()
    (OUT / "control_workspace_comparison.png").write_bytes(image_bytes)
    built = datetime.now(timezone.utc).isoformat()
    provenance = dict(built_utc=built, builder=str(Path(__file__).resolve()), builder_sha256=sha(Path(__file__)),
                      template_sha256=sha(OUT / "template.html"), core_directory=str(CORE),
                      corrected_aachen_target_directory=str(AA), sudmex_target_directory=str(SUD),
                      historical_aachen="Invalidated by anatomical over-extraction; no historical Aachen result rows imported.",
                      patient_records_included=False, analyses_recomputed=False,
                      image=dict(path=str(image_path), sha256=sha(image_path), copied_byte_identical=True),
                      source_tables=sources, verification_receipts=receipts,
                      limits=["Aachen n15 has no sham; target subset n8 has no sham.",
                              "Control overlap is descriptive and does not certify temporal or regional signal quality.",
                              "Method-by-time and HRF associations are not TMS treatment effects.",
                              "Fresh corrected target 24-case MATLAB check unavailable: license-server error 96; existing production-MATLAB summary agreement and independent numerical checks passed.",
                              "Families and variants overlap; displayed rows are not independent hypotheses."])
    # Explicit inventory of unbundled tables prevents silent omission of analysis families.
    allowed = {str(s["root"] / s["name"]) for s in specs}
    provenance["unbundled_core_csv"] = [dict(path=str(p), reason="Participant/scan-level input or intermediate/redundant extraction table; aggregate inference and diagnostics are bundled.")
                                        for p in sorted(CORE.rglob("*.csv")) if str(p) not in allowed]
    data = dict(built_utc=built, sources=sources, rows=rows, diagnostics=diagnostics,
                provenance=provenance, counts=dict(reported_rows=len(rows), diagnostic_rows=sum(len(d["rows"]) for d in diagnostics),
                                                 source_tables=len(sources), roles=dict(Counter(r["role"] for r in rows))))
    serialized = json.dumps(data, ensure_ascii=False, allow_nan=False, separators=(",", ":")).replace("</", "<\\/")
    template = (OUT / "template.html").read_text(encoding="utf-8")
    html = template.replace("__DASHBOARD_DATA__", serialized).replace("__CONTROL_IMAGE__", "data:image/png;base64," + base64.b64encode(image_bytes).decode())
    assert "__DASHBOARD_DATA__" not in html
    write_utf8(OUT / "index.html", html)
    write_utf8(OUT / "provenance.json", json.dumps(provenance, indent=2, ensure_ascii=False))
    write_utf8(OUT / "data.json", json.dumps(data, ensure_ascii=False, allow_nan=False))
    build_receipt = dict(built_utc=built, passed=True, counts=data["counts"],
                         source_hashes_unchanged=True, original_summary_values_recomputed=False,
                         historical_aachen_rows_included=0, patient_level_clinical_rows_included=0,
                         missing_values_serialized_as_null=True,
                         browser_qa="See browser_qa.json; build receipt alone does not claim browser QA.",
                         outputs={str(p.relative_to(OUT)):sha(p) for p in sorted([
                             OUT / name for name in ("index.html", "data.json", "provenance.json", "template.html", "README.md", "DESIGN.md")
                         ] + [OUT / s["download"] for s in sources])},
                         control_image_sha256=sha(OUT / "control_workspace_comparison.png"))
    write_utf8(OUT / "build_receipt.json", json.dumps(build_receipt, indent=2))
    print(json.dumps(dict(passed=True, counts=data["counts"], index_sha256=sha(OUT / "index.html"))))


if __name__ == "__main__":
    main()
