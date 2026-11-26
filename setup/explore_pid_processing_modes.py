"""Exploratory six-mode PhiID analysis using corrected TMS signals only.

Production preprocessing, estimator helpers, HCP and psychedelic inputs remain
unchanged. The six disjoint modes cover all 16 atoms (Mediano et al., 2025,
doi:10.1073/pnas.2423297122). They are predictive information descriptions;
the mode names do not establish neuronal causation.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "output/corrected_tms_analysis/20261126T112054_284547Z"
CORRECTED = ROOT / "reruns/aachen_skullstrip_skip_20261125"
sys.path[:0] = [str(ROOT / "pipeline/scripts/analysis"), str(ROOT / "setup")]
import tms_publication as tms
import tms_numit as tn
from analyse_tms_familywise import exact_max_t
from test_corrected_gwt_expectations import direct_hc3

MODES = {
    "storage": ["rtr", "xtx", "yty", "sts"],
    "transfer": ["xty", "ytx"],
    "copy": ["xtr", "ytr"],
    "erasure": ["rtx", "rty"],
    "joint_to_individual": ["str", "stx", "sty"],
    "individual_to_joint": ["rts", "xts", "yts"],
}
SOURCE_FILES = {}


def track(path):
    SOURCE_FILES[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def save_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, default=str) + "\n", encoding="utf-8")


def independent_atoms(sigma):
    """Independent logdet implementation and explicit triangular inversion."""
    def mi(a, b):
        def ld(ix):
            sign, value = np.linalg.slogdet(sigma[:, ix, :][:, :, ix])
            assert np.all(sign > 0)
            return value
        return .5 * (ld(a) + ld(b) - ld(sorted(a+b)))
    xa, xb, ya, yb = [mi([i], [j]) for i,j in [(0,2),(0,3),(1,2),(1,3)]]
    xya, xyb, xab, yab = mi([0,1],[2]), mi([0,1],[3]), mi([0],[2,3]), mi([1],[2,3])
    out = {"rtr": np.minimum.reduce([xa, xb, ya, yb])}
    def add(name, value, previous):
        out[name] = value - sum(out[key] for key in previous)
    add("rtx", np.minimum(xa, ya), ["rtr"])
    add("rty", np.minimum(xb, yb), ["rtr"])
    add("rts", np.minimum(xab, yab), ["rtr","rtx","rty"])
    add("xtr", np.minimum(xa, xb), ["rtr"])
    add("ytr", np.minimum(ya, yb), ["rtr"])
    add("str", np.minimum(xya, xyb), ["rtr","xtr","ytr"])
    add("xtx", xa, ["rtr","rtx","xtr"])
    add("xty", xb, ["rtr","rty","xtr"])
    add("ytx", ya, ["rtr","rtx","ytr"])
    add("yty", yb, ["rtr","rty","ytr"])
    add("stx", xya, ["rtr","rtx","xtr","xtx","ytr","ytx","str"])
    add("sty", xyb, ["rtr","rty","xtr","xty","ytr","yty","str"])
    add("xts", xab, ["rtr","rtx","rty","rts","xtr","xtx","xty"])
    add("yts", yab, ["rtr","rtx","rty","rts","ytr","ytx","yty"])
    add("sts", mi([0,1],[2,3]), list(out))
    return out


def modes(atoms):
    return {m: sum(atoms[k] for k in names) for m,names in MODES.items()}


def extraction(out):
    plan = pd.read_csv(track(CORE / "raw_paired_changes.csv"))
    scans = pd.read_csv(track(CORE / "raw_scan_metrics.csv"))
    ref = tms.load_reference(ROOT)
    for p in [*ref["input_files"], Path(__file__), Path(tn.pg.__file__), Path(tms.__file__),
              ROOT / "setup/analyse_tms_familywise.py", ROOT / "setup/test_corrected_gwt_expectations.py"]:
        track(p)
    release = track(ROOT / "output/aachen_alignment_review/corrected_cohort/cohort_release_receipt.json")
    assert hashlib.sha256(release.read_bytes()).hexdigest() == "d02401c42be2ea20c9e6de008a30c9665d92e7151467bb0aa4c873a8dfd81895"
    ids = set(zip(plan.subject_id, plan.baseline_session)) | set(zip(plan.subject_id, plan.followup_session))
    scans = scans.loc[[(r.subject_id,r.session) in ids for r in scans.itertuples()]]
    assert len(scans) == len(ids) == 140
    assert plan.groupby("dataset").size().to_dict() == {"aachen_itbs":30,"sudmex_tms":40}
    edges = np.array(np.triu_indices(454, k=1)).T
    i,j = edges.T
    ws, labels = ref["workspace"], ref["modules"]
    striatum = ref["labels"].Structure_Name.eq("Striatum").to_numpy()
    assert striatum.sum() == 20
    def between(a,b):
        return (a[i] & b[j]) | (b[i] & a[j])
    scopes = {"whole_brain": np.ones(len(i),bool), "within_workspace":ws[i]&ws[j],
              "outside_workspace":~ws[i]&~ws[j], "workspace_boundary":ws[i]!=ws[j],
              "FPN_striatum":between(labels==5,striatum),
              "FPN_DAN":between(labels==5,labels==2), "FPN_DMN":between(labels==5,labels==6)}
    assert np.all(sum(scopes[k].astype(int) for k in ["within_workspace","outside_workspace","workspace_boundary"]) == 1)
    save_json(out / "scope_pair_counts.json", {k:int(v.sum()) for k,v in scopes.items()})
    rows, checks = [], []
    for k,r in enumerate(scans.itertuples()):
        base = CORRECTED if r.dataset == "aachen_itbs" else ROOT
        path = track(base / "pipeline/results/timeseries" / f"{r.prefix}_timeseries_deconv.npy")
        ts = np.load(path)
        assert ts.shape[1] == 454 and np.isfinite(ts).all() and np.all(np.std(ts,axis=0) > 0)
        sigma = tn.pg.lagged_sigma_from_timeseries(ts, 1, edges)
        atoms = tn.pg.phiid_atoms(sigma)
        assert all(np.isfinite(v).all() for v in atoms.values())
        group = modes(atoms)
        np.testing.assert_allclose(sum(group.values()), atoms["tdmi"], atol=1e-10)
        sample = np.linspace(0,len(edges)-1,128,dtype=int)
        independent = independent_atoms(sigma[sample])
        err = max(float(np.max(abs(independent[a]-atoms[a][sample]))) for a in tn.pg.ATOM_NAMES)
        assert err < 2e-8, (r.prefix, "independent_atoms", err)
        sw = modes(tn.pg.phiid_atoms(sigma[sample][:,[1,0,3,2],:][:,:,[1,0,3,2]]))
        for name in MODES: np.testing.assert_allclose(group[name][sample],sw[name],atol=2e-8)
        errors = []
        for measure in tms.MATRIX_NAMES:
            old = np.load(track(base / "pipeline/results/phiid_matrices" / f"{r.prefix}_{measure}.npy"))[i,j]
            error = float(np.max(abs(atoms[measure]-old)))
            assert error < 2e-8, (r.prefix,measure,error)
            errors.append(error)
        checks.append(dict(prefix=r.prefix, independent_atom_error=err, max_saved_production_error=max(errors), n_timepoints=len(ts)))
        row = {n:getattr(r,n) for n in ["dataset","subject_id","session","arm","mean_fd"]}
        for scope, mask in scopes.items():
            for name, values in {**group, **{a:atoms[a] for a in tn.pg.ATOM_NAMES}, "tdmi":atoms["tdmi"]}.items():
                row[f"{scope}__{name}"] = float(values[mask].mean())
        rows.append(row)
        if (k+1)%10 == 0: print(f"Extracted {k+1}/140 scans; saved-output error {max(errors):.2g}",flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "scan_modes_and_atoms.csv",index=False)
    pd.DataFrame(checks).to_csv(out / "numerical_checks.csv",index=False)
    return frame,plan,list(scopes)


def inference(scans,plan,scopes,out):
    metrics = [f"{s}__{m}" for s in scopes for m in MODES]
    si = scans.set_index(["subject_id","session"])
    rows = []
    for r in plan.itertuples():
        a,b = si.loc[(r.subject_id,r.baseline_session)],si.loc[(r.subject_id,r.followup_session)]
        row = {n:getattr(r,n) for n in ["dataset","subject_id","arm","comparison"]}
        row["delta_fd"] = b.mean_fd-a.mean_fd
        for m in metrics:
            row[m+"_pre"],row[m+"_post"],row[m+"_delta"] = a[m],b[m],b[m]-a[m]
        rows.append(row)
    pairs = pd.DataFrame(rows)
    pairs.to_csv(out / "participant_changes.csv",index=False)
    aa = pairs.query("dataset == 'aachen_itbs'").groupby("subject_id").mean(numeric_only=True).sort_index()
    sud = pairs.query("dataset == 'sudmex_tms'").copy()
    sud["active"],sud["intercept"] = sud.arm.eq("active").astype(float),1.
    assert len(aa)==15 and len(sud)==40
    values = aa[[m+"_delta" for m in metrics]].to_numpy()
    pexact,pmax = exact_max_t(values)
    records = []
    for scope in scopes:
        columns = [m for m in metrics if m.startswith(scope+"__")]
        pscope = exact_max_t(aa[[m+"_delta" for m in columns]].to_numpy())[1]
        for j,m in enumerate(columns):
            k = metrics.index(m); x = aa[m+"_delta"]
            mean,se = x.mean(),stats.sem(x)
            ci = mean + np.array([-1,1])*stats.t.ppf(.975,14)*se
            records.append(dict(dataset="aachen_itbs",scope=scope,mode=m.split("__")[1],metric=m,model="paired_change",n=15,
                                estimate=mean,ci_low=ci[0],ci_high=ci[1],p=pexact[k],p_t=stats.ttest_1samp(x,0).pvalue,
                                p_maxT6=pscope[j],p_maxT42=pmax[k],pre_mean=aa[m+"_pre"].mean(),post_mean=aa[m+"_post"].mean()))
            for model, predictors in [("baseline_motion_adjusted",["intercept","active",m+"_pre","delta_fd"]),
                                      ("baseline_adjusted",["intercept","active",m+"_pre"])]:
                fit = direct_hc3(sud,m+"_post",predictors,"active")
                # Independently verify coefficients/uncertainty through statsmodels.
                check = sm.OLS(sud[m+"_post"],sud[predictors]).fit(cov_type="HC3",use_t=True)
                np.testing.assert_allclose([fit["estimate"],fit["se"],fit["p"]],
                                          [check.params.active,check.bse.active,check.pvalues.active],rtol=1e-7,atol=1e-10)
                records.append(dict(dataset="sudmex_tms",scope=scope,mode=m.split("__")[1],metric=m,model=model,
                                    active_change=sud.loc[sud.active.eq(1),m+"_delta"].mean(),
                                    sham_change=sud.loc[sud.active.eq(0),m+"_delta"].mean(),**fit))
    results = pd.DataFrame(records)
    results["p_holm6"] = results.groupby(["dataset","model","scope"]).p.transform(lambda p: multipletests(p,method="holm")[1])
    results["p_holm42"] = results.groupby(["dataset","model"]).p.transform(lambda p: multipletests(p,method="holm")[1])
    results.to_csv(out / "mode_treatment_results.csv",index=False)
    primary = results.loc[~results.model.eq("baseline_adjusted")]
    conjunction=[]
    for m in metrics:
        a = primary.loc[primary.dataset.eq("aachen_itbs") & primary.metric.eq(m)].iloc[0]
        b = primary.loc[primary.dataset.eq("sudmex_tms") & primary.metric.eq(m)].iloc[0]
        same = np.sign(a.estimate)==np.sign(b.estimate)
        conjunction.append(dict(metric=m,scope=a.scope,mode=a["mode"],same_direction=bool(same),
                                aachen_estimate=a.estimate,sudmex_estimate=b.estimate,p_aachen=a.p,p_sudmex=b.p,
                                p=max(a.p,b.p) if same else 1.))
    cc=pd.DataFrame(conjunction)
    cc["p_holm6"] = cc.groupby("scope").p.transform(lambda p:multipletests(p,method="holm")[1])
    cc["p_holm42"] = multipletests(cc.p,method="holm")[1]
    cc.to_csv(out / "cross_dataset_conjunction.csv",index=False)
    print("PRIMARY WORKSPACE BOUNDARY",flush=True)
    print(cc.loc[cc.scope.eq("workspace_boundary")].to_string(index=False),flush=True)
    print("COMPLETE 42-ENDPOINT VIEW sorted by conjunction p",flush=True)
    print(cc.sort_values("p").head(12).to_string(index=False),flush=True)
    return results,cc


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--reuse",type=Path);args=parser.parse_args()
    out=args.reuse or ROOT / "output/pid_processing_modes" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    out.mkdir(parents=True,exist_ok=True)
    print("OUTPUT:",out,flush=True)
    assert sorted(a for names in MODES.values() for a in names)==sorted(tn.pg.ATOM_NAMES)
    spec={"created_utc":datetime.now(timezone.utc).isoformat(),"status":"Exploratory follow-up; not preregistered. Specification saved before outcome extraction.",
          "question":"Does prefrontal stimulation alter specific modes of predictive information processing consistently across cohorts?",
          "modes":MODES,"reference":"Mediano et al. PNAS 2025; doi:10.1073/pnas.2423297122",
          "primary_scope":"Workspace-boundary pairs, motivated by the requested workspace versus outside question. All six modes form one family.",
          "additional_scopes":["whole_brain","within_workspace","outside_workspace","FPN_striatum","FPN_DAN","FPN_DMN"],
          "families":"All six modes per scope (Holm6; Aachen maxT6); all 42 outcomes (Holm42; Aachen maxT42). Conjunctions max(two two-sided p), p=1 if signs differ; Holm6 and Holm42. Baseline-only SUDMEX sensitivity retains identical families.",
          "contrasts":"Aachen 15 independent people: average two complete visit changes. SUDMEX 22 active/18 sham: followup ~ arm + baseline + deltaFD, HC3/t.",
          "normalization":"Raw nats. No division by total mutual information and no summation of conditional percentiles. Any future normalized-mode analysis requires its own matched null calibration.",
          "prohibitions":"No HRF association models, no inference of stability from nonsignificance, no raw/normalized GWT rerouting or preprocessing change, no HCP/psychedelic recomputation.",
          "validation":"Every pair's five saved production composites; sampled independent full 16-atom reconstruction with LAPACK determinants; source-swap invariance; six-mode sum equals TDMI; HC3 independent implementations; unchanged source hashes."}
    if args.reuse:
        frame=pd.read_csv(out/"scan_modes_and_atoms.csv");plan=pd.read_csv(CORE/"raw_paired_changes.csv")
        scopes=list(json.loads((out/"scope_pair_counts.json").read_text()))
        old=json.loads((out/"source_sha256.json").read_text());SOURCE_FILES.update(old)
    else:
        save_json(out/"analysis_specification.json",spec)
        frame,plan,scopes=extraction(out)
        save_json(out/"source_sha256.json",SOURCE_FILES)
    results,conjunction=inference(frame,plan,scopes,out)
    for relative,digest in SOURCE_FILES.items():
        assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==digest,relative
    save_json(out/"verification.json",{"passed":True,"scans":len(frame),"mode_results":len(results),"conjunctions":len(conjunction),
                                     "source_files_unchanged":len(SOURCE_FILES),"checks":"Full numerical check table and independent estimators passed.",
                                     "shared_holm6_results":int((conjunction.p_holm6<.05).sum()),"shared_holm42_results":int((conjunction.p_holm42<.05).sum())})
    (out.parent/"LATEST.txt").write_text(out.name+"\n",encoding="utf-8")
    print("COMPLETE",out,flush=True)


if __name__ == "__main__":
    with threadpool_limits(limits=2): main()
