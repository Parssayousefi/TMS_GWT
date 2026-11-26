# -*- coding: utf-8 -*-
"""Final gate: every machine-checkable number AND every named artefact in
MANUSCRIPT_STRUCTURE.md, against the freshly regenerated output folder."""
import re, sys
import numpy as np, pandas as pd
from pathlib import Path
P = Path(r"D:\projects\TMS_GWT\pipeline")
N = P / "results" / "group_analysis" / "manuscript_tables_uniform"
DOC = (P / "manuscript" / "MANUSCRIPT_STRUCTURE.md").read_text(encoding="utf-8")
def rd(f): return pd.read_csv(N / f, comment="#")
fails, checks = [], []
def chk(label, got, want, tol=5e-4):
    if isinstance(want, (str, list, tuple)):
        ok = got == want
    elif want != 0 and abs(want) < 1e-3:
        ok = abs(got - want) <= abs(want) * .05
    else:
        ok = abs(got - want) <= tol
    checks.append((label, ok, got, want)); fails.append(label) if not ok else None

PSY, STIM, PROP = "classic psychedelic", "stimulant/entactogen", "propofol"
cl = rd("supp_class_level_tests.csv")

# ---- Claim 2: workspace-mean -------------------------------------------
wm = cl[cl.block == "workspace-mean PC shift"]
def wq(k, m): return float(wm[(wm.drug_class == k) & (wm.measure == m)].q_class.iloc[0])
chk("WM psy red 4.6e-30", wq(PSY, "redundancy"), 4.64e-30)
chk("WM stim red 2.4e-7", wq(STIM, "redundancy"), 2.37e-7)
chk("WM prop red .0056",  wq(PROP, "redundancy"), .00563, 5e-5)
chk("WM psy syn 3.1e-9",  wq(PSY, "synergy"), 3.11e-9)
chk("WM stim syn .75",    wq(STIM, "synergy"), .748, 5e-3)
chk("WM prop syn 1.0",    wq(PROP, "synergy"), 1.0)
chk("WM red sig 3 classes", int((wm[wm.measure=="redundancy"].q_class < .05).sum()), 3)
chk("WM syn sig 1 class",   int((wm[wm.measure=="synergy"].q_class < .05).sum()), 1)
chk("WM prop syn sign-inconsistent",
    str(wm[(wm.drug_class==PROP)&(wm.measure=="synergy")].consistent_sign.iloc[0]), "False")

# ---- Claim 2: network ----------------------------------------------------
nw = cl[(cl.block=="network x role") & (cl.analysis_tier=="primary")]
fpn = nw[(nw.drug_class==PSY)&(nw.Network=="Cont")&(nw.Role=="broadcaster")]
dmn = nw[(nw.drug_class==PSY)&(nw.Network=="Default")&(nw.Role=="gateway")]
chk("primary families 6", len(nw), 6)
chk("FPN q 2.2e-11", float(fpn.q_class.iloc[0]), 2.25e-11)
chk("FPN delta -4.32", float(fpn.mean_delta.iloc[0]), -4.322, 5e-3)
chk("FPN 4 of 4", int(fpn.n_arms_sig_uncorrected.iloc[0]), 4)
chk("DMN q .018", float(dmn.q_class.iloc[0]), .0184, 5e-4)
chk("DMN delta +1.32", float(dmn.mean_delta.iloc[0]), 1.322, 5e-3)
chk("DMN 2 of 4", int(dmn.n_arms_sig_uncorrected.iloc[0]), 2)
chk("primary survivors 2", int((nw.q_class < .05).sum()), 2)
ex = cl[(cl.block=="network x role") & (cl.analysis_tier=="exploratory")]
chk("exploratory families 42", len(ex), 42)
chk("exploratory survivors 12", int((ex.q_class < .05).sum()), 12)
for k, n in ((PSY, 8), (STIM, 4), (PROP, 0)):
    chk(f"exploratory {k} {n}", int((ex[ex.drug_class==k].q_class < .05).sum()), n)
chk("psy FPN gateway 1.2e-10",
    float(ex[(ex.drug_class==PSY)&(ex.Network=="Cont")&(ex.Role=="gateway")].q_class.iloc[0]),
    1.24e-10)

# ---- Claim 1 -------------------------------------------------------------
t1 = rd("table1_cross_dataset_replication.csv")
chk("workspace 235", int(t1.N_workspace.iloc[0]), 235)
chk("HCP gw/bc", (int(t1.N_gateways.iloc[0]), int(t1.N_broadcasters.iloc[0])), (99, 136))
ctrl = t1[t1.Dataset != "hcp"]
chk("ctrl DMN gw", tuple(ctrl.DMN_GW_pct), (95.3, 95.3, 96.9))
chk("ctrl DAN bc", tuple(ctrl.DAN_BC_pct), (100.0, 100.0, 93.3))
t1c = rd("table1c_saphire_drug_arms.csv")
chk("drug DMN range", (float(t1c.DMN_GW_pct.min()), float(t1c.DMN_GW_pct.max())), (93.8, 100.0))

# ---- Claim 3 -------------------------------------------------------------
cs = rd("supp_rsa_class_separation.csv"); raw = cs[cs.scale=="raw"].set_index("measure")
chk("CSI syn .552", float(raw.loc["Synergy PC","csi"]), .5519, 1e-3)
chk("CSI red .555", float(raw.loc["Redundancy PC","csi"]), .5552, 1e-3)
chk("CSI fc .382",  float(raw.loc["Mean FC","csi"]), .3820, 1e-3)
po = rd("supp_rsa_class_separation_pooled.csv"); pr0 = po[po.scale=="raw"]
chk("pooled +.172", float(pr0.difference.iloc[0]), .1715, 1e-3)
chk("pooled 2/126", str(pr0.exact_two.iloc[0]), "2/126")

# ---- 13.6d study confound (new) -----------------------------------------
sc = rd("supp_rsa_study_confound.csv")
sub = sc[sc.block == "csi by pair subset"].set_index("subset")
chk("cross-study contrast +.078", float(sub.loc["cross-study only","pooled_contrast"]), .0783, 1e-3)
chk("cross-study exact 48/126", str(sub.loc["cross-study only","exact_two"]), "48/126")
chk("cross-study csi syn .422", float(sub.loc["cross-study only","csi_synergy"]), .4217, 1e-3)
chk("cross-study csi red .500", float(sub.loc["cross-study only","csi_redundancy"]), .4998, 1e-3)
chk("cross-study csi fc .382",  float(sub.loc["cross-study only","csi_meanfc"]), .3825, 1e-3)
con = sc[sc.block == "same-study contribution"].set_index("measure")
chk("same-study adds syn +.130", float(con.loc["Synergy PC","same_study_contribution"]), .1302, 1e-3)
chk("same-study adds red +.055", float(con.loc["Redundancy PC","same_study_contribution"]), .0554, 1e-3)
chk("same-study adds fc -.0005", round(float(con.loc["Mean FC","same_study_contribution"]), 4), -0.0005)
pur = sc[sc.block == "nearest-neighbour purity"].set_index("measure")
for m, v in (("Synergy PC","9/9"), ("Redundancy PC","9/9"), ("Mean FC","8/9")):
    chk(f"purity {m}", str(pur.loc[m,"purity"]), v)

# ---- coverage control: joint null (new, primary) -------------------------
cvj = rd("supp_rsa_receptor_coverage_control.csv")
jn = cvj[cvj.block == "affinity beyond coverage (joint null)"].set_index("measure")
chk("joint syn partial 46/720", str(jn.loc["Synergy PC","partial_exact"]), "46/720")
chk("joint syn partial p .064", float(jn.loc["Synergy PC","partial_p"]), .0639, 1e-3)
chk("joint syn MRM p .192",     float(jn.loc["Synergy PC","mrm_p"]), .1917, 1e-3)
chk("joint red partial p .650", float(jn.loc["Redundancy PC","partial_p"]), .6500, 1e-3)
rel = rd("supp_rsa_reliability.csv"); mm = rel.groupby("measure").reliability.mean()
chk("rel syn .370", float(mm["Synergy PC"]), .3701, 1e-3)
chk("rel red .520", float(mm["Redundancy PC"]), .5197, 1e-3)
chk("rel fc .639",  float(mm["Mean FC"]), .63904, 1e-4)
ARMS = ["basel_lam/lsd","basel_lpm/lsd","basel_lpm/psilocybin","basel_lpm/mescaline",
        "basel_lam/mdma","basel_lam/amph","cambridge_propofol_uniform/restdeep",
        "cambridge_propofol_uniform/restlight","cambridge_propofol_uniform/restrecovery"]
SERO = set(ARMS[:4])
for m, f, want in (("syn","supp_synergy_r_matrix.csv",9),
                   ("red","supp_redundancy_r_matrix.csv",9),
                   ("fc","supp_meanfc_cross_cohort_r_matrix.csv",8)):
    R = pd.read_csv(N/f, index_col=0).loc[ARMS, ARMS].values.copy()
    np.fill_diagonal(R, np.nan)
    chk(f"NN purity {m}", sum((ARMS[int(np.nanargmax(R[i]))] in SERO)==(ARMS[i] in SERO)
                              for i in range(9)), want)

# ---- Claim 4 -------------------------------------------------------------
m6 = rd("supp_rsa_mantel6_both.csv")
def mv(e,m,c): return m6[(m6.external==e)&(m6.measure==m)][c].iloc[0]
for e, m, r_, k in (("ASC (shape)","Synergy PC",.8955,"3/720"),
                    ("Receptor (shape)","Synergy PC",.8124,"12/720"),
                    ("ASC (shape)","Redundancy PC",.73047,"29/720"),
                    ("Receptor (shape)","Redundancy PC",.6733,"42/720"),
                    ("ASC (shape)","Mean FC",.4686,"54/720"),
                    ("Receptor (shape)","Mean FC",.2914,"218/720")):
    chk(f"{m} x {e[:3]} rho", float(mv(e,m,"mantel_r")), r_, 1e-3)
    chk(f"{m} x {e[:3]} exact", str(mv(e,m,"exact_two")), k)
mx = rd("supp_rsa_maxT.csv")
chk("maxT syn x ASC", str(mx[mx.test.str.contains("Synergy") & mx.test.str.contains("ASC")].exact_maxT.iloc[0]), "3/720")
ad = rd("supp_rsa_anchor_dependence.csv")
cb = ad[ad.block=="class-indicator anchor (ceiling-bound)"].set_index("quantity")
chk("class ceiling .859", float(cb.loc["ARITHMETIC CEILING","value"]), .8593, 1e-3)
chk("class anchor syn .824", float(cb.loc["Synergy PC x class","value"]), .8235, 1e-3)
chk("class anchor red = ceiling", float(cb.loc["Redundancy PC x class","value"]), .8593, 1e-3)
chk("class anchor fc .609", float(cb.loc["Mean FC x class","value"]), .6087, 1e-3)
cr = ad[ad.block=="class-restricted null"].set_index("quantity")
for m, k in (("Synergy PC","3/48"),("Redundancy PC","29/48"),("Mean FC","25/48")):
    chk(f"restricted {m} x ASC", str(cr.loc[f"{m} x ASC (shape)","exact_restricted"]), k)
dep = ad[ad.block=="anchor dependence"].set_index("quantity")
chk("rho(ASC,REC) .875", float(dep.loc["rho(ASC, receptor)","value"]), .8751, 1e-3)
chk("rho partial .515", float(dep.loc["rho(ASC, receptor | class)","value"]), .5151, 1e-3)
chk("rho within -.162",
    float(dep.loc["rho(ASC, receptor) within serotonergic class","value"]), -.1622, 1e-3)
cov = rd("supp_rsa_receptor_coverage_control.csv")
cp = cov[cov.block=="affinity beyond coverage"].set_index("measure")
chk("cov syn .613", float(cp.loc["Synergy PC","rho_partial"]), .613, 2e-3)
chk("cov red .174", float(cp.loc["Redundancy PC","rho_partial"]), .174, 2e-3)
chk("cov fc .035",  float(cp.loc["Mean FC","rho_partial"]), .035, 2e-3)
wc = rd("supp_rsa_within_class.csv")
ws = wc[(wc.design.str.startswith("serotonergic"))&(wc.measure=="Synergy PC")&
        (wc.external=="ASC (shape)")]
chk("within-class syn .718", float(ws.mantel_r.iloc[0]), .7182, 1e-3)

# ---- controls ------------------------------------------------------------
ac = rd("supp_autocorrelation_control.csv")
chk("dAC1 sig 0/9", int((ac[ac.block=="per-arm autocorrelation"].q_ac1 < .05).sum()), 0)
adj = ac[ac.block=="PC effect adjusted for dAC1"]
chk("dAC1 flips 0/18", int(adj["flips"].astype(str).str.lower().eq("true").sum()), 0)
mp = rd("supp_motion_pc_coupling.csv")
chk("dFD flips 0/18", int(mp["flips"].astype(str).str.lower().eq("true").sum()), 0)
t3 = rd("table3_significant_role_switches.csv")
chk("role switches 74", len(t3), 74)
chk("MDMA+propofol switches 0",
    int(t3.drug_task.isin(["mdma","restdeep","restlight","restrecovery"]).sum()), 0)

# ---- no Bayes factor survives anywhere -----------------------------------
bfcols = []
for f in N.glob("*.csv"):
    try: c = pd.read_csv(f, comment="#", nrows=1).columns.tolist()
    except Exception: continue
    if any("BF" in x or "Bayes" in x for x in c): bfcols.append(f.name)
chk("no BF column in any CSV", bfcols, [])
chk("no BF/Bayes in the document",
    [l for l in DOC.splitlines() if "BF" in l or "Bayes" in l], [])

# ---- every artefact the document names actually exists -------------------
named = set(re.findall(r"`(supp_[a-z0-9_]+\.csv|table\d[a-z]?_[a-z0-9_]+\.csv)`", DOC))
missing = sorted(f for f in named if not (N / f).exists())
chk("named CSVs all exist", missing, [])
figs = set(re.findall(r"`(fig[A-Za-z0-9_]+)`", DOC))
FF = P / "manuscript" / "figures_final"
man = rd("supp_figure_manifest.csv")
# 12 staged labels x 2 formats. Fig2A + Fig2B were merged into a single Fig2
# display item for the Nature Communications 10-item cap (2026-08-16).
chk("figure manifest 24 files", len(man), 24)
chk("figures_final has 24", len(list(FF.glob("*.p*g")) + list(FF.glob("*.pdf"))), 24)
chk("Fig2 staged as one item", sorted(man.label.unique().tolist()),
    ["Fig1", "Fig2", "Fig3", "Fig4", "FigS1", "FigS2", "FigS3", "FigS4",
     "FigS5", "FigS6", "FigS7", "FigS8"])

for label, ok, got, want in checks:
    print(f"  {'PASS' if ok else 'FAIL':4s}  {label:36s} got={str(got)[:46]}  want={str(want)[:32]}")
print(f"\n{len(checks)} checks, {len(fails)} failures")
if fails: print("FAILED:", fails); sys.exit(1)
