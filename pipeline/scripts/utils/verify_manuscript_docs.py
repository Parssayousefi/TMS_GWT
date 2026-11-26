# -*- coding: utf-8 -*-
"""Final gate on the BUILT manuscript: every reported number in
Methods_and_Results.docx and Supplementary_Material.docx, checked against the
CSVs the notebook regenerated in the same execution.

Reads the .docx (not the markdown), so pandoc mangling is caught too.
"""
import re, sys
# console may be cp1252; the report contains unicode (approx sign, en dash)
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass
import pandas as pd
from docx import Document
from pathlib import Path

P = Path(r"D:\projects\TMS_GWT\pipeline")
N = P / "results" / "group_analysis" / "manuscript_tables_uniform"
M = P / "manuscript"


def rd(f):
    return pd.read_csv(N / f, comment="#")


def text(docx):
    d = Document(str(docx))
    parts = [p.text for p in d.paragraphs]
    for t in d.tables:
        for row in t.rows:
            parts.extend(c.text for c in row.cells)
    # normalise the unicode the APA pipeline emits
    s = "\n".join(parts)
    for a, b in (("\u2212", "-"), ("\u2013", "-"), ("\u2014", "-"),
                 ("\u00a0", " "), ("\u2009", " "), ("\u202f", " ")):
        s = s.replace(a, b)
    return re.sub(r"[ \t]+", " ", s)


MAIN = text(M / "Methods_and_Results.docx")
SUPP = text(M / "Supplementary_Material.docx")
BOTH = MAIN + "\n" + SUPP

fails, checks = [], []


def has(label, needle, hay=None, where="main"):
    hay = {"main": MAIN, "supp": SUPP, "both": BOTH}[where] if hay is None else hay
    ok = needle in hay
    checks.append((label, ok, needle[:58], where))
    if not ok:
        fails.append(label)


def num(x, dp=3, sign=False):
    """.552 / -.003 / +.172 formatting used throughout the manuscript."""
    s = f"{abs(float(x)):.{dp}f}"
    s = s[1:] if s.startswith("0.") else s
    if float(x) < 0:
        return "-" + s
    return ("+" + s) if sign else s


# ------------------------------------------------------- Claim 1
t1 = rd("table1_cross_dataset_replication.csv")
has("workspace 235 of 454", f"{int(t1.N_workspace.iloc[0])} of 454")
has("gateways/broadcasters",
    f"{int(t1.N_gateways.iloc[0])} were gateways and {int(t1.N_broadcasters.iloc[0])} broadcasters")
ctrl = t1[t1.Dataset != "hcp"]
has("control DMN gateway %", ", ".join(f"{v}%" for v in ctrl.DMN_GW_pct).replace(
    f", {list(ctrl.DMN_GW_pct)[-1]}%", f" and {list(ctrl.DMN_GW_pct)[-1]}%"))
t1c = rd("table1c_saphire_drug_arms.csv")
_lo, _hi = float(t1c.DMN_GW_pct.min()), float(t1c.DMN_GW_pct.max())
# accept "100%" as well as "100.0%" -- both render the same value
_hi_s = f"{_hi:g}"
has("drug DMN range", f"{_lo:g}% to {_hi_s}%")

# ------------------------------------------------------- Claim 2
cl = rd("supp_class_level_tests.csv")
PSY, STIM, PROP = "classic psychedelic", "stimulant/entactogen", "propofol"
wm = cl[cl.block == "workspace-mean PC shift"]


def wq(k, m):
    return float(wm[(wm.drug_class == k) & (wm.measure == m)].q_class.iloc[0])


# Prose reports class-level q below 1e-3 as an order of magnitude; the exact
# value lives in Figure 2c and Supplementary Table S1a. Assert BOTH: that the
# rounded form the text uses matches the computed exponent, and that no
# mantissa survived in the main text.
# In the .docx text layer the italic q has no asterisks and the exponent is a
# superscript run, so `*q* ~ 10^-30^` in the markdown extracts as `q ~ 10-30`
# (the minus is normalised to a hyphen by text() above).
import math


def _oom(q):
    return f"10-{-int(math.floor(math.log10(q)))}"


for lbl, key, meas in (("WM psy red oom", PSY, "redundancy"),
                       ("WM psy syn oom", PSY, "synergy"),
                       ("WM stim red oom", STIM, "redundancy")):
    has(lbl, f"q \u2248 {_oom(wq(key, meas))}")
has("WM prop red .0056", ".0056")
has("WM stim syn .75", ".75")
_mant = re.findall(r"q = [0-9]\.[0-9] \u00d7 10", MAIN)
checks.append(("no mantissa left in main prose", not _mant, str(_mant), "main"))
if _mant:
    fails.append("no mantissa left in main prose")
nw = cl[(cl.block == "network x role") & (cl.analysis_tier == "primary")]
fpn = nw[(nw.drug_class == PSY) & (nw.Network == "Cont") & (nw.Role == "broadcaster")].iloc[0]
dmn = nw[(nw.drug_class == PSY) & (nw.Network == "Default") & (nw.Role == "gateway")].iloc[0]
has("FPN mean delta", f"{abs(float(fpn.mean_delta)):.2f} regions")
has("FPN arms 4 of 4", f"{int(fpn.n_arms_sig_uncorrected)} of 4 psychedelic arms")
has("DMN mean delta", f"{float(dmn.mean_delta):.2f} regions")
has("DMN q .018", num(float(dmn.q_class), 3))
has("DMN arms 2 of 4", f"{int(dmn.n_arms_sig_uncorrected)} of 4 arms")
ex = cl[(cl.block == "network x role") & (cl.analysis_tier == "exploratory")]
has("exploratory 42 families", f"{len(ex)} exploratory families")
has("exploratory 12 survive", f"{int((ex.q_class < .05).sum())} survived correction")
for k, lab in ((PSY, "psychedelic"), (STIM, "stimulant/entactogen"), (PROP, "propofol")):
    n_ = int((ex[ex.drug_class == k].q_class < .05).sum())
    has(f"exploratory {lab} {n_}/14", f"{n_} of 14 {lab}" if k != PROP else "14 propofol families")

# ------------------------------------------------------- Claim 3
cs = rd("supp_rsa_class_separation.csv")
raw = cs[cs.scale == "raw"].set_index("measure")
has("CSI triple", f"synergy {num(raw.loc['Synergy PC','csi'])}, redundancy "
                  f"{num(raw.loc['Redundancy PC','csi'])}, and mean functional connectivity "
                  f"{num(raw.loc['Mean FC','csi'])}")
po = rd("supp_rsa_class_separation_pooled.csv")
pr = po[po.scale == "raw"].iloc[0]
has("pooled contrast", f"{num(pr.difference, 3, sign=True)} (*p* = {pr.exact_two}"
    .replace("*p*", "p").replace("(p = ", "(p = "))
has("pooled contrast value", num(pr.difference, 3, sign=True))
has("pooled exact", str(pr.exact_two))
sc = rd("supp_rsa_study_confound.csv")
sub = sc[sc.block == "csi by pair subset"].set_index("subset")
cx = sub.loc["cross-study only"]
has("cross-study contrast", num(cx.pooled_contrast, 3, sign=True))
has("cross-study exact", str(cx.exact_two))
has("cross-study syn/red", f"({num(cx.csi_redundancy)}) in fact exceeds cross-study synergy "
                           f"({num(cx.csi_synergy)})")
con = sc[sc.block == "same-study contribution"].set_index("measure")
has("same-study syn", num(con.loc["Synergy PC", "same_study_contribution"], 3, sign=True))
has("same-study red", num(con.loc["Redundancy PC", "same_study_contribution"], 3, sign=True))
has("same-study fc", num(con.loc["Mean FC", "same_study_contribution"], 4))
rel = rd("supp_rsa_reliability.csv").groupby("measure").reliability.mean()
has("reliability triple", f"{num(rel['Synergy PC'])} for synergy and {num(rel['Redundancy PC'])} "
                          f"for redundancy against {num(rel['Mean FC'])} for connectivity")

# ------------------------------------------------------- Claim 4
m6 = rd("supp_rsa_mantel6_both.csv")


def mv(e, m, c):
    return m6[(m6.external == e) & (m6.measure == m)][c].iloc[0]


for e, m, lbl in (("ASC (shape)", "Synergy PC", "syn x ASC"),
                  ("Receptor (shape)", "Synergy PC", "syn x REC"),
                  ("ASC (shape)", "Redundancy PC", "red x ASC"),
                  ("Receptor (shape)", "Redundancy PC", "red x REC"),
                  ("ASC (shape)", "Mean FC", "fc x ASC"),
                  ("Receptor (shape)", "Mean FC", "fc x REC")):
    has(f"{lbl} rho", num(mv(e, m, "mantel_r")))
    has(f"{lbl} exact", str(mv(e, m, "exact_two")))
ad = rd("supp_rsa_anchor_dependence.csv")
dep = ad[ad.block == "anchor dependence"].set_index("quantity")
has("anchor rho .875", num(dep.loc["rho(ASC, receptor)", "value"]))
has("anchor partial .515", num(dep.loc["rho(ASC, receptor | class)", "value"]))
has("anchor within -.162",
    num(dep.loc["rho(ASC, receptor) within serotonergic class", "value"]))
cb = ad[ad.block == "class-indicator anchor (ceiling-bound)"].set_index("quantity")
has("class ceiling .859", num(cb.loc["ARITHMETIC CEILING", "value"]), where="both")
cr = ad[ad.block == "class-restricted null"].set_index("quantity")
has("restricted syn 3/48", str(cr.loc["Synergy PC x ASC (shape)", "exact_restricted"]),
    where="both")
cov = rd("supp_rsa_receptor_coverage_control.csv")
fx = cov[cov.block == "affinity beyond coverage"].set_index("measure")
jn = cov[cov.block == "affinity beyond coverage (joint null)"].set_index("measure")
has("coverage syn partial", num(fx.loc["Synergy PC", "rho_partial"]), where="both")
has("coverage red partial", num(fx.loc["Redundancy PC", "rho_partial"]), where="both")
has("coverage fc partial", num(fx.loc["Mean FC", "rho_partial"]), where="both")
has("joint syn partial exact", str(jn.loc["Synergy PC", "partial_exact"]), where="both")
has("joint syn mrm exact", str(jn.loc["Synergy PC", "mrm_exact"]), where="both")
mx = rd("supp_rsa_maxT.csv")
has("maxT syn x ASC",
    str(mx[mx.test.str.contains("Synergy") & mx.test.str.contains("ASC")].exact_maxT.iloc[0]))

# ------------------------------------------------------- controls
has("motion 0 of 18", "none** of the 18 verdicts".replace("**", "") if False else "18 verdicts")
has("autocorr 0 of 9", "0 of 9")
t3 = rd("table3_significant_role_switches.csv")
has("role switches 74", f"{len(t3)} significant", where="both")

# ------------------------------------------------ role switches in the narrative
t3 = t3.assign(arm=t3.dataset + "/" + t3.drug_task)
PSY_ARMS = {"basel_lam/lsd", "basel_lpm/lsd", "basel_lpm/psilocybin", "basel_lpm/mescaline"}
BG, GB = "Broadcaster -> Gateway", "Gateway -> Broadcaster"
cont = t3[t3.RSN_name == "Cont"]
vis = t3[t3.RSN_name == "Visual"]
has("switches total in 4.2", f"Of the 235 workspace regions, {len(t3)}")
has("switches psychedelic", f"{int(t3.arm.isin(PSY_ARMS).sum())} of them in the four classic")
has("switches none MDMA/propofol", "none at all under MDMA or at any depth of propofol")
has("Cont unidirectional",
    f"{int((cont.switch_dir == BG).sum())} of {len(cont)} broadcaster to gateway")
has("Visual direction",
    f"{int((vis.switch_dir == GB).sum())} of {len(vis)} gateway to broadcaster")

# ------------------------------------------ orphaned supplementary tables cited
for s in ("S9", "S13", "S14", "S16", "S19"):
    has(f"cites Supplementary Table {s}", f"Supplementary Table {s}")

# Main-text Table 2 was dropped in favour of Figure 2c; tables renumbered 1-3.
# docx text layer: the caption is bold-formatted, so there are no asterisks
_defs = re.findall(r"^Table (\d+)\.", MAIN, re.M)
checks.append(("main tables numbered 1-3", _defs == ["1", "2", "3"], str(_defs), "main"))
if _defs != ["1", "2", "3"]:
    fails.append("main tables numbered 1-3")
_refs = set(re.findall(r"Table (\d+)\b", MAIN))
checks.append(("no reference above Table 3", _refs <= {"1", "2", "3"}, str(sorted(_refs)), "main"))
if not _refs <= {"1", "2", "3"}:
    fails.append("no reference above Table 3")
has("Figure 2 has a panel c", "The class-level inference itself")

# Claim 2's headline used to report class q with no magnitude. The class q is a
# Fisher combination, so there is no single combined effect size; the member
# arms' g_z RANGE is reported instead. Hedges' g can exceed 1, so APA keeps the
# leading zero and num() must not be used -- it strips it.
_f2d = rd("supp_fig2_dissociation_data.csv")
_grp = {"Serotonergic psychedelic": "psychedelic", "Empathogen": "stimulant/entactogen",
        "Stimulant": "stimulant/entactogen", "Anaesthetic": "propofol"}
_f2d = _f2d.assign(grp=_f2d["class"].map(_grp))
for _col, _lbl in (("red_g", "redundancy"), ("syn_g", "synergy")):
    for _g, _s in _f2d.groupby("grp"):
        _lo, _hi = f"{_s[_col].min():.2f}", f"{_s[_col].max():.2f}"
        # text() normalises en dash and U+2212 to ASCII hyphen, so the needle
        # must use hyphens even though the manuscript renders typographic dashes.
        # Negative lower bounds are written "-0.05 to 0.47", positives "a-b".
        _needle = (f"{_lo} to {_hi}" if _s[_col].min() < 0 else f"{_lo}-{_hi}")
        has(f"Claim 2 {_lbl} g range, {_g}", _needle)

# Figure 3 panel b became the within-study matrices. Three things must hold
# together or the caption describes a figure we are not shipping: that the
# matrices are computed within study, that mean FC is on the same 235 workspace
# regions as the decomposed measures (the only reason a column is comparable),
# and the corroborative-not-inferential caveat, which is the honest limit of a
# three-arm block. Needles avoid italics, subscripts and en dashes, none of
# which survive the .docx text layer.
has("Figure 3b says within each study",
    "each study so that scanner, protocol and control condition")
has("Figure 3b says mean FC is workspace-scoped",
    "scored on the same 235 workspace regions")
has("Figure 3b keeps the corroborative caveat", "corroborates")
# Table S3 used to say propofol's cross-cohort similarity was "displayed in
# Figure 3", which was true only while panel b held the 9x9 cross-cohort
# matrices. Panel b is now within-study, so its propofol column contains only
# the three sedation levels and no pharmacological arm.
has("Table S3 repointed off the old panel b",
    "contains only the three sedation levels", where="supp")
has("Figure 3b says no exact null at three arms", "no exact null is available")

# The Table S19 note used to point at main-text Table 2, which was dropped in
# favour of Figure 2c. Assert the repointed wording, and that no supplementary
# cross-reference names a main-text table that no longer exists.
has("S19 note repointed", "loss reported in the main text", where="supp")

# The notebook stages all 12 figure labels; the Word supplement reproduces two
# and renumbers them, so caption number != staged filename. Assert the embedded
# images are the INTENDED sources, via the manifest -- this catches the failure
# mode where renumbering silently swaps in a figure that was cut.
_supp_md = (M / "Supplementary_Material_NatComm.md").read_text(encoding="utf-8")
_emb = re.findall(r"figures_final/(Fig[A-Za-z0-9]+)\.png", _supp_md)
_man = rd("supp_figure_manifest.csv")
_man = _man[_man.ext == "png"].set_index("label").source.to_dict()
_want = {"fig2_panelD_switcher_anatomy.png",
         "fig2_panelC_workspace_gateway_broadcaster_brains.png"}
_got = {_man.get(l, "MISSING") for l in _emb}
checks.append(("supp embeds the intended figures", _got == _want,
               str(sorted(_got)), "supp"))
if _got != _want:
    fails.append("supp embeds the intended figures")

# ---- the two label corrections (audit, 2026-08-16) -------------------------
# The 74 role switches come from a directional single-measure rule under
# workspace FDR, not from McNemar's test. Verify the claim against the data as
# well as the wording: min(syn_p, red_p) < .05 for every row, while the stored
# p_corrected column is the max and exceeds .05 for many.
_sw = rd("table3_significant_role_switches.csv")
_mn = _sw[["syn_p", "red_p"]].min(axis=1)
_mx = _sw[["syn_p", "red_p"]].max(axis=1)
checks.append(("role-switch rule is single-measure", bool((_mn < .05).all()),
               f"max of per-row min p = {_mn.max():.4f}", "main"))
if not (_mn < .05).all():
    fails.append("role-switch rule is single-measure")
checks.append(("p_corrected is the max, not the test", bool(_mx.round(9).equals(_sw.p_corrected.round(9))),
               f"{int((_sw.p_corrected > .05).sum())} of {len(_sw)} rows > .05", "main"))
for _bad, _lbl in (("McNemar", "no McNemar claim in main"), ("mcnemar", "no mcnemar claim in main")):
    _ok = _bad not in MAIN and _bad not in SUPP
    checks.append((_lbl, _ok, f"absent: {_bad}", "both"))
    if not _ok:
        fails.append(_lbl)
has("role-switch rule described", "directional single-measure rule")
has("p_corrected disambiguated", "it is not the test statistic", where="supp")

# Normality rejection rates: 4.8-19.4% overall, not the 4.0-13.2% printed before.
has("normality range corrected", "4.8-19.4%", where="supp")
_stale = "4.0\u201313.2" in BOTH
checks.append(("no stale normality range", not _stale, "absent: 4.0-13.2", "both"))
if _stale:
    fails.append("no stale normality range")
_supp_refs = set(re.findall(r"Table (\d+) of the main text", SUPP))
checks.append(("supp names no dropped main table", _supp_refs <= {"1", "2", "3"},
               str(sorted(_supp_refs)), "supp"))
if not _supp_refs <= {"1", "2", "3"}:
    fails.append("supp names no dropped main table")
has("panel b evidence-vs-effect caveat", "Length in (b) is evidence, not effect size")

# Figure S5's caption must point at the role-switch table, not the old main Table 3
ok = "summarised in Table S19" in SUPP and "summarised in Table 3," not in SUPP
checks.append(("FigS5 ref repaired", ok, "summarised in Table S19", "supp"))
if not ok:
    fails.append("FigS5 ref repaired")

# ------------------------------------------------------- hygiene
for bad, lbl in (("Bayes", "no Bayes factors"), ("BF10", "no BF10"),
                 ("Figure 2A", "no old Fig2A label"), ("Figure 2B", "no old Fig2B label"),
                 ("matrices with no projection", "old Fig3b description gone"),
                 ("pharmacological arms is displayed in Figure 3",
                  "stale Table S3 pointer gone"),
                 ("only analysis in which", "coverage overclaim retracted"),
                 ("the one analysis in the paper", "coverage overclaim retracted (2)")):
    ok = bad not in BOTH
    checks.append((lbl, ok, f"absent: {bad}", "both"))
    if not ok:
        fails.append(lbl)

man = rd("supp_figure_manifest.csv")
ok = f"all {len(man)} staged files" in SUPP
checks.append(("staged-file count", ok, f"all {len(man)} staged files", "supp"))
if not ok:
    fails.append("staged-file count")

# figures actually embedded
# The notebook still builds and stages all 12 figure labels; the Word
# supplement prints only the two that earn their place.
for d, n_exp, lbl in ((Document(str(M / "Methods_and_Results.docx")), 4, "main"),
                      (Document(str(M / "Supplementary_Material.docx")), 2, "supp")):
    n_img = len(d.inline_shapes)
    ok = n_img == n_exp
    checks.append((f"{lbl} embeds {n_exp} figures", ok, f"found {n_img}", lbl))
    if not ok:
        fails.append(f"{lbl} figure count")

for label, ok, needle, where in checks:
    print(f"  {'PASS' if ok else 'FAIL':4s}  [{where:4s}] {label:34s} {needle}")
print(f"\n{len(checks)} checks, {len(fails)} failures")
if fails:
    print("FAILED:", fails)
    sys.exit(1)
