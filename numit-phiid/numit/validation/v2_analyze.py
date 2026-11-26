# V2 analysis: paired drug-vs-control t on ws-mean NuMIT quantiles, per null
# configuration - do sign and significance survive the anchoring choice?
# TEMPORARY - tmp-only investigation.
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

OUT = r"D:\projects\TMS_GWT\numit-phiid\numit\validation\v2_results"
ARMS = {
    "basel_lpm": ("psilocybin", "placebo"),
    "berlin_ketamine": ("ketamine", "placebo"),
    "monash_psilo": ("psilocybin", "baseline"),
}
CONFIGS = ["prod", "vanilla", "armpool", "df5"]
MEASURES = ["synergy", "redundancy", "synergy_total", "phi_r"]


def cond_files(ds, cond):
    out = {}
    for f in glob.glob(os.path.join(OUT, f"dataset-{ds}_sub-*.npz")):
        if f"_ses-{cond}_" not in os.path.basename(f):
            continue
        out[re.search(r"sub-([^_]+)", os.path.basename(f)).group(1)] = f
    return out


rows = []
sat_rows = []
for ds, (drug, ctrl) in ARMS.items():
    fd, fc = cond_files(ds, drug), cond_files(ds, ctrl)
    subs = sorted(set(fd) & set(fc))
    dd = {s: dict(np.load(fd[s])) for s in subs}
    cc = {s: dict(np.load(fc[s])) for s in subs}
    for m in MEASURES:
        row = {"arm": f"{ds}/{drug}", "n": len(subs), "measure": m}
        for cfg in CONFIGS:
            dv = np.array([dd[s][f"{m}_{cfg}"] for s in subs])
            cv = np.array([cc[s][f"{m}_{cfg}"] for s in subs])
            t, p = stats.ttest_rel(dv, cv)
            row[f"{cfg}_t"] = round(float(t), 2)
            row[f"{cfg}_p"] = float(f"{p:.3f}")
        rows.append(row)
    srow = {"arm": f"{ds}/{drug}"}
    for cfg in CONFIGS:
        sats = [d[f"sat_{cfg}"] for d in list(dd.values()) + list(cc.values())]
        srow[f"sat_{cfg}"] = round(float(np.mean(sats)), 3)
    sat_rows.append(srow)

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "..", "v2_table.csv"), index=False)
pd.set_option("display.width", 220)
print(df.to_string(index=False))
print("\n=== saturated-quantile fraction (synergy) per config ===")
print(pd.DataFrame(sat_rows).to_string(index=False))
print("\nVERDICT criteria: same sign across configs everywhere; significance")
print("(p<.05) preserved for the arms/measures significant under 'prod'.")
