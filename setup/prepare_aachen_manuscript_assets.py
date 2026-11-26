"""Build manuscript assets from completed Aachen outputs without rerunning the pipeline."""
from pathlib import Path
import hashlib
import json
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'manuscript' / 'aachen_brief_report'
FIG = OUT / 'figures'
SUP = OUT / 'supplementary_data'
for folder in [OUT, FIG, SUP, OUT/'sources', OUT/'qa']:
    folder.mkdir(parents=True, exist_ok=True)
RAW = ROOT/'output/tms_publication/20261125T145928_921789Z'
NORM = ROOT/'output/tms_numit/20261125T155252_384543Z'
FW = ROOT/'output/tms_familywise/20261125T161051_102285Z'
ROB = ROOT/'output/tms_exploration/20261125T155441_285503Z'
inputs = {}


def read(path):
    inputs[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path)


def ci(x):
    x = np.asarray(x, dtype=float)
    mean = x.mean()
    margin = stats.t.ppf(.975, len(x)-1)*stats.sem(x)
    return [float(mean), float(mean-margin), float(mean+margin)]


main = read(FW/'main_comparisons_familywise.csv')
main = main.loc[main.dataset.eq('aachen_itbs')].copy()
main.to_csv(SUP/'Table_S1_all_25_neural_endpoints.csv', index=False)
global_results = main.loc[main.family.eq('global_phiid')].copy()
global_results.to_csv(OUT/'table_1_global_results.csv', index=False)
motion = read(FW/'aachen_motion_global_familywise.csv')
motion.to_csv(SUP/'Table_S2_motion_adjustment.csv', index=False)
extended = read(FW/'extended_exploration_familywise.csv')
extended.loc[extended.dataset.eq('aachen_itbs')].to_csv(SUP/'Table_S3_secondary_exploration.csv', index=False)
seeds = read(FW/'numit_global_family_seed_stability.csv')
seeds.to_csv(SUP/'Table_S4_null_reference_stability.csv', index=False)
robust = read(ROB/'aachen_selected_global_robustness.csv')
robust.to_csv(SUP/'Table_S5_additional_robustness.csv', index=False)
loo = read(ROB/'aachen_selected_global_leave_one_out.csv')
loo.to_csv(SUP/'Table_S6_leave_one_out.csv', index=False)
demographics = []
for line in (ROOT/'data/aachen_itbs/derivatives/README.txt').read_text().splitlines():
    match = re.match(r'^(\d+)\s+([\d/]+)\s+([fm])\s+(\d+)', line)
    if match and match[1] != '2':
        demographics.append({'subject_id':int(match[1]), 'order':match[2], 'sex':match[3], 'age':int(match[4])})
demo = pd.DataFrame(demographics)
assert len(demo) == 15

branches = {}
for branch, location in [('raw',RAW),('normalized',NORM)]:
    scans = read(location/'participant_scan_metrics.csv')
    scans = scans.loc[scans.dataset.eq('aachen_itbs')].copy()
    pairs = read(location/'participant_paired_changes.csv')
    pairs = pairs.loc[pairs.dataset.eq('aachen_itbs')].copy()
    scans.to_csv(SUP/f'{branch}_participant_scan_metrics.csv', index=False)
    pairs.to_csv(SUP/f'{branch}_participant_paired_changes.csv', index=False)
    branches[branch] = (scans,pairs)

raw_scans, raw_pairs = branches['raw']
pivot = raw_scans.pivot(index='subject_id',columns='session',values='global_phi')
norm_pivot = branches['normalized'][0].pivot(index='subject_id',columns='session',values='global_phi')
pre = pivot[['pre90','pre120']].mean(axis=1)
post = pivot[['post90','post120']].mean(axis=1)
d90 = pivot.post90-pivot.pre90
d120 = pivot.post120-pivot.pre120
delta = post-pre
np.testing.assert_allclose(delta, (d90+d120)/2, atol=1e-14)
row = global_results.loc[global_results.branch.eq('raw') & global_results.metric.eq('global_phi')].iloc[0]
np.testing.assert_allclose(ci(delta), [row.estimate,row.ci_low,row.ci_high], atol=1e-12)
np.testing.assert_allclose(stats.ttest_1samp(delta,0).pvalue, row.p, atol=1e-12)
assert len(pivot) == 15 and len(raw_scans) == 60
assert raw_scans.n_valid_pairs.eq(102831).all()
fd = raw_scans.pivot(index='subject_id',columns='session',values='mean_fd')
fd_pre = fd[['pre90','pre120']].mean(axis=1)
fd_post = fd[['post90','post120']].mean(axis=1)
dz = float(delta.mean()/delta.std(ddof=1))
summary = {
    'n':15,'scans':60,'women':int(demo.sex.eq('f').sum()),'men':int(demo.sex.eq('m').sum()),
    'age_mean':float(demo.age.mean()),'age_sd':float(demo.age.std(ddof=1)),
    'age_range':[int(demo.age.min()),int(demo.age.max())],
    'order_90_first':int(demo.order.eq('1/2/3/4').sum()),
    'order_120_first':int(demo.order.eq('3/4/1/2').sum()),
    'raw_pre_mean':float(pre.mean()),'raw_pre_sd':float(pre.std(ddof=1)),
    'raw_post_mean':float(post.mean()),'raw_post_sd':float(post.std(ddof=1)),
    'raw_delta_ci':ci(delta),'positive_changes':int(delta.gt(0).sum()),
    'raw_phi_t':float(stats.ttest_1samp(delta,0).statistic),
    'raw_phi_dz':dz,'raw_phi_gz_approx':dz*(1-3/(4*14-1)),
    'percent_change_of_group_mean':float(100*delta.mean()/pre.mean()),
    'intensity_90_ci':ci(d90),'intensity_120_ci':ci(d120),
    'intensity_difference_ci':ci(d120-d90),
    'intensity_change_correlation':float(d90.corr(d120)),
    'fd_pre_mean_mm':float(fd_pre.mean()),'fd_post_mean_mm':float(fd_post.mean()),
    'fd_delta_ci_mm':ci(fd_post-fd_pre),
    'normalized_delta_ci_percentile_points':ci(100*(norm_pivot[['post90','post120']].mean(axis=1)-norm_pivot[['pre90','pre120']].mean(axis=1))),
    'timepoint_counts':raw_scans.n_timepoints.value_counts().sort_index().to_dict(),
    'main_rows':global_results.to_dict(orient='records'),
    'motion_rows':motion.to_dict(orient='records'),
}

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,
                    'axes.spines.right':False,'axes.labelsize':9,'axes.titlesize':10,
                    'pdf.fonttype':42,'ps.fonttype':42})
blue, orange = '#236B8E','#BD5B2C'
fig, axes = plt.subplots(1,3,figsize=(10.8,3.8),layout='constrained')
for ax, before, after, ylabel, title in [
    (axes[0],pre,post,'Mean pairwise Φ (nats)','A  Raw Φ'),
    (axes[1],100*norm_pivot[['pre90','pre120']].mean(axis=1),
     100*norm_pivot[['post90','post120']].mean(axis=1),'Mean conditional percentile','B  Normalized Φ')]:
    for x,y in zip(before,after):
        ax.plot([0,1],[x,y],color='#B9BFC3',lw=.9,zorder=1)
    ax.scatter(np.zeros(len(before)),before,color=blue,s=23,zorder=2)
    ax.scatter(np.ones(len(after)),after,color=orange,s=23,zorder=2)
    ax.set_xticks([0,1],['Pre average','Post average'])
    ax.set_xlim(-.28,1.28)
    ax.set_ylabel(ylabel)
    ax.set_title(title,loc='left',fontweight='bold',pad=10)
ax=axes[2]
for x,y in zip(d90,d120):
    ax.plot([0,1],[x,y],color='#C4C8CB',lw=.8,zorder=1)
ax.scatter(np.zeros(15),d90,color=blue,s=23,zorder=2)
ax.scatter(np.ones(15),d120,color=orange,s=23,zorder=2)
for pos,values in [(0,d90),(1,d120)]:
    m,l,h=ci(values)
    ax.errorbar(pos+.17,m,yerr=[[m-l],[h-m]],fmt='D',color='black',capsize=3,ms=4,lw=1.3,zorder=3)
ax.axhline(0,color='#555555',ls='--',lw=.8)
ax.set_xticks([0,1],['90% rMT','120% rMT'])
ax.set_xlim(-.28,1.35)
ax.set_ylabel('Post minus pre Φ (nats)')
ax.set_title('C  Visit-specific changes',loc='left',fontweight='bold',pad=10)
fig.savefig(FIG/'Figure_1_participant_changes.png',dpi=300)
fig.savefig(FIG/'Figure_1_participant_changes.pdf')
plt.close(fig)

# Supplementary figure shows existing leave-one-out estimates only.
phi_loo=loo.loc[loo.metric.eq('global_phi')].copy()
estimate_col='estimate' if 'estimate' in phi_loo else 'mean_change'
assert len(phi_loo)==15
fig,ax=plt.subplots(figsize=(6.8,4.5),layout='constrained')
ax.scatter(phi_loo[estimate_col],np.arange(15),color=blue,s=24)
if {'ci_low','ci_high'}.issubset(phi_loo.columns):
    ax.hlines(np.arange(15),phi_loo.ci_low,phi_loo.ci_high,color=blue,lw=1)
ax.axvline(0,color='#555555',ls='--',lw=.8)
ax.axvline(delta.mean(),color=orange,lw=1,label='Full sample estimate')
subject_col=next(c for c in ['omitted_subject','excluded_subject','left_out','subject_id','excluded'] if c in phi_loo)
labels=phi_loo[subject_col].astype(str).tolist()
ax.set_yticks(np.arange(15),labels)
ax.set_ylabel('Omitted participant')
ax.set_xlabel('Mean post minus pre Φ (nats)')
ax.legend(frameon=False,loc='lower right')
fig.savefig(FIG/'Figure_S1_leave_one_out.png',dpi=300)
fig.savefig(FIG/'Figure_S1_leave_one_out.pdf')
plt.close(fig)
summary['leave_one_out_columns']=loo.columns.tolist()
summary['leave_one_out_range']=phi_loo[estimate_col].agg(['min','max']).to_dict()

prepared=json.loads((NORM/'cache/prepared.json').read_text())
summary['normalization_pool_size']=prepared['baseline_pool_sizes']['aachen_itbs']
summary['normalization_tdmi_range']=prepared['tdmi_ranges']['aachen_itbs']
cal=read(NORM/'synthetic_calibration_summary.csv')
cal.loc[cal.dataset.eq('aachen_itbs')].to_csv(SUP/'Table_S7_null_model_calibration.csv',index=False)
(OUT/'verified_results.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(OUT/'source_file_hashes.json').write_text(json.dumps(inputs,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in summary.items() if k not in ['main_rows','motion_rows']},indent=2))
