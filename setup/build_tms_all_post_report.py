"""Create the all-follow-up figure and add the completed analysis to the notebook."""
from pathlib import Path
import json
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/tms_all_post_average/20261125T162521_698094Z'


def main():
    receipt=json.loads((OUT/'independent_verification.json').read_text())
    assert receipt['status']=='passed'
    table=pd.read_csv(OUT/'neural_pre_post_results.csv')
    labels={'global_synergy':'Persistent synergy','global_redundancy':'Persistent redundancy',
            'global_synergy_total':'Total synergy','global_phi':'Whole-minus-sum Φ','global_phi_r':'Revised ΦR'}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'savefig.dpi':180,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
    for ax,branch,scale,title,unit in zip(axes,['raw','normalized'],[1000,100],
                                         ['Raw information','Conditional-null normalization'],
                                         ['Mean change (10⁻³ nats)','Mean change (percentile points)']):
        block=table.query("dataset=='sudmex_tms' and cohort=='all_eligible'")
        block=block.loc[block.branch.eq(branch)].set_index('metric').loc[list(labels)]
        y=np.arange(5)[::-1]
        for j,(m,row) in enumerate(block.iterrows()):
            color='#087f8c' if row.p_maxT_family<.05 else '#657280'
            ax.errorbar(row.estimate*scale,y[j],xerr=np.array([[row.estimate-row.ci_low],[row.ci_high-row.estimate]])*scale,
                        fmt='o',color=color,capsize=4,ms=6,lw=1.8)
            ax.text(1.02,y[j],f'{row.p_maxT_family:.3f}',transform=ax.get_yaxis_transform(),va='center',color=color)
        ax.axvline(0,color='#9ba4af',ls='--',lw=1)
        ax.set_yticks(y,[labels[m] for m in block.index]);ax.set_ylim(-.6,4.8)
        ax.set_xlabel(unit);ax.set_title(title,loc='left',fontweight='bold')
        ax.text(1.02,4.65,'FWER p',transform=ax.get_yaxis_transform(),fontsize=9)
        ax.grid(axis='x',alpha=.15)
    fig.suptitle('SUDMEX: baseline versus mean of all available follow-ups (41 participants)',fontsize=13)
    fig.get_layout_engine().set(rect=(0,.08,1,.88))
    fig.text(.02,.014,'Bars: marginal 95% intervals. Family max-T uses 200,000 joint sign configurations.\nBoth original arms are included; this is not a sham-controlled treatment contrast.',fontsize=9)
    for ext in ['png','pdf']:fig.savefig(OUT/('sudmex_all_post_global_family.'+ext),bbox_inches='tight')
    plt.close(fig)
    path=ROOT/'output/jupyter-notebook/TMS_Exploratory_Findings.ipynb'
    backup=OUT/'TMS_Exploratory_Findings_before_all_post.ipynb'
    if not backup.exists():shutil.copy2(path,backup)
    nb=nbformat.read(path,as_version=4)
    tag='all-post-average-20261125'
    nb.cells=[c for c in nb.cells if tag not in c.metadata.get('tags',[])]
    md=nbformat.v4.new_markdown_cell;code=nbformat.v4.new_code_cell
    cells=[md('## Baseline versus the average of every available post measurement\n\nSUDMEX includes 41 baseline scans and 88 follow-ups from 41 people, each contributing one equal-weight participant average. Visits range from t1 to t4 and include the original sham arm\'s open-label phase. This is a longitudinal temporal comparison, not a sustained randomized active–sham contrast. Aachen is unchanged because mean(post90, post120) − mean(pre90, pre120) equals the original average of matched changes.'),
           code("ALL_POST = ROOT/'output/tms_all_post_average/20261125T162521_698094Z'\nassert json.loads((ALL_POST/'independent_verification.json').read_text())['status'] == 'passed'\ncontributions = pd.read_csv(ALL_POST/'scan_contributions.csv')\ndisplay(contributions.groupby(['dataset','phase']).size().rename('scans').reset_index())\nall_post = pd.read_csv(ALL_POST/'neural_pre_post_results.csv')\ndisplay(all_post.query(\"cohort == 'all_eligible' and metric in ['global_synergy','global_synergy_total','global_phi','global_phi_r','mean_synergy_pc','gateways_FPN']\")[[ 'branch','dataset','metric','n','estimate','ci_low','ci_high','p_maxT_family','p_holm_family']])\ndisplay(Image(filename=str(ALL_POST/'sudmex_all_post_global_family.png')))"),
           md('Raw SUDMEX synergy and total synergy increase, and mean synergy participation decreases, after family-wise correction. These do not survive normalization in the full sample. Phi and Phi_R do not show supported changes. Restricting to the original 40-person sample preserves the raw findings and adds a borderline normalized total-synergy decrease; it fails in the full sample and under Holm.\n\nBaseline adjustment again removes the FPN original-arm contrast. All seven clinical scales improve in the overall pre/post comparison, but no original-arm difference passes clinical-family correction.'),
           code("arms = pd.read_csv(ALL_POST/'original_arm_contrasts.csv')\ndisplay(arms.query(\"cohort == 'all_eligible' and metric == 'gateways_FPN'\")[[ 'branch','model','estimate','ci_low','ci_high','p_holm_family']])\ndisplay(pd.read_csv(ALL_POST/'clinical_all_post_results.csv')[[ 'contrast','metric','n','estimate','p_holm_clinical7']])\ndisplay(all_post.query(\"dataset == 'sudmex_tms' and branch == 'normalized' and metric == 'global_synergy_total'\")[[ 'cohort','n','estimate','p_maxT_family','p_holm_family']])"),
           md('**Interpretation:** this is a compact summary of observed follow-ups, with unequal horizons and treatment histories. Keep the acute randomized comparison for intervention inference and use the all-post mean as a longitudinal sensitivity. The full findings and scan accounting are in `docs/TMS_ALL_POST_AVERAGE_20261125.md`.'),
           code("display(pd.read_csv(ALL_POST/'temporal_information_descriptive_summary.csv'))\nprint('All-post report:', ROOT/'docs/TMS_ALL_POST_AVERAGE_20261125.md')\nprint('Verified:', json.loads((ALL_POST/'independent_verification.json').read_text())['status'])")]
    for cell in cells:cell.metadata['tags']=[tag]
    nb.cells.extend(cells)
    nbformat.write(nb,path)
    print('Added all-post section:',path)
    print('Figure:',OUT/'sudmex_all_post_global_family.png')


if __name__=='__main__':main()
