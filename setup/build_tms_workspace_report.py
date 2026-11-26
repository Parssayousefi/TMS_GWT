"""Build reviewable workspace audit, plots, and a small executable notebook."""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'setup'))
import analyse_tms_workspace_scope as analysis

OUT=ROOT/'output/tms_workspace_scope'/(ROOT/'output/tms_workspace_scope/LATEST.txt').read_text().strip()
MAN=ROOT/'manuscript/aachen_brief_report'
FIG=MAN/'figures'
LABELS={'synergy':'Persistent synergy','redundancy':'Persistent redundancy','synergy_total':'Total source synergy','phi':'Φ','phi_r':'ΦR'}
SCOPE_LABELS={'whole_brain':'Whole brain','within_workspace':'Within workspace','outside_workspace':'Outside workspace','between':'Workspace–outside links'}


def table_md(frame):
    def cell(x):
        if isinstance(x,(int,float,np.number)):return '' if pd.isna(x) else f'{x:.6g}'
        return str(x).replace('|','/').replace('\n',' ')
    return '\n'.join(['| '+' | '.join(frame.columns)+' |','| '+' | '.join(['---']*len(frame.columns))+' |']+
                     ['| '+' | '.join(cell(x) for x in row)+' |' for row in frame.itertuples(index=False,name=None)])


def diagnostics():
    motion=[]
    for branch in ['raw','normalized']:
        frame=pd.read_csv(OUT/(branch+'_aachen_pooled_changes.csv'))
        # Expanded inventory includes all four scopes and all 20 original role/PC endpoints.
        result=analysis.test_aachen(frame,analysis.METRICS)
        result.to_csv(OUT/(branch+'_aachen_expanded40_sensitivity.csv'),index=False)
        for scope in analysis.SCOPES:
            for m in analysis.MEASURES:
                d=pd.DataFrame({'y':frame[scope+'__'+m+'_delta'],'intercept':1.,'fd':frame.delta_fd})
                motion.append({'branch':branch,'scope':scope,'metric':m,**analysis.hc3(d,'y',['intercept','fd'],'intercept')})
    motion=pd.DataFrame(motion)
    motion['p_holm_scope5']=motion.groupby(['branch','scope']).p.transform(analysis.holm)
    motion['p_holm_all20']=motion.groupby('branch').p.transform(analysis.holm)
    motion.to_csv(OUT/'aachen_all_spatial_motion.csv',index=False)


def figures():
    FIG.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(8.8,6.7),layout='constrained')
    for row,branch in enumerate(['raw','normalized']):
        frame=pd.read_csv(OUT/(branch+'_aachen_pooled_changes.csv'))
        factor=1 if branch=='raw' else 100
        for col,scope in enumerate(analysis.SCOPES[:2]):
            ax=axes[row,col];a=frame[scope+'__phi_baseline'].to_numpy()*factor;b=frame[scope+'__phi_followup'].to_numpy()*factor
            for x,y in zip(a,b):ax.plot([0,1],[x,y],color='#64748b',alpha=.58,lw=.9,marker='o',ms=3.5)
            ax.plot([0,1],[a.mean(),b.mean()],color='#9b382b',lw=2.4,marker='D',ms=6,label='Group mean')
            ax.set_xticks([0,1],['Pre average','Post average']);ax.set_xlim(-.2,1.2)
            ax.set_title(f"{'ABCD'[row*2+col]}  {SCOPE_LABELS[scope]}",loc='left',fontweight='bold')
            ax.set_ylabel('Mean Φ (nats)' if branch=='raw' else 'Mean normalized Φ\n(conditional percentile)')
            ax.grid(axis='y',alpha=.15)
            if row==0 and col==0:ax.legend(frameon=False,fontsize=9)
    for ext in ['png','pdf']:fig.savefig(FIG/('Figure_1_participant_changes.'+ext),dpi=300,bbox_inches='tight')
    plt.close(fig)
    results=pd.read_csv(OUT/'aachen_all_spatial20_results.csv')
    fig,axes=plt.subplots(2,2,figsize=(9.6,6.3),layout='constrained')
    colors=['#334155','#7336a2','#097d75','#c97c21']
    for row,branch in enumerate(['raw','normalized']):
        factor=1 if branch=='raw' else 100
        for col,m in enumerate(['phi','phi_r']):
            ax=axes[row,col]
            for i,scope in enumerate(analysis.SCOPES):
                r=results.loc[results.branch.eq(branch)&results.metric.eq(scope+'__'+m)].iloc[0]
                e=r.estimate*factor;lo=r.ci_low*factor;hi=r.ci_high*factor
                ax.errorbar(e,3-i,xerr=[[e-lo],[hi-e]],fmt='o',color=colors[i],capsize=3,ms=6,lw=1.8)
            ax.axvline(0,color='#555555',ls=':',lw=1)
            ax.set_yticks([3,2,1,0],[SCOPE_LABELS[s] for s in analysis.SCOPES] if col==0 else ['']*4)
            ax.set_ylim(-.5,3.6);ax.set_title(f"{'ABCD'[row*2+col]}  {LABELS[m]}"+(" normalized" if branch=='normalized' else ''),loc='left',fontweight='bold')
            ax.set_xlabel('Mean change (nats)' if branch=='raw' else 'Mean change (percentile points)')
            ax.ticklabel_format(axis='x',style='plain',useOffset=False);ax.grid(axis='x',alpha=.15)
    for ext in ['png','pdf']:fig.savefig(FIG/('Figure_2_spatial_scope.'+ext),dpi=300,bbox_inches='tight')
    plt.close(fig)
    leave_one_out_figure()


def leave_one_out_figure():
    from matplotlib.ticker import MaxNLocator
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42})
    # Leave-one-out intervals remain marginal, including when the full-sample CI excludes zero.
    loo=pd.read_csv(OUT/'aachen_leave_one_out.csv')
    fig,axes=plt.subplots(1,2,figsize=(8.8,5.2),layout='constrained')
    for ax,scope in zip(axes,analysis.SCOPES[:2]):
        d=loo.loc[loo.branch.eq('raw')&loo.metric.eq(scope+'__phi')].sort_values('omitted_subject')
        y=np.arange(len(d))
        ax.errorbar(d.estimate,y,xerr=np.array([d.estimate-d.ci_low,d.ci_high-d.estimate]),fmt='o',ms=4,capsize=2,color='#7336a2')
        ax.axvline(0,color='#555555',ls=':',lw=1)
        ax.set_yticks(y,[s.rsplit('_',1)[1] for s in d.omitted_subject]);ax.invert_yaxis()
        ax.set_title(SCOPE_LABELS[scope]);ax.set_xlabel('Mean Φ change (nats)');ax.set_ylabel('Participant omitted')
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    for ext in ['png','pdf']:fig.savefig(FIG/('Figure_S1_leave_one_out.'+ext),dpi=300,bbox_inches='tight')
    plt.close(fig)


def report():
    main=pd.read_csv(OUT/'aachen_main_results.csv')
    short=main.loc[main.metric.str.contains('__'),['branch','metric','estimate','ci_low','ci_high','p','p_maxT_family','p_maxT_whole_and_workspace10','p_maxT_all']].copy()
    for branch,filename in [('raw','table_1_scope_results.csv'),('normalized','table_2_scope_results.csv')]:
        frame=short[short.branch.eq(branch)].copy()
        if branch=='normalized':frame[['estimate','ci_low','ci_high']]*=100
        frame['scope']=frame.metric.str.split('__').str[0].map(SCOPE_LABELS)
        frame['measure']=frame.metric.str.split('__').str[1].map(LABELS)
        frame.to_csv(MAN/filename,index=False)
    text=f'''# Fixed HCP workspace audit and TMS reanalysis

The fixed HCP workspace was already used correctly for participation coefficients, gateways and broadcasters. The five global information averages used all 454 regions. This rerun adds means restricted to pairs whose **both endpoints** are in the unchanged 235-region HCP workspace (27,495 unique pairs), and compares them with the 102,831-pair whole-brain means.

The GWT publication notebook and Yousefi paper report participation and role outcomes, not a matching mean-Φ outcome. Luppi tested pairwise ΦR using network-based statistics and subsequently examined workspace overlap. This scalar workspace analysis is therefore an extension of the shared estimator, not an exact repetition of Luppi's inferential procedure.

## What was rerun

All 193 completed TMS scans were re-audited from their existing pairwise matrices. All 25 historical scalar outcomes were reproduced. Conditional normalization was reevaluated for the 189 scans contributing to either the acute comparisons or all-follow-up analysis, retaining the original shared reference distributions and both seeds. Preprocessing and ΦID estimation were unchanged. No HCP or psychedelic participant analysis was rerun.

Aachen retains 15 independent participants: each contributes the mean of their two matched post-minus-pre changes. SUDMEX retains the randomized two-week comparison (40 participants) separately from its all-post longitudinal sensitivity (41 participants). The latter averages each person's available follow-ups, including open-label treatment, and does not estimate a randomized active-versus-sham effect.

## Families and correction

The five information measures are **one family**, not five families. Exact joint sign flips and the maximum absolute t statistic preserve their dependence. Scope-specific five-test results, combined whole-brain/workspace ten-test results, and a 30-endpoint sensitivity including the original 20 role/PC endpoints are all retained. Four-scope diagnostics use a joint 20-test family, with all 40 neural endpoints also provided as a sensitivity. These nested analyses remain exploratory; they do not make the complete historical search confirmatory. Raw and normalized analyses are reported together as dependent method sensitivities.

## Aachen results

Raw changes are nats. Normalized entries below are conditional-CDF units; multiply estimates and intervals by 100 for percentile points. Confidence intervals are marginal.

{table_md(short)}

No within-workspace information endpoint survives its five-measure family correction. Workspace Φ rises by 0.001299 nats (95% CI −0.000293 to 0.002892; maxT p=0.2841). Workspace ΦR rises by 0.001981 nats (0.000173 to 0.003790; uncorrected p=0.0340, maxT p=0.1088). Normalization does not strengthen these findings: the respective five-measure maxT p values are 0.5195 and 0.1882. Motion-adjusted workspace Φ and ΦR intervals include zero; the nominal p values are 0.1828 and 0.1111.

The original whole-brain Φ finding is reproduced. Its raw maxT p is 0.0239 across five information measures and 0.0429 across the ten whole-brain/workspace outcomes. The normalized equivalents are 0.0483 and 0.0782. Neither survives the expanded 30-endpoint correction. Participation and role outcomes are unchanged by this audit; the mask had already been applied to them.

## Spatial diagnosis

The 219 non-workspace regions define 23,871 outside–outside pairs; 51,465 pairs connect workspace to outside regions. The raw Φ change outside the workspace is 0.002809 nats (95% CI 0.001198 to 0.004420), and normalized Φ rises by 1.151 percentile points (0.513 to 1.790). These survive correction across all 20 information/scope combinations (raw p=0.00946; normalized p=0.01270) and the 40-endpoint sensitivity (0.03851 and 0.03479). Within the five outside-workspace information measures, motion-adjusted Holm p values are 0.03983 and 0.02262; across all 20 scope/information measures, they are 0.15930 and 0.09048.

However, the **direct difference in change outside versus within the workspace** does not survive correction across the five information measures (raw Φ maxT p=0.1652; normalized Φ p=0.1681). A significant effect outside and a nonsignificant effect inside are not evidence of a significant difference between locations. The defensible statement is that the clearest corrected evidence concerns pairs outside the fixed workspace; a selective or preferential outside-workspace effect has not been established.

## SUDMEX and clinical outcomes

For the randomized comparison, all five within-workspace information measures have Holm-adjusted p=1.00 in both the baseline-only and baseline-plus-motion models, for raw and normalized data. Raw baseline-adjusted workspace Φ is −0.000125 nats (95% CI −0.004117 to 0.003866), and ΦR is +0.000007 (−0.003930 to 0.003945). There is no evidence here for a randomized workspace-integration increase.

The all-post sensitivity shows raw increases in workspace persistent synergy and total source synergy (five-measure maxT p=0.00116 and 0.00110), but normalized total source synergy decreases (p=0.00813); normalized persistent synergy is not supported. This is a metric- and normalization-dependent longitudinal pattern, with mixed treatment stages and unequal follow-up availability. It should not replace the randomized comparison or be called a replicated Aachen effect. No clinical–neural association survives its declared family correction. Aachen provides no matched craving or impulsivity endpoint.

## Implication for the manuscript

Retain a clearly exploratory paper about information dynamics after iTBS. Do not describe the Aachen global Φ increase as increased workspace synergy, a replicated Luppi/Yousefi result, or improved consciousness. The new fixed-mask analysis belongs prominently in the main results. Report the outside-workspace decomposition as a diagnostic extension, including the nonsignificant direct spatial contrast. The study remains small, lacks sham in Aachen, and is not evidence of clinical efficacy.

## Verification and reproducibility

The independent verifier reconstructed workspace means by indexing both matrix axes, reproduced participant averaging and all 32,768 workspace sign-flip statistics with SciPy, and checked graph-partition identities and input hashes. The direct-indexing discrepancy was below 1.1×10⁻¹⁶ nats. All frozen inputs remained unchanged.

- Run directory: `{OUT.relative_to(ROOT).as_posix()}`.
- Analysis: `setup/analyse_tms_workspace_scope.py`.
- Verification: `setup/verify_tms_workspace_scope.py`.
- This report, expanded diagnostics and figures: `setup/build_tms_workspace_report.py`.
- Executable review notebook: `output/jupyter-notebook/TMS_Fixed_Workspace_Reanalysis.ipynb`.
- Source methods: [Yousefi et al.](https://www.biorxiv.org/content/10.64898/2026.09.28.755172v1.full), [Luppi et al.](https://andrealuppi.github.io/publication/luppi2023elife/luppi2023elife.pdf).
'''
    (ROOT/'docs/TMS_FIXED_WORKSPACE_REANALYSIS_20261125.md').write_text(text,encoding='utf-8')
    (OUT/'README.md').write_text(text,encoding='utf-8')
    return short


def notebook():
    path=ROOT/'output/jupyter-notebook/TMS_Fixed_Workspace_Reanalysis.ipynb'
    nb=nbformat.read(path,as_version=4)
    md=nbformat.v4.new_markdown_cell;code=nbformat.v4.new_code_cell
    nb.cells=[md('# Fixed HCP workspace audit and TMS reanalysis\n\nRead the completed TMS rerun and its independent numerical verification. The frozen 235-region mask is unchanged. The original GWT role/PC definitions already used it; the new analysis adds information averages with both endpoints inside it. No preprocessing or frozen cohort analysis is launched.'),
      code("from pathlib import Path\nimport json\nimport numpy as np\nimport pandas as pd\nfrom IPython.display import display, Image\nROOT = Path(r'D:/projects/TMS_GWT')\nRUN = ROOT / 'output/tms_workspace_scope' / '"+OUT.name+"'\nverification = json.loads((RUN / 'independent_verification.json').read_text())\nassert verification['status'] == 'passed'\nprint('Verified:', verification['raw_scans'], 'scans; Aachen n=15; randomized SUDMEX n=40')\npd.set_option('display.max_colwidth', 45)"),
      md('## Fixed mask and estimand\n\nWithin workspace means both nodes belong to the fixed mask: 235 × 234 / 2 = 27,495 pairs. PCs still use all 454 destinations before workspace averaging and ranking. The five information measures form one family. The ten-test analysis jointly covers whole-brain and workspace summaries.'),
      code("mask = pd.read_csv(RUN / 'fixed_hcp_workspace.csv')\nassert mask.workspace.sum() == 235\ndisplay(mask.groupby('network').workspace.agg(['sum', 'count']).rename(columns={'sum':'Workspace regions','count':'Atlas regions'}))"),
      md('## Aachen averaged changes\n\nRaw changes are nats; normalized changes are shown in percentile points. Intervals are marginal, and raw/normalized results are dependent sensitivities.'),
      code("main = pd.read_csv(RUN / 'aachen_main_results.csv')\ncols = ['branch','metric','estimate','ci_low','ci_high','p','p_maxT_family','p_maxT_whole_and_workspace10']\nshow = main.loc[main.metric.str.contains('__'), cols].copy()\nshow.loc[show.branch.eq('normalized'), ['estimate','ci_low','ci_high']] *= 100\ndisplay(show.round(6))"),
      code("display(Image(filename=str(ROOT / 'manuscript/aachen_brief_report/figures/Figure_1_participant_changes.png'), width=850))"),
      md('No within-workspace information measure survives correction. Raw whole-brain Φ survives the joint ten-test correction; normalized whole-brain Φ does not. The unchanged role and PC results remain available below.'),
      code("display(main.loc[~main.metric.str.contains('__'), ['branch','metric','estimate','p','p_maxT_family']].round(5))"),
      md('## Motion and spatial diagnostics\n\nOutside-workspace Φ has stronger evidence, but the direct outside-minus-inside comparison does not survive correction. Significance in one location and nonsignificance in another is not a test of their difference.'),
      code("display(pd.read_csv(RUN / 'aachen_robustness.csv').query(\"metric.str.startswith('within_workspace') and check == 'motion_HC3'\", engine='python').round(6))\ndisplay(pd.read_csv(RUN / 'aachen_spatial_interactions.csv').round(6))"),
      code("display(Image(filename=str(ROOT / 'manuscript/aachen_brief_report/figures/Figure_2_spatial_scope.png'), width=900))"),
      md('## SUDMEX randomized comparison\n\nBaseline-adjusted post-treatment active-minus-sham effects, with motion sensitivity. This preserves the randomized two-week question.'),
      code("sud = pd.read_csv(RUN / 'sudmex_randomized_results.csv')\ndisplay(sud.loc[sud.metric.str.startswith('within_workspace') & sud.model.ne('unadjusted_change_difference'), ['branch','model','metric','estimate','ci_low','ci_high','p','p_holm_family']].round(6))"),
      md('## All-post longitudinal sensitivity and clinical associations\n\nThis mixes randomized and later open-label periods and unequal numbers of follow-ups. It is kept separate from intervention inference. Clinical associations use complete available cases and mask flagged source values.'),
      code("display(pd.read_csv(RUN / 'sudmex_all_post_workspace_results.csv')[['branch','metric','estimate','ci_low','ci_high','p_maxT_family']].round(6))\nclinical = pd.read_csv(RUN / 'clinical_neural_associations.csv')\nprint('Clinical associations surviving correction over 56 tests per branch:', int((clinical.p_holm56 < .05).sum()))"),
      md('## Reproduction and interpretation\n\n`setup/analyse_tms_workspace_scope.py` derives the results from existing TMS pairwise matrices; `setup/verify_tms_workspace_scope.py` independently checks means, intervals, permutation p values and immutable inputs. `setup/build_tms_workspace_report.py` produces the review artifacts.\n\nThe results support an exploratory description of information change after iTBS. They do not demonstrate increased workspace synergy, a selective spatial effect, improved consciousness, or clinical efficacy. See the complete audit in `docs/TMS_FIXED_WORKSPACE_REANALYSIS_20261125.md`.')]
    nb.metadata.kernelspec={'display_name':'Python (gwt)','language':'python','name':'gwt'}
    nbformat.validate(nb);nbformat.write(nb,path)


def main():
    with threadpool_limits(limits=1):diagnostics()
    figures();report();notebook()
    receipt={'created_utc':datetime.now(timezone.utc).isoformat(),'analysis_run':OUT.name,
             'source_script_sha256':analysis.tms.sha256(Path(__file__)),
             'supplementary_diagnostics':['raw_aachen_expanded40_sensitivity.csv','normalized_aachen_expanded40_sensitivity.csv','aachen_all_spatial_motion.csv'],
             'input_sha256':{p.name:analysis.tms.sha256(p) for p in OUT.glob('*.csv')},
             'note':'Added spatial diagnostics retained in complete inventory; no changes to primary estimates'}
    analysis.save_json(OUT/'report_provenance.json',receipt)
    print('Built audit, manuscript plots and notebook from',OUT)


if __name__=='__main__':main()
