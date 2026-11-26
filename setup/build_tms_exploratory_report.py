"""Build reviewable figures, a findings report and a lightweight results notebook."""
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat as nbf
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'output/tms_publication/20261125T145928_921789Z'
FWE=ROOT/'output/tms_familywise/20261125T161051_102285Z'
NUMIT=ROOT/'output/tms_numit/20261125T155252_384543Z'
ROUND1=ROOT/'output/tms_exploration/20261125T154020_735179Z'
FOLLOWUP=ROOT/'output/tms_exploration/20261125T155441_285503Z'
DIAG=ROOT/'output/tms_diagnostics/20261125T152110_809161Z'
OUT=ROOT/'output/tms_exploratory_report/20261125'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    table=pd.read_csv(FWE/'main_comparisons_familywise.csv')
    pairs=pd.read_csv(SOURCE/'participant_paired_changes.csv')
    scans=pd.read_csv(SOURCE/'participant_scan_metrics.csv')
    raw_extended=pd.read_csv(ROUND1/'neural_exploration.csv')
    clinical=pd.read_csv(FWE/'clinical_familywise.csv')
    verification=json.loads((FWE/'independent_verification.json').read_text())
    assert verification['status']=='passed_for_estimable_contrasts'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':180,'pdf.fonttype':42})
    labels={'global_synergy':'Persistent synergy','global_redundancy':'Persistent redundancy',
            'global_synergy_total':'Total synergy','global_phi':'Whole-minus-sum Φ','global_phi_r':'Revised ΦR'}
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    for ax,branch,scale,title,unit in zip(axes,['raw','numit'],[1000,100],
                                         ['Raw information','Conditional-null normalization'],
                                         ['Mean paired change (10⁻³ nats)','Mean paired change (percentile points)']):
        b=table.loc[table.dataset.eq('aachen_itbs')&table.branch.eq(branch)&table.family.eq('global_phiid')].set_index('metric').loc[list(labels)]
        y=np.arange(len(b))[::-1]
        for k,(_,row) in enumerate(b.iterrows()):
            color='#087f8c' if row.name=='global_phi' else '#657280'
            ax.errorbar(row.estimate*scale,y[k],xerr=np.array([[row.estimate-row.ci_low],[row.ci_high-row.estimate]])*scale,
                        fmt='o',color=color,capsize=4,ms=6,lw=1.8)
            ax.text(1.02,y[k],f'{row.p_maxT_family:.3f}',transform=ax.get_yaxis_transform(),va='center',color=color)
        ax.axvline(0,color='#9ba4af',ls='--',lw=1)
        ax.set_yticks(y,[labels[m] for m in b.index]); ax.set_ylim(-.6,4.8)
        ax.set_xlabel(unit); ax.set_title(title,loc='left',fontweight='bold')
        ax.text(1.02,4.65,'FWER p',transform=ax.get_yaxis_transform(),fontsize=9)
        ax.grid(axis='x',alpha=.15)
    fig.suptitle('Aachen: average change across two iTBS visits (15 participants)',fontsize=13)
    fig.get_layout_engine().set(rect=(0,.06,1,.90))
    fig.text(.02,.012,'Bars: marginal 95% confidence intervals. FWER p: exact max-|t| across all five global measures; 32,768 joint sign flips.',fontsize=9)
    for ext in ['png','pdf']: fig.savefig(OUT/('aachen_global_family.'+ext),bbox_inches='tight')
    plt.close(fig)

    a90=pairs.query("dataset=='aachen_itbs' and comparison=='post90_minus_pre90'").set_index('subject_id').sort_index()
    a120=pairs.query("dataset=='aachen_itbs' and comparison=='post120_minus_pre120'").set_index('subject_id').sort_index()
    values=np.column_stack([a90.global_phi_delta,a120.global_phi_delta,(a90.global_phi_delta+a120.global_phi_delta)/2])*1000
    fig,ax=plt.subplots(figsize=(8.4,4.7),layout='constrained')
    jitter=np.random.default_rng(2026112516).uniform(-.075,.075,len(values))
    for row,j in zip(values,jitter):
        ax.plot(np.array([0,1])+j,row[:2],color='#aab3bc',alpha=.5,lw=.8,zorder=1)
        ax.scatter(np.arange(3)+j,row,color=['#657280','#657280','#087f8c'],s=23,alpha=.7,zorder=2)
    for i,contrast in enumerate(['90_change','120_change','pooled_change']):
        row=raw_extended.loc[raw_extended.dataset.eq('aachen_itbs')&raw_extended.metric.eq('global_phi')&raw_extended.contrast.eq(contrast)].iloc[0]
        ax.errorbar(i+.20,row.estimate*1000,yerr=np.array([[row.estimate-row.ci_low],[row.ci_high-row.estimate]])*1000,
                    fmt='D',color='#202e3b',capsize=4,ms=7,lw=2,zorder=4)
    ax.axhline(0,color='#9ba4af',ls='--',lw=1)
    ax.set_xticks([0,1,2],['90% intensity\n15 paired visits','120% intensity\n15 paired visits','Within-person average\n15 participants'])
    ax.set_ylabel('Global Φ: post minus pre (10⁻³ nats)')
    ax.set_title('Averaging reduces variability; both visit means are positive',loc='left',fontweight='bold')
    ax.text(.01,.98,'Individual visit tests: p = .077 and .149 (uncorrected)\nAverage: family-wise max-T p = .024; intensity difference: p = .773',
            transform=ax.transAxes,va='top',fontsize=9,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
    ax.set_ylim(min(values.min(),-.006*1000)-1,max(values.max(),.008*1000)+5)
    ax.grid(axis='y',alpha=.15)
    for ext in ['png','pdf']: fig.savefig(OUT/('aachen_averaging.'+ext),bbox_inches='tight')
    plt.close(fig)
    if '--figures-only' in sys.argv:
        return

    # Preserve an explicit scan/contrast accounting and exact pooling check.
    scans.groupby(['dataset','session']).size().rename('usable_scans').to_csv(OUT/'scan_counts.csv')
    pairs.groupby(['dataset','comparison','arm']).size().rename('participants').to_csv(OUT/'main_contrast_counts.csv')
    average_error=0.
    for metric in [m for m in pairs if m.endswith('_delta') and m!='delta_fd']:
        name=metric[:-6]
        expected=((a90[name+'_followup']+a120[name+'_followup'])-(a90[name+'_baseline']+a120[name+'_baseline']))/2
        actual=(a90[metric]+a120[metric])/2
        average_error=max(average_error,float(abs(expected-actual).max()))
    assert average_error<1e-12
    # Complete audit inventory retains all earlier broad tests and their p values.
    prior=pd.read_csv(ROUND1/'complete_new_test_inventory.csv').assign(branch='raw_exploration')
    normalized=pd.read_csv(NUMIT/'neural_contrasts.csv')
    delayed=pd.read_csv(FOLLOWUP/'delayed_active_results.csv')
    inventory=pd.concat([prior[['branch','table','row','p']],
                         pd.DataFrame({'branch':'normalized','table':'neural_contrasts','row':normalized.index,'p':normalized.p}),
                         pd.DataFrame({'branch':'delayed_active','table':'delayed_active_results','row':delayed.index,'p':delayed.p})],ignore_index=True)
    invalid=pd.read_csv(FWE/'numerically_unestimable_contrasts.csv')
    inventory['numerically_valid']=True
    for r in invalid.loc[invalid.reference_seed.eq('raw')].itertuples():
        idx=raw_extended.index[raw_extended.dataset.eq(r.dataset)&raw_extended.contrast.eq(r.contrast)&raw_extended.metric.eq(r.metric)]
        inventory.loc[inventory.branch.eq('raw_exploration')&inventory.table.eq('neural')&inventory.row.isin(idx),'numerically_valid']=False
    inventory.to_csv(OUT/'complete_additional_exploration_inventory.csv',index=False)

    report='''# TMS exploratory findings — 25 November 2026

The clearest new neural finding is an increase in whole-brain mean pairwise
whole-minus-sum Phi (Φ) in Aachen, averaged across the two iTBS intensity visits.
It passes family-wise max-T correction across the five global information
measures in both the raw and conditional-null-normalized analyses. This is a
candidate exploratory paper result, not evidence that TMS increases consciousness
or reproduces the psychedelic workspace signature.

## What was averaged

| Dataset | Main participants | Pre scans | Post scans | Analysis |
|---|---:|---:|---:|---|
| SUDMEX | 40: 22 active, 18 sham | 40 t0 | 40 t1 | One pre/post pair per person; active–sham comparison adjusted for that outcome's baseline |
| Aachen | 15 | 15 pre90 + 15 pre120 | 15 post90 + 15 post120 | Compute both matched changes, then average within each participant |

Aachen's main response is `((post90-pre90)+(post120-pre120))/2`. Each post scan
is paired with the baseline from its own visit. The sample size is 15 people,
not 30 independent observations. No fMRI time series, covariance matrices, or
PhiID matrices were averaged before participant metrics were computed. Averaging
the two metric changes equals subtracting the average of the two baseline
metrics from the average of the two post metrics, because all 15 people have
complete observations. The repeated-measures time effect gives the identical
parametric test: F(1,14)=9.0556, p=.009378.

For Phi, mean changes were +0.002043 nats at 90% (95% CI −0.000254 to 0.004341,
uncorrected p=.0771) and +0.001541 at 120% (−0.000624 to 0.003706, p=.1492).
The intensity difference was −0.000503 (−0.004164 to 0.003159, p=.7728).
The two visit changes correlated −.346 across participants; their average is
less variable than either visit alone. This supports reporting an average
across visits, not claiming independent replication at both intensities or
established equivalence of intensities. Eleven of 15 average changes were positive.

All usable completed scans were: SUDMEX t0=41, t1=43, t14=12, t2=19, t3=13,
t4=5 (133 total), and Aachen 15 at each of four sessions (60 total).
The SUDMEX main comparison requires both t0 and t1 and therefore uses 80 scans.
Later visits were not averaged into it. Eleven originally sham participants had
usable t1/t14 pairs and complete three-visit trajectories; delayed-active
comparisons produced no supported neural result. This open-label phase is
confounded with period, expectancy, and retention.

## Correction and surviving findings

All analyses are exploratory. Family-wise error is controlled within scientific
outcome families, separately for each dataset's main comparison: five global
information measures, two mean participation coefficients, the original FPN
gateway endpoint, seven other network gateway counts, eight broadcaster counts,
and two total-role counts. This controls error within each family, not across
the entire manuscript. The original FPN focus comes from the initial TMS
proposal. Secondary intensity discoveries are corrected across the 90%, 120%
and intensity-difference comparisons together within their outcome family.

Aachen uses two-sided, exact single-step max-|t| correction over all 32,768 joint
participant sign flips. The same sign applies to all of a person's outcomes,
preserving their dependence. Validity requires appropriate joint symmetry of
the paired errors; the sign flips do not supply a randomized sham comparison.
Holm correction of parametric p values is shown alongside max-T. SUDMEX uses
baseline-adjusted HC3 regression with within-family Holm correction. The seven
clinical scales form one family; the 35 brain–clinical associations form another.
We do not select the smaller correction separately for each endpoint.

| Aachen global measure | Average change (95% marginal CI) | Family max-T p | Family Holm p |
|---|---|---:|---:|
| Raw Φ | +0.001792 nats (0.000515 to 0.003070) | **.02393** | **.04689** |
| Raw ΦR | +0.002220 nats (0.000356 to 0.004085) | .06580 | .09169 |
| Normalized Φ | +0.7123 percentile points (0.1630 to 1.2617) | **.04828** | .07360 |
| Normalized ΦR | +0.7980 percentile points (0.1274 to 1.4687) | .07690 | .09206 |

Only Phi survives the global family max-T correction. No other main neural
family contains a surviving effect, including gateway counts and participation
coefficients. No extended network, subcortical, clinical-moderation, or multivariate
finding survives its family correction. With all 25 Aachen neural metrics in
one family, raw Phi max-T p=.1434 and normalized Phi p=.2246; family scope must
therefore be clear. The raw Phi paired Hedges' g_z is 0.735.

For baseline-adjusted SUDMEX clinical outcomes, cognitive impulsivity decreased
by 4.776 points (Holm p=.03560) and motor impulsivity by 4.931 (p=.04198).
Total BIS impulsivity decreased by 11.564 (p=.05290) and craving VAS by 1.582
(p=.07489); neither passes family-wise correction. These clinical data and
benefits were already reported by the original trial authors. Their reanalysis
provides clinical context, not the new paper's principal novelty. Brain–clinical
associations remain unsupported.

## Baseline, motion, and robustness

The SUDMEX FPN change contrast was +3.141 gateways before adjustment. Active
participants started with 13.955 FPN gateways versus 16.556 in sham. Adjustment
for baseline alone reduced the contrast to +1.207 (95% CI −1.990 to 4.404,
p=.449); motion alone left it at +3.068. Joint adjustment gave +1.085
(−2.130 to 4.300, p=.498). Baseline imbalance accounts algebraically for 94.5%
of the fitted attenuation, versus 5.5% for measured mean-FD change. Those are
regression contributions, not causal percentages. Regression to the mean is
compatible with this pattern but is not proved. Mean FD cannot exclude every
form of residual motion or physiological noise.

Aachen's raw Phi average stayed positive in every leave-one-person-out analysis
(+0.001392 to +0.002030 nats). Each uncorrected leave-one-out p was below .021;
the two-sided sign-flip p was .00452 and the Wilcoxon p was .00537. These are
dependent robustness checks, not additional independent replications.

Motion adjustment barely changed the raw Phi estimate (+0.001787), but its
HC3 uncertainty was larger: nominal p=.03329, five-global-measure Holm p=.16647.
Normalized Phi behaved similarly (estimate +0.7191 percentile points, nominal
p=.03417, Holm p=.17087). Thus there is little measured-motion attenuation of
the effect size, but corrected support is not retained by every model. No sham
condition in Aachen distinguishes stimulation effects from common pre/post
effects such as time, arousal, scanner drift, or repeated measurement.

## What normalization adds

The [NuMIT paper](https://doi.org/10.1371/journal.pcbi.1013629) compares information
components with a null distribution matched for total mutual information.
It demonstrates PID; our pairwise PhiID application is a methodological
extension. It belongs after PhiID estimation and before participation/role
calculation, not as a rerun of MRI preprocessing or spatial normalization.

We used a shared baseline-derived VAR coefficient pool within each dataset,
the same reference for both conditions, and two independent 32,768-draw null
references per TDMI grid point. All 140 main scans and all 25 endpoints were
processed. Interpolated CDFs differed from direct held-out CDFs by at most
0.000845. Independent scipy Lyapunov solutions matched, and 1,536 finite-record
simulations assessed calibration at retained lengths of 251 and 295 samples.
Finite-record bias remains: for example, average synergy percentile shifts
relative to population values were about −.020 to −.026. This validates useful
numerical properties, not the full generative adequacy of a VAR(1) null for
filtered, deconvolved fMRI. Null-pool uncertainty is not propagated.

Normalized Phi's max-T p values were .04852 and .04822 in the two independent
reference seeds, versus .04828 combined. Its paired-change Monte Carlo error
was small relative to between-person variation. Some discrete role counts
were less stable. Normalization therefore supplies a consistent, borderline
sensitivity for Phi, not an independently replicated discovery. Quantiles are
not nats, and separately normalized atoms/composites are not additive.

## Comparison and proposed manuscript

[Luppi et al.](https://elifesciences.org/articles/88173/figures) emphasized
redundancy-corrected Phi_R and its loss in workspace connections during
unconsciousness. The original whole-minus-sum Phi could instead increase during
unconsciousness. Our whole-brain Phi increase is consequently not evidence of
enhanced consciousness. The closer Phi_R comparison is positive but does not
pass the chosen family correction, and our spatial averaging also differs.

[Yousefi et al.](https://www.biorxiv.org/content/10.64898/2026.09.28.755172v1.full)
reported FPN/DAN gateway recruitment and changes in the distribution of
synergistic and redundant participation. We have not established that pattern
in these TMS samples. A difference in statistical significance does not prove
a difference between interventions; the data do not establish a neural–clinical
dissociation or a shared mechanism. We use published comparisons only, without
private psychedelic data or maps.

**Working title:** Temporal information balance following prefrontal stimulation:
an exploratory PhiID study.

**Central claim:** Across two iTBS visits in healthy volunteers, pre/post changes
in mean pairwise whole-minus-sum Phi provide exploratory evidence of a shift in
temporal information balance, including under conditional-null normalization.
Parallel analysis of a randomized clinical rTMS dataset does not establish
frontoparietal gateway recruitment. These results motivate separating global
information quantities from the regional organization of the synergistic
workspace when studying stimulation.

This supports a cautious exploratory report. It is not yet a strong causal
mechanism paper: the clearest signal is in the small uncontrolled dataset,
depends on averaging visits and the correction family, is weaker for Phi_R,
and loses corrected support in the HC3 motion sensitivity. A sham-controlled
replication of the stated global family and both intensities would be the
decisive next study. Further searching this dataset for a preferred significant
endpoint is not needed to make the present result interpretable.

## Review artifacts

- Original completed publication run: `output/tms_publication/20261125T145928_921789Z`.
- Baseline/motion diagnostics: `output/tms_diagnostics/20261125T152110_809161Z`.
- Normalized branch: `output/tms_numit/20261125T155252_384543Z`.
- Final family-wise tables and numerical verification: `output/tms_familywise/20261125T161051_102285Z`.
- Figures and scan accounting: `output/tms_exploratory_report/20261125`.
- Results notebook: `output/jupyter-notebook/TMS_Exploratory_Findings.ipynb`.

Independent checks verified 977 estimable neural contrasts, 100 delayed-active
contrasts and seven baseline-adjusted clinical models. One secondary amygdala
intensity model had undefined HC3 covariance due to unit leverage and is marked
unestimable in the final table. None of the main findings depends on that row.
Production processing, original outputs, frozen cohorts and D:/projects/GWT
were not modified.
'''
    report_path=ROOT/'docs/TMS_EXPLORATORY_FINDINGS_20261125.md'
    report_path.write_text(report,encoding='utf-8')

    notebook=ROOT/'output/jupyter-notebook/TMS_Exploratory_Findings.ipynb'
    nb=nbf.read(notebook,as_version=4)
    markdown=nbf.v4.new_markdown_cell; code=nbf.v4.new_code_cell
    nb.cells=[
        markdown('# TMS exploratory findings\n\nThis notebook reviews completed TMS-only analyses. It does not rerun preprocessing or analyze private psychedelic data. The full report is `docs/TMS_EXPLORATORY_FINDINGS_20261125.md`.\n\nQuestion: does prefrontal stimulation change global information measures or the organization of the synergistic workspace?'),
        code("from pathlib import Path\nimport json\nimport pandas as pd\nfrom IPython.display import display, Image\nROOT = Path(r'D:/projects/TMS_GWT')\nSOURCE = ROOT/'output/tms_publication/20261125T145928_921789Z'\nFWE = ROOT/'output/tms_familywise/20261125T161051_102285Z'\nNUMIT = ROOT/'output/tms_numit/20261125T155252_384543Z'\nFIGURES = ROOT/'output/tms_exploratory_report/20261125'\nverification = json.loads((FWE/'independent_verification.json').read_text())\nassert verification['status'] == 'passed_for_estimable_contrasts'\nprint('Verified neural contrasts:', verification['independently_verified_neural_contrasts'])"),
        markdown('## Scan counts and averaging\n\nSUDMEX: 40 baseline/two-week pairs, 22 active and 18 sham, with no averaging across follow-up visits. Aachen: 15 participants, each with two separate pre/post visits. Its main test uses the mean of the 90% and 120% changes within each person. The sample size remains 15.'),
        code("display(pd.read_csv(FIGURES/'scan_counts.csv'))\ndisplay(pd.read_csv(FIGURES/'main_contrast_counts.csv'))\ndisplay(Image(filename=str(FIGURES/'aachen_averaging.png')))"),
        markdown('## Global information family\n\nExact two-sided max-T controls family-wise error across all five global measures, preserving their dependence through joint participant sign flips. Holm is shown alongside it. Error is controlled within each family, not across all families. The intervals shown are marginal 95% intervals. Normalized values are conditional-null percentiles, not nats.'),
        code("results = pd.read_csv(FWE/'main_comparisons_familywise.csv')\nglobal_results = results.query(\"dataset == 'aachen_itbs' and family == 'global_phiid'\")\ndisplay(global_results[['branch','metric','estimate','ci_low','ci_high','p_maxT_family','p_holm_family','p_maxT_all25']])\ndisplay(Image(filename=str(FIGURES/'aachen_global_family.png')))"),
        markdown('## Other neural and clinical outcomes\n\nOnly global Phi passes a main neural family correction. SUDMEX FPN evidence weakens after baseline adjustment. Cognitive and motor impulsivity pass Holm across the seven clinical scales; total impulsivity and craving do not. Clinical benefits are already known from the source trial, and brain–clinical associations are unsupported.'),
        code("display(results.query(\"metric == 'gateways_FPN'\")[[ 'branch','dataset','estimate','ci_low','ci_high','p_holm_family','p_maxT_family']])\nclinical = pd.read_csv(FWE/'clinical_familywise.csv')\ndisplay(clinical[['metric','n','estimate','ci_low','ci_high','p_holm_clinical7']])"),
        markdown('## Baseline and motion checks\n\nSUDMEX FPN attenuation is principally associated with baseline imbalance. In Aachen, motion adjustment barely changes the Phi estimate, but wider HC3 uncertainty does not retain five-measure Holm support. The source diagnostics explain the attenuation algebra and its causal limits.'),
        code("display(pd.read_csv(ROOT/'output/tms_diagnostics/20261125T152110_809161Z/fpn_model_comparison.csv')[[ 'model','estimate','ci_low','ci_high','p']])\ndisplay(pd.read_csv(FWE/'aachen_motion_global_familywise.csv').query(\"metric in ['global_phi','global_phi_r']\")[[ 'branch','metric','estimate','ci_low','ci_high','p','p_holm_global5']])"),
        markdown('## Normalization quality and repeatability\n\nThis extends a published PID null-model approach to PhiID. Finite-record bias and uncertainty in the baseline-derived null remain limitations. The normalized result is a sensitivity analysis of a different estimand. It is not a second sample or an independent replication.'),
        code("display(json.loads((NUMIT/'validation_before_effects.json').read_text()))\ndisplay(pd.read_csv(FWE/'numit_global_family_seed_stability.csv').query(\"metric in ['global_phi','global_phi_r']\"))"),
        markdown('## Interpretation\n\n**Candidate narrative:** prefrontal iTBS is followed by a change in temporal information balance averaged across two visits, with uncertain evidence for the regional gateway recruitment reported for psychedelics.\n\nPhi and revised Phi_R must remain distinct: Luppi\'s principal loss-of-consciousness result concerns Phi_R, while the stronger present result concerns Phi. This does not measure an increase in consciousness. Aachen has no sham condition, each intensity is individually uncertain, and the clinical trial does not provide a confirming neural effect.\n\nAll analyses are exploratory. Full result tables, numerical limitations and analysis decisions remain in the project record.'),
        code("print('Complete findings report:', ROOT/'docs/TMS_EXPLORATORY_FINDINGS_20261125.md')\nprint('Additional exploratory test rows:', len(pd.read_csv(FIGURES/'complete_additional_exploration_inventory.csv')))\nprint('Independent checks:', verification['status'])")]
    nb.metadata['kernelspec']={'display_name':'Python (gwt)','language':'python','name':'python3'}
    nbf.write(nb,notebook)
    (OUT/'build_receipt.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),
                                                   'max_pooling_identity_error':average_error,
                                                   'report':str(report_path),'notebook':str(notebook),
                                                   'additional_inventory_rows':len(inventory)},indent=2)+'\n')
    print('REPORT:',report_path)
    print('NOTEBOOK:',notebook)
    print('FIGURES:',OUT)


if __name__=='__main__': main()
