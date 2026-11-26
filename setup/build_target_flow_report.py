"""Readable result tables, effect plot and executed review notebook."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import nbformat
from nbclient import NotebookClient

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/targeted_information_flow'/(ROOT/'output/targeted_information_flow/LATEST.txt').read_text().strip()
LABELS={'net_te':'Net outward transfer entropy','net_unique':'Net unique-source contribution',
        'net_synergistic':'Net synergistic contribution','persistent_synergy':'Persistent synergy',
        'persistent_redundancy':'Persistent redundancy'}

def mdtable(d,columns):
    rows=[list(columns)]+[['---']*len(columns)]+[[str(v).replace('|','\\|') for v in row] for row in d[columns].itertuples(index=False,name=None)]
    return '\n'.join('| '+' | '.join(row)+' |' for row in rows)

def main():
    d=pd.read_csv(OUT/'inference.csv')
    main=d.loc[d.definition.eq('common_baseline_target')&d.variant.eq('lag1_full')&d.model.eq('baseline')].copy()
    show=main.copy();show['Measure']=show.metric.map(LABELS)
    show['Scope']=show.scope.map({'left':'Left target','left_minus_right':'Left minus right'})
    show['Active−sham (nats)']=show.estimate.map(lambda x:f'{x:+.6f}')
    show['95% CI']=show.apply(lambda r:f'[{r.ci_low:+.6f}, {r.ci_high:+.6f}]',axis=1)
    show['Uncorrected permutation p']=show.p_permutation.round(4)
    show['FWER p (10 tests)']=show.p_maxT10.round(4)
    cols=['Scope','Measure','Active−sham (nats)','95% CI','Uncorrected permutation p','FWER p (10 tests)']
    show[cols].to_csv(OUT/'sudmex_readable_primary_results.csv',index=False)
    fig,axes=plt.subplots(5,1,figsize=(8.5,10),layout='constrained')
    for ax,(metric,label) in zip(axes,LABELS.items()):
        block=main.loc[main.metric.eq(metric)].set_index('scope')
        for j,scope in enumerate(['left','left_minus_right']):
            r=block.loc[scope]
            ax.errorbar(r.estimate,1-j,xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],fmt='o',capsize=4,
                        color=['#245c85','#ae5b38'][j],markersize=6)
        ax.axvline(0,color='#777777',lw=.8,ls='--');ax.set_yticks([1,0],['Left target','Left minus right'])
        ax.set_ylim(-.6,1.6);ax.set_title(label,loc='left',fontsize=11,fontweight='bold')
        ax.xaxis.set_major_locator(MaxNLocator(4));ax.ticklabel_format(axis='x',style='sci',scilimits=(-3,3))
        ax.spines[['top','right']].set_visible(False)
        ax.set_xlabel('Baseline-adjusted active − sham difference (nats)',fontsize=9)
    fig.suptitle('SUDMEX: no corrected evidence of altered target-centred flow\n40 participants; marginal 95% HC3 intervals; all maxT10 p ≥ 0.406',fontsize=13)
    fig.savefig(OUT/'sudmex_primary_effects.png',dpi=150);fig.savefig(OUT/'sudmex_primary_effects.pdf');plt.close(fig)
    sens=d.loc[d.dataset.eq('sudmex_tms')].groupby(['definition','variant','model']).agg(n=('n','first'),minimum_maxT10_p=('p_maxT10','min')).reset_index()
    sens.minimum_maxT10_p=sens.minimum_maxT10_p.round(4)
    checked=json.loads((OUT/'independent_verification.json').read_text())
    matlab=json.loads((OUT/'matlab_verification.json').read_text())
    report=f'''# Target-based TMS information-flow analysis

The sham-controlled SUDMEX comparison does not support a stimulation-related change in net outward information flow, persistent synergy or persistent redundancy at the DLPFC target neighbourhood. All ten primary tests have family-wise p ≥ {main.p_maxT10.min():.3f}. This is insufficient evidence for a change, not proof of no effect. All original Aachen results are **suspended pending corrected anatomical preprocessing and downstream reanalysis**; their present numbers must not be used as evidence of absence, of incorrect stimulation placement, or of a validated global effect.

## What was tested

SUDMEX contributed 40 complete acute-phase MRI pairs: 22 active and 18 sham, baseline versus the two-week visit only. The fixed left target is a 10-mm neighbourhood around the mean of all 21 published baseline MNI localizations: (−33.14, 43.71, 37.14) mm. It was defined before extracting targeted information outcomes. The mirrored right region is a spatial comparator, not an untreated control. Individual baseline coordinates were available for 19 usable participants (12 active, 7 sham); these form a separate localization sensitivity. The source specifies MNI but not its precise template variant, so localization is approximate. [Data descriptor](https://www.nature.com/articles/s41597-024-03242-y); [public localization workbook](https://zenodo.org/records/10409461).

The unchanged deconvolved 454-region time series were analysed at one TR (2 s). For target parcel X and remote parcel Y, outward transfer entropy is I(Xpast;Yfuture | Ypast), and inward transfer entropy reverses X and Y. Their difference describes net directional predictability. In the ΦID lattice:

- TEout = xtr + xty + str + sty; TEin = ytr + ytx + str + stx.
- Net unique-source contribution = (xtr + xty) − (ytr + ytx).
- Net synergistic contribution = sty − stx.
- Net TE equals those two net contributions summed, at the scan level.
- Persistent synergy is sts and persistent redundancy is rtr. These last two are symmetric pairwise measures; they do not encode a direction of flow.

The synergistic TE contribution requires the joint past of both regions; it cannot be interpreted as information originating solely in the stimulated region. This follows the [ΦID framework](https://doi.org/10.1073/pnas.2423297122). Seed parcels were weighted by anatomical overlap; all remote destinations received equal weight, excluding both seed sets. Both sides used the same destination pool. Effects are participant-level averages, not independent tests of thousands of edges. We did not search destination parcels or alter regions according to their effects.

## Primary randomized comparison

Post-treatment outcomes were regressed on baseline and active/sham assignment, with HC3 uncertainty. Ten tests combine five measures with two questions: change at the left target, and change at the left target relative to the right comparator. Shared treatment-label permutations (32,768 draws, preserving 22/18 counts) give single-step max-|t| family-wise p values. Label exchangeability among imaging completers is assumed; this is not a reconstruction of undisclosed trial randomization blocks. Marginal confidence intervals are not multiplicity-adjusted.

{mdtable(show,cols)}

Even before correction, none of the ten primary tests has p < .05. The smallest p concerns lower persistent redundancy in the left-minus-right comparison: p={main.p_permutation.min():.3f} uncorrected and {main.p_maxT10.min():.3f} corrected. It is not a supported positive finding.

## Robustness checks

No test survives ten-test correction after adding change in mean framewise displacement, changing the lag to 4 s, matching pre/post record lengths within participant, or using individual baseline locations. Equal-length checks reuse the already deconvolved signals; they do not refit HRFs. Motion adjustment is a diagnostic sensitivity, because post-treatment motion could itself be influenced by treatment.

{mdtable(sens,list(sens.columns))}

Leave-one-participant-out coefficients are saved separately as influence diagnostics; no omission was selected as the main analysis. All five measures, both scopes, both localization definitions and the three signal/lag variants remain available in the inference table. These bounded correction families do not correct the entire historical project search.

## Aachen alignment status

Public masks were found for eight usable participants. The independent audit established that the deposited Aachen T1s were already skull-stripped but underwent forced repeat extraction. The resulting masks removed substantial outer cortex and compromised anatomy-dependent processing. The deposited target masks align with the full native anatomy; their appearance outside the clipped standardized T1 did not establish incorrect stimulation placement. Numerical reproduction of a saved image using its own transform was insufficient to validate registration. The corrected subject-06 pilot preserves 99.8% of native support and passes visual native/standard-space and atlas-placement checks. The remaining cohort is being rerun in an isolated branch, followed by CONN, parcellation, deconvolution and ΦID. All old Aachen estimates, including whole-brain Φ, remain provisional until recomputed. See `PROVENANCE/AACHEN_ANATOMICAL_REPROCESSING_20261125.md` and `output/aachen_alignment_review/pilot_release_receipt.json`. The deposited seed is an intended analysis ROI; exact coil or peak-field location is not independently verified.

The separate SUDMEX audit covered 44 anatomical images and 142 functional runs. It found no comparable T1 over-extraction, so a blanket anatomical rerun is not indicated. Inferior functional coverage and some local alignment limitations remain; these were not resolved by the numerical checks on information estimates. See `output/sudmex_registration_audit/SUDMEX_FULL_PREPROCESSING_QC.md`.

## What can be called causal

Randomized active-versus-sham assignment provides a basis for estimating a treatment effect on these measured information summaries, subject to complete-case selection and other trial limitations. This analysis found no convincing such effect. Pairwise lagged fMRI predictability does not identify direct neuronal causal connections: unmeasured common inputs, indirect paths, haemodynamic modelling and slow sampling remain relevant. A net-flow average can also hide offsetting changes across destinations or simultaneous changes in both directions. The results therefore neither establish nor rule out every possible circuit effect.

## Verification and reproducibility

Original MATLAB PhiIDFull agreed for {matlab['cases']} fixed validation pairs across both datasets and both lags, all 16 atoms per pair: maximum error {matlab['max_absolute_error']:.2g} nats. Reconstructed five original matrix measures agreed within {checked['max_saved_matrix_error']:.2g} nats. Transfer entropy independently matched Gaussian conditional mutual information within {checked['max_independent_CMI_error']:.2g} nats. Synthetic unidirectional data and source swapping verified the direction convention. Independent full-design HC3/permutation calculations reproduced inference, and {checked['input_hashes_checked']} input hashes were checked.

Run: `{OUT.relative_to(ROOT).as_posix()}`. Reproduction scripts: `setup/analyse_targeted_information_flow.py`, `setup/verify_targeted_information_flow.py`, and `setup/verify_target_flow_matlab.m`. All calculations read frozen TMS signals only; HCP/psychedelic cohorts and production processing code were not rerun or changed. New directed measures were not normalized using the old symmetric-measure NuMIT references.
'''
    control_base=ROOT/'output/control_workspace_qc'
    aachen_qc_latest=control_base/'LATEST_aachen_corrected.txt'
    aachen_qc_run=None
    if aachen_qc_latest.exists():
        aachen_qc_run=control_base/aachen_qc_latest.read_text().strip()
        aachen_qc=json.loads((aachen_qc_run/'summary.json').read_text())
        corrected_status=json.loads((ROOT/'logs/aachen_corrected_20261125/downstream/status.json').read_text())
        assert corrected_status['phase']=='complete' and corrected_status['phiid_scans']==60
        report=report.replace('suspended pending corrected anatomical preprocessing and downstream reanalysis',
                              'suspended pending statistical reanalysis of the corrected inputs')
        previous='The corrected subject-06 pilot preserves 99.8% of native support and passes visual native/standard-space and atlas-placement checks. The remaining cohort is being rerun in an isolated branch, followed by CONN, parcellation, deconvolution and ΦID.'
        assert previous in report
        report=report.replace(previous,'All 15 participants and 60 runs have now passed independent anatomical review. The isolated corrected branch has completed CONN, parcellation, deconvolution and ΦID with output validation. Corrected native support retention is 99.63–100%; regional functional-coverage limitations remain documented.')
        report=report.replace('output/aachen_alignment_review/pilot_release_receipt.json',
                              'output/aachen_alignment_review/corrected_cohort/cohort_release_receipt.json')
    control_run=control_base/(control_base/'LATEST.txt').read_text().strip()
    control=json.loads((control_run/'summary.json').read_text())
    control_table=pd.DataFrame([{'Condition':name,'N':s['n'],'Workspace regions':s['n_workspace'],
                                'Shared with HCP':s['shared_regions'],'Dice':round(s['dice_hcp'],3),
                                'Gradient rho':round(s['gradient_rho_hcp'],3)}
                               for name,s in control['conditions'].items()])
    control_note='The workspace was derived independently from the control scans with unchanged GWT core functions. The HCP mask and gradient are byte-identical to the original GWT project reference. The matched sham sample has 18 participants at both visits; later open-label treatment is excluded. This is descriptive pattern QC, not proof of valid registration at every parcel or a treatment-effect result. Aachen workspace QC awaits corrected signals.'
    if aachen_qc_run:
        control_note=control_note.replace('Aachen workspace QC awaits corrected signals.',
                                         'Corrected Aachen baseline-workspace results are reported separately below.')
    report+='\n## Control workspace quality check\n\n'+control_note+'\n\n'+mdtable(control_table,list(control_table.columns))+'\n'
    if aachen_qc_run:
        aachen_note='Aachen has no sham arm. These workspaces use the two prestimulation visits, separately and averaged within each of the same 15 participants. The HCP mask remains fixed for treatment inference. Bootstraps, split halves and uniform leave-one-participant-out outputs describe sampling stability; no omission selects the primary map.'
        aachen_table=pd.DataFrame([{'Condition':name,'N':s['n'],'Workspace regions':s['n_workspace'],
                                    'Shared with HCP':s['shared_regions'],'Dice':round(s['dice_hcp'],3),
                                    'Gradient rho':round(s['gradient_rho_hcp'],3)}
                                   for name,s in aachen_qc['conditions'].items()])
        report+='\n## Corrected Aachen control workspace\n\n'+aachen_note+'\n\n'+mdtable(aachen_table,list(aachen_table.columns))+'\n'
    (ROOT/'docs/TMS_TARGETED_INFORMATION_FLOW_20261125.md').write_text(report,encoding='utf-8')
    (OUT/'RESULTS.md').write_text(report,encoding='utf-8')
    nbpath=ROOT/'output/jupyter-notebook/TMS_Targeted_Information_Flow.ipynb'
    nb=nbformat.read(nbpath,as_version=4)
    cells=[
      nbformat.v4.new_markdown_cell('# Targeted TMS information flow\n\n'+report.split('## What was tested')[0].split('\n\n',1)[1]),
      nbformat.v4.new_markdown_cell('## Definitions\n\n'+report.split('## What was tested\n\n',1)[1].split('## Primary randomized comparison',1)[0]),
      nbformat.v4.new_code_cell(f"from pathlib import Path\nimport json\nimport pandas as pd\nfrom IPython.display import display, Image\nROOT = Path(r'{ROOT}')\nRUN = ROOT / r'{OUT.relative_to(ROOT)}'\nresults = pd.read_csv(RUN / 'inference.csv')\nverification = json.loads((RUN / 'independent_verification.json').read_text())\nassert verification['status'] == 'passed'\nprint('Verified run:', RUN.name)"),
      nbformat.v4.new_markdown_cell('## Primary sham-controlled comparison\n\nFive measures × two spatial questions form a single ten-test family. All marginal intervals include zero. Full outputs retain raw and corrected p values.'),
      nbformat.v4.new_code_cell("display(pd.read_csv(RUN / 'sudmex_readable_primary_results.csv'))\ndisplay(Image(filename=str(RUN / 'sudmex_primary_effects.png')))"),
      nbformat.v4.new_markdown_cell('## Sensitivities\n\nMotion, lag, paired record length and individual baseline targets do not yield corrected evidence. The localization subset is smaller and is not an independent replication.'),
      nbformat.v4.new_code_cell("sud = results.query(\"dataset == 'sudmex_tms'\")\ndisplay(sud.groupby(['definition','variant','model']).agg(n=('n','first'),minimum_maxT10_p=('p_maxT10','min')).round(4))"),
      nbformat.v4.new_markdown_cell('## Aachen\n\n'+report.split('## Aachen alignment status\n\n',1)[1].split('## What can be called causal',1)[0]),
      nbformat.v4.new_markdown_cell('## Verification and interpretation\n\n'+report.split('## What can be called causal\n\n',1)[1].split('## Verification and reproducibility',1)[0]),
      nbformat.v4.new_code_cell("print('Historical numerical verification: its anatomical commentary predates and is superseded by the corrected-cohort audit above.')\ndisplay(verification)\ndisplay(json.loads((RUN / 'matlab_verification.json').read_text()))"),
      nbformat.v4.new_markdown_cell('## Control workspace quality check\n\n'+control_note),
      nbformat.v4.new_code_cell(f"CONTROL_RUN = ROOT / r'{control_run.relative_to(ROOT)}'\ncontrol = json.loads((CONTROL_RUN / 'summary.json').read_text())\ndisplay(pd.DataFrame(control['conditions']).T[['n','n_workspace','shared_regions','dice_hcp','gradient_rho_hcp']])\ndisplay(control['sham_pre_post'])\ndisplay(Image(filename=str(CONTROL_RUN / 'control_workspace_comparison.png')))" )]
    if aachen_qc_run:
        cells.extend([
            nbformat.v4.new_markdown_cell('## Corrected Aachen control workspace\n\n'+aachen_note),
            nbformat.v4.new_code_cell(f"AACHEN_CONTROL_RUN = ROOT / r'{aachen_qc_run.relative_to(ROOT)}'\naachen_control = json.loads((AACHEN_CONTROL_RUN / 'summary.json').read_text())\ndisplay(pd.DataFrame(aachen_control['conditions']).T[['n','n_workspace','shared_regions','dice_hcp','gradient_rho_hcp','leave_one_out_dice_hcp_range']])\ndisplay(aachen_control['baseline_visit_similarity'])\ndisplay(Image(filename=str(AACHEN_CONTROL_RUN / 'control_workspace_comparison.png')))"),
        ])
    nb.cells=cells
    nb.metadata.kernelspec={'display_name':'Python (gwt)','language':'python','name':'gwt'}
    # Preserve the scaffold, but execute this compact review directly with the current interpreter kernel.
    nbformat.write(nb,nbpath)
    NotebookClient(nb,timeout=180,kernel_name='gwt',resources={'metadata':{'path':str(ROOT)}}).execute()
    nbformat.write(nb,nbpath)
    assert not any(o.output_type=='error' for c in nb.cells if c.cell_type=='code' for o in c.outputs)
    (OUT/'report_verification.json').write_text(json.dumps({'notebook':str(nbpath),'code_cells':sum(c.cell_type=='code' for c in nb.cells),'execution_errors':0},indent=2),encoding='utf-8')
    print('Report and executed notebook:',nbpath)

if __name__=='__main__':main()
