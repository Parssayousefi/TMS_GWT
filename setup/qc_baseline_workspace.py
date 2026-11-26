"""TMS-only baseline workspace QC using the unchanged GWT implementation.

Frozen HCP/GWT masks and gradients are read for comparison, never recomputed.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,importlib.util,json,warnings
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--dataset',choices=['sudmex_tms','aachen_itbs'],default='sudmex_tms')
    parser.add_argument('--task',default='t0');parser.add_argument('--permit-provisional-aachen',action='store_true')
    args=parser.parse_args()
    if args.dataset=='aachen_itbs' and not args.permit_provisional_aachen:
        raise RuntimeError('Aachen must wait for corrected preprocessing; do not use old alignment as a validated QC result.')
    script=ROOT/'pipeline/scripts/analysis/05_group_average_workspace.py'
    spec=importlib.util.spec_from_file_location('gwt_workspace_frozen',script);gwt=importlib.util.module_from_spec(spec);spec.loader.exec_module(gwt)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out=ROOT/'output/baseline_workspace_qc'/stamp;out.mkdir(parents=True)
    # Production function operates only on the explicitly selected TMS baseline.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='Mean of empty slice')
        result=gwt.compute_workspace_for_dataset(args.dataset,project_dir=ROOT/'pipeline',task=args.task,synergy_suffix='synergy')
    saved=gwt.save_single_dataset_results(result,out)
    report=gwt.print_single_dataset_report(result)
    (out/'production_workspace_report.txt').write_text(report,encoding='utf-8')
    reference=ROOT/'pipeline/results/group_analysis/hcp/rest'
    original=Path('D:/projects/GWT/pipeline/results/group_analysis/hcp/rest')
    files=['hcp_rest_workspace_mask.npy','hcp_rest_synergy_redundancy_gradient.npy',
           'hcp_rest_synergy_node_strength.npy','hcp_rest_redundancy_node_strength.npy']
    checks=[]
    for name in files:
        local=reference/name;source=original/name
        checks.append({'file':name,'copy_sha256':hashlib.sha256(local.read_bytes()).hexdigest(),
                       'original_sha256':hashlib.sha256(source.read_bytes()).hexdigest() if source.exists() else None})
    pd.DataFrame(checks).to_csv(out/'reference_copy_verification.csv',index=False)
    assert all(r['original_sha256'] is None or r['copy_sha256']==r['original_sha256'] for r in checks)
    hcp=np.load(reference/files[0]).astype(bool);grad=np.load(reference/files[1]);mask=result['workspace_mask']
    assert mask.shape==hcp.shape==(454,)
    overlap=int((mask&hcp).sum());union=int((mask|hcp).sum())
    syn=result['all_matrices']['synergy'];red=result['all_matrices']['redundancy'];n=len(syn)
    offdiag=~np.eye(454,dtype=bool)
    assert np.isfinite(syn[:,offdiag]).all() and np.isfinite(red[:,offdiag]).all()
    strengths_syn=np.nansum(syn,axis=2);strengths_red=np.nansum(red,axis=2)
    np.testing.assert_allclose(strengths_syn.mean(0),result['synergy_strength'],atol=1e-10)
    np.testing.assert_allclose(strengths_red.mean(0),result['redundancy_strength'],atol=1e-10)
    # Independently reconstruct ranks (pandas rather than production scipy rankdata).
    independent=pd.Series(strengths_syn.mean(0)).rank()-pd.Series(strengths_red.mean(0)).rank()
    np.testing.assert_array_equal(independent.to_numpy()>0,mask)
    np.testing.assert_array_equal(independent.to_numpy(),result['gradient'])
    rng=np.random.default_rng(2026112531)
    boot_masks=[];bootstrap=[];split=[]
    for draw in range(2000):
        idx=rng.integers(n,size=n)
        m,g,_,_=gwt.identify_workspace(strengths_syn[idx].mean(0),strengths_red[idx].mean(0))
        boot_masks.append(m)
        bootstrap.append({'draw':draw,'n_workspace':int(m.sum()),'dice_hcp':gwt.dice_coefficient(m,hcp),
                           'dice_full_sample':gwt.dice_coefficient(m,mask),'rho_hcp':stats.spearmanr(g,grad).statistic})
    for draw in range(1000):
        order=rng.permutation(n);a=order[:n//2];b=order[n//2:]
        ma,ga,_,_=gwt.identify_workspace(strengths_syn[a].mean(0),strengths_red[a].mean(0))
        mb,gb,_,_=gwt.identify_workspace(strengths_syn[b].mean(0),strengths_red[b].mean(0))
        split.append({'draw':draw,'dice_halves':gwt.dice_coefficient(ma,mb),'rho_halves':stats.spearmanr(ga,gb).statistic})
    bootstrap=pd.DataFrame(bootstrap);split=pd.DataFrame(split)
    bootstrap.to_csv(out/'participant_bootstrap.csv',index=False);split.to_csv(out/'split_half_stability.csv',index=False)
    labels=np.load(reference/'hcp_rest_rsn_labels.npy',allow_pickle=True)
    nodes=pd.DataFrame({'roi':np.arange(1,455),'network':labels,'baseline_workspace':mask,'hcp_workspace':hcp,
                        'baseline_gradient':result['gradient'],'hcp_gradient':grad,'bootstrap_membership_probability':np.mean(boot_masks,axis=0)})
    nodes.to_csv(out/'regional_comparison.csv',index=False)
    networks=nodes.groupby('network').agg(n_regions=('roi','size'),baseline_workspace=('baseline_workspace','sum'),hcp_workspace=('hcp_workspace','sum'))
    networks.to_csv(out/'network_comparison.csv')
    summary={'dataset':args.dataset,'baseline_session':args.task,'n_baseline_subjects':n,'n_regions':454,
      'baseline_workspace_size':int(mask.sum()),'hcp_workspace_size':int(hcp.sum()),'overlap_regions':overlap,
      'union_regions':union,'dice':gwt.dice_coefficient(mask,hcp),'jaccard':overlap/union,
      'fraction_hcp_retained':overlap/hcp.sum(),'fraction_baseline_in_hcp':overlap/mask.sum(),
      'gradient_spearman_rho':stats.spearmanr(result['gradient'],grad).statistic,
      'bootstrap_dice_hcp_95_percentiles':bootstrap.dice_hcp.quantile([.025,.975]).tolist(),
      'bootstrap_gradient_rho_95_percentiles':bootstrap.rho_hcp.quantile([.025,.975]).tolist(),
      'split_half_dice_median':split.dice_halves.median(),'split_half_dice_95_percentiles':split.dice_halves.quantile([.025,.975]).tolist(),
      'baseline_original_GWT_reference_mask_identical':checks[0]['copy_sha256']==checks[0]['original_sha256'],
      'interpretation':'Descriptive anatomical-pattern QC. No spatial significance claim or predefined pass/fail threshold; does not validate registration.',
      'inference':'Bootstrap resamples participants in TMS only; HCP reference fixed. No naive independent-parcel p values.',
      'subjects':result['subjects'],'workspace_script_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),
      'outputs':str(saved),'reference_checks':checks}
    (out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    ax=axes[0];ax.scatter(grad,result['gradient'],s=12,alpha=.65,c=np.where(mask==hcp,'#2a7196','#c56b40'))
    ax.axhline(0,color='gray',lw=.7);ax.axvline(0,color='gray',lw=.7)
    ax.set_xlabel('Frozen HCP synergy − redundancy rank');ax.set_ylabel('TMS baseline synergy − redundancy rank')
    ax.set_title(f'Gradient correlation: ρ = {summary["gradient_spearman_rho"]:.3f}')
    axes[1].hist(split.dice_halves,bins=25,color='#9cc3d7',edgecolor='white')
    axes[1].axvline(summary['dice'],color='#ac522f',lw=2,label='Baseline versus HCP')
    axes[1].set_xlabel('Dice overlap');axes[1].set_ylabel('Random balanced splits');axes[1].legend(frameon=False,fontsize=9)
    axes[1].set_title('Within-dataset split-half stability')
    fig.suptitle(f'{args.dataset}: independent baseline workspace QC, n={n}')
    fig.savefig(out/'baseline_workspace_qc.png',dpi=160);fig.savefig(out/'baseline_workspace_qc.pdf');plt.close(fig)
    text=f'''# Baseline workspace quality control

The unchanged GWT workspace algorithm was run on **{n} usable {args.dataset} baseline scans ({args.task})**, using persistent synergy and redundancy, arithmetic group means, node strengths, and the sign of synergy-minus-redundancy ranks. Post-treatment scans were not used. The primary intervention analysis has 40 complete pairs; baseline QC may also include a participant whose follow-up was unavailable.

The independently derived baseline workspace contains **{int(mask.sum())} regions**, versus **{int(hcp.sum())} in the frozen HCP reference**. They share **{overlap} regions**: Dice **{summary['dice']:.3f}**, Jaccard **{summary['jaccard']:.3f}**, and gradient Spearman correlation **{summary['gradient_spearman_rho']:.3f}**. The baseline workspace includes **{summary['fraction_hcp_retained']:.1%}** of HCP workspace regions. This is a descriptive comparison, not a test that assumes spatially independent parcels.

The HCP mask copied into this project is byte-identical to the one in the original GWT project, as are the reference gradient and node-strength files. Therefore “our mask” and that GWT HCP mask are the same reference, not two independent replications. The original GWT files were read only; no HCP or psychedelic computation was rerun.

In 2,000 participant bootstraps, HCP Dice ranged from {summary['bootstrap_dice_hcp_95_percentiles'][0]:.3f} to {summary['bootstrap_dice_hcp_95_percentiles'][1]:.3f} between the 2.5th and 97.5th percentiles. Across 1,000 balanced random split halves, median Dice was {summary['split_half_dice_median']:.3f} (2.5th–97.5th percentiles {summary['split_half_dice_95_percentiles'][0]:.3f}–{summary['split_half_dice_95_percentiles'][1]:.3f}). These describe sample stability conditional on the processing and fixed HCP reference, not external validation. Clinical population, scan length and acquisition differences can affect similarity.

Aachen baseline QC will be repeated after correcting its anatomical preprocessing. A plausible workspace from the old Aachen signals would not repair or validate the identified brain-extraction error. No settings or thresholds are tuned to maximize reference overlap.

Outputs: `{out.relative_to(ROOT).as_posix()}`. The production workspace function and export function were called directly; an independent pandas ranking reproduced the mask exactly. New QC exports are separate from production and frozen reference outputs.
'''
    (out/'README.md').write_text(text,encoding='utf-8')
    (ROOT/'docs/TMS_BASELINE_WORKSPACE_QC_20261125.md').write_text(text,encoding='utf-8')
    (out.parent/f'LATEST_{args.dataset}.txt').write_text(out.name,encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['subjects','reference_checks']},indent=2))

if __name__=='__main__':main()
