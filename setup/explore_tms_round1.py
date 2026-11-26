"""Transparent exploratory TMS screen, preserving the complete test inventory."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'output/tms_publication/20261125T145928_921789Z'
sys.path.insert(0,str(ROOT/'pipeline/scripts/analysis'))
import tms_publication as tms
KEY = ['dataset','subject_id','session']


def regression(frame, outcome, predictors, term):
    data = frame[[outcome,*predictors]].replace([np.inf,-np.inf],np.nan).dropna()
    if data[outcome].std()<1e-13 or data[term].nunique()<2 and term!='intercept':
        return {'n':len(data),'estimate':np.nan,'p':np.nan,'ci_low':np.nan,'ci_high':np.nan,'status':'constant'}
    keep = [c for c in predictors if c==term or data[c].nunique()>1]
    x = data[keep].astype(float)
    if 'intercept' not in x:
        x = sm.add_constant(x,has_constant='add')
    if np.linalg.matrix_rank(x)<x.shape[1] or len(x)<=x.shape[1]+2:
        return {'n':len(data),'estimate':np.nan,'p':np.nan,'ci_low':np.nan,'ci_high':np.nan,'status':'rank_or_sample_limit'}
    model = sm.OLS(data[outcome].astype(float),x).fit(cov_type='HC3',use_t=True)
    low,high = model.conf_int().loc[term]
    return {'n':len(data),'estimate':model.params[term],'ci_low':low,'ci_high':high,
            'p':model.pvalues[term],'se':model.bse[term],'status':'ok'}


def one_sample(values):
    values = np.asarray(values,float)
    assert np.isfinite(values).all()
    if values.std(ddof=1)<1e-13:
        return {'n':len(values),'estimate':values.mean(),'ci_low':np.nan,'ci_high':np.nan,'p':np.nan,'status':'constant'}
    result = tms.paired_summary(values)
    return result


def orthogonal_columns(matrix, center=False):
    if center: matrix = matrix-matrix.mean(0)
    u,s,v = np.linalg.svd(matrix,full_matrices=False)
    rank = int((s>max(s)*1e-10).sum())
    return u[:,:rank],v[:rank].T


def residual_permutation_test(outcome,baseline,arm,rng,n_perm=9999):
    # Freedman-Lane approximation: all outcomes for a participant move together.
    scale = baseline.std(0,ddof=1)
    b = (baseline[:,scale>1e-13]-baseline[:,scale>1e-13].mean(0))/scale[scale>1e-13]
    q0,_ = orthogonal_columns(np.column_stack([np.ones(len(b)),b]))
    m0 = np.eye(len(b))-q0@q0.T
    a = m0@arm
    aa = a@a
    assert aa>1e-8
    m = m0-np.outer(a,a)/aa
    _, basis = orthogonal_columns(outcome,center=True)
    y = outcome@basis
    null_fit = q0@q0.T@y
    null_residual = m0@y

    def statistic(batch):
        beta = np.einsum('i,bij->bj',a,batch)/aa
        residual = np.einsum('ij,bjk->bik',m,batch)
        sse = np.einsum('bni,bnj->bij',residual,residual)
        inv = np.linalg.pinv(sse,rcond=1e-10)
        return aa*np.einsum('bi,bij,bj->b',beta,inv,beta)

    observed = float(statistic(y[None])[0])
    extreme = 0
    for first in range(0,n_perm,256):
        count = min(256,n_perm-first)
        orders = np.array([rng.permutation(len(y)) for _ in range(count)])
        simulated = null_fit[None]+null_residual[orders]
        extreme += int((statistic(simulated)>=observed-1e-12).sum())
    return {'statistic':observed,'p':(extreme+1)/(n_perm+1),'n':len(y),
            'n_permutations':n_perm,'outcome_rank':basis.shape[1],'nuisance_rank':q0.shape[1],
            'method':'Freedman_Lane_baseline_profile_adjusted'}


def sign_flip_test(values):
    _,basis = orthogonal_columns(values,center=True)
    y = values@basis
    n = len(y)
    yty = y.T@y
    def statistic(means):
        covariance = (yty[None]-n*np.einsum('bi,bj->bij',means,means))/(n-1)
        inv = np.linalg.pinv(covariance,rcond=1e-10)
        return n*np.einsum('bi,bij,bj->b',means,inv,means)
    observed = float(statistic(y.mean(0)[None])[0])
    extreme = 0
    for first in range(0,2**n,1024):
        numbers = np.arange(first,min(first+1024,2**n))
        signs = 2*((numbers[:,None]>>np.arange(n))&1)-1
        extreme += int((statistic(signs@y/n)>=observed-1e-10).sum())
    return {'statistic':observed,'p':extreme/2**n,'n':n,'n_permutations':2**n,
            'outcome_rank':basis.shape[1],'method':'exact_paired_vector_sign_flip_symmetry_test'}


def main():
    out = ROOT/'output/tms_exploration'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True)
    paths = [SOURCE/'participant_scan_metrics.csv',SOURCE/'participant_regional_metrics.csv',
             SOURCE/'participant_paired_changes.csv',SOURCE/'clinical_neural_participant_changes.csv',
             ROOT/'pipeline/atlases/roi_to_network_mapping.csv',ROOT/'data/aachen_itbs/derivatives/README.txt',
             ROOT/'PROVENANCE/TMS_EXPLORATORY_ROUND1_20261125.json']
    hashes = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    scan = pd.read_csv(paths[0]).set_index(KEY)
    regions = pd.read_csv(paths[1])
    original_pairs = pd.read_csv(paths[2])
    atlas = pd.read_csv(paths[4])
    atlas = atlas.loc[atlas.ROI>0].copy()
    atlas['anatomy'] = atlas.Substructure.replace({'vTHA':'THA','dTHA':'THA'})
    regions = regions.merge(atlas[['ROI','anatomy']],left_on='roi',right_on='ROI',validate='many_to_one')
    groups = regions.groupby(KEY+['network'])[['synergy_pc','redundancy_pc','rank_score']].mean().unstack('network')
    groups.columns = ['network_'+measure+'_'+network for measure,network in groups.columns]
    anatomies = ['CAU','PUT','NAc','HIP','AMY','GP','THA']
    sub = regions.loc[regions.anatomy.isin(anatomies)].groupby(KEY+['anatomy']).agg(
        gateways=('gateway','sum'),rank_score=('rank_score','mean')).unstack('anatomy')
    sub.columns = ['anatomy_'+measure+'_'+region for measure,region in sub.columns]
    scan = scan.join(groups,validate='one_to_one').join(sub,validate='one_to_one')
    metrics = [*tms.METRICS,*groups.columns,*sub.columns]
    assert len(metrics)==63 and scan[metrics].notna().all().all()
    scan.reset_index()[[*KEY,'arm',*metrics]].to_csv(out/'scan_metrics_extended.csv',index=False)
    domain = {m:('GWT' if m in tms.METRICS else 'network_continuous' if m.startswith('network_') else 'subcortical_anatomy') for m in metrics}
    pd.DataFrame([{'metric':m,'domain':domain[m]} for m in metrics]).to_csv(out/'metric_inventory.csv',index=False)
    counts = regions.groupby(['dataset','subject_id','session','anatomy']).size().groupby('anatomy').first()
    counts.to_csv(out/'anatomical_workspace_sizes.csv')
    pairs = []
    for row in original_pairs.itertuples():
        a,b = [scan.loc[(row.dataset,row.subject_id,s)] for s in [row.baseline_session,row.followup_session]]
        item = {'dataset':row.dataset,'subject_id':row.subject_id,'arm':row.arm,'comparison':row.comparison,
                'active':int(row.arm=='active'),'delta_fd':row.delta_fd}
        for metric in metrics:
            item[metric+'_baseline'] = a[metric]
            item[metric+'_followup'] = b[metric]
            item[metric+'_delta'] = b[metric]-a[metric]
        pairs.append(item)
    pairs = pd.DataFrame(pairs)
    assert len(pairs)==70 and not pairs.duplicated(['dataset','subject_id','comparison']).any()
    pairs.to_csv(out/'paired_changes_extended.csv',index=False)
    sud = pairs.query("dataset=='sudmex_tms'").copy()
    aachen = pairs.query("dataset=='aachen_itbs'").copy()
    a90 = aachen.query("comparison=='post90_minus_pre90'").set_index('subject_id').sort_index()
    a120 = aachen.query("comparison=='post120_minus_pre120'").set_index('subject_id').sort_index()
    assert list(a90.index)==list(a120.index) and len(a90)==15
    order = {}
    for line in paths[5].read_text().splitlines():
        match = re.match(r'^\s*(\d+)\s+([1342/]+)\s+',line)
        if match: order['aachen_itbs_'+match.group(1).zfill(2)] = int(match.group(2).startswith('3'))
    order_values = np.array([order[s] for s in a90.index],float)
    rows = []
    for metric in metrics:
        result = regression(sud,metric+'_followup',['active',metric+'_baseline'],'active')
        rows.append({'dataset':'sudmex_tms','contrast':'baseline_adjusted_active_sham',
                     'metric':metric,'domain':domain[metric],'reference_only':False,**result})
        d90,d120 = a90[metric+'_delta'],a120[metric+'_delta']
        for name,values in [('90_change',d90),('120_change',d120),('pooled_change',(d90+d120)/2),('intensity_120_minus_90',d120-d90)]:
            rows.append({'dataset':'aachen_itbs','contrast':name,'metric':metric,'domain':domain[metric],
                         'reference_only':name=='pooled_change' and metric in tms.METRICS,**one_sample(values)})
        # The intercept estimates the intensity contrast at equal baseline/FD
        # differences and balanced order. Centering the observed differences
        # would force the intercept back to the raw mean and defeat adjustment.
        covariates = pd.DataFrame({'outcome':d120-d90,'intercept':1.,'order_centered':order_values-.5,
                                   'baseline_difference':a120[metric+'_baseline']-a90[metric+'_baseline'],
                                   'fd_difference':a120.delta_fd-a90.delta_fd})
        result = regression(covariates,'outcome',['intercept','order_centered','baseline_difference','fd_difference'],'intercept')
        rows.append({'dataset':'aachen_itbs','contrast':'intensity_adjusted_order_baseline_motion',
                     'metric':metric,'domain':domain[metric],'reference_only':False,**result})
    neural = pd.DataFrame(rows)
    neural['family'] = neural.dataset+'__'+neural.contrast
    neural['q_block'] = np.nan
    for _,b in neural.groupby('family'):
        good = b.p.notna()
        neural.loc[b.index[good],'q_block'] = multipletests(b.loc[good,'p'],method='fdr_bh')[1]
    neural.to_csv(out/'neural_exploration.csv',index=False)
    print('Computed',len(neural),'neural comparisons including reference rows.',flush=True)

    clinical = pd.read_csv(paths[3]); clinical['active'] = clinical.arm.eq('active').astype(int)
    moderators = ['gateways_FPN','gateways_DAN','mean_synergy_pc','mean_redundancy_pc','global_phi_r']
    moderation = []
    for outcome in ['craving_vas','impulsivity_bis_total']:
        for metric in moderators:
            d = clinical[[outcome+'_delta',outcome+'_t0',metric+'_baseline','active']].dropna().copy()
            d['neural_z'] = (d[metric+'_baseline']-d[metric+'_baseline'].mean())/d[metric+'_baseline'].std()
            d['interaction'] = d.active*d.neural_z
            result = regression(d,outcome+'_delta',['active',outcome+'_t0','neural_z','interaction'],'interaction')
            moderation.append({'clinical_metric':outcome,'neural_metric':metric,'family':'clinical_moderation',**result})
    moderation = pd.DataFrame(moderation)
    moderation['q_block'] = multipletests(moderation.p,method='fdr_bh')[1]
    moderation.to_csv(out/'clinical_baseline_moderation.csv',index=False)

    rng = np.random.default_rng(2026112504)
    omnibus = []
    for measure in ['synergy_pc','redundancy_pc','rank_score']:
        profile = ['network_'+measure+'_'+net for net in tms.NETWORKS]
        y = sud[[m+'_followup' for m in profile]].to_numpy()
        b = sud[[m+'_baseline' for m in profile]].to_numpy()
        result = residual_permutation_test(y,b,sud.active.to_numpy(),rng)
        omnibus.append({'dataset':'sudmex_tms','contrast':'baseline_profile_adjusted','profile':measure,'family':'multivariate_profiles',**result})
        d90 = a90[[m+'_delta' for m in profile]].to_numpy()
        d120 = a120[[m+'_delta' for m in profile]].to_numpy()
        for name,y in [('pooled_change',(d120+d90)/2),('intensity_120_minus_90',d120-d90)]:
            result = sign_flip_test(y)
            omnibus.append({'dataset':'aachen_itbs','contrast':name,'profile':measure,'family':'multivariate_profiles',**result})
        print('Multivariate profile completed:',measure,flush=True)
    omnibus = pd.DataFrame(omnibus)
    omnibus['q_block'] = multipletests(omnibus.p,method='fdr_bh')[1]
    omnibus.to_csv(out/'multivariate_profiles.csv',index=False)

    inventory = []
    for name,table in [('neural',neural),('moderation',moderation),('omnibus',omnibus)]:
        for index,row in table.iterrows():
            if name=='neural' and row.reference_only: continue
            inventory.append({'table':name,'row':index,'p':row.p})
    inventory = pd.DataFrame(inventory)
    inventory['q_all_new_tests'] = np.nan
    good = inventory.p.notna()
    inventory.loc[good,'q_all_new_tests'] = multipletests(inventory.loc[good,'p'],method='fdr_bh')[1]
    inventory.to_csv(out/'complete_new_test_inventory.csv',index=False)
    for name,table,filename in [('neural',neural,'neural_exploration'),('moderation',moderation,'clinical_baseline_moderation'),('omnibus',omnibus,'multivariate_profiles')]:
        q = inventory.loc[inventory.table.eq(name)].set_index('row').q_all_new_tests
        table['q_all_new_tests'] = q.reindex(table.index)
        table.to_csv(out/(filename+'.csv'),index=False)
    for name,digest in hashes.items(): assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    record = {'created_utc':datetime.now(timezone.utc).isoformat(),'status':'exploratory_round_complete',
              'input_sha256':hashes,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'metrics':len(metrics),'neural_table_rows':len(neural),'previously_reported_reference_rows':int(neural.reference_only.sum()),
              'moderation_tests':len(moderation),'multivariate_tests':len(omnibus),
              'new_test_rows':len(inventory),'finite_new_p_values':int(good.sum()),
              'aachen_order_120_first':int(order_values.sum()),'aachen_order_90_first':int((1-order_values).sum()),
              'participants':{'sudmex_tms':40,'aachen_itbs':15},'no_new_exclusions':True}
    (out/'provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    print('OUTPUT:',out)
    print('\nTop neural rows (complete table retained):')
    print(neural.loc[~neural.reference_only].sort_values('p')[['dataset','contrast','metric','estimate','ci_low','ci_high','p','q_block','q_all_new_tests']].head(18).to_string(index=False))
    print('\nMultivariate:\n',omnibus[['dataset','contrast','profile','p','q_block','q_all_new_tests']].to_string(index=False))
    print('\nModeration:\n',moderation[['clinical_metric','neural_metric','estimate','p','q_block','q_all_new_tests']].to_string(index=False))


if __name__=='__main__':
    with threadpool_limits(limits=1): main()
