"""Independent numerical checks for the completed TMS exploratory branches."""
from datetime import datetime, timezone
import itertools
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
ROUND1 = ROOT/'output/tms_exploration/20261125T154020_735179Z'
FOLLOWUP = ROOT/'output/tms_exploration/20261125T155441_285503Z'
NUMIT = ROOT/'output/tms_numit/20261125T155252_384543Z'
FWE = sorted(p for p in (ROOT/'output/tms_familywise').iterdir() if p.is_dir())[-1]
sys.path.insert(0,str(ROOT/'pipeline/scripts/analysis'))
sys.path.insert(0,str(ROOT/'setup'))
import tms_numit as tn
import explore_tms_round1 as ex
import analyse_tms_familywise as fw
UNESTIMABLE = []


def close(a,b):
    np.testing.assert_allclose(a,b,rtol=2e-6,atol=2e-9,equal_nan=True)


def hc3(y,x,index):
    """Independent HC3 covariance using the design pseudoinverse."""
    x,y = np.asarray(x,float),np.asarray(y,float)
    pinv = np.linalg.pinv(x)
    beta = pinv@y
    residual = y-x@beta
    hat = np.sum(x*pinv.T,axis=1)
    if np.max(hat)>=1-1e-10:
        raise ArithmeticError('Unit leverage makes HC3 covariance undefined')
    covariance = (pinv*(residual/(1-hat))**2)@pinv.T
    se = np.sqrt(covariance[index,index])
    df = len(y)-x.shape[1]
    half = stats.t.ppf(.975,df)*se
    return [beta[index],beta[index]-half,beta[index]+half,2*stats.t.sf(abs(beta[index]/se),df)]


def verify_contrasts(results,pairs):
    checked = 0
    sud = pairs.loc[pairs.dataset.eq('sudmex_tms')].copy()
    active = sud.arm.eq('active').astype(float).to_numpy()
    a90 = pairs.loc[pairs.comparison.eq('post90_minus_pre90')].set_index('subject_id').sort_index()
    a120 = pairs.loc[pairs.comparison.eq('post120_minus_pre120')].set_index('subject_id').sort_index()
    order = {}
    for line in (ROOT/'data/aachen_itbs/derivatives/README.txt').read_text().splitlines():
        match = re.match(r'^\s*(\d+)\s+([1342/]+)\s+',line)
        if match: order['aachen_itbs_'+match.group(1).zfill(2)] = int(match.group(2).startswith('3'))
    for row in results.itertuples():
        if not np.isfinite(row.p): continue
        m = row.metric
        if row.dataset=='sudmex_tms':
            if row.contrast=='unadjusted_change_active_sham':
                x = sud.loc[sud.arm.eq('active'),m+'_delta'].to_numpy()
                y = sud.loc[sud.arm.eq('sham'),m+'_delta'].to_numpy()
                test = stats.ttest_ind(x,y,equal_var=False)
                interval = test.confidence_interval()
                expected = [x.mean()-y.mean(),interval.low,interval.high,test.pvalue]
            else:
                columns = [np.ones(len(sud)),active,sud[m+'_baseline'].to_numpy()]
                if 'motion' in row.contrast: columns.append(sud.delta_fd.to_numpy())
                expected = hc3(sud[m+'_followup'],np.column_stack(columns),1)
        else:
            d90,d120 = a90[m+'_delta'].to_numpy(),a120[m+'_delta'].to_numpy()
            if row.contrast=='intensity_adjusted_order_baseline_motion':
                covariates = [np.ones(len(d90)),np.array([order[s]-.5 for s in a90.index]),
                              (a120[m+'_baseline']-a90[m+'_baseline']).to_numpy(),
                              (a120.delta_fd-a90.delta_fd).to_numpy()]
                keep = [v for i,v in enumerate(covariates) if i==0 or np.std(v)>1e-13]
                try:
                    expected = hc3(d120-d90,np.column_stack(keep),0)
                except ArithmeticError as error:
                    UNESTIMABLE.append({'dataset':row.dataset,'contrast':row.contrast,'metric':row.metric,
                                        'reference_seed':getattr(row,'reference_seed','raw'),
                                        'original_p':row.p,'reason':str(error)})
                    continue
            else:
                x = {'90_change':d90,'120_change':d120,'pooled_change':(d90+d120)/2,
                     'intensity_120_minus_90':d120-d90}[row.contrast]
                test = stats.ttest_1samp(x,0)
                interval = test.confidence_interval()
                expected = [x.mean(),interval.low,interval.high,test.pvalue]
        close([row.estimate,row.ci_low,row.ci_high,row.p],expected)
        checked += 1
    return checked


def verify_multivariate():
    rng = np.random.default_rng(123)
    x = rng.normal(size=(9,2))+.2
    observed = ex.sign_flip_test(x)
    def statistic(y):
        mean = y.mean(0)
        return len(y)*mean@np.linalg.inv(np.cov(y,rowvar=False))@mean
    target = statistic(x)
    simulated = [statistic(x*np.array(s)[:,None]) for s in itertools.product([-1,1],repeat=len(x))]
    close(observed['statistic'],target)
    close(observed['p'],np.mean(np.array(simulated)>=target-1e-10))
    # Direct one-outcome Freedman-Lane reference with the same permutation RNG.
    y = rng.normal(size=(20,1))
    baseline = rng.normal(size=(20,2))
    arm = np.r_[np.zeros(10),np.ones(10)]
    got = ex.residual_permutation_test(y,baseline,arm,np.random.default_rng(77),n_perm=199)
    x0 = np.column_stack([np.ones(20),baseline])
    design = np.column_stack([x0,arm])
    null_fit = x0@np.linalg.lstsq(x0,y,rcond=None)[0]
    residual = y-null_fit
    aa = np.sum((arm-x0@np.linalg.lstsq(x0,arm,rcond=None)[0])**2)
    def statistic_one(z):
        beta = np.linalg.lstsq(design,z,rcond=None)[0]
        e = z-design@beta
        return float(aa*beta[-1,0]**2/np.sum(e**2))
    target = statistic_one(y)
    permutation_rng = np.random.default_rng(77)
    sims = [statistic_one(null_fit+residual[permutation_rng.permutation(20)]) for _ in range(199)]
    close(got['statistic'],target)
    close(got['p'],(1+np.sum(np.array(sims)>=target-1e-12))/200)


def main():
    raw = pd.read_csv(ROUND1/'neural_exploration.csv')
    raw_pairs = pd.read_csv(ROUND1/'paired_changes_extended.csv')
    count = verify_contrasts(raw,raw_pairs)
    normalized = pd.read_csv(NUMIT/'neural_contrasts_all_seeds.csv',dtype={'reference_seed':str})
    pairs = pd.read_csv(NUMIT/'paired_changes_all_seeds.csv',dtype={'reference_seed':str})
    for seed,frame in normalized.groupby('reference_seed'):
        count += verify_contrasts(frame,pairs.loc[pairs.reference_seed.eq(seed)])
    delayed_results=pd.read_csv(FOLLOWUP/'delayed_active_results.csv')
    delayed_data=pd.read_csv(FOLLOWUP/'delayed_active_participant_changes.csv')
    delayed_checked=0
    for row in delayed_results.itertuples():
        if not np.isfinite(row.p): continue
        frame=delayed_data.loc[delayed_data.contrast.eq(row.contrast)]
        if row.model=='motion_adjusted':
            expected=hc3(frame[row.metric],np.column_stack([np.ones(len(frame)),frame.mean_fd]),0)
        else:
            test=stats.ttest_1samp(frame[row.metric],0); ci=test.confidence_interval()
            expected=[frame[row.metric].mean(),ci.low,ci.high,test.pvalue]
        close([row.estimate,row.ci_low,row.ci_high,row.p],expected)
        delayed_checked+=1
    clinical=pd.read_csv(ROOT/'output/tms_publication/20261125T145928_921789Z/clinical_neural_participant_changes.csv')
    clinical_results=pd.read_csv(FWE/'clinical_familywise.csv')
    for row in clinical_results.itertuples():
        frame=clinical[[row.metric+'_delta',row.metric+'_t0','arm']].dropna()
        expected=hc3(frame[row.metric+'_delta'],np.column_stack([np.ones(len(frame)),frame.arm.eq('active'),frame[row.metric+'_t0']]),1)
        close([row.estimate,row.ci_low,row.ci_high,row.p],expected)
    motion_path=FWE/'aachen_motion_global_familywise.csv'
    if motion_path.exists():
        for row in pd.read_csv(motion_path).itertuples():
            location=ROOT/'output/tms_publication/20261125T145928_921789Z' if row.branch=='raw' else NUMIT
            p=pd.read_csv(location/'participant_paired_changes.csv')
            data=p.loc[p.dataset.eq('aachen_itbs')].groupby('subject_id').mean(numeric_only=True)
            expected=hc3(data[row.metric+'_delta'],np.column_stack([np.ones(len(data)),data.delta_fd]),0)
            close([row.estimate,row.ci_low,row.ci_high,row.p],expected)
    verify_multivariate()
    fw.check_max_t()
    # Numerical reference ordering, chunking, and independent Gaussian solver.
    grid=np.array([.5,1.]); reference=np.array([[0.,1.,1.,3.],[0.,2.,2.,4.]])
    tdmi=np.array([.5,1.,.75,.5,1.]); obs=np.array([1.,2.,1.,-1.,5.])
    q=tn.conditional_cdf(grid,reference,tdmi,obs)
    close(q,[.5,.5,.375,0.,1.])
    perm=np.array([4,2,0,3,1])
    close(tn.conditional_cdf(grid,reference,tdmi[perm],obs[perm]),q[perm])
    close(np.r_[tn.conditional_cdf(grid,reference,tdmi[:2],obs[:2]),tn.conditional_cdf(grid,reference,tdmi[2:],obs[2:])],q)
    rng=np.random.default_rng(601)
    sigma=rng.normal(size=(100,4,4)); sigma=sigma@np.swapaxes(sigma,1,2)+np.eye(4)
    stable=tn.stable_phiid_atoms(sigma); analytic=tn.pg.phiid_atoms(sigma)
    for measure in stable: close(stable[measure],analytic[measure])
    # Report seed stability of the global family using exactly the same max-T.
    seed_results=[]
    metrics=['global_'+m for m in tn.nn.MEASURES]
    for seed,frame in pairs.groupby('reference_seed'):
        pooled=frame.loc[frame.dataset.eq('aachen_itbs')].groupby('subject_id').mean(numeric_only=True)
        raw_p,adjusted=fw.exact_max_t(pooled[[m+'_delta' for m in metrics]].to_numpy())
        for m,p,p_fwe in zip(metrics,raw_p,adjusted):
            seed_results.append({'reference_seed':seed,'metric':m,'estimate':pooled[m+'_delta'].mean(),'p_sign_flip':p,'p_maxT_global5':p_fwe})
    pd.DataFrame(seed_results).to_csv(FWE/'numit_global_family_seed_stability.csv',index=False)
    main_table=pd.read_csv(FWE/'main_comparisons_familywise.csv')
    for _,group in main_table.groupby(['branch','dataset','family']):
        ok=group.p.notna()
        close(group.loc[ok,'p_holm_family'],multipletests(group.loc[ok,'p'],method='holm')[1])
    pd.DataFrame(UNESTIMABLE).to_csv(FWE/'numerically_unestimable_contrasts.csv',index=False)
    record={'status':'passed_for_estimable_contrasts','created_utc':datetime.now(timezone.utc).isoformat(),
            'independently_verified_neural_contrasts':count,
            'independently_verified_delayed_active_contrasts':delayed_checked,
            'independently_verified_clinical_baseline_models':len(clinical_results),
            'unestimable_hc3_rows':UNESTIMABLE,
            'checks':['Independent SciPy t tests and intervals','Independent pseudoinverse HC3 covariance',
                      'Brute-force multivariate sign flips','Independent univariate Freedman-Lane calculation',
                      'Brute-force SciPy max-T check','Holm family corrections','CDF ties, boundaries, ordering and chunk invariance',
                      'Independent LAPACK Gaussian information implementation parity','Normalization global-family seed sensitivity']}
    (FWE/'independent_verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
    print(pd.DataFrame(seed_results).to_string(index=False))


if __name__=='__main__':
    with threadpool_limits(limits=1): main()
