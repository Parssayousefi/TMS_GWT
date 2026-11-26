"""TMS-only temporal-information audit and bounded NuMIT method pilot.

Does not run a normalized treatment-effect analysis or change production files.
Uses pure mathematical helpers only; old cohort drivers are never imported.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import solve_discrete_lyapunov
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'output/tms_publication/20261125T145928_921789Z'
sys.path.insert(0, str(ROOT/'numit-phiid/numit'))
sys.path.insert(0, str(ROOT/'pipeline/scripts/analysis'))
import phiid_gauss as pg
import numit_null as nn
import tms_publication as tms


def percentile(observed, null, tdmi_error):
    valid = np.isfinite(null) & np.isfinite(tdmi_error) & (tdmi_error <= 1e-6)
    denominator = valid.sum(axis=1)
    numerator = ((null < observed[:, None]) & valid).sum(axis=1)
    numerator = numerator + .5*((null == observed[:, None]) & valid).sum(axis=1)
    q = np.divide(numerator, denominator, out=np.full(len(observed), np.nan), where=denominator>0)
    q[~np.isfinite(observed)] = np.nan
    return q, valid.mean()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('diagnostic_directory', type=Path)
    args = parser.parse_args()
    parent = args.diagnostic_directory.resolve()
    if ROOT not in parent.parents or not parent.is_dir():
        raise ValueError('Use an existing TMS_GWT diagnostic run directory')
    out = parent/'numit_feasibility'
    out.mkdir(exist_ok=False)
    seed = 2026112503
    rng = np.random.default_rng(seed)
    pairs = pd.read_csv(SOURCE/'participant_paired_changes.csv')
    scans = pd.read_csv(SOURCE/'participant_scan_metrics.csv')
    keys = set()
    for row in pairs.itertuples():
        for session in [row.baseline_session, row.followup_session]:
            keys.add((row.dataset, row.subject_id, session))
    scans = scans[[tuple(row) in keys for row in scans[['dataset','subject_id','session']].to_numpy()]].copy()
    assert len(scans) == 140 and set(scans.dataset) == {'sudmex_tms', 'aachen_itbs'}
    scans = scans.sort_values(['dataset','subject_id','session'])
    edges = np.column_stack(np.triu_indices(454, 1))
    selected = rng.choice(len(edges), 1024, replace=False)
    parity_idx = selected[:256]
    pilot_idx = selected[:256]
    baseline_sessions = {'sudmex_tms': {'t0'}, 'aachen_itbs': {'pre90','pre120'}}
    pilot_subjects = set(scans.loc[scans.dataset.eq('sudmex_tms')].groupby('arm').subject_id.min())
    pilot_subjects.add(scans.loc[scans.dataset.eq('aachen_itbs'), 'subject_id'].min())
    pools = {name: [] for name in baseline_sessions}
    pilot_scans, rows, input_hashes = [], [], {}
    max_parity = {m: 0. for m in nn.MEASURES}
    max_phi_identity = 0.
    start = time.perf_counter()
    for index, scan in enumerate(scans.itertuples()):
        path = ROOT/'pipeline/results/timeseries'/(scan.prefix+'_timeseries_deconv.npy')
        assert path.name.startswith('dataset-'+scan.dataset+'_')
        input_hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        ts = np.load(path)
        assert ts.shape == (scan.n_timepoints, 454) and np.isfinite(ts).all()
        sigma = pg.lagged_sigma_from_timeseries(ts, 1, edges)
        tdmi = pg.tdmi_from_sigma4(sigma)
        assert np.isfinite(tdmi).all() and (tdmi>0).all()
        obs = pg.phiid_atoms(sigma[parity_idx])
        for measure in nn.MEASURES:
            p = ROOT/'pipeline/results/phiid_matrices'/(scan.prefix+'_'+measure+'.npy')
            production = np.load(p, mmap_mode='r')[edges[parity_idx,0], edges[parity_idx,1]]
            error = float(np.max(np.abs(obs[measure]-production)))
            max_parity[measure] = max(max_parity[measure], error)
            np.testing.assert_allclose(obs[measure], production, rtol=1e-7, atol=2e-9)
        past, future = ts[:-1], ts[1:]
        past = (past-past.mean(0))/past.std(0)
        future = (future-future.mean(0))/future.std(0)
        ac = (past*future).mean(0)
        auto_mi = -.5*np.log1p(-ac**2)
        identity = float(tdmi.mean()-2*auto_mi.mean())
        max_phi_identity = max(max_phi_identity, abs(identity-scan.global_phi))
        np.testing.assert_allclose(identity, scan.global_phi, rtol=1e-7, atol=2e-9)
        a, v, _ = pg.var1_from_sigma(sigma[selected])
        stable = np.max(np.abs(np.linalg.eigvals(a)), axis=1)<1
        assert np.isfinite(a).all() and np.isfinite(v).all()
        rows.append({'dataset':scan.dataset, 'subject_id':scan.subject_id, 'session':scan.session,
                     'arm':scan.arm, 'mean_fd':scan.mean_fd, 'n_timepoints':scan.n_timepoints,
                     'mean_tdmi':tdmi.mean(), 'median_tdmi':np.median(tdmi),
                     'mean_autocorrelation':ac.mean(), 'mean_self_mi':auto_mi.mean(),
                     'sampled_stable_var_fraction':stable.mean(),
                     **{name:getattr(scan,'global_'+name) for name in nn.MEASURES}})
        if scan.session in baseline_sessions[scan.dataset]:
            pools[scan.dataset].append(a[stable])
        if scan.subject_id in pilot_subjects:
            pilot_scans.append({'dataset':scan.dataset, 'subject_id':scan.subject_id,
                                'session':scan.session, 'arm':scan.arm,
                                'tdmi':tdmi[pilot_idx], 'obs':obs})
        if (index+1)%20 == 0:
            print(f'Temporal-information and solver audit: {index+1}/140 scans', flush=True)
    scan_table = pd.DataFrame(rows)
    scan_table.to_csv(out/'temporal_information_scans.csv',index=False)
    scan_table.groupby(['dataset','arm','session'])[['mean_tdmi','mean_autocorrelation','mean_self_mi']].mean().to_csv(out/'temporal_information_arm_means.csv')
    pools = {name:np.concatenate(values) for name,values in pools.items()}
    temporal_pairs = []
    lookup = scan_table.set_index(['dataset','subject_id','session'])
    for row in pairs.itertuples():
        a = lookup.loc[(row.dataset,row.subject_id,row.baseline_session)]
        b = lookup.loc[(row.dataset,row.subject_id,row.followup_session)]
        item = {'dataset':row.dataset,'subject_id':row.subject_id,'arm':row.arm,'comparison':row.comparison}
        for metric in ['mean_tdmi','mean_autocorrelation','mean_self_mi']:
            item[metric+'_baseline'] = a[metric]
            item[metric+'_delta'] = b[metric]-a[metric]
        temporal_pairs.append(item)
    temporal_pairs = pd.DataFrame(temporal_pairs)
    temporal_pairs.to_csv(out/'temporal_information_changes.csv',index=False)
    results = []
    for metric in ['mean_tdmi','mean_autocorrelation','mean_self_mi']:
        sud = temporal_pairs.query("dataset == 'sudmex_tms'")
        summary = tms.welch_summary(sud.loc[sud.arm.eq('active'),metric+'_delta'],sud.loc[sud.arm.eq('sham'),metric+'_delta'])
        results.append({'dataset':'sudmex_tms','metric':metric,**summary})
        pooled = temporal_pairs.query("dataset == 'aachen_itbs'").groupby('subject_id')[metric+'_delta'].mean()
        results.append({'dataset':'aachen_itbs','metric':metric,**tms.paired_summary(pooled)})
    pd.DataFrame(results).to_csv(out/'temporal_information_descriptive_contrasts.csv',index=False)

    # Pure scale invariance of participation/roles on one TMS scan, not NuMIT.
    reference = tms.load_reference(ROOT)
    first = scans.iloc[0]
    matrices = {m:np.load(ROOT/'pipeline/results/phiid_matrices'/(first.prefix+'_'+m+'.npy')) for m in nn.MEASURES}
    before, _ = tms.metric_values(matrices,reference)
    after, _ = tms.metric_values({m:value*(100 if m=='synergy' else .01) for m,value in matrices.items()},reference)
    pc_error = max(abs(before[m]-after[m]) for m in ['mean_synergy_pc','mean_redundancy_pc'])
    assert pc_error < 1e-12
    assert all(before[m]==after[m] for m in before if m.startswith(('gateways_','broadcasters_')))

    # Check the stationary covariance helper independently on random stable VARs.
    covariance_error = 0.
    for _ in range(100):
        a = rng.normal(size=(2,2)); a *= .75/max(abs(np.linalg.eigvals(a)))
        c = rng.normal(size=(2,2)); v = c@c.T + .1*np.eye(2)
        expected = solve_discrete_lyapunov(a,v)
        actual = pg.stationary_sigma0(a,v)
        covariance_error = max(covariance_error,float(np.max(abs(expected-actual))))
        np.testing.assert_allclose(expected,actual,rtol=1e-10,atol=1e-10)

    # Pilot only: same shared baseline pool for both visits/arms within each dataset.
    # The Gaussian-A comparison approximates the default prior with a large fixed pool.
    gaussian_pool = rng.normal(size=(20000,2,2))
    pilot_rows = []
    main_q = {}
    pilot_start = time.perf_counter()
    for scan_index, scan in enumerate(pilot_scans):
        for name, pool in [('shared_baseline_pool',pools[scan['dataset']]), ('gaussian_A_pool',gaussian_pool)]:
            null = nn.null_atoms_for_targets(pool,scan['tdmi'],512,np.random.default_rng(seed+1000+scan_index))
            for measure in nn.MEASURES:
                q, valid_fraction = percentile(scan['obs'][measure],null[measure],null['_tdmi_err'])
                assert np.isfinite(q).all() and ((q>=0)&(q<=1)).all()
                pilot_rows.append({k:scan[k] for k in ['dataset','subject_id','session','arm']} | {
                    'null':name,'measure':measure,'edges':len(q),'nulls_per_edge':512,
                    'mean_quantile':q.mean(),'endpoint_saturation_fraction':((q==0)|(q==1)).mean(),
                    'valid_null_fraction':valid_fraction,'max_tmi_error':np.nanmax(null['_tdmi_err'])})
                if name=='shared_baseline_pool': main_q[(scan_index,measure)] = q
        print(f'NuMIT method pilot: {scan_index+1}/{len(pilot_scans)} scans',flush=True)
    pd.DataFrame(pilot_rows).to_csv(out/'pilot_quantile_quality.csv',index=False)

    repeat_rows = []
    for scan_index,scan in enumerate(pilot_scans):
        null = nn.null_atoms_for_targets(pools[scan['dataset']],scan['tdmi'],512,np.random.default_rng(seed+2000+scan_index))
        for measure in nn.MEASURES:
            q,_ = percentile(scan['obs'][measure],null[measure],null['_tdmi_err'])
            diff = abs(q-main_q[(scan_index,measure)])
            repeat_rows.append({'dataset':scan['dataset'],'subject_id':scan['subject_id'],'session':scan['session'],
                                'measure':measure,'mean_absolute_quantile_difference':diff.mean(),
                                'max_absolute_quantile_difference':diff.max()})
    pd.DataFrame(repeat_rows).to_csv(out/'pilot_monte_carlo_repeatability.csv',index=False)

    calibration = []
    for dataset,pool in pools.items():
        targets = scan_table.loc[scan_table.dataset.eq(dataset),'mean_tdmi'].quantile([.1,.5,.9]).to_numpy()
        for target in targets:
            observations = nn.null_atoms_for_targets(pool,np.array([target]),128,rng)
            reference_null = nn.null_atoms_for_targets(pool,np.full(128,target),512,rng)
            assert np.isfinite(observations['_tdmi_err']).all()
            assert np.max(observations['_tdmi_err']) <= 1e-6
            for measure in nn.MEASURES:
                q, valid = percentile(observations[measure][0],reference_null[measure],reference_null['_tdmi_err'])
                assert np.isfinite(q).all()
                ks = stats.kstest(q,'uniform')
                calibration.append({'dataset':dataset,'target_tdmi':target,'measure':measure,
                                    'observations':128,'reference_nulls':512,'mean_quantile':q.mean(),
                                    'ks_uniform_D':ks.statistic,'ks_p_diagnostic_only':ks.pvalue,
                                    'valid_reference_fraction':valid})
    pd.DataFrame(calibration).to_csv(out/'population_null_calibration.csv',index=False)

    for name,digest in input_hashes.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    helper_paths = [ROOT/'numit-phiid/numit/phiid_gauss.py',ROOT/'numit-phiid/numit/numit_null.py']
    record = {'created_utc':datetime.now(timezone.utc).isoformat(),'status':'pilot_complete',
              'seed':seed,'primary_scans_audited':len(scans),'pilot_scans':len(pilot_scans),
              'pilot_edges_per_scan':256,'pilot_nulls_per_edge':512,
              'max_solver_error_vs_production_by_measure':max_parity,
              'max_global_phi_tdmi_identity_error':max_phi_identity,
              'uniform_weight_rescaling_max_pc_difference':pc_error,
              'uniform_weight_rescaling_role_counts_unchanged':True,
              'max_lyapunov_error_vs_scipy':covariance_error,
              'baseline_pool_sizes':{k:len(v) for k,v in pools.items()},
              'total_elapsed_seconds':time.perf_counter()-start,
              'pilot_elapsed_seconds':time.perf_counter()-pilot_start,
              'limits':['No normalized treatment-effect tests', 'No full-edge gateway stability assessment',
                        'Synthetic calibration is at population covariance level, not finite-sample time series',
                        'PhiID application is an extension of the published PID method'],
              'input_sha256':input_hashes,
              'helper_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in helper_paths}}
    (out/'provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    print('OUTPUT:',out,flush=True)
    print(pd.DataFrame(results)[['dataset','metric','estimate','ci_low','ci_high','p']].to_string(index=False))
    print('Pilot summary:',pd.DataFrame(pilot_rows).groupby(['null','measure'])[['endpoint_saturation_fraction','valid_null_fraction']].mean().to_string())


if __name__=='__main__':
    with threadpool_limits(limits=1):
        main()
