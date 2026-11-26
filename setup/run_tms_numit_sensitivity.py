"""Validated, isolated TMS-only NuMIT-inspired PhiID sensitivity analysis.

This is an exploratory method extension, not a preregistered efficacy analysis.
Run in the gwt environment. Existing inputs and production outputs are read-only.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import solve_discrete_lyapunov
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'output/tms_publication/20261125T145928_921789Z'
sys.path.insert(0, str(ROOT/'pipeline/scripts/analysis'))
sys.path.insert(0, str(ROOT/'setup'))
import tms_publication as tms
import tms_numit as tn
import phiid_gauss as pg
import explore_tms_round1 as exploratory

MEASURES = tms.MATRIX_NAMES
SEEDS = [2026112511, 2026112512]
COUNT = 32768
STEP = .005
KEYS = ['dataset', 'subject_id', 'session']


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str)+'\n', encoding='utf-8')


def load_scans():
    pairs = pd.read_csv(SOURCE/'participant_paired_changes.csv')
    scans = pd.read_csv(SOURCE/'participant_scan_metrics.csv')
    wanted = {(r.dataset, r.subject_id, s) for r in pairs.itertuples()
              for s in [r.baseline_session, r.followup_session]}
    scans = scans[[tuple(row) in wanted for row in scans[KEYS].to_numpy()]]
    scans = scans.sort_values(KEYS).reset_index(drop=True)
    assert len(scans) == 140 and len(pairs) == 70
    assert not scans.duplicated(KEYS).any()
    return scans, pairs


def prepare(out, scans):
    cache = out/'cache'
    cache.mkdir(exist_ok=True)
    if (cache/'prepared.json').exists():
        return json.loads((cache/'prepared.json').read_text())
    edges = np.column_stack(np.triu_indices(454, 1))
    selected = np.random.default_rng(2026112503).choice(len(edges), 1024, replace=False)
    np.save(cache/'edges.npy', edges)
    np.save(cache/'pool_edge_indices.npy', selected)
    pools = {dataset: [] for dataset in tms.DATASETS}
    ranges = {dataset: [np.inf, -np.inf] for dataset in tms.DATASETS}
    hashes, sample_targets = {}, {dataset: [] for dataset in tms.DATASETS}
    for i, scan in enumerate(scans.itertuples()):
        path = ROOT/'pipeline/results/timeseries'/(scan.prefix+'_timeseries_deconv.npy')
        hashes[str(path.relative_to(ROOT))] = tms.sha256(path)
        ts = np.load(path)
        assert ts.shape == (scan.n_timepoints, 454) and np.isfinite(ts).all()
        sigma = pg.lagged_sigma_from_timeseries(ts, 1, edges)
        tdmi = pg.tdmi_from_sigma4(sigma)
        assert np.isfinite(tdmi).all() and (tdmi>0).all()
        np.save(cache/(scan.prefix+'_tdmi.npy'), tdmi)
        ranges[scan.dataset][0] = min(ranges[scan.dataset][0], float(tdmi.min()))
        ranges[scan.dataset][1] = max(ranges[scan.dataset][1], float(tdmi.max()))
        sample_targets[scan.dataset].extend(tdmi[selected].tolist())
        if scan.session in ['t0', 'pre90', 'pre120']:
            a, _, _ = pg.var1_from_sigma(sigma[selected])
            stable = np.max(abs(np.linalg.eigvals(a)), axis=1) < 1
            assert np.isfinite(a).all()
            pools[scan.dataset].append(a[stable])
        atoms = pg.phiid_atoms(sigma[selected[:64]])
        for m in MEASURES:
            path = ROOT/'pipeline/results/phiid_matrices'/(scan.prefix+'_'+m+'.npy')
            hashes[str(path.relative_to(ROOT))] = tms.sha256(path)
            matrix = np.load(path, mmap_mode='r')
            np.testing.assert_allclose(atoms[m], matrix[edges[selected[:64],0],edges[selected[:64],1]], rtol=1e-7, atol=2e-9)
        if (i+1)%20 == 0: print(f'Input preparation: {i+1}/140', flush=True)
    for dataset, pieces in pools.items():
        pool = np.concatenate(pieces)
        np.save(cache/(dataset+'_baseline_A_pool.npy'), pool)
        pools[dataset] = len(pool)
    record = {'input_sha256': hashes, 'tdmi_ranges': ranges, 'baseline_pool_sizes': pools,
              'calibration_targets': {d: np.quantile(v,[.1,.5,.9]).tolist() for d,v in sample_targets.items()}}
    save_json(cache/'prepared.json', record)
    return record


def build_references(out, prepared):
    refs = {}
    for dataset in tms.DATASETS:
        pool = np.load(out/'cache'/(dataset+'_baseline_A_pool.npy'))
        lower, upper = prepared['tdmi_ranges'][dataset]
        # Extra coverage for estimated TDMI in finite-record calibration.
        lower = max(STEP, np.floor((lower-.30)/STEP)*STEP)
        upper = np.ceil((upper+.30)/STEP)*STEP
        grid = np.arange(round(lower/STEP), round(upper/STEP)+1)*STEP
        for seed in SEEDS:
            folder = out/'references'/dataset/str(seed)
            folder.mkdir(parents=True, exist_ok=True)
            paths = {m: folder/(m+'.npy') for m in MEASURES}
            if not (folder/'complete.json').exists():
                np.save(folder/'tdmi_grid.npy', grid)
                arrays = {m: np.lib.format.open_memmap(p, mode='w+', dtype='float32', shape=(len(grid), COUNT)) for m,p in paths.items()}
                records = []
                for i, target in enumerate(grid):
                    draws, qc = tn.valid_reference(pool, target, COUNT, seed)
                    for m in MEASURES: arrays[m][i] = np.sort(draws[m])
                    records.append({'target_tdmi': target, **qc})
                    if (i+1)%50 == 0 or i+1 == len(grid):
                        print(f'Reference {dataset}, seed {seed}: {i+1}/{len(grid)} grid points', flush=True)
                for value in arrays.values(): value.flush()
                pd.DataFrame(records).to_csv(folder/'draw_quality.csv', index=False)
                save_json(folder/'complete.json', {'seed': seed, 'draws_per_target': COUNT, 'grid_step': STEP, 'points': len(grid)})
            refs[(dataset,seed)] = (np.load(folder/'tdmi_grid.npy'), {m: np.load(p,mmap_mode='r') for m,p in paths.items()})
    return refs


def validate_interpolation(out, scans, refs):
    edges = np.load(out/'cache/edges.npy')
    rng = np.random.default_rng(2026112513)
    records = []
    for dataset in tms.DATASETS:
        pool = np.load(out/'cache'/(dataset+'_baseline_A_pool.npy'))
        block = scans.loc[scans.dataset.eq(dataset)]
        # Held-out targets across the actual range and participants, chosen
        # without reference to treatment effects. Includes extreme TDMI pairs.
        queries = []
        for scan in block.iloc[np.linspace(0,len(block)-1,8,dtype=int)].itertuples():
            tdmi = np.load(out/'cache'/(scan.prefix+'_tdmi.npy'))
            indices = [int(tdmi.argmin()), int(tdmi.argmax()), *rng.choice(len(tdmi),2,replace=False)]
            matrices = {m: np.load(ROOT/'pipeline/results/phiid_matrices'/(scan.prefix+'_'+m+'.npy'),mmap_mode='r') for m in MEASURES}
            for k in indices:
                queries.append((float(tdmi[k]), {m: float(v[tuple(edges[k])]) for m,v in matrices.items()}))
        for seed in SEEDS:
            grid, tables = refs[(dataset,seed)]
            for i,(target,observed) in enumerate(queries):
                direct,_ = tn.valid_reference(pool,target,COUNT,seed)
                for m in MEASURES:
                    expected = float(tn.empirical_cdf(np.sort(direct[m]),observed[m]))
                    actual = float(tn.conditional_cdf(grid,tables[m],np.array([target]),np.array([observed[m]]))[0])
                    records.append({'dataset': dataset, 'seed': seed, 'query': i, 'measure': m,
                                    'target_tdmi': target, 'direct': expected, 'interpolated': actual,
                                    'absolute_error': abs(actual-expected)})
            print(f'Grid interpolation checked: {dataset}, seed {seed}',flush=True)
    table = pd.DataFrame(records)
    table.to_csv(out/'interpolation_validation.csv',index=False)
    maximum = float(table.absolute_error.max())
    if maximum > .002:
        raise RuntimeError(f'Grid interpolation error {maximum:.6f} exceeds 0.002; refine for numerical accuracy before effects')
    return maximum


def independent_var(pool, target, rng):
    """Independent scipy root/Lyapunov construction of a draw from the prior."""
    a = pool[rng.integers(len(pool))].copy()
    a /= max(abs(np.linalg.eigvals(a)))
    c11,c22 = np.sqrt(rng.chisquare(2)), np.sqrt(rng.chisquare(1))
    c = np.array([[c11,0.],[rng.normal(),c22]])
    v = c@c.T
    if np.linalg.det(v)<1e-8: v += .05*np.eye(2)
    def objective(radius):
        s = solve_discrete_lyapunov(radius*a,v)
        return .5*(np.linalg.slogdet(s)[1]-np.linalg.slogdet(v)[1])-target
    radius = brentq(objective,1e-6,1-1e-6,xtol=1e-12)
    a *= radius
    s = solve_discrete_lyapunov(a,v)
    sigma = np.block([[s,s@a.T],[a@s,s]])
    independent_tdmi = .5*(2*np.linalg.slogdet(s)[1]-np.linalg.slogdet(sigma)[1])
    assert abs(independent_tdmi-target)<1e-7
    # MI is invariant to diagonal scaling. Normalize before the analytic
    # cofactor helper, exactly as for reference nulls, to avoid cancellation
    # when the Wishart draw has a large dynamic range.
    sigma = pg.normalize_sigma(sigma)
    assert abs(float(tn.stable_phiid_atoms(sigma)['tdmi'])-target)<1e-6
    return a,v,s,sigma


def calibrate_finite_records(out, prepared, refs):
    rng = np.random.default_rng(2026112514)
    records, summaries = [], []
    max_model_error, max_cofactor_error, imprecise_cofactor_draws = 0., 0., 0
    for dataset in tms.DATASETS:
        pool = np.load(out/'cache'/(dataset+'_baseline_A_pool.npy'))
        for target in prepared['calibration_targets'][dataset]:
            for length in [251,295]:
                rows = []
                for replicate in range(128):
                    a,v,s,sigma = independent_var(pool,target,rng)
                    max_model_error = max(max_model_error,float(np.max(abs(pg.stationary_sigma0(a,v)-s))))
                    ts = np.empty((length,2))
                    ts[0] = rng.multivariate_normal(np.zeros(2),s)
                    innovations = rng.multivariate_normal(np.zeros(2),v,size=length-1)
                    for t in range(1,length): ts[t] = a@ts[t-1]+innovations[t-1]
                    empirical = np.corrcoef(np.column_stack([ts[:-1],ts[1:]]),rowvar=False)
                    for mode, covariance in [('population',sigma),('finite_record',empirical)]:
                        atoms = tn.stable_phiid_atoms(covariance)
                        fast_atoms = pg.phiid_atoms(covariance)
                        discrepancy = max(float(abs(atoms[m]-fast_atoms[m])) for m in ['tdmi',*MEASURES])
                        max_cofactor_error = max(max_cofactor_error,discrepancy)
                        imprecise_cofactor_draws += int(discrepancy>1e-6)
                        tdmi = float(atoms['tdmi'])
                        for m in MEASURES:
                            quantiles = []
                            covered = True
                            for seed in SEEDS:
                                grid,tables = refs[(dataset,seed)]
                                if tdmi < grid[0] or tdmi > grid[-1]:
                                    covered = False
                                    break
                                quantiles.append(float(tn.conditional_cdf(grid,tables[m],np.array([tdmi]),np.array([atoms[m]]))[0]))
                            item = {'dataset':dataset,'target_tdmi':target,'length':length,'replicate':replicate,
                                    'mode':mode,'measure':m,'estimated_tdmi':tdmi,'covered':covered,
                                    'quantile':np.mean(quantiles) if covered else np.nan}
                            records.append(item)
                            rows.append(item)
                frame = pd.DataFrame(rows)
                for (mode,m),group in frame.groupby(['mode','measure']):
                    q = group['quantile'].dropna().to_numpy()
                    ks = stats.kstest(q,'uniform')
                    summaries.append({'dataset':dataset,'target_tdmi':target,'length':length,'mode':mode,'measure':m,
                                      'n':len(q),'uncovered':int((~group.covered).sum()),'mean_quantile':q.mean(),
                                      'se_mean':q.std(ddof=1)/np.sqrt(len(q)),'ks_D':ks.statistic,'ks_p_diagnostic_only':ks.pvalue,
                                      'lower_tail_fraction':float((q<.05).mean()),'upper_tail_fraction':float((q>.95).mean())})
                print(f'Independent calibration: {dataset}, TDMI {target:.3f}, T={length}',flush=True)
    pd.DataFrame(records).to_csv(out/'synthetic_calibration_records.csv',index=False)
    table = pd.DataFrame(summaries)
    table.to_csv(out/'synthetic_calibration_summary.csv',index=False)
    return {'max_independent_lyapunov_error':max_model_error,
            'max_analytic_cofactor_error_in_calibration':max_cofactor_error,
            'calibration_covariances_with_cofactor_error_over_1e_6':imprecise_cofactor_draws,
            'max_finite_record_mean_bias':float(abs(table.loc[table['mode'].eq('finite_record'),'mean_quantile']-.5).max()),
            'max_population_mean_bias':float(abs(table.loc[table['mode'].eq('population'),'mean_quantile']-.5).max()),
            'uncovered_rows':int(table.uncovered.sum()),
            'interpretation':'Diagnostic calibration under the specified VAR prior only; not validation that real filtered BOLD is a VAR(1) process.'}


def extract_metrics(out, scans, refs):
    reference = tms.load_reference(ROOT)
    edges = np.load(out/'cache/edges.npy')
    tables, regional, saturation = [], [], []
    for i,scan in enumerate(scans.itertuples()):
        tdmi = np.load(out/'cache'/(scan.prefix+'_tdmi.npy'))
        q = {str(seed): {} for seed in SEEDS}
        q['combined'] = {}
        for m in MEASURES:
            raw = np.load(ROOT/'pipeline/results/phiid_matrices'/(scan.prefix+'_'+m+'.npy'))
            observed = raw[edges[:,0],edges[:,1]]
            for seed in SEEDS:
                grid,values = refs[(scan.dataset,seed)]
                q[str(seed)][m] = tn.conditional_cdf(grid,values[m],tdmi,observed)
            q['combined'][m] = (q[str(SEEDS[0])][m]+q[str(SEEDS[1])][m])/2
            saturation.append({**{k:getattr(scan,k) for k in KEYS},'measure':m,
                               'endpoint_fraction':float(((q['combined'][m]==0)|(q['combined'][m]==1)).mean())})
        for seed,values in q.items():
            matrices = {}
            for m,values_m in values.items():
                matrix = np.full((454,454),np.nan)
                matrix[edges[:,0],edges[:,1]] = matrix[edges[:,1],edges[:,0]] = values_m
                matrices[m] = matrix
            metrics,profile = tms.metric_values(matrices,reference)
            item = {k:getattr(scan,k) for k in [*KEYS,'arm','mean_fd','prefix','n_timepoints','n_valid_rois','valid_roi_signature']}
            tables.append({**item,'reference_seed':seed,**metrics})
            if seed == 'combined':
                for k in [*KEYS,'arm']: profile[k] = getattr(scan,k)
                regional.append(profile)
        if (i+1)%20 == 0: print(f'Full-graph normalized metrics: {i+1}/140',flush=True)
    table = pd.DataFrame(tables)
    table.to_csv(out/'scan_metrics_all_seeds.csv',index=False)
    table.loc[table.reference_seed.eq('combined')].to_csv(out/'participant_scan_metrics.csv',index=False)
    pd.concat(regional,ignore_index=True).to_csv(out/'participant_region_metrics.csv',index=False)
    pd.DataFrame(saturation).to_csv(out/'endpoint_saturation.csv',index=False)
    inventory = scans[KEYS+['arm']].assign(eligible=True)
    pairs = []
    for seed,block in table.groupby('reference_seed'):
        p = tms.paired_changes(inventory,block)
        p['reference_seed'] = seed
        pairs.append(p)
    pairs = pd.concat(pairs,ignore_index=True)
    pairs.to_csv(out/'paired_changes_all_seeds.csv',index=False)
    pairs.loc[pairs.reference_seed.eq('combined')].to_csv(out/'participant_paired_changes.csv',index=False)
    return table,pairs


def monte_carlo_stability(out, table, pairs):
    records = []
    for kind,frame,keys,metrics in [('scans',table,KEYS,tms.METRICS),
                                   ('paired_changes',pairs,['dataset','subject_id','comparison'],[m+'_delta' for m in tms.METRICS])]:
        for dataset in tms.DATASETS:
            a = frame.loc[frame.dataset.eq(dataset)&frame.reference_seed.eq(str(SEEDS[0]))].set_index(keys).sort_index()
            b = frame.loc[frame.dataset.eq(dataset)&frame.reference_seed.eq(str(SEEDS[1]))].set_index(keys).sort_index()
            assert list(a.index)==list(b.index)
            for metric in metrics:
                difference = a[metric]-b[metric]
                scale = ((a[metric]+b[metric])/2).std(ddof=1)
                # SD of one reference estimate = RMS difference / sqrt(2).
                noise = np.sqrt(np.mean(difference**2)/2)
                records.append({'level':kind,'dataset':dataset,'metric':metric,'n':len(a),
                                'rms_seed_difference':np.sqrt(np.mean(difference**2)),
                                'max_absolute_seed_difference':abs(difference).max(),
                                'data_sd':scale,'estimated_single_reference_mc_sd':noise,
                                'mc_sd_over_data_sd':noise/scale if scale>0 else np.nan})
    pd.DataFrame(records).to_csv(out/'monte_carlo_metric_stability.csv',index=False)


def statistical_tables(out,pairs):
    order = {}
    for line in (ROOT/'data/aachen_itbs/derivatives/README.txt').read_text().splitlines():
        match = re.match(r'^\s*(\d+)\s+([1342/]+)\s+',line)
        if match: order['aachen_itbs_'+match.group(1).zfill(2)] = int(match.group(2).startswith('3'))
    records = []
    for seed,p in pairs.groupby('reference_seed'):
        sud = p.query("dataset=='sudmex_tms'").copy()
        sud['active'] = sud.arm.eq('active').astype(int)
        a90 = p.query("dataset=='aachen_itbs' and comparison=='post90_minus_pre90'").set_index('subject_id').sort_index()
        a120 = p.query("dataset=='aachen_itbs' and comparison=='post120_minus_pre120'").set_index('subject_id').sort_index()
        assert list(a90.index)==list(a120.index) and len(a90)==15 and len(sud)==40
        for metric in tms.METRICS:
            for contrast,covariates in [('baseline_adjusted_active_sham',['active',metric+'_baseline']),
                                         ('baseline_motion_adjusted_active_sham',['active',metric+'_baseline','delta_fd'])]:
                result = exploratory.regression(sud,metric+'_followup',covariates,'active')
                records.append({'dataset':'sudmex_tms','contrast':contrast,'metric':metric,'reference_seed':seed,**result})
            result = tms.welch_summary(sud.loc[sud.active.eq(1),metric+'_delta'],sud.loc[sud.active.eq(0),metric+'_delta'])
            records.append({'dataset':'sudmex_tms','contrast':'unadjusted_change_active_sham','metric':metric,'reference_seed':seed,**result})
            d90,d120 = a90[metric+'_delta'],a120[metric+'_delta']
            for contrast,values in [('90_change',d90),('120_change',d120),('pooled_change',(d90+d120)/2),('intensity_120_minus_90',d120-d90)]:
                records.append({'dataset':'aachen_itbs','contrast':contrast,'metric':metric,'reference_seed':seed,**exploratory.one_sample(values)})
            cov = pd.DataFrame({'outcome':d120-d90,'intercept':1.,'order_centered':[order[s]-.5 for s in a90.index],
                                'baseline_difference':a120[metric+'_baseline']-a90[metric+'_baseline'],
                                'fd_difference':a120.delta_fd-a90.delta_fd})
            result = exploratory.regression(cov,'outcome',['intercept','order_centered','baseline_difference','fd_difference'],'intercept')
            records.append({'dataset':'aachen_itbs','contrast':'intensity_adjusted_order_baseline_motion','metric':metric,'reference_seed':seed,**result})
    table = pd.DataFrame(records)
    table['q_block'] = table.groupby(['reference_seed','dataset','contrast']).p.transform(tms.bh)
    table['q_normalized_branch'] = table.groupby('reference_seed').p.transform(tms.bh)
    table['original_metric_family'] = table.metric.map(lambda m:tms.METRIC_INFO[m][0])
    table['q_original_family'] = table.groupby(['reference_seed','dataset','contrast','original_metric_family']).p.transform(tms.bh)
    table.to_csv(out/'neural_contrasts_all_seeds.csv',index=False)
    main = table.loc[table.reference_seed.eq('combined')].copy()
    # Union of new tests in both rounds, in addition to the complete branch.
    previous = pd.read_csv(ROOT/'output/tms_exploration/20261125T154020_735179Z/complete_new_test_inventory.csv')
    inventory = pd.concat([previous[['table','row','p']].assign(branch='raw_exploration'),
                           pd.DataFrame({'table':'normalized_contrasts','row':main.index,'p':main.p,'branch':'normalized'})],ignore_index=True)
    inventory['q_all_exploratory_rounds'] = tms.bh(inventory.p)
    main['q_all_exploratory_rounds'] = inventory.loc[inventory.branch.eq('normalized')].set_index('row').q_all_exploratory_rounds.reindex(main.index)
    main.to_csv(out/'neural_contrasts.csv',index=False)
    inventory.to_csv(out/'all_exploratory_tests.csv',index=False)
    print(main.sort_values('p')[['dataset','contrast','metric','estimate','ci_low','ci_high','p','q_block','q_normalized_branch','q_all_exploratory_rounds']].head(15).to_string(index=False),flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--resume',type=Path)
    args = parser.parse_args()
    out = args.resume.resolve() if args.resume else ROOT/'output/tms_numit'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    if ROOT not in out.parents or 'tms_numit' not in out.parts: raise ValueError('Use isolated TMS output')
    out.mkdir(parents=True,exist_ok=bool(args.resume))
    print('OUTPUT:',out,flush=True)
    start = time.perf_counter()
    scans,_ = load_scans()
    prepared = prepare(out,scans)
    refs = build_references(out,prepared)
    interpolation = validate_interpolation(out,scans,refs)
    calibration = calibrate_finite_records(out,prepared,refs)
    save_json(out/'validation_before_effects.json',{'max_interpolation_error':interpolation,**calibration})
    table,pairs = extract_metrics(out,scans,refs)
    monte_carlo_stability(out,table,pairs)
    statistical_tables(out,pairs)
    for relative,digest in prepared['input_sha256'].items():
        assert tms.sha256(ROOT/relative)==digest
    sources = [Path(__file__),ROOT/'pipeline/scripts/analysis/tms_numit.py',
               ROOT/'pipeline/scripts/analysis/tms_publication.py',ROOT/'setup/explore_tms_round1.py',
               ROOT/'numit-phiid/numit/phiid_gauss.py',ROOT/'numit-phiid/numit/numit_null.py',
               ROOT/'PROVENANCE/TMS_NUMIT_SENSITIVITY_PLAN_20261125.json',
               SOURCE/'participant_scan_metrics.csv',SOURCE/'participant_paired_changes.csv',
               ROOT/'data/aachen_itbs/derivatives/README.txt',*tms.load_reference(ROOT)['input_files']]
    save_json(out/'provenance.json',{'status':'exploratory_sensitivity_complete','created_utc':datetime.now(timezone.utc).isoformat(),
                                  'seeds':SEEDS,'nulls_per_grid_point_per_seed':COUNT,'grid_step_nats':STEP,
                                  'primary_scans':140,'participants':{'sudmex_tms':40,'aachen_itbs':15},
                                  'elapsed_seconds':time.perf_counter()-start,'validation':calibration,
                                  'max_interpolation_error':interpolation,'input_sha256':prepared['input_sha256'],
                                  'source_sha256':{str(p.relative_to(ROOT)):tms.sha256(p) for p in sources},
                                  'normalization_units':'Conditional-null percentiles, not nats; each composite normalized separately',
                                  'limitations':['PhiID extension of a PID method','Baseline-derived null prior treated as fixed; uncertainty in learned prior not propagated',
                                                 'Synthetic calibration tests stationary VAR(1), not filtered/deconvolved BOLD generative validity',
                                                 'Post-hoc exploratory sensitivity; all contrast rows retained; seed analyses are numerical repeats, not independent replications']})
    print('COMPLETE:',out,flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1): main()
