"""Independent numerical and design checks for the fixed-workspace rerun."""
from pathlib import Path
import json
import sys
import unittest

import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'setup'))
import analyse_tms_workspace_scope as run


class SelectionChecks(unittest.TestCase):
    def test_both_endpoints_required_and_diagonal_excluded(self):
        edges,masks=run.edge_scopes([True,False,True,False])
        np.testing.assert_array_equal(edges[masks['within_workspace']],[[0,2]])
        np.testing.assert_array_equal(edges[masks['outside_workspace']],[[1,3]])
        self.assertEqual(masks['between'].sum(),4)
        self.assertTrue((edges[:,0]<edges[:,1]).all())

    def test_known_block_values_and_weighted_identity(self):
        _,masks=run.edge_scopes([True,True,False,False,False])
        values=np.zeros(10)
        for k,s in enumerate(run.SCOPES[1:],1):values[masks[s]]=k
        means=run.scoped_means(values,masks)
        self.assertEqual(means['within_workspace'],1.)
        self.assertEqual(means['outside_workspace'],2.)
        self.assertEqual(means['between'],3.)
        self.assertAlmostEqual(means['whole_brain'],(1+3*2+6*3)/10)

    def test_missing_edges_fail_instead_of_changing_denominator(self):
        _,masks=run.edge_scopes([True,True,False])
        with self.assertRaises(ValueError):run.scoped_means([np.nan,1,2],masks)


def verify(out):
    results=pd.read_csv(out/'aachen_main_results.csv')
    ws=np.load(ROOT/'pipeline/results/group_analysis/hcp/rest/hcp_rest_workspace_mask.npy').astype(bool)
    raw=pd.read_csv(out/'raw_scan_metrics.csv')
    maximum=0.
    # Direct indexing of both axes is independent of the production upper-edge mask.
    for scan in raw.itertuples():
        for m in ['phi','phi_r']:
            a=np.load(ROOT/'pipeline/results/phiid_matrices'/(scan.prefix+'_'+m+'.npy'))
            actual=float(np.nanmean(a[np.ix_(ws,ws)]))
            expected=raw.loc[scan.Index,'within_workspace__'+m]
            maximum=max(maximum,abs(actual-expected))
            np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=1e-12)
    checks=[]
    for branch in ['raw','normalized']:
        scans=pd.read_csv(out/(branch+'_scan_metrics.csv'))
        paired=pd.read_csv(out/(branch+'_paired_changes.csv'))
        pooled=pd.read_csv(out/(branch+'_aachen_pooled_changes.csv')).set_index('subject_id')
        assert len(pooled)==15 and len(paired.loc[paired.dataset.eq('sudmex_tms')])==40
        for m in run.METRICS:
            a=scans.loc[scans.dataset.eq('aachen_itbs')].pivot(index='subject_id',columns='session',values=m).sort_index()
            direct=a[['post90','post120']].mean(1)-a[['pre90','pre120']].mean(1)
            np.testing.assert_allclose(direct,pooled[m+'_delta'],atol=1e-12,rtol=1e-12)
            reconstructed=sum(pooled[s+'__'+m.split('__')[-1]+'_delta']*weight for s,weight in
                              [('within_workspace',27495/102831),('outside_workspace',23871/102831),('between',51465/102831)]) if m.startswith('whole_brain__') else None
            if reconstructed is not None:np.testing.assert_allclose(reconstructed,direct,atol=1e-12,rtol=1e-12)
        cols=['within_workspace__'+m+'_delta' for m in run.MEASURES]
        x=pooled[cols].to_numpy()
        # scipy independently computes every randomized t, rather than using the sum-of-squares shortcut.
        labels=np.arange(2**15)
        signs=2*((labels[:,None]>>np.arange(15))&1)-1
        permuted=stats.ttest_1samp(signs[:,:,None]*x[None,:,:],0,axis=1).statistic
        observed=abs(stats.ttest_1samp(x,0,axis=0).statistic)
        corrected=(abs(permuted).max(1)[:,None]>=observed[None,:]-1e-10).mean(0)
        for k,m in enumerate(run.MEASURES):
            row=results.loc[results.branch.eq(branch)&results.metric.eq('within_workspace__'+m)].iloc[0]
            independent=stats.ttest_1samp(x[:,k],0)
            ci=independent.confidence_interval()
            np.testing.assert_allclose([row.estimate,row.p,row.ci_low,row.ci_high,row.p_maxT_family],
                                       [x[:,k].mean(),independent.pvalue,ci.low,ci.high,corrected[k]],rtol=1e-10,atol=1e-12)
        # The three disjoint graph blocks must reconstruct the whole graph per scan and metric.
        for m in run.MEASURES:
            reconstructed=sum(scans[s+'__'+m]*count/102831 for s,count in
                              [('within_workspace',27495),('outside_workspace',23871),('between',51465)])
            np.testing.assert_allclose(reconstructed,scans['whole_brain__'+m],atol=1e-12,rtol=1e-12)
        checks.append(branch+' pooled means, scipy intervals, exhaustive maxT and graph partition identities passed')
    allpost=pd.read_csv(out/'sudmex_all_post_participant_changes.csv')
    assert allpost.groupby('branch').size().eq(41).all()
    assert allpost.n_pre.eq(1).all() and allpost.n_post.between(1,5).all()
    # Document that fewer effective tests are achieved by dependence, not by dropping measures.
    for branch,block in results.groupby('branch'):
        assert len(block)==30
        assert block.metric.str.startswith('within_workspace__').sum()==5
        good=block.p_maxT_family.notna()
        assert (block.loc[good,'p_maxT_family']>=block.loc[good,'p_sign_flip']-1e-12).all()
        assert (block.loc[good,'p_maxT_all']>=block.loc[good,'p_maxT_family']-1e-12).all()
    audit=json.loads((out/'provenance.json').read_text())
    for relative,digest in audit['input_sha256'].items():assert run.tms.sha256(ROOT/relative)==digest
    record={'status':'passed','direct_workspace_indexing_max_error':maximum,'checks':checks,
            'fixed_mask_regions':int(ws.sum()),'raw_scans':len(raw),'all_inputs_unchanged':True,
            'aachen_independent_participants':15,'sudmex_acute_participants':40,'sudmex_all_followup_participants':41}
    run.save_json(out/'independent_verification.json',record)
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(SelectionChecks)
    if not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():raise SystemExit(1)
    if '--unit-only' not in sys.argv:
        out=ROOT/'output/tms_workspace_scope'/(ROOT/'output/tms_workspace_scope/LATEST.txt').read_text().strip()
        with threadpool_limits(limits=1):verify(out)
