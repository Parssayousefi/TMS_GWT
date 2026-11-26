"""Repeat the existing bounded spatial/representation/HRF checks on corrected inputs."""
from pathlib import Path
import argparse
import json
import sys
import types
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'setup'),str(ROOT/'pipeline/scripts/analysis')]
import analyse_corrected_tms as core
import analyse_tms_workspace_scope as ws


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--analysis',required=True,type=Path);args=parser.parse_args()
    analysis=args.analysis.resolve();out=analysis/'sensitivities';out.mkdir(exist_ok=True)
    assert (analysis/'normalization/references/aachen_itbs/2026112512/complete.json').exists()
    assert core.tms.sha256(core.RELEASE)==core.RELEASE_SHA
    originals=[ROOT/'setup/deepen_aachen_spatial_effect.py',ROOT/'setup/verify_aachen_spatial_followup.py']
    def load(original,name):
        source=original.read_text(encoding='utf-8')
        source=source.replace("ROOT=Path(__file__).resolve().parents[1]",'ROOT=Path('+repr(str(ROOT))+')')
        source=source.replace("ROOT/'pipeline/results/timeseries'","CORRECTED/'pipeline/results/timeseries'")
        source=source.replace("ROOT/'pipeline/results/phiid_matrices'","CORRECTED/'pipeline/results/phiid_matrices'")
        if name=='corrected_follow':
            source=source.replace("SOURCE=ROOT/'output/tms_workspace_scope'/(ROOT/'output/tms_workspace_scope/LATEST.txt').read_text().strip()",'SOURCE=Path('+repr(str(analysis))+')')
            source=source.replace("out=ROOT/'output/aachen_spatial_followup'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')",'out=Path('+repr(str(out))+')')
            source=source.replace('out.mkdir(parents=True);','out.mkdir(parents=True,exist_ok=True);')
            source=source.replace("(ROOT/'output/aachen_spatial_followup/LATEST.txt').write_text(out.name+'\\n')",'pass # Historical latest pointer preserved')
            source=source.replace('What supports or limits the apparent outside-workspace Phi increase?','Repeat the original bounded spatial checks after anatomical correction, irrespective of significance.')
            source=source.replace("'normalized_reference_unchanged':True","'normalized_reference_rebuilt_from_corrected_baselines':True")
        else:
            source=source.replace('import deepen_aachen_spatial_effect as follow','import corrected_follow as follow')
            source=source.replace("out=ROOT/'output/aachen_spatial_followup'/(ROOT/'output/aachen_spatial_followup/LATEST.txt').read_text().strip()",'out=Path('+repr(str(out))+')')
        generated=out/(name+'.py');generated.write_text(source,encoding='utf-8')
        module=types.ModuleType(name);module.__file__=str(generated);sys.modules[name]=module
        exec(compile(source,str(generated),'exec'),module.__dict__);module.CORRECTED=core.CORRECTED
        return module,generated
    follow,fpath=load(originals[0],'corrected_follow')
    follow.base=types.SimpleNamespace(**{n:getattr(ws,n) for n in dir(ws) if not n.startswith('__')})
    follow.base.NUMIT=analysis/'normalization'
    follow.main()
    verifier,vpath=load(originals[1],'corrected_follow_verify');verifier.main()
    core.save_json(out/'runtime_routing_provenance.json',{'corrected_root':str(core.CORRECTED),'analysis':str(analysis),'release_sha256':core.RELEASE_SHA,
        'changes':'Explicit corrected signal/matrix/source/null-reference/output paths only; old statistical checks unchanged; old latest pointer preserved.',
        'source_sha256':{str(p.relative_to(ROOT)):core.tms.sha256(p) for p in [Path(__file__),*originals,fpath,vpath]}})


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
