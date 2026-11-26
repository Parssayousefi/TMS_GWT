"""Isolate an Aachen-only rerun; preserve production code and prior outputs.

Only absolute run paths and the explicit MATLAB DATASET setting are changed
in copied scripts. Numerical processing is byte-for-byte unchanged otherwise.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,shutil
import yaml

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'pipeline'
DEST=ROOT/'reruns/aachen_skullstrip_skip_20261125/pipeline'

def main():
    assert DEST.resolve().is_relative_to(ROOT.resolve())
    if DEST.exists():raise RuntimeError('Stage already exists; verify/reuse it rather than overwrite.')
    DEST.mkdir(parents=True)
    for name in ['scripts','configs','atlases','tools/spm12','tools/rsHRF','tools/npy-matlab','docs/luppi_supps/software/MMI_PhiID']:
        shutil.copytree(SOURCE/name,DEST/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for name in ['results/timeseries','results/phiid_matrices','results/qc','data/conn_input','results/timeseries_intermediate']:
        (DEST/name).mkdir(parents=True,exist_ok=True)
    changes=[]
    names=['scripts/preprocessing/step01b_conn_denoise.m','scripts/preprocessing/step02_hrf_deconvolution_matlab.m','scripts/analysis/step04_compute_phiid_matlab.m']
    for name in names:
        source=(SOURCE/name).read_bytes();replacement=source
        old=str(SOURCE).encode();new=str(DEST).encode()
        assert old in source
        replacement=replacement.replace(old,new)
        if name.endswith('step01b_conn_denoise.m'):
            old_dataset=b"DATASET   = 'missouri_psilo'";new_dataset=b"DATASET   = 'aachen_itbs'"
            assert source.count(old_dataset)==1
            replacement=replacement.replace(old_dataset,new_dataset)
            reconstructed=replacement.replace(new_dataset,old_dataset).replace(new,old)
        else:reconstructed=replacement.replace(new,old)
        assert reconstructed==source,'Unapproved processing change'
        (DEST/name).write_bytes(replacement)
        changes.append({'file':name,'source_sha256':hashlib.sha256(source).hexdigest(),
                        'staged_sha256':hashlib.sha256(replacement).hexdigest(),'change':'Absolute pipeline run directory; CONN DATASET=aachen_itbs only where applicable.'})
    # Point only the Aachen dataset at newly generated anatomy-dependent outputs.
    path=DEST/'configs/datasets/aachen_itbs.yaml'
    config=yaml.safe_load(path.read_text(encoding='utf-8'))
    new_derivatives=(ROOT/'data/aachen_itbs/derivatives/fmriprep_skullstrip_skip_20261125').as_posix()
    config['paths']['local']['fmriprep_dir']=new_derivatives;config['paths']['fmriprep_dir']=new_derivatives
    config['fmriprep_extra_args']='--skull-strip-t1w skip'
    path.write_text(yaml.safe_dump(config,sort_keys=False),encoding='utf-8')
    # Keep registry limited to Aachen to make accidental all-dataset processing impossible.
    registry_path=DEST/'configs/dataset_configs.yaml'
    registry=yaml.safe_load(registry_path.read_text(encoding='utf-8'))
    registry['datasets']={'aachen_itbs':registry['datasets']['aachen_itbs']};registry['default_dataset']='aachen_itbs'
    registry['global_paths']['local']['project_dir']=DEST.as_posix()
    registry['global_paths']['local']['work_dir']=(ROOT/'tmp/fmriprep_aachen_skullstrip_skip_20261125').as_posix()
    registry['global_paths']['local']['log_dir']=(DEST.parent/'logs').as_posix()
    registry_path.write_text(yaml.safe_dump(registry,sort_keys=False),encoding='utf-8')
    unchanged=[]
    for p in (DEST/'scripts').rglob('*'):
        if p.is_file() and p.suffix in ['.py','.m','.ps1','.sh']:
            rel=p.relative_to(DEST).as_posix()
            if rel not in names:
                assert p.read_bytes()==(SOURCE/rel).read_bytes(),rel
                unchanged.append(rel)
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'source':str(SOURCE),'stage':str(DEST),
      'changes':changes,'unchanged_processing_scripts_verified':unchanged,
      'config_changes':'Aachen-only registry; corrected derivative directory; permitted --skull-strip-t1w skip; stage output/log paths.',
      'raw_bids_root':config['paths']['local']['bids_root'],
      'safety':'No original processing file, existing result, frozen HCP/psychedelic file or GWT project file modified. Only Aachen signals can enter the new result folders.',
      'status':'prepared, not yet processed; pilot anatomical validation required'}
    (DEST.parent/'stage_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    (ROOT/'PROVENANCE/AACHEN_ISOLATED_PIPELINE_20261125.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('Prepared isolated pipeline:',DEST)
    print('Byte-identical scripts:',len(unchanged),'path/dataset-only copies:',len(changes))

if __name__=='__main__':main()
