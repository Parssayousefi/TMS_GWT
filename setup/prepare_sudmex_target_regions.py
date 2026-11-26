"""Define DLPFC target neighbourhoods from published baseline localization only."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import nibabel as nib
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/sudmex_target_anatomy'


def main():
    OUT.mkdir(exist_ok=True)
    path=ROOT/'data/sudmex_target_rois/14_DLPFC localization information.xlsx'
    coords=pd.read_excel(path,sheet_name='data')
    baseline=coords.loc[coords.Ses.eq('T0')].copy()
    center=baseline[['MNIx','MNIy','MNIz']].mean().to_numpy()
    source=ROOT/'output/tms_workspace_scope'/(ROOT/'output/tms_workspace_scope/LATEST.txt').read_text().strip()
    pairs=pd.read_csv(source/'raw_paired_changes.csv').query("dataset=='sudmex_tms'")
    atlas_path=ROOT/'pipeline/atlases/schaefer_400_tian_54_combined.nii.gz'
    atlas=nib.load(atlas_path);labels=np.asarray(atlas.dataobj).astype(int)
    xyz=nib.affines.apply_affine(atlas.affine,np.indices(labels.shape).reshape(3,-1).T)
    rows=[];qc=[]
    definitions=[('common_baseline_target',subject,center) for subject in pairs.subject_id]
    for r in baseline.itertuples():
        subject='sudmex_tms_'+str(r.ID).zfill(3)
        if subject in set(pairs.subject_id):definitions.append(('individual_baseline_target',subject,np.array([r.MNIx,r.MNIy,r.MNIz])))
    for definition,subject,left in definitions:
        for side in ['l','r']:
            point=left.copy();point[0]=abs(point[0]) if side=='r' else -abs(point[0])
            mask=(np.linalg.norm(xyz-point,axis=1)<=10).reshape(labels.shape)
            count=np.bincount(labels[mask],minlength=455).astype(float)[1:455];count[400:]=0
            coverage=count.sum()/mask.sum()
            qc.append({'definition':definition,'subject_id':subject,'side':side,'x':point[0],'y':point[1],'z':point[2],
                       'radius_mm':10,'cortical_fraction':coverage,'eligible':coverage>.1})
            if coverage<=.1:continue
            count/=count.sum()
            for roi,weight in enumerate(count,1):rows.append({'definition':definition,'subject_id':subject,'side':side,'roi':roi,'weight':weight})
            if definition=='common_baseline_target' and subject==pairs.subject_id.iloc[0]:
                nib.save(nib.Nifti1Image(mask.astype('uint8'),atlas.affine),OUT/f'common_{side}_target_neighbourhood.nii.gz')
    pd.DataFrame(rows).to_csv(OUT/'subject_region_weights.csv',index=False)
    pd.DataFrame(qc).to_csv(OUT/'target_region_qc.csv',index=False)
    baseline.to_csv(OUT/'published_baseline_localizations.csv',index=False)
    (OUT/'anatomy_provenance.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),
      'definition':'10-mm neighbourhood around pooled mean of all published T0 MNI coordinates; mirror x for contralateral comparator.',
      'center_mni_mm':center.tolist(),'n_baseline_localizations':len(baseline),
      'individual_sensitivity':'Only own T0 coordinates; later-coordinate rows are not used.',
      'limitations':'Approximate target neighbourhood in common MNI coordinates; anatomical template variant is not specified by the source dictionary. Not a voxelwise or electric-field analysis.',
      'functional_outcomes_read':False,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
      'atlas_sha256':hashlib.sha256(atlas_path.read_bytes()).hexdigest()},indent=2),encoding='utf-8')
    print('Baseline coordinate count:',len(baseline),'Common center:',center)
    print(pd.DataFrame(qc).groupby(['definition','side']).eligible.agg(['size','sum']).to_string())


if __name__=='__main__':main()
