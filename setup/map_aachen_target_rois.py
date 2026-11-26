"""Map deposited anatomical DLPFC masks to the existing analysis atlas.

No functional outcome is read. Original masks and fMRIPrep outputs are read-only.
The deposited masks have SPM-recentered headers but are resliced on the original
functional voxel grid. Check grid geometry before restoring the raw grid header.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import sys
import numpy as np
import pandas as pd
import nibabel as nib

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tmp/target_roi_deps'))
import SimpleITK as sitk
sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(4)

SESSIONS={1:'pre90',2:'post90',3:'pre120',4:'post120'}
SOURCE=ROOT/'data/aachen_target_rois'
OUT=ROOT/'output/aachen_target_anatomy'


def main():
    OUT.mkdir(exist_ok=True)
    atlas_path=ROOT/'pipeline/atlases/schaefer_400_tian_54_combined.nii.gz'
    atlas=nib.load(atlas_path); labels=np.asarray(atlas.dataobj).astype(int)
    reference=sitk.ReadImage(str(atlas_path))
    files=sorted(SOURCE.glob('*rFC*.nii'))
    inventory=[]
    hashes={}
    def track(p):
        hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
        return p
    track(atlas_path)
    for p in files:
        m=re.search(r'sub-(\d+)_ses-(\d+)_rFC([lr])',p.name)
        inventory.append({'subject':m[1],'session':int(m[2]),'side':m[3],'path':p})
    inventory=pd.DataFrame(inventory)
    summaries=[];weights=[];transforms=[];coverage_rows=[]
    for subject,block in inventory.groupby('subject'):
        if subject=='02':continue
        available=set(block.session)
        chosen=sorted(available&{1,3}) or [min(available)]
        fp=ROOT/'data/aachen_itbs/derivatives/fmriprep'/('sub-'+subject)
        anat=fp/'ses-anat/anat'
        a2m_path=anat/f'sub-{subject}_ses-anat_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5'
        m2a_path=anat/f'sub-{subject}_ses-anat_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5'
        a2m=sitk.ReadTransform(str(track(a2m_path)));m2a=sitk.ReadTransform(str(track(m2a_path)))
        masks={'l':[],'r':[]}
        for session in chosen:
            condition=SESSIONS[session];stem=f'sub-{subject}_ses-{condition}_task-rest'
            func=fp/('ses-'+condition)/'func'
            native_path=func/(stem+'_desc-coreg_boldref.nii.gz')
            native=nib.load(native_path); native_sitk=sitk.ReadImage(str(native_path))
            raw=ROOT/'data/aachen_itbs'/('sub-'+subject)/('ses-'+condition)/'func'/(stem+'_bold.nii.gz')
            raw=nib.load(raw)
            np.testing.assert_allclose(native.affine,raw.affine,atol=1e-4)
            coreg_path=func/(stem+'_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt')
            coreg=sitk.ReadTransform(str(track(coreg_path)))
            composite=sitk.CompositeTransform([coreg,a2m])
            # Validate pull-transform order against the saved standardized BOLD reference.
            standard_path=func/(stem+'_space-MNI152NLin2009cAsym_res-2_boldref.nii.gz')
            standard=sitk.ReadImage(str(standard_path))
            reconstructed=sitk.Resample(native_sitk,standard,composite,sitk.sitkLinear,0,sitk.sitkFloat32)
            x=sitk.GetArrayFromImage(reconstructed).ravel();y=sitk.GetArrayFromImage(standard).ravel()
            sel=(x>np.percentile(x[x>0],15))&(y>np.percentile(y[y>0],15))
            agreement=float(np.corrcoef(x[sel],y[sel])[0,1])
            if agreement<.93:raise ValueError(f'Failed transform order validation: {subject} {condition} {agreement}')
            transforms.append({'subject':subject,'session':condition,'bold_reference_correlation':agreement})
            for side in ['l','r']:
                p=block.loc[block.session.eq(session)&block.side.eq(side),'path'].item()
                image=nib.load(track(p));arr=np.asarray(image.dataobj)
                assert image.shape==native.shape
                corners=np.array(np.meshgrid(*[[-(d-1)/2,(d-1)/2] for d in image.shape])).reshape(3,-1).T
                grid_error=float(np.linalg.norm(corners@(image.affine[:3,:3]-native.affine[:3,:3]).T,axis=1).max())
                assert grid_error<.5, f'Functional voxel-grid discrepancy {subject} {condition} {side}: {grid_error} mm'
                native_mask=sitk.GetImageFromArray((arr>0).astype('uint8').transpose(2,1,0))
                native_mask.CopyInformation(native_sitk)
                mapped=sitk.Resample(native_mask,reference,composite,sitk.sitkNearestNeighbor,0,sitk.sitkUInt8)
                data=sitk.GetArrayFromImage(mapped).transpose(2,1,0).astype(float)
                assert data.sum()>0
                point=nib.affines.apply_affine(native.affine,np.argwhere(arr>0).mean(0))*[-1,-1,1]
                t1point=coreg.GetInverse().TransformPoint(point.tolist())
                mnipoint=np.asarray(m2a.TransformPoint(t1point))*[-1,-1,1]
                back=np.asarray(a2m.TransformPoint((mnipoint*[-1,-1,1]).tolist()))
                residual=float(np.linalg.norm(back-t1point))
                center=nib.affines.apply_affine(atlas.affine,np.argwhere(data>0).mean(0))
                assert residual<.5 and np.linalg.norm(center-mnipoint)<4
                assert (center[0]<0 if side=='l' else center[0]>0)
                masks[side].append(data)
                sitk.WriteImage(mapped,str(OUT/f'sub-{subject}_{condition}_{side}_mapped.nii.gz'))
                summaries.append({'subject':subject,'session':condition,'side':side,
                  'native_voxels':int((arr>0).sum()),'mni_voxels':int(data.sum()),
                  'mni_x':center[0],'mni_y':center[1],'mni_z':center[2],
                  'header_dx':float(image.affine[0,3]-native.affine[0,3]),
                  'header_dy':float(image.affine[1,3]-native.affine[1,3]),
                  'header_dz':float(image.affine[2,3]-native.affine[2,3]),
                  'inverse_transform_error_mm':residual,'centroid_agreement_mm':float(np.linalg.norm(center-mnipoint)),
                  'linear_grid_corner_error_mm':grid_error,
                  'baseline_mask':session in [1,3]})
        for side,data_list in masks.items():
            mean=np.mean(data_list,axis=0)
            cortical=(labels>=1)&(labels<=400)
            coverage=float(mean[cortical].sum()/mean.sum())
            coverage_rows.append({'subject_id':'aachen_itbs_'+subject,'side':side,'cortical_fraction':coverage,'eligible':coverage>.1})
            nib.save(nib.Nifti1Image(mean.astype('float32'),atlas.affine),OUT/f'sub-{subject}_{side}_fixed_mask.nii.gz')
            if coverage<=.1:
                print('Low atlas coverage:',subject,side,coverage,flush=True)
                continue
            counts=np.bincount(labels.ravel(),weights=mean.ravel(),minlength=455)[1:455]
            counts[400:]=0
            counts/=counts.sum()
            for roi,weight in enumerate(counts,1):
                weights.append({'subject_id':'aachen_itbs_'+subject,'side':side,'roi':roi,'weight':weight,
                                'cortical_fraction':coverage,'source_sessions':','.join(map(str,chosen))})
        print('Mapped',subject,'from sessions',chosen,flush=True)
    pd.DataFrame(summaries).to_csv(OUT/'mask_mapping_qc.csv',index=False)
    pd.DataFrame(transforms).to_csv(OUT/'transform_validation.csv',index=False)
    pd.DataFrame(weights).to_csv(OUT/'subject_region_weights.csv',index=False)
    pd.DataFrame(coverage_rows).to_csv(OUT/'parcel_coverage.csv',index=False)
    (OUT/'anatomy_provenance.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),
      'software':'SimpleITK '+sitk.Version_VersionString(),'input_sha256':hashes,
      'selection':'Available baseline masks, averaged within participant; earliest available mask only if no baseline mask exists.',
      'parcel_eligibility':'More than 10 percent of mapped mask mass must overlap the 400 cortical atlas parcels; lower coverage is recorded, not replaced by outcome-selected parcels.',
      'space':'Native functional voxel grids restored after checking equal shape/linear affine; fMRIPrep transforms project to the existing atlas.',
      'limitation':'Mean-functional masks and fMRIPrep reference can differ by residual within-run motion. These are ROI masks, not electric-field models.',
      'functional_outcomes_read':False},indent=2),encoding='utf-8')
    print('Finished:',OUT)


if __name__=='__main__':main()
