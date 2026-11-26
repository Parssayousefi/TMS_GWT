"""Audit source mask alignment before interpreting any localized result."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import nibabel as nib
from nibabel.processing import resample_from_to
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/aachen_target_alignment_audit'

def show(ax,data,mask,title):
    center=np.rint(np.argwhere(mask>0).mean(0)).astype(int)
    x=center[0]
    ax.imshow(data[x].T,origin='lower',cmap='gray',vmin=0,vmax=np.percentile(data[data>0],98))
    ax.contour(mask[x].T,levels=[.5],colors=['#ff6535'],linewidths=1.6)
    ax.set_title(title,fontsize=10);ax.axis('off')

def main():
    OUT.mkdir(exist_ok=True)
    rows=[];header=[]
    for subject in ['01','03','04','06','13']:
        stem=f'sub-{subject}_ses-pre90_task-rest'
        func=ROOT/'data/aachen_itbs/derivatives/fmriprep'/f'sub-{subject}'/'ses-pre90/func'
        source=next((ROOT/'data/aachen_target_rois').glob(f'sub-{subject}_ses-1_rFCl*.nii'))
        mask=nib.load(source);m=np.asarray(mask.dataobj)>0
        native=nib.load(func/(stem+'_desc-coreg_boldref.nii.gz'));n=np.asarray(native.dataobj)
        raw=nib.load(ROOT/'data/aachen_itbs'/f'sub-{subject}'/'ses-pre90/func'/(stem+'_bold.nii.gz'))
        std=nib.load(func/(stem+'_space-MNI152NLin2009cAsym_res-2_boldref.nii.gz'))
        mapped=nib.load(ROOT/'output/aachen_target_anatomy'/f'sub-{subject}_pre90_l_mapped.nii.gz')
        target=np.asarray(resample_from_to(mapped,std,order=0).dataobj)
        anatpath=next((ROOT/'data/aachen_itbs/derivatives/fmriprep'/f'sub-{subject}'/'ses-anat/anat').glob('*space-MNI*desc-preproc_T1w.nii.gz'))
        t1=np.asarray(resample_from_to(nib.load(anatpath),std,order=1).dataobj)
        fig,axes=plt.subplots(1,3,figsize=(13,5),layout='constrained')
        show(axes[0],n,m,'Source mask voxels on native BOLD reference')
        show(axes[1],np.asarray(std.dataobj),target,'Mapped mask on saved MNI BOLD reference')
        show(axes[2],t1,target,'Same mask on MNI T1 anatomy')
        fig.suptitle(f'Aachen subject {subject}: left deposited mask; sagittal slice through its centroid')
        fig.savefig(OUT/f'sub-{subject}_native_to_mni.png',dpi=125);plt.close(fig)
        center=np.argwhere(m).mean(0)
        record={'subject':subject,'center_voxel':center.tolist(),'mask_shape':mask.shape,'raw_shape':raw.shape,'native_shape':native.shape}
        for name,img in [('mask',mask),('raw',raw),('native',native)]:
            record[name+'_affine']=img.affine.tolist()
            record[name+'_qform']=img.get_qform().tolist();record[name+'_sform']=img.get_sform().tolist()
            record[name+'_codes']=[int(img.header['qform_code']),int(img.header['sform_code'])]
            record[name+'_world_center']=nib.affines.apply_affine(img.affine,center).tolist()
        header.append(record)
        brainmask=nib.load(func/(stem+'_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz'))
        bm=np.asarray(brainmask.dataobj)>0
        rows.append({'subject':subject,'mask_native_voxels':int(m.sum()),'native_bold_median_in_mask':float(np.median(n[m])),
                     'native_bold_p50_brain':float(np.median(n[n>np.percentile(n,60)])),
                     'fraction_mapped_in_fmriprep_bold_brainmask':float((target*bm).sum()/target.sum())})
    pd.DataFrame(rows).to_csv(OUT/'source_mask_checks.csv',index=False)
    (OUT/'headers.json').write_text(json.dumps(header,indent=2),encoding='utf-8')
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=='__main__':main()
