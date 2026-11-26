"""Render anatomical target checks before reading information outcomes."""
from pathlib import Path
import numpy as np
import pandas as pd
import nibabel as nib
from nibabel.processing import resample_from_to
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]

def main():
    out=ROOT/'output/aachen_target_anatomy'
    atlas=nib.load(ROOT/'pipeline/atlases/schaefer_400_tian_54_combined.nii.gz')
    labels=np.asarray(atlas.dataobj)
    qc=pd.read_csv(out/'parcel_coverage.csv')
    subjects=sorted(qc.subject_id.unique())
    for page in range(2):
        fig,axes=plt.subplots(4,2,figsize=(11,13),layout='constrained')
        for row,subject in enumerate(subjects[page*4:(page+1)*4]):
            sid=subject.rsplit('_',1)[1]
            path=next((ROOT/'data/aachen_itbs/derivatives/fmriprep'/f'sub-{sid}'/'ses-anat/anat').glob('*space-MNI*desc-preproc_T1w.nii.gz'))
            t1=np.asarray(resample_from_to(nib.load(path),atlas,order=1).dataobj)
            for col,side in enumerate(['l','r']):
                mask=np.asarray(nib.load(out/f'sub-{sid}_{side}_fixed_mask.nii.gz').dataobj)
                # Sagittal view at the mask's maximum-mass x slice.
                x=int(np.argmax(mask.sum(axis=(1,2))))
                ax=axes[row,col]
                ax.imshow(t1[x].T,origin='lower',cmap='gray',vmin=0,vmax=np.percentile(t1[t1>0],99))
                ax.contour((labels[x]>0).T,levels=[.5],colors=['#58b9df'],linewidths=.45,alpha=.55)
                ax.contour(mask[x].T,levels=[.25],colors=['#ff6d35'],linewidths=1.7)
                r=qc.loc[qc.subject_id.eq(subject)&qc.side.eq(side)].iloc[0]
                ax.set_title(f'{sid} {"LEFT target" if side=="l" else "RIGHT comparator"} | cortical overlap {r.cortical_fraction:.0%}'+(' | excluded' if not r.eligible else ''),fontsize=10)
                ax.set_xlim(10,110);ax.set_ylim(15,92);ax.axis('off')
        fig.suptitle('Deposited targets mapped to existing analysis space\nOrange: mapped mask; blue: atlas boundary. Sagittal slice through each mask.',fontsize=13)
        fig.savefig(out/f'anatomy_qc_{page+1}.png',dpi=130);plt.close(fig)
    sud=ROOT/'output/sudmex_target_anatomy'
    fig,axes=plt.subplots(1,2,figsize=(10,5),layout='constrained')
    for ax,side in zip(axes,['l','r']):
        mask=np.asarray(nib.load(sud/f'common_{side}_target_neighbourhood.nii.gz').dataobj)
        x=int(np.argmax(mask.sum(axis=(1,2))))
        ax.imshow(labels[x].T>0,origin='lower',cmap='gray',alpha=.5)
        ax.contour(mask[x].T,levels=[.5],colors=['#d74b2a'],linewidths=2)
        ax.set_title('SUDMEX '+('left target' if side=='l' else 'right comparator'))
        ax.set_xlim(10,110);ax.set_ylim(15,92);ax.axis('off')
    fig.suptitle('Common 10-mm neighbourhood around the mean of 21 published baseline locations')
    fig.savefig(sud/'anatomy_qc.png',dpi=130);plt.close(fig)

if __name__=='__main__':main()
