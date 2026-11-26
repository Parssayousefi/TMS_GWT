"""Anatomical display of independently derived baseline/reference membership."""
from pathlib import Path
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from nilearn import plotting

ROOT=Path(__file__).resolve().parents[1]

def main():
    base=ROOT/'output/baseline_workspace_qc';out=base/(base/'LATEST_sudmex_tms.txt').read_text().strip()
    d=pd.read_csv(out/'regional_comparison.csv')
    atlas=nib.load(ROOT/'pipeline/atlases/schaefer_400_tian_54_combined.nii.gz');a=np.asarray(atlas.dataobj).astype(int)
    s=d.baseline_workspace.to_numpy(bool);h=d.hcp_workspace.to_numpy(bool)
    values=np.zeros(455,dtype=np.int16)
    values[1:][s&h]=1;values[1:][s&~h]=2;values[1:][~s&h]=3
    img=nib.Nifti1Image(values[a],atlas.affine)
    nib.save(img,out/'baseline_hcp_membership_comparison.nii.gz')
    for name,mask in [('baseline',s),('hcp_reference',h)]:
        v=np.r_[0,mask.astype('uint8')]
        nib.save(nib.Nifti1Image(v[a],atlas.affine),out/f'{name}_workspace.nii.gz')
    bg=ROOT/'tmp/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz'
    colors=['#24799e','#dd8745','#9a67ae']
    fig=plt.figure(figsize=(13,4.1))
    plotting.plot_roi(img,bg_img=bg,cmap=ListedColormap(colors),vmin=1,vmax=3,
      display_mode='z',cut_coords=[-20,0,20,40,60],colorbar=False,alpha=.8,
      figure=fig,axes=[.01,.18,.98,.7],annotate=True,draw_cross=False,black_bg=False)
    fig.suptitle('SUDMEX baseline workspace versus the frozen HCP/GWT reference',fontsize=14,y=.98)
    fig.legend(handles=[Patch(facecolor=colors[0],label='Shared: 179 regions'),
                        Patch(facecolor=colors[1],label='SUDMEX only: 48'),
                        Patch(facecolor=colors[2],label='HCP only: 56')],
               loc='lower center',bbox_to_anchor=(.5,.02),ncol=3,frameon=False)
    fig.savefig(out/'workspace_overlap_brain.png',dpi=170);fig.savefig(out/'workspace_overlap_brain.pdf');plt.close(fig)
    print(out/'workspace_overlap_brain.png')

if __name__=='__main__':main()
