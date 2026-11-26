# Publication analysis notebook

`Global_Workspace_Complete_Pipeline_publication.ipynb` contains the downstream analyses and figure generation for *How Do Psychoactive Drugs Alter Consciousness? Measuring the Effects of Classic Psychedelics, MDMA, and Amphetamine on the Global Workspace*.

The notebook is intended to be distributed with the associated input data and upstream pipeline sources. Its saved outputs can be read without executing it. The accompanying HTML is a static rendering of the executed notebook.

## Reproduction

Install `requirements_publication.txt` in a compatible Python environment. Run from the `pipeline` directory:

```text
python -m pip install -r requirements_publication.txt
python -m ipykernel install --user --name python3
jupyter nbconvert --to html --execute --ExecutePreprocessor.timeout=-1 Global_Workspace_Complete_Pipeline_publication.ipynb
```

The verified execution environment is recorded in `output/publication-notebook-2026-09-26/execution_environment.json`. First-time surface annotation downloads require internet access. The public Schaefer annotation files used for the checked run are also included under that output directory's `plot_assets/schaefer_annot` folder.

The data root is detected from the working directory or its parents, or can be supplied through the `GWT_ROOT` environment variable. These paths are relative to the data root:

```text
pipeline/
  Global_Workspace_Complete_Pipeline_publication.ipynb
  atlases/
  results/phiid_matrices/
  scripts/
data/
  hcp/derivatives/fmriprep/
  basel_lam/                         # raw questionnaire/study sheets and fMRIPrep confounds
  basel_lpm/                         # raw questionnaire/study sheets and fMRIPrep confounds
output/publication-notebook-2026-09-26/
```

The notebook generates new output files under `output/publication-notebook-2026-09-26/`. It does not overwrite the original notebook, manuscript, or previous verified run. Running it again replaces its own generated outputs. Optional archived group/QC summaries are used only for diagnostic comparisons when present.

## Analysis scope

The notebook starts from existing per-scan ΦID matrices and their retained-volume summaries. It recomputes the workspace, participation coefficients, motion screen and retained-window FD covariates, role assignments, questionnaire scoring, samples, statistical tests, tables, and figures. It does not execute upstream fMRIPrep, CONN denoising, HRF deconvolution, or ΦID estimation. Their source scripts and parameters are documented in Part I.

The analysis includes the six active conditions: LSD (LAM), LSD (LPM), psilocybin, mescaline, MDMA, and amphetamine. Propofol is not included. Exact test families, randomization counts, seeds, and numerical definitions are retained from the verified corrected notebook.

The generated manuscript figures are Figures 2–5 and Supplementary Figures 1–8. Figure 1 is a schematic and is not generated here. Supplementary Figure 1 is drawn after Supplementary Figure 2 so both use the same validated surface annotations and subcortical mapping. Full-precision tables are saved as CSV; figure labels round only for display.

## Figure rendering

For Supplementary Figures 1–2, each cortical hemisphere receives its own 200 Schaefer parcels: left/lh uses ROIs 1–200, right/rh uses ROIs 201–400. Executable assertions verify both aliases and reject unknown hemisphere labels. Grey regions are outside the workspace; positive within-workspace PC rank differences are violet (gateway) and negative differences orange (broadcaster). Tied ranks belong to neither role.

The S1 colorbar has violet/“More Gateway” on the left and orange/“More Broadcaster” on the right. S1 retains the manuscript's continuous cortical/subcortical map and two percentage panels, with VAN used consistently. Its numerical source is `suppfig1_hcp_workspace_by_network.csv`.

S2 retains the six-condition 2 × 3 layout and its existing normalization and surface geometry. The publication version corrects a hemisphere-alias mismatch in the previous rendering. The earlier S2 image must therefore be replaced by the publication output even if it has already been inserted into the manuscript. This change does not alter the regional data or any statistical computation.

S4 and Figure 5 use three-decimal q-values (or < .001) and two-decimal correlations in their annotations. Significance decisions use full-precision values.

The manuscript-editing prompt in the output directory is a separate working artifact; it is not supplementary material for publication.
