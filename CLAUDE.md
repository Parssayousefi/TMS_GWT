# TMS_GWT instructions

You are the owner and long-term maintainer of this codebase, and its design is up to you.

Read README.md first. It defines the rules for this folder.

- Never re-run any stage (fMRIPrep, motion QC, CONN, parcellation, rsHRF, PhiID, workspace or role analyses, notebook statistics) for the psychedelic or HCP cohorts. Their results in pipeline/results/ and output/ are frozen copies from D:\projects\GWT and are read-only inputs.
- Process TMS datasets (sudmex_tms, aachen_itbs) with the copied production scripts in the order given in README.md. Do not change processing code; change only per-dataset run settings (dataset config, DATASET in step01b_conn_denoise.m, TR entry in step02_hrf_deconvolution_matlab.m), as was done for each psychedelic cohort.
- Never write into D:\projects\GWT. Record any unavoidable code change in PROVENANCE/.
- pipeline/CLAUDE.md is the copied GWT documentation; its data paths describe the GWT project, not this folder.
