# Fixed-workspace audit and bounded Aachen follow-up

The user requested checking whether the fixed HCP workspace had been applied as in the GWT paper, rerunning affected analyses, and investigating the outside-workspace Phi finding.

The production estimator, preprocessing scripts, atlas and HCP mask were not modified. No HCP or psychedelic stage was run. The original 25 TMS scalar metrics were reproduced from the existing TMS matrices. New mean information summaries restrict both pair endpoints to the existing 235-region HCP mask. PCs and roles retain their original all-destination/within-workspace definitions.

New analysis scripts are `setup/analyse_tms_workspace_scope.py`, `setup/verify_tms_workspace_scope.py`, `setup/deepen_aachen_spatial_effect.py`, and `setup/verify_aachen_spatial_followup.py`. They write separate timestamped outputs under `output/tms_workspace_scope/` and `output/aachen_spatial_followup/`, each with immutable input hashes and verification receipts.

The bounded follow-up computes Gaussian temporal-information sensitivities from TMS signals only: all existing deconvolved samples, first 294 samples of those signals, and the already denoised signals before deconvolution. It does not rerun or overwrite HRF estimation or production PhiID. The independent covariance calculation checks equivalence to saved MATLAB Phi/PhiR in the original branch. All alternative results, direct processing interactions, and numerical checks are retained.

The manuscript was revised to distinguish whole-brain and workspace summaries, include the nonsignificant direct spatial contrast, and report failure to reproduce the increase in denoised BOLD. Preliminary source text is archived in the manuscript folder. The raw and normalized reference analyses remain dependent exploratory analyses, with explicit family-wise correction sets.

The review notebook was executed using the gwt environment with a local Jupyter kernel. Manuscript documents were generated using the bundled python-docx runtime. The packaged document renderer lacked LibreOffice on Windows; a hidden Word instance exported only the project documents read-only to QA PDFs, and bundled Poppler produced the review images. Documents have not been submitted or publicly uploaded.
