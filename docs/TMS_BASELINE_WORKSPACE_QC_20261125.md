# Baseline workspace quality control

The unchanged GWT workspace algorithm was run on **41 usable sudmex_tms baseline scans (t0)**, using persistent synergy and redundancy, arithmetic group means, node strengths, and the sign of synergy-minus-redundancy ranks. Post-treatment scans were not used. The primary intervention analysis has 40 complete pairs; baseline QC may also include a participant whose follow-up was unavailable.

The independently derived baseline workspace contains **227 regions**, versus **235 in the frozen HCP reference**. They share **179 regions**: Dice **0.775**, Jaccard **0.633**, and gradient Spearman correlation **0.676**. The baseline workspace includes **76.2%** of HCP workspace regions. This is a descriptive comparison, not a test that assumes spatially independent parcels.

The HCP mask copied into this project is byte-identical to the one in the original GWT project, as are the reference gradient and node-strength files. Therefore “our mask” and that GWT HCP mask are the same reference, not two independent replications. The original GWT files were read only; no HCP or psychedelic computation was rerun.

In 2,000 participant bootstraps, HCP Dice ranged from 0.737 to 0.797 between the 2.5th and 97.5th percentiles. Across 1,000 balanced random split halves, median Dice was 0.804 (2.5th–97.5th percentiles 0.774–0.835). These describe sample stability conditional on the processing and fixed HCP reference, not external validation. Clinical population, scan length and acquisition differences can affect similarity.

Aachen baseline QC will be repeated after correcting its anatomical preprocessing. A plausible workspace from the old Aachen signals would not repair or validate the identified brain-extraction error. No settings or thresholds are tuned to maximize reference overlap.

Outputs: `output/baseline_workspace_qc/20261125T195614_096696Z`. The production workspace function and export function were called directly; an independent pandas ranking reproduced the mask exactly. New QC exports are separate from production and frozen reference outputs.
