"""Refresh and execute the current results notebook without preprocessing."""
from pathlib import Path
import json
import nbformat as nbf
from nbclient import NotebookClient
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/jupyter-notebook/TMS_Corrected_Analyses.ipynb'
RUN='output/corrected_tms_analysis/20261126T112054_284547Z'
TARGET='output/corrected_aachen_target_flow/20261126T113156_656159Z'
nb=nbf.read(OUT,as_version=4)
md=nbf.v4.new_markdown_cell;code=nbf.v4.new_code_cell
nb.cells=[
    md('# Corrected TMS information and workspace analyses\n\nCurrent results, 26 November 2026. Aachen uses the corrected anatomical branch; historical Aachen treatment effects are superseded. SUDMEX inputs are unchanged. All tests are exploratory and retain their full stated family. No preprocessing or psychedelic calculation runs in this notebook.'),
    code("from pathlib import Path\nimport json\nimport pandas as pd\nfrom IPython.display import display, Image\nROOT=Path(r'D:\\projects\\TMS_GWT')\nRUN=ROOT/'"+RUN+"'\nTARGET=ROOT/'"+TARGET+"'\nassert json.loads((RUN/'status.json').read_text())['normalized']=='complete'\npd.set_option('display.max_rows',100)\nread=lambda name: pd.read_csv(RUN/name)"),
    md('## Whole brain and fixed workspace\nFive information measures per scope; the complete forty-endpoint correction is also visible. Confidence intervals are marginal, not simultaneous. Normalized information values are CDF units; multiply by 100 for percentile points.'),
    code("for branch in ['raw','normalized']:\n    print(branch)\n    d=read(branch+'_aachen_all40.csv')\n    display(d[d.metric.str.startswith(('whole_brain__','within_workspace__'))][['metric','estimate','ci_low','ci_high','p','p_maxT_family','p_maxT_whole_and_workspace10','p_maxT_all']])"),
    md('## Four expectations motivated by Yousefi\nFPN and DAN gateways, redundancy participation, synergy participation. Directional expectation: increase; tests remain two-sided. A coefficient describes network distribution, not the amount of redundant information.'),
    code("for branch in ['raw','normalized']:\n    display(read('gwt_expectations_'+branch+'.csv')[['dataset','metric','n','estimate','ci_low','ci_high','p','p_holm4','p_maxT4']])"),
    md('## All information and role endpoints\nThe full forty-outcome inventory prevents a selected nominal result from standing in for the entire analysis.'),
    code("display(read('raw_aachen_all40.csv')[['metric','family','estimate','p','p_maxT_family','p_maxT_all']])"),
    md('## Spatial scope and robustness'),
    code("display(read('aachen_spatial_interactions.csv'))\ndisplay(read('aachen_all_spatial_motion.csv')[['branch','metric','estimate','ci_low','ci_high','p_holm_scope5','p_holm20']])\ndisplay(read('aachen_intensity_results.csv')[['branch','metric','contrast','estimate','ci_low','ci_high','p_holm_family']])"),
    md('## Signal representation and haemodynamics\nMethod interactions and HRF associations are diagnostics on the same signals; they do not identify a stimulation effect.'),
    code("display(read('sensitivities/estimator_sensitivity_results.csv'))\ndisplay(read('sensitivities/method_interactions.csv'))\ndisplay(read('sensitivities/hrf_phi_associations.csv'))"),
    md('### Whole-brain haemodynamic follow-up\nThree whole-brain tests extend the six original associations; Holm9 covers all nine. Width relationships survive rank sensitivity but are not independent of the deconvolution model. Participant omissions test influence, not replication.'),
    code("display(read('haemodynamic_followup/association_results.csv'))\ndisplay(read('haemodynamic_followup/rank_sensitivity.csv'))\ndisplay(read('haemodynamic_followup/influence_summary.csv'))"),
    md('## Corrected target flow\nEight people with usable bilateral masks. Primary tests, adjusted models and intensity sensitivities remain separate. Fresh MATLAB sentinel validation is unavailable because of license-server error96; independent Python and existing MATLAB-summary checks passed.'),
    code("display(pd.read_csv(TARGET/'inference.csv'))\ndisplay(pd.read_csv(TARGET/'aachen_baseline_motion_sensitivity.csv'))\ndisplay(pd.read_csv(TARGET/'aachen_intensity_descriptive.csv'))"),
    md('## SUDMEX acute randomized comparison and later follow-ups\nThe two-week baseline-adjusted arm contrast estimates the randomized comparison. All-post averaging is a longitudinal sensitivity with later open-label treatment.'),
    code("display(read('raw_sudmex_acute.csv').query(\"metric == 'gateways_FPN'\"))\ndisplay(read('sudmex_all_post_all40_results.csv')[['branch','metric','n','estimate','ci_low','ci_high','p_maxT_family','p_maxT_all_metrics']])"),
    md('## Craving and impulsivity\nComplete seven-by-eight outcome family; source-flagged changes remain masked.'),
    code("display(read('clinical_neural_associations.csv')[['branch','clinical_metric','neural_metric','n','estimate','ci_low','ci_high','p','p_holm56']])"),
    md('## Recovery of the reference workspace\nThe control maps overlap the canonical workspace but are not identical. They do not replace the frozen HCP mask in treatment analyses.'),
    code("display(Image(filename=str(ROOT/'output/control_workspace_qc/combined_20261125T234309_541220Z/control_workspace_comparison.png'),width=1100))"),
    md('## Provenance and interpretation\nThe corrected global Φ estimate and redundancy-participation signal do not survive their relevant families. The study recovers recognizable baseline architecture without establishing the predicted acute workspace reorganization. Full sensitivities and localized model-dependent findings remain in the dashboard and supplement.'),
    code("display(json.loads((RUN/'analysis_specification.json').read_text()))\nprint('Manuscript:',ROOT/'manuscript/aachen_brief_report/corrected_20261126')\nprint('Dashboard:',ROOT/'dashboard/index.html')")]
nb.metadata.kernelspec={'display_name':'Python (gwt)','language':'python','name':'gwt'}
nbf.write(nb,OUT)
NotebookClient(nb,timeout=180,kernel_name='gwt',resources={'metadata':{'path':str(ROOT)}}).execute()
nbf.write(nb,OUT)
receipt={'notebook':str(OUT),'executed_code_cells':sum(c.cell_type=='code' for c in nb.cells),'errors':sum(o.output_type=='error' for c in nb.cells if c.cell_type=='code' for o in c.outputs)}
assert receipt['errors']==0
(OUT.parent/'TMS_Corrected_Analyses.execution.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt))
