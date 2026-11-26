"""Populate the skill-scaffolded TMS publication notebook with runnable analyses."""
from pathlib import Path
from textwrap import dedent
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'output/jupyter-notebook/TMS_Global_Workspace_Publication.ipynb'
notebook = nbf.read(DEST, as_version=4)
cells = []


def md(text):
    cells.append(nbf.v4.new_markdown_cell(dedent(text).strip()))


def code(text):
    cells.append(nbf.v4.new_code_cell(dedent(text).strip()))


md('''
# Effects of DLPFC stimulation on the global workspace

**TMS-only publication notebook · analysis specification v0.1 · 25 November 2026**

Does prefrontal stimulation change the brain's synergistic integration and redundant dissemination, and do those changes relate to craving and impulsivity?

This notebook follows the structure and metric definitions of `Global_Workspace_Complete_Pipeline_publication.ipynb`, keeping the sections needed for this paper. It reads the completed TMS ΦID matrices and clinical sheet; it does not launch preprocessing or read psychedelic participant data. The HCP workspace is a **fixed, previously computed reference**.

**Status:** an executable analysis draft, not a registered analysis plan. Output timestamps record when results are first inspected. No claims of preregistration should be made until an actual registration exists. Missing processing outputs are reported as pending, never replaced with simulated results or zeros.

**References:** [Yousefi et al., Psychedelics Restructure the Global Workspace of Consciousness](https://www.biorxiv.org/content/10.64898/2026.09.28.755172v1.full); [Luppi et al., eLife](https://doi.org/10.7554/eLife.88173); [SUDMEX clinical release](https://zenodo.org/records/10409461).
''')
md('''
## Reproduction map

| Part | Publication output |
|---|---|
| I | Analysis settings, fixed workspace, processing status and cohort flow |
| II | Participant-level global ΦID, participation coefficients and regional roles |
| III | SUDMEX active-versus-sham changes; primary FPN gateway result |
| IV | Aachen paired iTBS changes, with intensities averaged within participant |
| V | Craving/impulsivity changes and adjusted brain–clinical associations |
| VI | Motion sensitivity and descriptive open-label follow-up |
| VII | Comparison with published GWT findings; exports and provenance |

Run all cells in order using the **gwt** Python environment. Each execution writes to a new `output/tms_publication/<UTC timestamp>/` directory. A dataset's inferential sections wait for all its planned, motion-QC-eligible primary pairs; other completed datasets can still be analysed.
''')
md('''
## Part I — Settings and planned contrasts

**Proposed primary endpoint:** SUDMEX change in **FPN gateway count**, active versus sham, from T0 to T1. The primary analysis uses a two-sided Welch test of participant change scores, with an unstandardized 95% confidence interval. Baseline- and motion-adjusted HC3 regression is a sensitivity analysis.

**Secondary neural families:** other seven network gateway counts; eight broadcaster counts; two mean workspace PCs; five global ΦID summaries; and two total role counts. Benjamini–Hochberg correction is applied separately within these explicitly listed families. Broadcaster and gateway counts are related outcomes, not independent corroborations. FPN is the only proposed primary endpoint here; the original paper's FPN/DMN testing family is not silently treated as this study's preregistration.

**Aachen:** pair post90 with pre90 and post120 with pre120; average the two changes within each participant before testing. There are 15 independent participants, not 30 independent observations. Intensity-specific changes are descriptive. Without sham, this supports a pre/post observation, not the same causal contrast as SUDMEX.

**Clinical analyses:** exploratory; seven score outcomes form one FDR family. The configured neural measures × seven clinical measures form one association family. Associations adjust for assigned arm, both baseline scores and change in mean FD; clinical change is the outcome. Craving and impulsivity are not direct measures of conscious experience.
''')
code('''
from pathlib import Path
from datetime import datetime, timezone
import sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('module://matplotlib_inline.backend_inline')
import matplotlib.pyplot as plt
from IPython.display import display, Markdown

# Locate TMS_GWT from the notebook's working directory.
ROOT = next((p for p in [Path.cwd().resolve(), *Path.cwd().resolve().parents]
             if (p / 'data/sudmex_tms').is_dir() and (p / 'pipeline/results/phiid_matrices').is_dir()), None)
if ROOT is None:
    raise FileNotFoundError('Open this notebook from inside D:/projects/TMS_GWT.')
sys.path.insert(0, str(ROOT / 'pipeline/scripts/analysis'))
sys.path.insert(0, str(ROOT / 'pipeline/scripts/utils'))
import tms_publication as tms
from clinical_data import METRICS as CLINICAL_METRICS, load_clinical_effects

CONFIG = {
    'specification': 'v0.1', 'preregistration': 'not_registered', 'seed': 20261125,
    'primary_endpoint': tms.PRIMARY, 'alpha': 0.05,
    'clinical_source': 'outputs/tms_clinical_20261125/TMS_clinical_outcomes.xlsx',
    'include_flagged_clinical': False,
    'clinical_neural_metrics': ['gateways_FPN', 'gateways_DAN', 'mean_synergy_pc',
                               'mean_redundancy_pc', 'global_phi_r'],
    'clinical_covariates': ['assigned_arm', 'baseline_clinical', 'baseline_neural', 'delta_fd'],
    'inference_requires_complete_planned_pairs': True,
    'global_change_requires_matching_valid_roi_sets': True,
}
CLINICAL_PATH = ROOT / CONFIG['clinical_source']  # A master CSV may be supplied instead.
STAMP = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
OUT = ROOT / 'output/tms_publication' / STAMP
FIG = OUT / 'figures'
FIG.mkdir(parents=True, exist_ok=True)
INPUTS = [ROOT / 'pipeline/scripts/analysis/tms_publication.py', CLINICAL_PATH]
TABLES = {}
COLORS = {'active': '#6A3D9A', 'sham': '#7F8791', 'iTBS': '#207F89'}
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300, 'figure.dpi': 110})
pd.set_option('display.max_columns', 14)
pd.set_option('display.max_rows', 30)

def table(name, frame, columns=None, limit=25):
    TABLES[name] = frame.copy()
    frame.to_csv(OUT / (name + '.csv'), index=False)
    display((frame[columns] if columns else frame).head(limit).round(5))

def pending(text):
    display(Markdown('**Pending:** ' + text))

def savefig(fig, name):
    fig.savefig(FIG / (name + '.png'), bbox_inches='tight')
    fig.savefig(FIG / (name + '.pdf'), bbox_inches='tight')
    plt.show()
    plt.close(fig)

def forest(results, metrics, title, xlabel):
    frame = results.set_index('metric').reindex(metrics)
    fig, ax = plt.subplots(figsize=(7.2, max(2.6, 0.43*len(metrics))))
    x, low, high = [frame[c].to_numpy(float) for c in ['estimate', 'ci_low', 'ci_high']]
    ok = np.isfinite(x) & np.isfinite(low) & np.isfinite(high)
    ax.errorbar(x[ok], np.flatnonzero(ok), xerr=np.vstack([x[ok]-low[ok], high[ok]-x[ok]]),
                fmt='o', color='#6A3D9A', capsize=3)
    ax.set_yticks(range(len(metrics)))
    ax.set_yticklabels([tms.METRIC_INFO.get(m, (None, m))[1] for m in metrics])
    ax.invert_yaxis(); ax.axvline(0, color='#9AA2AD', lw=1)
    ax.set(xlabel=xlabel, title=title)
    fig.tight_layout()
    return fig

print('Root:', ROOT)
print('Outputs:', OUT.relative_to(ROOT))
display(pd.Series(CONFIG, name='Analysis specification').to_frame())
''')
md('''
### Fixed workspace and exact metric definitions

The 454-region atlas contains Schaefer-400 cortex and Tian-54 subcortex. Participation uses seven Yeo networks plus a single subcortical module. The fixed HCP mask contains 235 regions and is loaded without recomputing any HCP analysis.

- **Global ΦID:** mean of the unique off-diagonal pairs for persistent synergy (`sts`), persistent redundancy (`rtr`), total source synergy, Φ and ΦR. These global summaries are additional TMS outcomes, not a claim that all were primary outcomes in Yousefi et al. The production Gaussian estimator uses natural logarithms (nats).
- **Participation coefficient:** `1 − Σm(kim/ki)²`, calculated over all 454 possible destinations, then averaged over the fixed workspace.
- **Roles:** rank synergy PC and redundancy PC separately **within the workspace**. A positive rank difference is a gateway; negative is a broadcaster; tied ranks are neither. Counts are calculated per participant, not inferred from group-average matrices.
- Missing workspace PCs produce missing role-count outcomes, not artificial zero counts. Global change scores require the same valid regional set at both visits. Coverage is exported for inspection. The workspace size is fixed by design and cannot itself change here.
''')
code('''
reference = tms.load_reference(ROOT)
INPUTS += reference['input_files']
workspace_table = pd.DataFrame([
    {'network': name, 'atlas_regions': int((reference['modules'] == i).sum()),
     'workspace_regions': int(((reference['modules'] == i) & reference['workspace']).sum())}
    for i, name in enumerate(tms.NETWORKS)
])
table('workspace_definition', workspace_table)
assert workspace_table.workspace_regions.sum() == 235
''')
md('''
### Processing status and cohort flow

Existing motion QC defines the expected scans. A ΦID scan is ready only when all five matrices and the completed summary sidecar exist. The notebook checks symmetry, dimensions, missing-edge structure and `ΦR = Φ + redundancy` before analysis. This avoids interpreting a partially written result as a completed scan.

The clinical workbook retains all source visits, but the acute imaging comparison uses only the existing MRI cohort and its eligible T0/T1 pairs. T1-4 is **t14**, the open-label crossover visit; it is never substituted for T1.
''')
code('''
scan_inventory, files = tms.inventory(ROOT)
INPUTS += files
table('scan_inventory', scan_inventory, ['dataset', 'subject_id', 'session', 'qc_status', 'regional_ready', 'hrf_ready', 'phiid_ready'], limit=8)
processing = scan_inventory.groupby('dataset').agg(
    configured_scans=('prefix', 'size'), eligible_scans=('eligible', 'sum'),
    regional_ready=('regional_ready', 'sum'), hrf_ready=('hrf_ready', 'sum'), phiid_ready=('phiid_ready', 'sum')).reset_index()
table('processing_status', processing)
clinical = load_clinical_effects(CLINICAL_PATH, include_flagged=CONFIG['include_flagged_clinical'])
table('clinical_availability', pd.DataFrame([
    {'outcome': m, 'MRI_eligible_pairs': len(clinical), 'usable_changes': int(clinical[m + '_delta'].notna().sum())}
    for m in CLINICAL_METRICS
]))
''')
md('''
## Part II — Participant-level measures and paired changes

Only complete TMS matrices are loaded. The expected pair table is built from motion QC, so an unfinished participant remains visibly pending. There is no inference on an opportunistic subset of whichever scans finish first.
''')
code('''
scans, regions, files = tms.extract_scan_metrics(ROOT, scan_inventory, reference)
INPUTS += files
scans.to_csv(OUT / 'participant_scan_metrics.csv', index=False)
regions.to_csv(OUT / 'participant_regional_metrics.csv', index=False)
pairs = tms.paired_changes(scan_inventory, scans)
pairs.to_csv(OUT / 'participant_paired_changes.csv', index=False)
sudmex_pairs = pairs[pairs.dataset == 'sudmex_tms'].copy()
aachen_pairs = pairs[pairs.dataset == 'aachen_itbs'].copy()
aachen_pooled = tms.pool_aachen(pairs)
READY = {
    'sudmex_tms': not sudmex_pairs.empty and bool(sudmex_pairs.pair_ready.all()),
    'aachen_itbs': not aachen_pooled.empty and bool(aachen_pooled.pair_ready.all()),
}
flow = pairs.groupby(['dataset', 'comparison']).agg(
    planned_pairs=('subject_id', 'size'), completed_pairs=('pair_ready', 'sum'),
    matching_roi_coverage=('same_valid_rois', 'sum')).reset_index()
table('pair_flow', flow)
print('Primary cohorts ready:', READY)
if not scans.empty:
    table('scan_coverage', scans[['dataset', 'subject_id', 'session', 'n_timepoints', 'n_valid_rois', 'n_valid_workspace', 'full_workspace']], limit=8)
''')
md('''
## Part III — SUDMEX: active versus sham

The estimand is **mean(T1−T0)active − mean(T1−T0)sham**. Positive neural effects mean a greater increase under active stimulation. This is the between-arm contrast of within-participant changes, not a comparison of active and sham follow-up values alone.

Unstandardized effects and 95% confidence intervals are the main summaries. The accompanying Hedges' g standardizes the difference between arm-specific change scores; it is **not the paired g_z** used in the earlier crossover study. Their magnitudes should not be equated directly.
''')
code('''
sudmex_results = tms.neural_tests(sudmex_pairs, between_arms=True)
RESULT_COLUMNS = ['metric', 'family', 'estimate', 'ci_low', 'ci_high', 'p', 'q_family', 'effect_size']
if READY['sudmex_tms']:
    table('sudmex_neural_results', sudmex_results, ['n_active', 'n_sham'] + RESULT_COLUMNS)
    table('primary_FPN_result', sudmex_results[sudmex_results.metric == tms.PRIMARY])
else:
    pending('SUDMEX primary and secondary neural estimates require all planned acute MRI pairs.')
''')
code('''
if READY['sudmex_tms']:
    savefig(forest(sudmex_results, ['gateways_' + n for n in tms.NETWORKS],
                   'SUDMEX: network gateway changes', 'Active − sham difference in change (regions; 95% CI)'), 'sudmex_gateways')
    savefig(forest(sudmex_results, ['broadcasters_' + n for n in tms.NETWORKS],
                   'SUDMEX: network broadcaster changes', 'Active − sham difference in change (regions; 95% CI)'), 'sudmex_broadcasters')
    savefig(forest(sudmex_results, ['mean_synergy_pc', 'mean_redundancy_pc'],
                   'SUDMEX: mean workspace participation', 'Active − sham difference in PC change (95% CI)'), 'sudmex_participation')
    savefig(forest(sudmex_results, ['global_' + m for m in tms.MATRIX_NAMES],
                   'SUDMEX: global ΦID', 'Active − sham difference in mean pairwise information change (nats; 95% CI)'), 'sudmex_global_phiid')
''')
md('''
### Regional profiles

This section exports TMS-only participant regional changes and displays average active-minus-sham changes in PC and role rank. It is descriptive: regions are not independent participants, and the rank-score panel is not itself a regional significance map. No psychedelic maps are loaded or correlated with these profiles.
''')
code('''
sudmex_regional_changes = pd.DataFrame()
profile_arms = set(sudmex_pairs.loc[sudmex_pairs.same_valid_rois, 'arm'])
if READY['sudmex_tms'] and profile_arms == {'active', 'sham'}:
    sudmex_regional_changes = tms.regional_changes(regions, sudmex_pairs)
    sudmex_regional_changes.to_csv(OUT / 'sudmex_regional_changes.csv', index=False)
    columns = ['synergy_pc', 'redundancy_pc', 'rank_score']
    means = sudmex_regional_changes.groupby(['arm', 'roi'])[columns].mean()
    contrast = means.loc['active'] - means.loc['sham']
    contrast.reset_index().to_csv(OUT / 'sudmex_regional_effect_profile.csv', index=False)
    fig, axes = plt.subplots(3, 1, figsize=(11, 5.5), sharex=True)
    for ax, name in zip(axes, columns):
        ax.plot(contrast.index, contrast[name], lw=1, color='#6A3D9A')
        ax.axhline(0, color='#9AA2AD', lw=.8); ax.set_ylabel(name.replace('_', ' '))
    axes[0].set_title('SUDMEX regional profiles: active − sham change (descriptive)')
    axes[-1].set_xlabel('Atlas ROI index; fixed workspace regions only')
    fig.tight_layout(); savefig(fig, 'sudmex_regional_profiles')
elif READY['sudmex_tms']:
    pending('Regional comparison needs matching ROI coverage in participants from both arms.')
''')
md('''
## Part IV — Aachen: acute pre/post iTBS changes

The inferential analysis uses each participant's mean of the two paired intensity-specific changes. One-sample t-tests estimate the average change, with 95% confidence intervals and paired Hedges' g_z. FDR families match those above, but this is a separate secondary study. Similar direction across the two datasets would be convergence across different protocols and populations, not an exact replication.
''')
code('''
aachen_results = tms.neural_tests(aachen_pooled, between_arms=False)
if READY['aachen_itbs']:
    table('aachen_neural_results', aachen_results, ['n'] + RESULT_COLUMNS)
    aachen_pooled.to_csv(OUT / 'aachen_participant_pooled_changes.csv', index=False)
    intensity_descriptive = aachen_pairs.groupby('comparison').agg(
        participants=('subject_id', 'nunique'), FPN_change=('gateways_FPN_delta', 'mean'),
        FPN_change_sd=('gateways_FPN_delta', 'std'), synergy_PC_change=('mean_synergy_pc_delta', 'mean'),
        redundancy_PC_change=('mean_redundancy_pc_delta', 'mean')).reset_index()
    table('aachen_intensity_descriptive', intensity_descriptive)
    savefig(forest(aachen_results, ['gateways_' + n for n in tms.NETWORKS],
                   'Aachen: average paired gateway change', 'Mean within-participant change (regions; 95% CI)'), 'aachen_gateways')
    savefig(forest(aachen_results, ['mean_synergy_pc', 'mean_redundancy_pc'],
                   'Aachen: average paired participation change', 'Mean within-participant PC change (95% CI)'), 'aachen_participation')
    savefig(forest(aachen_results, ['global_' + m for m in tms.MATRIX_NAMES],
                   'Aachen: global ΦID', 'Mean within-participant information change (nats; 95% CI)'), 'aachen_global_phiid')
else:
    pending('Aachen estimates require both intensity pairs for every planned participant.')
''')
code('''
if READY['aachen_itbs']:
    values = aachen_pairs.pivot(index='subject_id', columns='comparison', values='gateways_FPN_delta')
    values = values[['post90_minus_pre90', 'post120_minus_pre120']]
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for _, row in values.iterrows():
        ax.plot([0, 1], row, marker='o', color='#207F89', alpha=.45, lw=.8)
    ax.axhline(0, color='#9AA2AD', lw=1)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['90% rMT', '120% rMT'])
    ax.set(ylabel='Post − pre FPN gateway count', title='Aachen intensity visits: one line per participant')
    fig.tight_layout(); savefig(fig, 'aachen_intensity_pairs')
''')
md('''
## Part V — Craving and impulsivity

The workbook supplies craving VAS, CCQ-General, CCQ-Now, BIS-11 total and its three subscales. Original scale units and valid zeros are preserved. Clinical change is T1−T0, so **negative means improvement**. Source-quality flags mask only the affected instrument; no participant-wide deletion or imputation is applied.

Treatment-effect tables concern the imaging cohort. Association models use complete cases per outcome and report their actual sample size. They estimate clinical change per one sample SD of neural change, adjusting for arm, baseline clinical score, baseline neural measure and change in mean FD (OLS with HC3 uncertainty). These are exploratory associations, not mediation, causal symptom mechanisms, or measures of altered consciousness. The five configured neural measures × seven clinical outcomes yield 35 planned tests in one BH family.
''')
code('''
clinical_joined, clinical_results, associations = tms.clinical_analyses(
    ROOT, sudmex_pairs, CLINICAL_PATH, CONFIG['clinical_neural_metrics'],
    include_flagged=CONFIG['include_flagged_clinical'])
if READY['sudmex_tms']:
    clinical_joined.to_csv(OUT / 'clinical_neural_participant_changes.csv', index=False)
    table('clinical_treatment_changes', clinical_results,
          ['metric', 'n_active', 'n_sham', 'estimate', 'ci_low', 'ci_high', 'p', 'q_family'])
    table('clinical_neural_associations', associations,
          ['neural_metric', 'clinical_metric', 'n', 'estimate', 'ci_low', 'ci_high', 'p', 'q_family', 'status'], limit=35)
else:
    pending('Clinical effect and association models will run with the complete SUDMEX imaging cohort.')
''')
code('''
if READY['sudmex_tms']:
    main_scores = ['craving_vas', 'craving_ccq_general', 'craving_ccq_now', 'impulsivity_bis_total']
    rng = np.random.default_rng(CONFIG['seed'])
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    for ax, metric in zip(axes.flat, main_scores):
        for x, arm in enumerate(['sham', 'active']):
            values = clinical_joined.loc[clinical_joined.arm == arm, metric + '_delta'].dropna()
            ax.scatter(x + rng.uniform(-.09, .09, len(values)), values, color=COLORS[arm], alpha=.75, s=24)
            if len(values):
                ax.plot([x-.18, x+.18], [values.mean()]*2, color='black', lw=2)
        ax.axhline(0, color='#9AA2AD', lw=.8)
        ax.set_xticks([0, 1]); ax.set_xticklabels(['Sham', 'Active'])
        ax.set(title=metric.replace('_', ' '), ylabel='T1 − T0 (original units)')
    fig.suptitle('Clinical changes in the imaging cohort; bars show means', y=1.01)
    fig.tight_layout(); savefig(fig, 'clinical_change_distributions')

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    for ax, metric in zip(axes, ['craving_vas', 'impulsivity_bis_total']):
        for arm in ['sham', 'active']:
            block = clinical_joined[clinical_joined.arm == arm]
            ax.scatter(block.gateways_FPN_delta, block[metric + '_delta'], color=COLORS[arm], label=arm, alpha=.75)
        ax.set(xlabel='Change in FPN gateways', ylabel='Change in ' + metric.replace('_', ' '))
        ax.axhline(0, color='#9AA2AD', lw=.7); ax.axvline(0, color='#9AA2AD', lw=.7)
    axes[0].legend(frameon=False)
    fig.suptitle('Unadjusted observations; adjusted inference is reported in the association table')
    fig.tight_layout(); savefig(fig, 'FPN_clinical_scatter')
''')
md('''
## Part VI — Motion sensitivity and later follow-up

The SUDMEX sensitivity model is `neural change ~ active arm + baseline neural measure + change in usable mean FD`, with HC3 uncertainty and the same endpoint families. Mean FD comes from the existing motion-QC records for retained data; QC thresholds and preprocessing are not rerun here. This is a robustness check, not a replacement primary selected after seeing significance.

HRF execution and numerical replay are documented in `logs/tms_pipeline_20261125/hrf_verification/`; the complete 193-scan audit is in `logs/tms_pipeline_20261125/final_validation/`. All outputs passed parameter, array-integrity and unchanged-input checks. Four scans each contain one fitted regional HRF with a 16-second peak (two Aachen and two SUDMEX); these descriptive flags remain recorded without introducing new exclusions. Numerical integrity does not establish perfect recovery of neuronal activity. Additional no-deconvolution or estimator-normalization analyses would need their own validated TMS inputs; they are not claimed to have been completed here.
''')
code('''
motion_results = tms.motion_sensitivity(sudmex_pairs)
if READY['sudmex_tms']:
    table('motion_baseline_adjusted_neural_effects', motion_results,
          ['metric', 'n', 'estimate', 'ci_low', 'ci_high', 'p', 'q_family', 'status'])
    main = sudmex_results[['metric', 'estimate', 'p', 'q_family']].merge(
        motion_results[['metric', 'estimate', 'p', 'q_family']], on='metric', suffixes=('_primary', '_adjusted'))
    table('motion_sensitivity_comparison', main)
else:
    pending('Motion sensitivity needs complete SUDMEX primary pairs.')
''')
md('''
### Later SUDMEX visits — descriptive only

T14 is the sham arm's open-label active-treatment crossover; T2/T3/T4 are later follow-ups. Small, changing samples and open-label exposure prevent interpreting later original-arm differences as the acute randomized contrast. Counts and observed means are reported only once all eligible SUDMEX scans are processed.
''')
code('''
all_sudmex_ready = bool(scan_inventory.loc[(scan_inventory.dataset == 'sudmex_tms') & scan_inventory.eligible, 'phiid_ready'].all())
if all_sudmex_ready:
    followup = scans[scans.dataset == 'sudmex_tms'].groupby(['arm', 'session']).agg(
        participants=('subject_id', 'nunique'), valid_FPN=('gateways_FPN', 'count'),
        mean_FPN_gateways=('gateways_FPN', 'mean'), mean_synergy_PC=('mean_synergy_pc', 'mean'),
        mean_redundancy_PC=('mean_redundancy_pc', 'mean')).reset_index()
    table('open_label_followup_descriptive', followup)
else:
    pending('Full follow-up summary awaits all eligible SUDMEX sessions; the acute contrast does not depend on this summary.')
''')
md('''
## Part VII — Comparison with Yousefi et al.

Compare the direction and anatomical specificity of the TMS findings with the **published** results. A shared FPN gateway effect would support a common network-level signature; it would not establish equivalent subjective states or a shared therapeutic mechanism. A non-significant TMS effect alone does not demonstrate that the interventions differ.

The prior study's drug-profile classification and cross-drug phenomenology/Mantel analyses are not transferred: this paper has a different design and no private psychedelic profiles. Direct between-study effect-size tests or spatial-map correlations would require appropriate public estimates/data and an explicit comparison model. The table below is a qualitative, uncertainty-aware comparison only.
''')
code('''
comparison = tms.comparison_to_paper(sudmex_results, aachen_results, alpha=CONFIG['alpha'])
table('comparison_with_published_GWT', comparison)
''')
md('''
### Export checks and provenance

Each run retains input hashes, settings, package versions, actual sample sizes and processing readiness. The notebook never writes into the original GWT project or the frozen publication-output directory. Empty inferential sections mean inputs are pending, not that the effect is zero.
''')
code('''
assert not pairs.duplicated(tms.KEYS).any()
assert len(aachen_pooled) == aachen_pairs.subject_id.nunique()
assert set(scans.dataset.unique()).issubset(tms.DATASETS)
assert not set(tms.METRICS) & set(CLINICAL_METRICS)
if not scans.empty:
    complete = scans[scans.full_workspace]
    assert np.allclose(complete.gateways_total + complete.broadcasters_total + complete.ties_total, 235)
if READY['sudmex_tms']:
    assert len(associations) == len(CONFIG['clinical_neural_metrics']) * len(CLINICAL_METRICS)
    assert not clinical_joined.duplicated(tms.KEYS).any()

INPUTS += [ROOT / 'pipeline/scripts/utils/clinical_data.py',
           ROOT / 'pipeline/configs/clinical_outcomes.yaml',
           ROOT / 'pipeline/Global_Workspace_Complete_Pipeline_publication.ipynb',
           ROOT / 'output/jupyter-notebook/TMS_Global_Workspace_Publication.ipynb']
provenance = tms.write_provenance(ROOT, OUT, INPUTS, CONFIG)
run_status = {'cohort_readiness': {k: bool(v) for k, v in READY.items()},
              'validated_scans': len(scans), 'planned_primary_pairs': len(pairs),
              'inferential_tables_created': [n for n in TABLES if 'results' in n or 'associations' in n],
              'preregistration': CONFIG['preregistration']}
(OUT / 'run_status.json').write_text(json.dumps(run_status, indent=2), encoding='utf-8')
display(Markdown(f'**Run saved:** `{OUT}`'))
display(pd.Series(run_status, name='Run status').to_frame())
''')

notebook.cells = cells
notebook.metadata.kernelspec = {'display_name': 'Python (gwt)', 'language': 'python', 'name': 'gwt'}
notebook.metadata.language_info = {'name': 'python', 'version': '3.9'}
notebook.metadata['tms_publication'] = {'specification': 'v0.1', 'TMS_only': True,
    'reference_notebook': 'pipeline/Global_Workspace_Complete_Pipeline_publication.ipynb'}
for index, cell in enumerate(notebook.cells):
    cell['id'] = f'tms-publication-{index:02d}'
nbf.validate(notebook)
nbf.write(notebook, DEST)
print(f'Created {DEST}: {len(cells)} cells, {sum(c.cell_type == "code" for c in cells)} code cells.')
