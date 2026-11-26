"""Prepare a single cross-dataset report from the released exploratory follow-up."""
from pathlib import Path
import hashlib
import json
import shutil
import zipfile
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / 'output/cross_dataset_spatial_hrf/20261126T121731_876575Z'
OUT = ROOT / 'manuscript/tms_cross_dataset_report/20261126'
OLD = ROOT / 'manuscript/aachen_brief_report/corrected_20261126'
TITLE = 'Haemodynamic covariation of information integration across two prefrontal stimulation datasets'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'figures').mkdir(exist_ok=True)
    (OUT / 'tables').mkdir(exist_ok=True)
    a = pd.read_csv(RUN / 'hrf_association_results.csv')
    spatial = pd.read_csv(RUN / 'spatial_results.csv')
    conjunction = pd.read_csv(RUN / 'hrf_cross_dataset_conjunction.csv')
    observations = pd.read_csv(RUN / 'hrf_observations.csv')
    lookup = lambda dataset, scope, branch='raw', model='motion_arm_adjusted': a.loc[a.dataset.eq(dataset) & a.outcome_scope.eq(scope) & a.predictor_scope.eq(scope) & a.parameter.eq('width_seconds') & a.branch.eq(branch) & a.model.eq(model)].iloc[0]
    aa, su = lookup('aachen_itbs', 'outside_workspace'), lookup('sudmex_tms', 'outside_workspace')
    aan, sun = lookup('aachen_itbs', 'outside_workspace', 'normalized'), lookup('sudmex_tms', 'outside_workspace', 'normalized')
    joint = lambda scope, branch: conjunction.loc[conjunction.branch.eq(branch) & conjunction.outcome_scope.eq(scope) & conjunction.predictor_scope.eq(scope) & conjunction.parameter.eq('width_seconds')].iloc[0].p_holm_family
    assert aa.p_holm15 < .05 and su.p_holm15 < .05 and aan.p_holm15 < .05 and sun.p_holm15 < .05
    names = {'aachen_itbs': 'Aachen', 'sudmex_tms': 'SUDMEX'}
    scopes = {'whole_brain': 'Whole brain', 'within_workspace': 'Within workspace', 'outside_workspace': 'Outside workspace'}
    colors = {'aachen_itbs': '#245a73', 'sudmex_tms': '#a34b3c'}
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2), layout='constrained')
    for ax, dataset in zip(axes, names):
        d = observations.query("branch == 'raw' and dataset == @dataset and outcome_scope == 'outside_workspace' and predictor_scope == 'outside_workspace' and parameter == 'width_seconds'")
        design = np.column_stack([np.ones(len(d)), d.delta_fd])
        if dataset == 'sudmex_tms': design = np.column_stack([design, d.active])
        x = (d.hrf_change.to_numpy() - d.hrf_change.mean()) / d.hrf_change.std(ddof=1)
        y = d.phi_change.to_numpy()
        x = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
        y = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
        for group, marker, label in [(1, 'o', 'Active'), (0, '^', 'Sham')] if dataset == 'sudmex_tms' else [(1, 'o', 'iTBS')]:
            selected = d.active.eq(group).to_numpy() if dataset == 'sudmex_tms' else np.ones(len(d), bool)
            ax.scatter(x[selected], y[selected], marker=marker, s=31, color=colors[dataset], alpha=.78, label=label)
        grid = np.array([x.min(), x.max()]); ax.plot(grid, lookup(dataset, 'outside_workspace').estimate * grid, color=colors[dataset], lw=1.7)
        ax.axhline(0, color='.8', lw=.7); ax.axvline(0, color='.8', lw=.7)
        ax.set_title(f'{names[dataset]} (n={len(d)})'); ax.set_xlabel('HRF width change\n(cohort SD; adjusted residual)'); ax.set_ylabel('Outside-workspace Φ change\n(nats; adjusted residual)'); ax.legend(frameon=False)
    fig.savefig(OUT / 'figures/Figure_1_hrf_width.png', dpi=220); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.3), layout='constrained')
    for ax, branch in zip(axes, ['raw', 'normalized']):
        scale = 1 if branch == 'raw' else 100
        for idx, scope in enumerate(scopes):
            for dataset, offset in [('aachen_itbs', .12), ('sudmex_tms', -.12)]:
                r = lookup(dataset, scope, branch)
                ax.errorbar(r.estimate * scale, 2 - idx + offset, xerr=[[scale * (r.estimate-r.ci_low)], [scale * (r.ci_high-r.estimate)]], color=colors[dataset], fmt='o', capsize=3, label=names[dataset] if idx == 0 else None)
        ax.axvline(0, color='.55', ls='--', lw=1); ax.set_yticks([2, 1, 0], scopes.values()); ax.set_ylim(-.5, 2.5)
        ax.set_xticks([-.012, -.008, -.004, 0] if branch == 'raw' else [-4, -2, 0])
        ax.set_title('Raw information' if branch == 'raw' else 'Conditional-null normalized'); ax.set_xlabel('Φ change per SD of HRF-width change\n' + ('(nats)' if branch == 'raw' else '(percentile points)')); ax.legend(frameon=False)
    fig.savefig(OUT / 'figures/Figure_2_scope_associations.png', dpi=220); plt.close(fig)
    shutil.copyfile(OLD / 'figures/Figure_S1_control_workspace_comparison.png', OUT / 'figures/Figure_S1_control_workspace_comparison.png')
    for f in ['spatial_results.csv', 'hrf_association_results.csv', 'hrf_cross_dataset_conjunction.csv', 'hrf_leave_one_out.csv', 'hrf_influence_summary.csv', 'analysis_specification.json', 'verification.json']:
        shutil.copyfile(RUN / f, OUT / 'tables' / f)
    shutil.copyfile(OLD / 'reference_audit.json', OUT / 'reference_audit.json')
    shutil.copyfile(OLD / 'references.bib', OUT / 'references.bib')
    main_text = f'''# {TITLE}

Article type: Brief Research Report

Keywords: transcranial magnetic stimulation; resting state fMRI; integrated information decomposition; haemodynamic response; global workspace

## Abstract

Information decomposition characterizes brain dynamics from fMRI signals whose temporal structure depends on haemodynamics. We examined whether changes in pairwise information integration differed inside and outside an independent workspace, and whether these changes covaried with the estimated haemodynamic response. Two datasets provided 15 healthy participants undergoing acute intermittent theta burst stimulation and 40 participants with cocaine use disorder in a two-week active-versus-sham trial. Both were analysed with the same parcellation, deconvolution and information estimator. No workspace-versus-outside treatment contrast survived correction across five information measures. Narrower estimated haemodynamic responses accompanied greater whole-brain Φ changes in both datasets. The outside-workspace association survived correction across fifteen haemodynamic associations in both raw and conditional-null-normalized analyses; cross-dataset conjunction p values were {joint('outside_workspace', 'raw'):.3f} and {joint('outside_workspace', 'normalized'):.3f}. Baseline adjustment retained the raw association, but some normalized and participant-omission checks lost corrected support. Neither replicated spatial selectivity nor stimulation specificity was established. These exploratory findings identify haemodynamic covariation in information-integration estimates across two datasets, while leaving its physiological and algorithmic contributions unresolved.

## Introduction

Transcranial magnetic stimulation can affect brain networks beyond its cortical target [@alkhasli2019; @garza2021]. Information decomposition offers a way to describe such distributed dynamics, distinguishing information shared by regions from information available only when they are considered jointly [@mediano2025]. Previous work identified a synergistic workspace [@luppi2024] and changes in its network participation during psychedelic conditions [@yousefi2026]. Whether prefrontal stimulation preferentially changes information processing within or beyond this workspace remains uncertain.

This question also depends on measurement. fMRI reflects haemodynamic responses, and blind deconvolution estimates regional response functions before inferring latent signals [@wu2013]. Information measures applied to those signals may therefore covary with fitted haemodynamic parameters. We examined spatial contrasts and this covariation in two prefrontal-stimulation datasets, retaining separate estimates for their different populations and treatment designs.

## Methods

### Datasets and common processing

Aachen provided resting-state fMRI from 15 healthy adults before and after two left-dorsolateral-prefrontal iTBS visits, at 90% and 120% of resting motor threshold (OpenNeuro ds001832, version 1.0.1) [@alkhasli2019]. Each visit used 600 pulses. One deposited participant was excluded because slice-timing metadata did not match the acquired slices. The 60 retained scans supplied two complete within-visit changes per participant, which were averaged. There was no sham condition.

SUDMEX contributed 40 participants with complete baseline and two-week scans from a randomized trial of 5-Hz left-prefrontal rTMS for cocaine use disorder: 22 assigned active stimulation and 18 sham (OpenNeuro ds003037, version 1.0.1) [@garza2021]. Later open-label visits were excluded from the present treatment comparison. The datasets were not pooled.

Both datasets used the established GWT processing framework [@yousefi2026]: fMRIPrep registration, CONN denoising, 454 Schaefer/Tian regional signals and rsHRF blind deconvolution [@wu2013]. Aachen's already stripped anatomical images were processed without a second brain extraction; corrected alignment and downstream outputs were validated. Gaussian minimum-mutual-information ΦID used a two-second lag. Supplementary Methods detail processing and remaining coverage limitations.

### Spatial information and haemodynamic summaries

Φ was joint past–future mutual information minus the sum of each region's own past–future mutual information [@barrett2011]. It is a pairwise whole-minus-sum measure, can be negative, and is not synergy alone. We averaged it across 102,831 unique whole-brain pairs, 27,495 pairs within the frozen 235-region HCP workspace and 23,871 pairs entirely outside it. The workspace was inherited from the previous pipeline [@yousefi2026] and was not selected from treatment results. Redundancy-corrected ΦR, persistent synergy, persistent redundancy and source synergy completed a five-measure spatial family.

Estimated HRF height, peak latency and full width at half maximum were averaged over the corresponding regions. Within-visit changes were calculated before averaging Aachen's two visits. SUDMEX used the baseline-to-two-week change. A previously computed conditional-null sensitivity placed each information measure relative to a model distribution matched on joint time-delayed mutual information [@liardi2025; @blackburne2026]. These normalized values are conditional percentiles, not additive information atoms.

### Statistical analysis

The direct spatial contrast was the outside-minus-inside difference in information change. Aachen used paired tests with all 32,768 joint sign patterns; SUDMEX estimated the active-minus-sham difference in follow-up spatial contrasts, adjusting for their baseline values and motion change. The five information measures received family-wise correction in each dataset and representation. Separate significant and nonsignificant compartment results were not treated as evidence of a spatial difference.

Exploratory haemodynamic models related Φ change to standardized HRF-parameter change, adjusting for mean framewise-displacement change and, in SUDMEX, treatment arm. HC3 uncertainty and two-sided t tests were used. One complete fifteen-test family comprised three HRF parameters at four scopes (whole brain, inside, outside and their direct difference), plus three tests relating the direct Φ difference to a common whole-brain HRF predictor. Holm correction covered all fifteen associations [@holm1979]. A same-direction conjunction used the larger of the two cohort p values, followed by Holm correction across the fifteen matching endpoints; it did not pool treatment effects. These analyses extend earlier exploratory Aachen results.

Sensitivity models added baseline Φ and HRF values, examined treatment-arm interactions and arm-specific slopes in SUDMEX, repeated raw models after each participant omission, and used rank correlations. Raw and normalized analyses and all nonsignificant outcomes were retained. Confidence intervals are marginal 95% intervals.

## Results

### Spatial contrasts did not establish a preferential treatment response

Control workspaces showed broad correspondence with the reference literature [@luppi2024; @yousefi2026]: Dice overlap was 0.743 in Aachen and 0.776/0.741 in matched SUDMEX sham scans (Supplementary Figure S1). This was descriptive spatial correspondence, not proof of identical masks.

No direct outside-versus-inside treatment contrast survived its five-measure family in either dataset, before or after normalization (Table 2). For raw Φ, Aachen's difference in change was −0.000771 nats (95% CI −0.003290 to 0.001747; maximum-t p=0.935). SUDMEX's baseline-and-motion-adjusted active-minus-sham contrast was −0.000472 nats (−0.007191 to 0.006247; Holm p=1.000). These intervals do not establish equivalence.

### Haemodynamic width covaried with Φ in both datasets

Outside the workspace, each one-SD increase in fitted HRF-width change was associated with a Φ-change difference of {aa.estimate:.6f} nats in Aachen (95% CI {aa.ci_low:.6f} to {aa.ci_high:.6f}; Holm15 p={aa.p_holm15:.6f}) and {su.estimate:.6f} nats in SUDMEX ({su.ci_low:.6f} to {su.ci_high:.6f}; p={su.p_holm15:.4f}). Thus narrowing accompanied larger Φ changes (Figure 1). Whole-brain width associations also survived within-dataset correction (p=0.0149 and 0.0213). Same-direction cross-dataset conjunction p values were {joint('whole_brain', 'raw'):.4f} for whole brain and {joint('outside_workspace', 'raw'):.4f} outside the workspace.

The outside-workspace association remained in the normalized analysis: slopes were {100*aan.estimate:.3f} and {100*sun.estimate:.3f} percentile points per SD in Aachen and SUDMEX, respectively (Holm15 p={aan.p_holm15:.4f} and {sun.p_holm15:.4f}; conjunction p={joint('outside_workspace', 'normalized'):.4f}; Figure 2). Whole-brain normalized associations did not survive correction in both datasets. Within-workspace associations were not supported in both datasets, and direct spatial tests did not establish that the width relationship was consistently stronger outside the workspace.

Baseline adjustment retained the raw outside-workspace association in both cohorts (Holm15 p=0.0271 and 0.0309), but the normalized Aachen association lost corrected support. Raw width slopes remained negative in every participant omission; corrected outside-workspace support persisted in 14/15 Aachen and 26/40 SUDMEX omissions. In SUDMEX, active and sham slopes were both negative, but neither arm-specific raw association survived its full family. The arm-by-width interaction was not supported (uncorrected p=0.338; Holm15 p=1.000). This does not prove equal slopes or a treatment-independent mechanism.

## Discussion

The shared finding across these datasets was covariation between fitted haemodynamic response width and estimated information integration. Narrower responses accompanied larger whole-brain Φ changes, and the outside-workspace association was also present after conditional-null normalization. This provides a common measurement-focused result across healthy participants undergoing acute iTBS and a clinical active-versus-sham rTMS trial.

The spatial findings require a narrower interpretation than preferential modulation outside a global workspace. Neither dataset established a direct outside-versus-inside treatment contrast, and the haemodynamic slope difference did not generalize across cohorts. A supported outside-workspace association alongside weaker within-workspace evidence does not itself demonstrate spatial selectivity. The fixed independent mask and full contrast inventory make this distinction explicit.

The biological interpretation remains unresolved. HRF parameters and deconvolved signals are estimated from the same BOLD data. Their association can therefore reflect physiological variation, properties of the estimator or both; it is not independent vascular validation or mediation of a stimulation effect. SUDMEX did not establish treatment specificity, and Aachen lacked sham. The cohorts differ in diagnosis, stimulation protocol and interval, supporting cross-dataset comparison of associations rather than a pooled treatment estimate. Small samples, regional signal coverage, exploratory analysis selection and sensitivity to some baseline and omission checks further limit certainty.

For this pipeline, apparent changes in Φ should be interpreted together with their estimated haemodynamic context. The results support a reproducible association across two datasets, while independent physiological measurements and controlled replication are needed to establish its source.

## Data availability

The source imaging is publicly available in OpenNeuro ds001832, version 1.0.1 (doi:10.18112/openneuro.ds001832.v1.0.1), and ds003037, version 1.0.1 (doi:10.18112/openneuro.ds003037.v1.0.1). Analysis code, complete result families, the dashboard and manuscript are available at https://github.com/Parssayousefi/TMS_GWT. No private psychedelic participant data or effect maps were used.

## Ethics statement

This secondary analysis used publicly released, deidentified data and involved no new recruitment or participant contact. Original acquisition approvals and consent procedures are described in the source studies [@alkhasli2019; @garza2021].

## Acknowledgments

We thank the original investigators for making these datasets available and the developers of the analysis tools.

## Generative AI statement

OpenAI Codex assisted with code preparation and manuscript drafting. Numerical claims were checked against saved outputs and independent matrix calculations. Responsibility for the final content remains with the investigators.

## Figures and tables

{{{{FIGURE1}}}}

Figure 1. Outside-workspace Φ changes and estimated HRF-width changes in both datasets. Points represent independent participants. Both axes are residualized for motion change and, in SUDMEX, treatment arm; HRF-width changes are standardized within cohort before residualization. Lines show the corresponding adjusted slopes. Shapes distinguish SUDMEX treatment arms. The models do not identify vascular causality.

{{{{FIGURE2}}}}

Figure 2. HRF-width associations across spatial scopes. Points and marginal 95% intervals show separately fitted Aachen and SUDMEX slopes. Raw estimates are nats per cohort SD of width change; normalized estimates are percentile points per SD. Spatial differences require direct tests, not comparison of individual p values.

Table 1. HRF-width associations in both cohorts.

{{{{WIDTH_TABLE}}}}

Raw estimates are nats per SD; normalized estimates are percentile points per SD. Holm correction covers the complete fifteen-association family in each dataset and representation. All height, latency, direct spatial and sensitivity results are supplied in Supplementary Data.

Table 2. Direct outside-minus-inside treatment contrasts.

{{{{SPATIAL_TABLE}}}}

Aachen estimates average its two within-visit changes. SUDMEX estimates active-minus-sham differences adjusted for baseline spatial contrast and motion change. They are different treatment estimands and are not pooled. All raw estimates are nats; intervals are marginal 95% intervals. The reported correction is exact maximum-t5 for Aachen and Holm5 for SUDMEX.

## References

{{{{REFERENCES}}}}
'''
    # Keep the existing, audited acquisition/processing account in the supplement.
    previous = (OLD / 'supplement_source.md').read_text(encoding='utf-8')
    acquisition = previous.split('### Acquisition and processing details', 1)[1].split('### Fixed mask and information aggregation', 1)[0]
    sudmethods = previous.split('### Contextual SUDMEX analyses', 1)[1].split('## Supplementary results', 1)[0].split('A longitudinal sensitivity', 1)[0]
    supplement = f'''# Supplementary material

{TITLE}

## Acquisition and processing

### Aachen

{acquisition.strip()}

### SUDMEX

{sudmethods.strip()}

The present cross-dataset analysis uses only the randomized baseline-to-two-week comparison. For both datasets, saved HRF parameters were read from the production MATLAB outputs, with 454 regions, three finite parameters and TR=2 seconds verified in every analysed scan. The 140 relevant scans comprise 60 Aachen scans and 80 SUDMEX scans. No processing was rerun for this follow-up.

## Spatial estimands and models

The independent HCP mask contains 235 regions and its complement 219. Information summaries average unique off-diagonal pairs; regional HRF parameters average regions, so the pair and region weights differ. The outside-minus-inside contrast is computed within each participant. Cross-boundary edges contribute to the whole-brain summary but belong to neither within-compartment mean. The direct contrast is not inferred by comparing significance in separate compartments.

Aachen's two complete within-visit contrasts are averaged before inference. Five direct spatial information contrasts use 32,768 exact joint sign flips and Holm correction. SUDMEX models the follow-up spatial difference as a function of treatment arm, its baseline value and motion change, using HC3 uncertainty. Unadjusted, baseline-only and arm-specific descriptive contrasts remain available. Every five-measure family is retained in both raw and normalized representations.

## Complete haemodynamic family

Fifteen associations are retained for each dataset, representation and model: height, peak latency and width at whole-brain, within-workspace, outside-workspace and outside-minus-inside scopes (twelve), plus those three whole-brain HRF predictors for the outside-minus-inside Φ outcome. The last three models test spatial slope differences using a common predictor. They differ from models that relate the HRF spatial difference to the Φ spatial difference. Holm15 controls each complete family, including nonsignificant tests.

Primary association models include motion change and, for SUDMEX, randomized arm. HRF change is standardized within dataset. Baseline sensitivities additionally include baseline Φ and the corresponding baseline HRF parameter. These are association models and do not adjust a treatment effect for an established mediator. Arm-by-HRF interaction models include both main effects and their product. Arm-specific slopes retain the full-cohort SD scale.

Same-direction cross-dataset conjunctions use the maximum of the two cohort p values for each endpoint, setting the conjunction p to one when slope directions disagree. Holm correction then covers all fifteen endpoints within each representation. This intersection-union procedure asks whether both cohorts support a same-direction association; it is not a meta-analysis or a test of equal effect sizes. The follow-up was motivated by already observed Aachen findings and remains exploratory.

## Robustness and numerical verification

The complete raw fifteen-model family was refitted after each participant omission: 225 Aachen and 600 SUDMEX fits. Correction was repeated within every omission. All width slopes at whole-brain, within-workspace and outside-workspace scopes retained their negative direction. Corrected whole-brain support remained in 12/15 Aachen and 38/40 SUDMEX omissions; outside-workspace support remained in 14/15 and 26/40. These dependent checks do not provide additional independent replications.

Rank sensitivities use twelve-decimal rounding to resolve numerical ties far below HRF measurement resolution. Whole-brain raw width correlations survive Holm15 in both datasets; outside-workspace raw rank support survives in Aachen but not SUDMEX. Normalized outside-workspace rank associations survive in both. Baseline adjustment preserves both raw outside-workspace associations but not the normalized Aachen association. The full tables retain these qualifications.

All 140 saved HRF parameter arrays were checked. The earlier 60 Aachen HRF summaries, ten spatial contrasts and nine haemodynamic associations were reproduced. Each newly fitted HC3 result was checked against a separate pseudoinverse-based sandwich-covariance calculation. The largest absolute discrepancy was below 2.4×10⁻¹². Input hashes were verified after analysis. No original matrices, masks or processing files were changed.

## Scope of interpretation

No direct spatial treatment contrast survived correction in either dataset. The shared HRF-width association does not establish preferential modulation outside the workspace, stimulation specificity, a physiological mechanism or a causal neuronal pathway. Nonsignificant interactions are not equivalence tests. The fitted response and information estimates share input signals and a deconvolution procedure. Existing denoised-versus-deconvolved Aachen checks are available in the full repository; they are method sensitivities, not independent physiological validation.

The repository retains the full earlier PID, participation, gateway/broadcaster, clinical, target-flow and all-post analyses. Later SUDMEX open-label visits are excluded from the present randomized comparison. Historical Aachen outputs preceding anatomical correction are clearly marked as superseded.

## Supplementary figure

{{{{SUPPLEMENT_FIGURE1}}}}

Figure S1. Control workspaces compared with the frozen HCP reference. The supplied image is unchanged. SUDMEX sham maps before and after two weeks use the same eighteen participants; Aachen maps average pretreatment visits within fifteen participants. Shared regions and Dice overlap document broad correspondence with previous work [@luppi2024; @yousefi2026], with incomplete overlap. These descriptive maps do not redefine the mask used in the statistical analyses.

## Supplementary data

The accompanying ZIP retains all seventy spatial result rows, 270 haemodynamic model/rank rows, thirty cross-dataset conjunction rows, 825 participant-omission fits and the influence summaries, together with the analysis specification and numerical verification. Rows from overlapping models and representations are not independent tests. Original coefficient and p-value precision is retained. No missing value is replaced by zero.

## References

{{{{REFERENCES}}}}
'''
    assert len(spatial) == 70 and len(a) == 270
    main_text = re.sub(r'(?<![A-Za-z0-9])-(?=\d)', '−', main_text)
    (OUT / 'manuscript_source.md').write_text(main_text, encoding='utf-8')
    (OUT / 'supplement_source.md').write_text(supplement, encoding='utf-8')
    with zipfile.ZipFile(OUT / 'TMS_Cross_Dataset_Supplementary_Data.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for p in sorted((OUT / 'tables').iterdir()): archive.write(p, 'TMS_Cross_Dataset_Supplementary_Data/' + p.name)
    receipt = {'analysis_directory': str(RUN), 'title': TITLE, 'author_information_included': False, 'abstract_words': len(main_text.split('## Abstract')[1].split('## Introduction')[0].split()), 'body_words': len(main_text.split('## Introduction')[1].split('## Data availability')[0].split()), 'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__), *RUN.glob('*.csv'), *OUT.glob('*source.md')]}}
    (OUT / 'preparation_receipt.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in receipt.items() if k != 'sha256'}))


if __name__ == '__main__':
    main()
