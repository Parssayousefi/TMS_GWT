"""Write the corrected report from saved, independently reviewable TMS tables."""
from pathlib import Path
import json
import shutil
import sys
import re
import zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'pipeline/scripts/analysis'))
import tms_publication as tms
RUN=ROOT/'output/corrected_tms_analysis/20261126T112054_284547Z'
OLD=ROOT/'manuscript/aachen_brief_report'
OUT=OLD/'corrected_20261126'
TARGET=ROOT/'output/corrected_aachen_target_flow/20261126T113156_656159Z'
QC=ROOT/'output/control_workspace_qc/combined_20261125T234309_541220Z'


def main():
    OUT.mkdir(exist_ok=True);(OUT/'figures').mkdir(exist_ok=True)
    data=OUT/'submission_data';data.mkdir(exist_ok=True)
    raw=pd.read_csv(RUN/'raw_aachen_all40.csv').set_index('metric')
    norm=pd.read_csv(RUN/'normalized_aachen_all40.csv').set_index('metric')
    exp=pd.read_csv(RUN/'gwt_expectations_raw.csv').query("dataset == 'aachen_itbs'").set_index('metric')
    phi=raw.loc['whole_brain__phi'];wphi=raw.loc['within_workspace__phi'];phir=raw.loc['whole_brain__phi_r'];nphi=norm.loc['whole_brain__phi']
    redundancy=raw.loc['whole_brain__redundancy'];rpc=exp.loc['mean_redundancy_pc'];fpn=exp.loc['gateways_FPN']
    assert not (raw.p_maxT_all<.05).any() and not (norm.p_maxT_all<.05).any()
    for branch,frame in [('raw',raw),('normalized',norm)]:
        selected=frame.loc[[s+'__'+m for s in ['whole_brain','within_workspace'] for m in tms.MATRIX_NAMES]].reset_index()
        selected['scope']=selected.metric.str.split('__').str[0].map({'whole_brain':'Whole brain','within_workspace':'Workspace'})
        selected['measure']=selected.metric.str.split('__').str[1].map({'synergy':'Persistent synergy','redundancy':'Persistent redundancy','synergy_total':'Source synergy','phi':'Φ','phi_r':'ΦR'})
        if branch=='normalized':selected[['estimate','ci_low','ci_high']]*=100
        selected.to_csv(OUT/('table_1_scope_results.csv' if branch=='raw' else 'table_2_scope_results.csv'),index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(9.6,7.2),layout='constrained')
    for row,branch in enumerate(['raw','normalized']):
        p=pd.read_csv(RUN/(branch+'_aachen_participant_changes.csv'))
        for col,scope in enumerate(['whole_brain','within_workspace']):
            ax=axes[row,col];m=scope+'__phi';scale=100 if branch=='normalized' else 1
            values=p[[m+'_baseline',m+'_followup']].to_numpy()*scale
            for a,b in values:ax.plot([0,1],[a,b],color='#7d9eac',alpha=.65,lw=1,marker='o',ms=3)
            ax.plot([0,1],values.mean(0),color='#1c4f65',marker='D',lw=2.5,ms=6,label='Mean')
            ax.set_xticks([0,1],['Pre mean','Post mean']);ax.set_xlim(-.16,1.16)
            ax.set_title(('Whole brain' if col==0 else 'Fixed workspace')+(' • raw' if row==0 else ' • normalized'))
            ax.set_ylabel('Mean pairwise Φ (nats)' if row==0 else 'Conditional percentile')
    fig.savefig(OUT/'figures/Figure_1_participant_changes.png',dpi=220);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(9.6,5.8),layout='constrained')
    for ax,m,title in zip(axes.flat,['gateways_FPN','gateways_DAN','mean_redundancy_pc','mean_synergy_pc'],['FPN gateways','DAN gateways','Redundancy participation','Synergy participation']):
        for y,branch,color in [(1,'raw','#1c4f65'),(0,'normalized','#8c5b9d')]:
            rr=pd.read_csv(RUN/('gwt_expectations_'+branch+'.csv')).query("dataset == 'aachen_itbs'").set_index('metric').loc[m]
            ax.errorbar(rr.estimate,y,xerr=[[rr.estimate-rr.ci_low],[rr.ci_high-rr.estimate]],fmt='o',color=color,capsize=4)
        ax.axvline(0,color='.55',lw=1,ls='--');ax.set_yticks([0,1],['Normalized','Raw']);ax.set_ylim(-.6,1.6)
        ax.set_title(title);ax.set_xlabel('Post minus pre (regions)' if m.startswith('gateways') else 'Post minus pre (coefficient)')
    fig.savefig(OUT/'figures/Figure_2_gwt_expectations.png',dpi=220);plt.close(fig)
    shutil.copyfile(QC/'control_workspace_comparison.png',OUT/'figures/Figure_S1_control_workspace_comparison.png')
    assert tms.sha256(QC/'control_workspace_comparison.png')==tms.sha256(OUT/'figures/Figure_S1_control_workspace_comparison.png')
    loo=pd.read_csv(RUN/'aachen_leave_one_out.csv').query("branch == 'raw'")
    fig,axes=plt.subplots(1,2,figsize=(9.5,6),layout='constrained')
    for ax,m,title in zip(axes,['whole_brain__phi','within_workspace__phi'],['Whole brain','Fixed workspace']):
        d=loo.loc[loo.metric.eq(m)].copy();y=np.arange(len(d))
        ax.errorbar(d.estimate,y,xerr=[d.estimate-d.ci_low,d.ci_high-d.estimate],fmt='o',capsize=3,color='#1c4f65')
        ax.axvline(0,color='.5',ls='--',lw=1);ax.set_yticks(y,d.omitted_subject.str.replace('aachen_itbs_',''))
        ax.set_title(title);ax.set_ylabel('Omitted participant');ax.set_xlabel('Mean Φ change (nats), 95% CI')
    fig.savefig(OUT/'figures/Figure_S2_leave_one_out.png',dpi=220);plt.close(fig)
    hemo=pd.read_csv(RUN/'haemodynamic_followup/participant_associations.csv')
    fig,axes=plt.subplots(1,2,figsize=(9.5,4.3),layout='constrained')
    for ax,scope,title in zip(axes,['whole_brain','outside_workspace'],['Whole brain','Outside workspace']):
        d=hemo.loc[hemo.scope.eq(scope)&hemo.parameter.eq('width_seconds')]
        design=np.column_stack([np.ones(len(d)),d.delta_fd])
        x=d.hrf_change.to_numpy();y=d.phi_change.to_numpy()
        x=x-design@np.linalg.lstsq(design,x,rcond=None)[0];y=y-design@np.linalg.lstsq(design,y,rcond=None)[0]
        ax.scatter(x,y,color='#1c4f65',s=28)
        for xx,yy,sid in zip(x,y,d.subject_id):ax.annotate(sid.split('_')[-1],(xx,yy),xytext=(4,4),textcoords='offset points',fontsize=7)
        grid=np.array([min(x),max(x)]);ax.plot(grid,np.dot(x,y)/np.dot(x,x)*grid,color='#8c5b9d',lw=1.5)
        ax.axhline(0,color='.75',lw=.7);ax.axvline(0,color='.75',lw=.7)
        ax.set_title(title);ax.set_xlabel('HRF width change residual (seconds)');ax.set_ylabel('Φ change residual (nats)')
    fig.savefig(OUT/'figures/Figure_S3_hrf_width_associations.png',dpi=220);plt.close(fig)
    original=(OLD/'manuscript_source.md').read_text(encoding='utf-8')
    front=original.split('## Abstract')[0]
    front=front.replace('Whole brain and workspace information dynamics following prefrontal theta burst stimulation','Information dynamics following prefrontal theta burst stimulation')
    front=front.split('\n\n',1)[0]+'\n\nArticle type:'+front.split('Article type:',1)[1]
    tail=original.split('## Data availability')[1].split('## Figures and tables')[0]
    tail=re.sub(r'## Author contributions\n.*?(?=## Acknowledgments\n)','',tail,flags=re.S)
    tail=tail.replace('[TODO: Insert the persistent public repository URL or DOI for the analysis code and submission materials before submission.]','Analysis code and submission materials are available at https://github.com/Parssayousefi/TMS_GWT.')
    pcts=lambda x:f'{100*x:.3f}'
    manuscript=front+f'''## Abstract

Prefrontal intermittent theta burst stimulation (iTBS) may influence distributed brain information processing, but a workspace response has not been established. We applied integrated information decomposition to resting-state fMRI from 15 healthy adults scanned before and after two active iTBS intensities. Five information measures and 20 participation and role outcomes were evaluated using 454 haemodynamically deconvolved regional signals and an independent 235-region workspace. Matched visit changes were averaged within participant. Control maps showed broad correspondence with the reference workspace. Whole-brain mean pairwise Φ changed by {phi.estimate:.5f} nats (95% confidence interval {phi.ci_low:.5f} to {phi.ci_high:.5f}; five-measure maximum-t p={phi.p_maxT_family:.3f}). No information measure survived correction within the whole-brain or workspace family, and conditional-null normalization did not establish an effect. Four outcomes motivated by psychedelic workspace findings also failed family-wise correction, including frontoparietal gateway recruitment. These exploratory results recover the expected spatial architecture without establishing a consistent acute reorganization after iTBS. Small sample size, absent sham and measurement dependence limit inference; nonsignificance does not demonstrate equivalence.

## Introduction

Prefrontal intermittent theta burst stimulation (iTBS) can affect connectivity beyond the stimulation target [@alkhasli2019]. Information decomposition offers a complementary description of these distributed dynamics by separating information shared redundantly by regions from information available only through their joint activity [@mediano2025]. Whether prefrontal stimulation measurably changes this information architecture remains unclear.

Previous work identified a synergistic workspace spanning association and subcortical regions [@luppi2024]. Yousefi and colleagues subsequently reported altered cross-network participation and recruitment of frontoparietal and dorsal-attention gateways during psychedelic conditions [@yousefi2026]. These findings motivate comparison across interventions, but do not establish that stimulating the dorsolateral prefrontal cortex (DLPFC) should reproduce a psychedelic response or increase mean Φ.

We examined information measures and workspace roles before and after two active iTBS visits in healthy adults. We first assessed whether control data recovered the previously described spatial pattern, then tested within-participant changes using the fixed independent workspace. We retained all planned metric families and examined normalization, motion, spatial scope and target-centred information flow.

## Methods

### Dataset and processing

We reanalysed OpenNeuro ds001832, version 1.0.1 [@alkhasli2019]. Participants received 600-pulse left-DLPFC iTBS at 90% and 120% of resting motor threshold on separate counterbalanced visits. Ten-minute resting scans preceded stimulation and began approximately seven minutes afterward; stimulation occurred outside the scanner. One of 16 participants was excluded because slice-timing metadata did not match acquisition dimensions. Fifteen participants contributed both complete pre/post pairs, yielding 60 scans. There was no sham condition.

Processing followed the established GWT framework [@yousefi2026]: fMRIPrep, CONN denoising, Schaefer-400/Tian-54 parcellation and blind haemodynamic deconvolution with rsHRF [@wu2013]. Already skull-stripped anatomical inputs were processed without repeat brain extraction. Anatomical alignment, retained signal coverage and all downstream outputs were checked. Gaussian minimum-mutual-information ΦID used a lag of one repetition time, two seconds [@luppi2024; @mediano2025]. Supplementary Methods specify processing and remaining regional coverage limitations.

### Information measures and workspace roles

Whole-minus-sum information, Φ, is joint past–future mutual information minus the sum of the two regions' own past–future mutual information [@barrett2011]. It is not synergy alone and may be negative. Redundancy-corrected ΦR adds persistent redundancy. Persistent synergy, persistent redundancy and total source-side synergy complete the five-measure information family.

We averaged each measure over 102,831 unique whole-brain pairs and over the 27,495 pairs wholly within the fixed 235-region HCP workspace used by Yousefi and colleagues [@yousefi2026]. Outside-workspace and cross-boundary summaries were retained as spatial sensitivities. These are pairwise averages, not integrated information of the brain considered as a single system. They also differ from Luppi's edgewise network-based inference [@luppi2024].

Weighted participation coefficients describe how each node's synergistic or redundant connections spread across eight modules. They use all 454 destinations; averages and role classification use the fixed workspace. A node is a gateway when its synergy-participation rank exceeds its redundancy-participation rank, and a broadcaster for the reverse ordering; ties are neither. We retained two mean coefficients, two total role counts and eight network counts for each role. Four literature-motivated expectations were evaluated together: increased frontoparietal (FPN) gateways, dorsal-attention (DAN) gateways, redundancy participation and synergy participation [@yousefi2026]. The evidence for synergy participation was less consistent across drug conditions; these comparisons are exploratory extensions to iTBS.

### Statistical analysis and sensitivities

Post-minus-pre differences were calculated within each visit and averaged within participant, retaining 15 independent observations. Time series were not pooled. Two-sided maximum-|t| tests used all 32,768 shared participant sign patterns. The five information measures formed a family at each spatial scope; broader corrections covered ten whole-brain/workspace measures, twenty information-by-scope measures and all forty information/role endpoints. The four literature-motivated outcomes received Holm and maximum-t correction as one family. These procedures assume jointly sign-symmetric null errors and are not sham-randomization tests. Confidence intervals are marginal 95% intervals. Complete families and all results are supplied.

A NuMIT-inspired sensitivity compared each measure with a bivariate autoregressive model null matched to joint time-delayed mutual information [@liardi2025; @blackburne2026]. Two independently generated references used the corrected baseline signals; normalized results are conditional percentiles and are not additive information atoms. Motion-adjusted HC3 models, intensity contrasts, participant omission, equal record lengths, spatial contrasts and signal/HRF diagnostics were retained. Target-centred flow was assessed in eight participants with usable deposited bilateral masks. SUDMEX supplied separate sham-workspace QC and contextual randomized and longitudinal analyses, detailed in the supplement.

## Results

### Recovery of the reference workspace

The independently derived Aachen baseline workspace shared 176 regions with the fixed HCP workspace (Dice=0.743). Matched SUDMEX sham maps showed Dice values of 0.776 before treatment and 0.741 after two weeks (Supplementary Figure S1). These findings show broad spatial correspondence with the architecture described by Luppi and Yousefi [@luppi2024; @yousefi2026]. The derived maps were not identical, and the independent HCP mask remained unchanged for intervention tests.

### Information measures and GWT expectations

Whole-brain Φ changed by {phi.estimate:.6f} nats (95% CI {phi.ci_low:.6f} to {phi.ci_high:.6f}; maximum-t p5={phi.p_maxT_family:.3f}, p10={phi.p_maxT_whole_and_workspace10:.3f}; Figure 1, Table 1). Within-workspace Φ changed by {wphi.estimate:.6f} nats ({wphi.ci_low:.6f} to {wphi.ci_high:.6f}; p5={wphi.p_maxT_family:.3f}). No information outcome survived its five-measure correction at either scope. Normalized whole-brain Φ changed by {pcts(nphi.estimate)} percentile points ({pcts(nphi.ci_low)} to {pcts(nphi.ci_high)}; p5={nphi.p_maxT_family:.3f}; Table 2). No endpoint survived the broader forty-test family in either representation.

FPN gateway count changed by {fpn.estimate:.2f} regions (95% CI {fpn.ci_low:.2f} to {fpn.ci_high:.2f}; Holm4 p={fpn.p_holm4:.3f}). None of the four GWT expectations survived correction (Figure 2). Redundancy participation increased by {rpc.estimate:.5f}, with an uncorrected p={rpc.p:.4f}, but Holm4 p={rpc.p_holm4:.3f}. This coefficient describes cross-network distribution; whole-brain persistent redundancy itself changed by only {redundancy.estimate:.6f} nats (uncorrected p={redundancy.p:.3f}).

### Spatial and sensitivity findings

No spatial information family established an outside-workspace increase, and the direct outside-minus-inside Φ contrast was uncertain (−0.000771 nats; five-measure maximum-t p=0.935). Equal-length estimation and participant omission did not establish a corrected global effect. Target-centred tests likewise yielded no corrected primary finding in the eight-person subset.

Diagnostics identified dependence on the signal model: within-workspace Φ changes differed between deconvolved and denoised representations (eight-test maximum-t p=0.00250), and outside-workspace Φ changes covaried with estimated HRF width and peak latency (six-test Holm p=0.000318 and 0.0331). These are method and association results, not evidence of a stimulation-induced information increase. Full diagnostics, the separate SUDMEX analyses and clinical associations are reported in the supplement.

## Discussion

Control scans recovered a broadly similar workspace to the reference literature, but the corrected iTBS analyses did not establish increased global pairwise Φ, increased workspace integration or recruitment of the gateway pattern reported under psychedelics. Spatial agreement and intervention sensitivity address different questions: recovering the expected architecture supports implementation consistency without guaranteeing detectable treatment effects.

Redundancy participation provided a nominal positive signal, but it did not survive the four-outcome correction and was not accompanied by greater whole-brain redundant information. A participation coefficient measures the distribution of weighted connections across modules, whereas persistent redundancy measures shared temporal information. Conflating these quantities would turn a tentative topological observation into an unsupported claim about information magnitude.

Several features constrain interpretation. Fifteen participants provide limited precision; the confidence intervals admit both decreases and increases in Φ. Absence of sham prevents separating stimulation from repeated scanning, arousal, expectancy and physiological changes. The two active visits are repeated observations of the same participants, not independent replications. Target masks do not measure the induced electric field, and the smaller localization subset further limits mechanistic inference. No subjective-experience or therapeutic outcome in Aachen supports claims about consciousness or clinical benefit.

The Gaussian estimator, atlas, regional coverage and estimated haemodynamic model also matter. Deconvolution was part of the intended processing chain, but its outputs remain model-dependent estimates. The observed method interaction and HRF associations constrain a direct neuronal interpretation; they neither prove a vascular artifact nor identify a biological mechanism. Conditional-null normalization addresses a specified model comparison, not all estimator bias or uncertainty in the learned reference.

The supported contribution is a reproducible assessment of information architecture after acute prefrontal iTBS: recognizable control-workspace organization, with uncertain intervention effects under the tested summaries. These estimates define the precision and measurement assumptions that a controlled replication must address.

## Data availability'''+tail+'''## Figures and tables

{{FIGURE1}}

Figure 1. Mean pairwise Φ before and after iTBS. Each line links a participant's average of the two pretreatment visits to the average of the two post-treatment visits; diamonds show group means. Raw and normalized whole-brain and fixed-workspace summaries use the same 15 people. Conditional percentiles are not nats. No displayed Φ outcome survives the five-information-measure correction.

{{FIGURE2}}

Figure 2. Four GWT-related outcomes after iTBS. Points and marginal 95% intervals show participant-averaged changes in FPN gateways, DAN gateways, redundancy participation and synergy participation. Raw and conditional-null-normalized graphs use the same participants and definitions. The literature-motivated direction is positive, with weaker prior consistency for synergy participation. No outcome survives correction across the four expectations. The coefficient axes differ from the count axes.

Table 1. Raw information changes across the whole brain and fixed workspace.

{{TABLE1}}

Changes and marginal 95% confidence intervals are in nats. Maximum-t p5 corrects the five information measures within each scope; p10 corrects both scopes jointly. All rows use 15 participant-averaged changes. Wider thirty- and forty-endpoint corrections are provided in Supplementary Data.

Table 2. Conditional-null-normalized information changes.

{{TABLE2}}

Changes and marginal 95% confidence intervals are in percentile points. Correction families match Table 1. Normalization is a dependent model sensitivity and does not preserve additive information identities.

## References

{{REFERENCES}}
'''
    manuscript=manuscript.replace('The original imaging dataset is publicly available under CC0','The Aachen imaging dataset is publicly available under CC0')
    manuscript=manuscript.replace('Participant-level scalar outputs, complete endpoint tables and numerical sensitivity results accompany this report as Supplementary Data.','Participant-level scalar outputs, complete endpoint tables and numerical sensitivity results accompany this report as Supplementary Data. Supplementary SUDMEX analyses use OpenNeuro ds003037, version 1.0.1, doi:10.18112/openneuro.ds003037.v1.0.1.')
    manuscript=manuscript.replace('These exploratory results recover the expected spatial architecture','Individual Φ changes covaried with fitted haemodynamic-response width in exploratory analyses. These results recover the expected spatial architecture')
    manuscript=manuscript.replace('Motion-adjusted HC3 models, intensity contrasts, participant omission, equal record lengths, spatial contrasts and signal/HRF diagnostics were retained.','Motion-adjusted HC3 models, intensity contrasts, participant omission, equal record lengths, spatial contrasts and signal/HRF diagnostics were retained. An expanded exploratory family related Φ change to three fitted HRF parameters at three scopes (nine Holm-corrected associations).')
    manuscript=manuscript.replace('Diagnostics identified dependence on the signal model: within-workspace Φ changes differed between deconvolved and denoised representations (eight-test maximum-t p=0.00250), and outside-workspace Φ changes covaried with estimated HRF width and peak latency (six-test Holm p=0.000318 and 0.0331). These are method and association results, not evidence of a stimulation-induced information increase. Full diagnostics, the separate SUDMEX analyses and clinical associations are reported in the supplement.','Diagnostics identified dependence on the signal model: within-workspace Φ changes differed between deconvolved and denoised representations (eight-test maximum-t p=0.00250). In the expanded nine-association family, greater Φ increases accompanied narrower fitted HRFs across the whole brain (−0.00459 nats per SD of width change, 95% CI −0.00693 to −0.00226; Holm p=0.00849) and outside the workspace (p=0.000477). Outside-workspace peak-latency association also survived (p=0.0463). Width slopes retained their direction in every participant omission, although corrected support did not persist in every omission (Supplementary Figure S3). These are same-signal model associations, not evidence that stimulation altered vascular physiology.')
    manuscript=manuscript.replace('The Gaussian estimator, atlas, regional coverage and estimated haemodynamic model also matter.','The association between narrower fitted HRFs and greater Φ changes provides a possible measurement-focused interpretation. The Gaussian estimator, atlas and regional coverage also matter.')
    manuscript=manuscript.replace('The observed method interaction and HRF associations constrain a direct neuronal interpretation; they neither prove a vascular artifact nor identify a biological mechanism.','Because the HRF is fitted from the same signals subsequently deconvolved, parameter–Φ coupling may arise partly from the estimation procedure. The observed associations neither prove a vascular artifact nor establish mediation of an iTBS effect. Independent haemodynamic measurements would be needed to separate physiological variation from modelling dependence.')
    (OUT/'manuscript_source.md').write_text(manuscript,encoding='utf-8')
    old_supp=(OLD/'supplement_source.md').read_text(encoding='utf-8')
    methods=old_supp.split('## Additional methods')[1].split('## Supplementary results')[0]
    methods=methods.replace('Susceptibility-distortion correction and FreeSurfer surface reconstruction were not performed.','Because the deposited anatomical images were already skull-stripped, fMRIPrep used --skull-strip-t1w skip. Susceptibility-distortion correction and FreeSurfer surface reconstruction were not performed.')
    methods=methods.replace('The historical whole-brain bootstrap used seed 2026112515; the workspace audit used seed 2026112526.','Current whole-brain and workspace checks used seed 2026112526; spatial follow-up checks used seed 2026112528.')
    methods=methods.replace('Independent log-determinant calculations agreed with the denoised-BOLD temporal information calculation to within 1.3×10⁻¹⁴ nats.','Independent log-determinant calculations agreed with the denoised-BOLD temporal information calculation to within 2.8×10⁻¹⁴ nats.')
    methods=methods.replace('After examining the fixed-workspace results, we evaluated two targets:','The retained spatial follow-up evaluated two targets:')
    methods=methods.replace('Finite-record calibration is retained in Table S7','Finite-record calibration is retained in Table S9')
    methods=methods.replace('For each node, the weighted participation coefficient is PCi=1−Σm(sim/si)², where sim is the sum of edge weights to module m and si is their sum across modules.','For each node, participation is one minus the sum, across modules, of the squared fractions of its total edge weight assigned to each module.')
    supplement='''# Supplementary material

Information dynamics following prefrontal theta burst stimulation

## Additional methods'''+methods+'''### Anatomical and signal quality control

All 15 corrected anatomical reconstructions retained 99.63–100% of the supplied native anatomical support; BOLD/T1 mask Dice ranged from 0.881 to 0.956. All 454 regional signals were finite and nonconstant in all 60 retained scans. Some inferior temporal and orbitofrontal parcels nevertheless had low functional brain-mask coverage. Nonzero regional signals and coherent registration do not guarantee full regional signal quality. We retained the predefined cohort and fixed reference rather than selecting participants or parcels according to treatment effects. Regional coverage tables accompany the reproducibility records.

### Recovery of the previously described workspace

We independently derived control workspaces with the same synergy-minus-redundancy strength-rank method used in the earlier pipeline [@luppi2024; @yousefi2026]. Persistent-synergy and persistent-redundancy matrices were arithmetically averaged, node strengths were ranked separately across all 454 regions, and a positive synergy-minus-redundancy rank identified workspace membership. For Aachen, the two pretreatment matrices were first averaged within participant and then across the same 15 participants. The SUDMEX comparison used the same 18 sham participants at baseline and after the two-week randomized phase [@garza2021]. Post-stimulation Aachen scans and later open-label SUDMEX visits did not define these control maps.

The maps recovered a broadly similar spatial architecture to the canonical workspace described by Luppi and Yousefi [@luppi2024; @yousefi2026] (Figure S1). The SUDMEX sham-baseline map contained 234 regions, sharing 182 with the 235-region HCP reference (Dice=0.776); the two-week map contained 232 regions, sharing 173 (Dice=0.741). The Aachen paired-baseline map contained 239 regions, sharing 176 (Dice=0.743). Corresponding strength-gradient Spearman correlations with HCP were 0.652, 0.615 and 0.596. Aachen FPN overlap was particularly high (43 shared regions; network Dice=0.905), with larger differences in some sensory/attention networks.

This is descriptive spatial correspondence, not proof of identical workspaces or a spatial-null significance test. The copied HCP mask and the HCP mask used in the GWT project are byte-identical and constitute one comparator, not two independent replications. Independently derived TMS maps served quality control only; intervention analyses retained the frozen HCP mask. Aachen baseline-visit maps had Dice=0.782 with one another, and leave-one-person-out maps had Dice=0.960–0.981 with the full baseline map. These checks do not establish a treatment effect or validate every regional signal.

### Literature motivated expectations

Four outcomes form one complete exploratory family: FPN gateways, DAN gateways, mean redundancy participation and mean synergy participation. The directional comparison is an increase, but all tests remain two-sided. The first three changes were reported across the classic-psychedelic conditions in the cited paper; synergy participation increased in three of four conditions [@yousefi2026]. Redundancy participation also increased with non-classic drugs, and DAN changes were not exclusive to classic psychedelics. The paper provides no matching directional prediction for global mean Φ or ΦR. We used published findings only, without private psychedelic participant data or maps. Broadcaster changes are complementary role summaries, not independent confirmation of gateway effects.

### Target centred information flow

Corrected individual deposited masks supplied cortical weights for eight Aachen participants and both hemispheres. All sixteen sides met the unchanged cortical-overlap criterion. Five quantities were retained: outgoing-minus-incoming transfer entropy, its unique and synergistic contributions, persistent synergy and persistent redundancy. Weighted target-to-other-region pair summaries excluded self-pairs. The left target and left-minus-right contrast were each tested as a five-measure family; Holm correction across ten tests supplied a wider sensitivity. Exact tests used all 256 participant sign patterns after averaging the two visit changes. Lag-two, equal-length, motion/baseline, intensity and participant-omission checks were retained. Right-sided masks were contralateral comparators, not evidence of right-sided stimulation. Pairwise predictive information does not identify a direct neuronal causal pathway or eliminate common causes.

### Contextual SUDMEX analyses

SUDMEX is a double-blind randomized trial of 5-Hz left-DLPFC rTMS in cocaine use disorder [@garza2021]. Our acute MRI comparison retained 40 complete pairs, 22 active and 18 sham, at baseline and two weeks. Models compared follow-up outcomes by arm with baseline adjustment and then baseline plus change in mean framewise displacement, using HC3 uncertainty. Sham maps in Figure S1 use the same 18 complete sham pairs. SUDMEX anatomy was audited separately; no Aachen-like repeat-extraction failure was found, although inferior coverage limitations remained.

A longitudinal sensitivity compared baseline with the within-person average of every available post measurement in 41 participants. Later visits included open-label treatment, so this analysis cannot estimate the randomized two-week treatment effect. Its complete forty-endpoint family used 200,000 joint random sign vectors, including the identity pattern. Clinical associations related seven published craving/impulsivity outcomes to eight neural summaries, adjusting for arm, both baseline values and motion change (56-test Holm family). Source-flagged clinical changes were masked; no missing outcome was imputed.

## Supplementary results

### Information and workspace outcomes

No Aachen information family or GWT expectation established a corrected intervention-associated change, before or after conditional normalization. The direct outside-minus-inside Φ difference was −0.000771 nats (95% CI −0.003290 to 0.001747; maximum-t p5=0.935), giving no support to preferential change outside the workspace. All leave-one-participant-out whole-brain Φ estimates remained positive (0.000399–0.002029 nats), but this sign consistency did not establish a reliable group effect (Figure S2). Equal-length whole-brain Φ was +0.001415 nats (95% CI −0.001457 to 0.004287). Complete intensity, motion, influence and reference-seed tables are supplied.

Raw redundancy participation increased by 0.00466 (95% CI 0.00054–0.00879), nominal p=0.0295, but four-expectation Holm p=0.118 and maximum-t p=0.114. Normalized redundancy participation increased by 0.00119, Holm4 p=0.153. Neither establishes the predicted GWT pattern. Raw whole-brain persistent redundancy was essentially unchanged (+0.000047 nats, p=0.954). Normalized DAN gateways decreased nominally (−0.47 regions, p=0.0382), opposite to the expected direction, and failed four-outcome correction (Holm p=0.153).

### Signal representation and estimated haemodynamics

Deconvolution is part of the intended processing [@yousefi2026; @wu2013]. Comparison with the denoised representation is a model sensitivity, not a requirement that preprocessing preserve every effect. The direct deconvolved-minus-denoised within-workspace Φ change was +0.002011 nats (eight-test maximum-t p=0.00250). Whole-brain and outside-workspace method interactions did not survive that family. No Φ change in the three-representation, four-scope family survived correction.

None of nine mean HRF-parameter changes survived correction. In the six motion-adjusted association models, however, outside-workspace Φ change was associated negatively with HRF width change (coefficient −0.00755 nats per SD, Holm p=0.000318) and positively with peak-latency change (+0.00713 nats per SD, Holm p=0.0331). These same-signal associations do not establish an independent vascular mechanism or a stimulation effect. They show that estimated temporal-information changes covary with the fitted haemodynamic model.

The fresh Aachen conditional reference used 30,720 baseline coefficient matrices, 369 grid points per seed and 32,768 retained draws per point. Maximum interpolation error was 0.000761 CDF units. Simulated finite records showed a maximum mean-percentile departure from 0.5 of 0.0631; normalization is therefore a model sensitivity, not a claim of universally unbiased estimates. The shared fitted prior's uncertainty is not propagated. No historical Aachen reference distribution was reused.

### Target centred and contextual results

No corrected primary target-centred Aachen effect was established in eight participants (minimum within-scope maximum-t p=0.180), and lag-two and equal-length families also had no corrected effect. Three left-minus-right estimates in the baseline-plus-motion HC3 sensitivity survived Holm correction across ten tests: net transfer entropy and its unique and synergistic contributions (each adjusted p=0.0476). These small-sample models have five residual degrees of freedom and maximum leverage up to 0.896. A separate intensity-descriptive family identified decreased left-target persistent synergy at 120% rMT (−0.04120 nats, 95% CI −0.06088 to −0.02152; Holm20 p=0.0331). These model- and intensity-specific findings use the same eight participants and do not establish a robust primary flow effect or causal propagation. The full tables and model diagnostics are retained. The SUDMEX randomized target analysis likewise had no corrected primary result.

SUDMEX's unadjusted acute FPN change contrast was +3.14 gateways (p=0.0486), but the baseline-adjusted estimate was +1.21 (p=0.449) and baseline-plus-motion-adjusted estimate +1.09 (95% CI −2.13 to 4.30; p=0.498). No four-expectation result survived Holm correction in either representation, and no clinical–neural association survived the 56-test family.

The longitudinal SUDMEX sensitivity retained changes that cannot be attributed specifically to randomized treatment: raw whole-brain persistent synergy rose by 0.01815 nats (forty-test maximum-t p=0.00525), while normalized within-workspace source-side synergy fell by 1.664 percentile points (p40=0.0431). Raw whole-brain persistent synergy did not increase after normalization. These results use different time points, estimands and representations from the Aachen acute comparison and do not form a coherent replicated TMS signature.

## Numerical audit

All sixty corrected Aachen scans retained 102,831 finite unique regional pairs. The identity ΦR=Φ+persistent redundancy, symmetry, scan identities and source hashes were checked. Direct covariance calculations reproduced saved MATLAB Φ/ΦR outputs; independent sign-pattern enumeration verified inference. Target-flow numerical identities and weighted summaries agreed with existing production MATLAB outputs. A fresh selected-pair MATLAB validation could not run because the license server was unavailable; this limitation and the retained verification inputs are documented with the target results. The control-workspace maps, overlap statistics and reference hashes passed a separate independent review. Production processing code and frozen HCP/psychedelic computations were unchanged.

## Supplementary tables and participant outputs

{{SUPPLEMENT_TABLE}}

Raw information estimates are in nats; normalized information CSV estimates use CDF units, so multiplication by 100 gives percentile points. Participation coefficients and region counts retain their own units. P values are two-sided; confidence intervals are marginal. Each table identifies its family. Missing values indicate unestimable models or untested diagnostics and must not be read as zero. Full test inventories are retained, including nonsignificant outcomes.

## Supplementary figures

{{SUPPLEMENT_FIGURE1}}

Figure S1. Recovery of the reference workspace in TMS control scans. Rows show SUDMEX sham participants before treatment and after the two-week randomized phase (the same 18 people), and Aachen pretreatment visits averaged within participant (15 people). Blue indicates regions shared with the fixed HCP workspace, orange regions present only in the TMS-derived workspace and purple regions present only in HCP. Axial coordinates are in millimetres. The figure documents broad spatial correspondence with prior workspace work [@luppi2024; @yousefi2026], with incomplete overlap. The maps do not replace the independent HCP mask used in intervention analyses.

{{SUPPLEMENT_FIGURE2}}

Figure S2. Leave-one-participant-out raw Φ changes. Each point uses the remaining fourteen participants' mean visit changes; bars are marginal 95% Student-t intervals. Labels identify the omitted participant. These dependent influence checks are not independent replications.

## References

References use the numbering of the main manuscript.

{{REFERENCES}}
'''
    supplement=supplement.replace('The fresh Aachen conditional reference used 30,720 baseline coefficient matrices,','A bounded follow-up added three whole-brain associations to the six existing models, retaining all nine under Holm correction. Whole-brain HRF parameters were averages across all 454 regions, calculated before within-visit differencing and participant averaging. Whole-brain Φ change was inversely associated with width change (−0.00459 nats per SD, 95% CI −0.00693 to −0.00226; Holm9 p=0.00849). In this expanded family, outside-width and outside-peak p values were 0.000477 and 0.0463. Width slopes remained negative in all fifteen participant omissions; Holm9 support remained in twelve omissions for whole-brain and fourteen for outside-workspace Φ. Unadjusted rank sensitivities gave whole-brain width ρ=−0.695 (Holm9 p=0.0257) and outside width ρ=−0.854 (p=0.000466). The direct spatial-contrast width association did not survive the motion-adjusted nine-test family. These are dependent robustness checks, and the HRF estimates and Φ outputs share the same input signals (Figure S3).\n\nThe fresh Aachen conditional reference used 30,720 baseline coefficient matrices,')
    supplement=supplement.replace('## References\n\nReferences use the numbering','{{SUPPLEMENT_FIGURE3}}\n\nFigure S3. Association between estimated HRF width and Φ changes. Points represent fifteen participants; labels show dataset identifiers. The two visit changes are averaged within participant. Axes show residuals after linear adjustment for change in mean framewise displacement, and lines show the corresponding partial slopes. Inferential coefficients in Table S11 use standardized HRF changes and HC3 uncertainty. Fitted HRF parameters and information estimates originate from the same signals, so the association does not distinguish physiological effects from mathematical coupling introduced by deconvolution.\n\n## References\n\nReferences use the numbering')
    (OUT/'supplement_source.md').write_text(supplement,encoding='utf-8')
    mapping=[]
    for number,(label,files) in enumerate([
        ('Main information and role outcomes',['raw_aachen_main30.csv','normalized_aachen_main30.csv','raw_aachen_all40.csv','normalized_aachen_all40.csv']),
        ('Four GWT expectations',['gwt_expectations_raw.csv','gwt_expectations_normalized.csv']),
        ('Spatial information and direct interaction',['raw_aachen_information20.csv','normalized_aachen_information20.csv','aachen_spatial_interactions.csv']),
        ('Motion and intensity',['aachen_all_spatial_motion.csv','aachen_intensity_results.csv','aachen_robustness.csv']),
        ('Influence and reference seeds',['aachen_leave_one_out.csv','normalization_seed_tests.csv']),
        ('SUDMEX acute and longitudinal',['raw_sudmex_acute.csv','normalized_sudmex_acute.csv','sudmex_all_post_all40_results.csv','sudmex_all_post_arm_results.csv']),
        ('Clinical associations',['clinical_neural_associations.csv'])],1):
        for f in files:shutil.copyfile(RUN/f,data/f)
        mapping.append([f'Table S{number}',label,'; '.join(files)])
    for sub,label in [('sensitivities','Signal representations networks visits and haemodynamics'),('normalization','Null calibration and reference stability')]:
        subout=data/sub;subout.mkdir(exist_ok=True)
        names=[]
        for f in (RUN/sub).glob('*.csv'):
            if sub=='normalization' and f.name in ['synthetic_calibration_records.csv','scan_metrics_all_seeds_scoped.csv','scan_metrics_all_seeds.csv','paired_changes_all_seeds.csv']:continue
            shutil.copyfile(f,subout/f.name);names.append(f.name)
        mapping.append([f'Table S{len(mapping)+1}',label,sub+'/*.csv'])
    targetout=data/'target_flow';targetout.mkdir(exist_ok=True)
    for f in TARGET.glob('*.csv'):shutil.copyfile(f,targetout/f.name)
    for f in TARGET.glob('*.json'):shutil.copyfile(f,targetout/f.name)
    mapping.append(['Table S10','Corrected Aachen target centred flow','target_flow/*.csv'])
    hdir=data/'haemodynamic_followup';hdir.mkdir(exist_ok=True)
    for f in (RUN/'haemodynamic_followup').iterdir():
        if f.is_file():shutil.copyfile(f,hdir/f.name)
    mapping.append(['Table S11','Whole brain haemodynamic association and influence follow-up','haemodynamic_followup/*.csv'])
    for f in ['raw_scan_metrics.csv','normalized_scan_metrics.csv','raw_paired_changes.csv','normalized_paired_changes.csv','raw_aachen_participant_changes.csv','normalized_aachen_participant_changes.csv','analysis_specification.json','raw_input_sha256.json','statistics_source_sha256.json','extended_checks_receipt.json']:
        shutil.copyfile(RUN/f,data/f)
    for f in ['summary.json','review.json','network_overlap_and_coverage.csv']:
        source=ROOT/'output/control_workspace_qc/aachen_corrected_20261125T234156_788094Z'/f
        shutil.copyfile(source,data/('aachen_control_'+f))
    shutil.copyfile(QC/'provenance.json',data/'control_figure_provenance.json')
    (OUT/'table_inventory.json').write_text(json.dumps(mapping,indent=2),encoding='utf-8')
    (data/'README.txt').write_text('Corrected TMS analyses, 26 November 2026.\nAachen input: reruns/aachen_skullstrip_skip_20261125/pipeline/results. SUDMEX original valid inputs retained.\nAachen participant estimates replace all historical Aachen effects.\nRaw information nats; normalized information CDF units unless explicitly percentile points; counts and participation coefficients keep their native units.\nMain clinical/imaging cohorts must not be pooled. All-post SUDMEX is a longitudinal sensitivity with later open-label visits.\nFamily-wise correction columns name their test family.\n',encoding='utf-8')
    claim_changes={'global_phi':{'historical':{'estimate':.00179217,'p_maxT5':.02393},'corrected':phi.to_dict()},
        'outside_workspace_narrative':'Historical preferential-outside/increase narrative withdrawn; corrected direct contrast uncertain.',
        'control_workspace':'Added unchanged user-supplied figure with measured partial overlap and Luppi/Yousefi citations.',
        'gwt_expectations':'Added four public-literature-motivated tests; none survived.',
        'hrf':'Replaced historical nonsignificant associations with corrected width/latency associations, explicitly diagnostic.',
        'retained':'Acquisition protocol, intended metric definitions and scientific results.',
        'author_information':'Omitted from the current draft at the user request; declarations will be completed later.'}
    (OUT/'revision_claim_audit.json').write_text(json.dumps(claim_changes,indent=2,default=str),encoding='utf-8')
    (OUT/'source_sha256.json').write_text(json.dumps({str(p.relative_to(ROOT)):tms.sha256(p) for p in [Path(__file__),QC/'control_workspace_comparison.png',RUN/'raw_aachen_all40.csv',RUN/'normalized_aachen_all40.csv',RUN/'gwt_expectations_raw.csv',RUN/'gwt_expectations_normalized.csv',RUN/'sensitivities/hrf_phi_associations.csv',*list(OUT.glob('*source.md'))]},indent=2),encoding='utf-8')
    with zipfile.ZipFile(OUT/'TMS_Corrected_Supplementary_Data.zip','w',zipfile.ZIP_DEFLATED) as z:
        for f in sorted(data.rglob('*')):
            if f.is_file():z.write(f,'TMS_Corrected_Supplementary_Data/'+f.relative_to(data).as_posix())
    counts={'abstract_words':len(manuscript.split('## Abstract')[1].split('## Introduction')[0].split()),'body_words':len(manuscript.split('## Introduction')[1].split('## Data availability')[0].split()),'data_files':sum(p.is_file() for p in data.rglob('*'))}
    (OUT/'preparation_receipt.json').write_text(json.dumps(counts,indent=2),encoding='utf-8');print(json.dumps(counts))


if __name__=='__main__':main()
