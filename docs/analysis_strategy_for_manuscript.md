# MMMSociability Manuscript Analysis Strategy

## Central Biological Hypothesis

Early behavioral adaptation after the first social instability exposure predicts later stress burden and aligns with hippocampal molecular adaptation programs.

The primary biological claim is not that resilience is simply more or less movement. The working model is that resilience reflects organized adaptation after social perturbation: controlled locomotor output, altered temporal structure, preserved or restored active/inactive organization, social reorganization, and recovery or stabilization dynamics.

## Core Manuscript Domains

1. Early behavioral adaptation
2. Temporal instability and organization
3. Ethological active/inactive phase organization
4. Social reorganization after regrouping
5. Behavioral-proteomic systems alignment

HMM, state-space, manifold, recurrence, energy-landscape and other nonlinear analyses are secondary or exploratory unless they directly support one of these domains.

## Primary Prediction Window

The primary prospective behavioral window is:

- first cage change
- first active phase
- first 12 h
- 5-min bins

The primary endpoint is `CombZ` or an equivalent later composite stress-burden metric.

The primary feature family is:

- `Movement_mean`
- `Movement_rmssd`
- `Entropy_acf1`

This family intentionally separates locomotor magnitude, short-timescale volatility, and temporal persistence of spatial organization.

## Prediction Strategy

`Analysis/08b_early_prediction_model_ladder.R` is the primary prediction script. The manuscript emphasis should be behavior-only prediction of later stress burden.

Group-aware models are useful sensitivity analyses, but if RES/SUS labels are derived from `CombZ`, they are not independent prospective predictors. They should not be used as the primary evidence for early behavior predicting later stress burden.

Primary reporting should include:

- LOOCV model ladder
- repeated grouped k-fold CV
- full-pipeline permutation testing
- bootstrap confidence intervals
- behavior-only model ladder
- behavior + group/sex sensitivity ladder
- duration sensitivity excluding short-duration epochs

## Duration / CC4 Robustness

Cage change 4 has shorter observation duration than other cage changes. The goal is not to remove CC4 biologically, but to prevent shorter observation duration from being interpreted as a behavioral phenotype.

All duration-sensitive analyses should report full-data and excluding-short-duration outputs. Main-text claims should be considered stable only when:

- effect direction is unchanged
- `abs(delta_cohen_d) < 0.30`

Raw counts, cumulative AUCs, path lengths and transition counts should not be compared across unequal duration without per-hour or per-epoch normalization.

## Phase Organization Interpretation

Mice are nocturnal, so active/inactive organization is ethologically meaningful. However, RFID home-cage metrics do not independently validate circadian disruption or sleep disruption.

Preferred language:

- active/inactive organization
- day/night behavioral structure
- phase-specific adaptation
- dark-phase adaptation
- light-phase rest-like organization

Avoid:

- circadian disruption
- sleep disruption

unless independently validated.

## Social Instability Interpretation

Social analyses should separate:

- proximity/co-occupancy proxies
- true dyadic networks
- partner instability
- preferred-partner persistence
- social fragmentation
- phase-specific social organization

If dyadic identity is unavailable, avoid strong graph-theory claims and frame outputs as animal-level social engagement or withdrawal proxies.

## Inactivity Caveat

`Analysis/16_sleep_like_inactivity_metrics.R` estimates sleep-like or rest-like inactivity structure from movement continuity and bouts.

Do not claim EEG-validated sleep. Use:

- sleep-like inactivity
- rest-like inactivity
- quiescence fragmentation
- inactivity fragmentation

## Proteomics Integration Logic

`Analysis/12_behavior_proteomics_integration.R` should avoid feature-explosion correlations as the primary story. The preferred analysis compares low-dimensional behavioral adaptation axes with curated hippocampal proteomic adaptation modules:

- locomotor/adaptation axis
- temporal organization axis
- phase organization axis
- social organization axis
- inactivity/quiescence axis

Proteomic axes:

- RNA/RNP/splicing
- translation/ribosome
- mitochondrial/OXPHOS
- proteostasis/endolysosomal
- synaptic/plasticity

Small-n behavior-proteomics findings should be labeled exploratory and effect-size focused.

## Evidence Tier Framework

Tier 1:
Primary prediction/adaptation findings from the first active phase after the first cage change.

Tier 2:
Mechanistic decomposition including instability, HMM/state switching, phase organization, social organization, adaptation kinetics and inactivity structure.

Tier 3:
Exploratory nonlinear/manifold analyses including recurrence maps, attractor depth, energy landscapes and heavy manifold optimization.

## Recommended Run Order

1. `Analysis/03_build_multiscale_behavior_metrics.R`
2. `Analysis/06_burstiness_temporal_instability.R`
3. `Analysis/08b_early_prediction_model_ladder.R`
4. `Analysis/11_gamm_trajectory_features.R`
5. `Analysis/15_behavioral_adaptation_kinetics.R`
6. `Analysis/17_ethological_phase_organization.R`
7. `Analysis/09_dynamic_social_networks.R`
8. `Analysis/10_hmm_behavioral_states.R`
9. `Analysis/07_behavioral_state_space.R`
10. `Analysis/16_sleep_like_inactivity_metrics.R`
11. `Analysis/12_behavior_proteomics_integration.R`
12. `Analysis/12_systems_neuroscience_summary.R`
13. `Analysis/14_nextgen_behavioral_phenotyping.R`
