# Multiscale behavior metrics: legacy-definition validation notes

This note documents how `Analysis/03_build_multiscale_behavior_metrics.R` maps onto the older MMMSociability metric definitions.

## Checked legacy sources

- `Functions/E9_SIS_AnimalPos-functions.R`
- `Analysis/E9_SIS_AnimalPos-analyzing.R`
- `Testing/E9_SIS_AnimalPos-analyzing v.2.0.0.R`
- `Analysis/05_build_dyadic_rfid_contacts.R`
- `Analysis/E9_SIS_AnimalPos-analyzing-shannon v.1.0.1.r`

## Position mapping

The multiscale script uses the same 2 x 4 RFID position layout implied by the legacy `find_id()` mapping:

| PositionID | x-grid | y-grid |
|---:|---:|---:|
| 1 | 1 | 1 |
| 2 | 2 | 1 |
| 3 | 3 | 1 |
| 4 | 4 | 1 |
| 5 | 1 | 2 |
| 6 | 2 | 2 |
| 7 | 3 | 2 |
| 8 | 4 | 2 |

This preserves the old eight-position cage representation while allowing Manhattan grid distance to be computed for movement distance and dyadic distance.

## Phase and time handling

The multiscale script assumes that `preprocessed_data/*_preprocessed.csv` already contains the legacy phase logic:

- active phase: 18:30–06:30
- inactive phase: 06:30–18:30
- `ConsecActive` and `ConsecInactive`
- synthetic phase-transition and half-hour rows added by preprocessing

The script therefore does not recompute phase structure from raw files; it consumes the already preprocessed event stream.

## Movement

Legacy AnimalPos scripts calculated movement from position changes between successive RFID position states. The multiscale script preserves this transition-count interpretation and additionally exports a grid-distance variant:

- `Movement`: number of RFID position transitions per animal and bin
- `MovementDistance`: summed Manhattan grid distance moved per animal and bin

`MovementDistance` is a new derived measure and should not be treated as identical to the legacy `Movement` column.

## Proximity

The old AnimalPos proximity output is based on time spent in the same RFID position with cage mates. The dyadic script formalized this as same-position dyadic contact.

For safety, downstream tables should ideally keep both forms:

- legacy-scale proximity seconds: summed same-position seconds across dyads for a focal animal
- normalized proximity fraction: same-position seconds divided by dyadic observation seconds

Recommended naming:

- `Proximity`: legacy-compatible same-position seconds, useful for matching old outputs
- `ProximityFraction`: normalized same-position fraction, preferred for comparing different bin sizes
- `AdjacentProximity`: adjacent-position seconds
- `AdjacentProximityFraction`: normalized adjacent-position fraction
- `MeanGridDistanceToOthers`: duration-weighted mean RFID grid distance to cage mates

## Entropy

The multiscale script calculates animal-level Shannon entropy from position-occupancy seconds within each bin. This corresponds to the intent of the Shannon entropy scripts, but at arbitrary bin sizes.

Formula:

```r
p <- seconds_at_position / sum(seconds_at_position)
Entropy <- -sum(p * log2(p))
```

For the 8-position RFID grid, the theoretical maximum is `log2(8) = 3` if an animal spends equal time in all positions in a bin.

## Current architectural implication

`03_build_multiscale_behavior_metrics.R` should be treated as the producer of canonical downstream files:

```text
analysis_ready/03_derived_metrics/10sec_based/all_behavior_metrics.csv
analysis_ready/03_derived_metrics/1min_based/all_behavior_metrics.csv
analysis_ready/03_derived_metrics/5min_based/all_behavior_metrics.csv
analysis_ready/03_derived_metrics/10min_based/all_behavior_metrics.csv
analysis_ready/03_derived_metrics/30min_based/all_behavior_metrics.csv
analysis_ready/03_derived_metrics/phase_based/all_behavior_metrics.csv
```

For analyses comparing across different temporal resolutions, use normalized columns such as `ProximityFraction`, not raw `Proximity` seconds.

For analyses intended to reproduce earlier phase/half-hour AnimalPos outputs, use `Movement`, `Proximity`, and `Entropy` at the corresponding scale.
