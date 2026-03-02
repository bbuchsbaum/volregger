# QC Guidance

## Per-frame checks

1. Flag frames with high `fd` spikes relative to run median.
2. Cross-check with `dvars` and `robust_dvars` peaks.
3. Inspect `fit_failure` and `censor_suggest` as conservative screening metrics.

## Engine-specific checks

For `rigid_spline`, inspect `packet_correction_mag`:

- near-zero values suggest volume-level rigid poses were already sufficient.
- persistent high values suggest meaningful intra-TR variability and justify packet-time application.

## Spin-history repair checks

`repair_spin_history()` returns:

- `spike_mask`: detected high-confidence slice/time events.
- `gains`: multiplicative correction factors.

Review gain magnitude distribution before committing corrections into downstream modeling.
