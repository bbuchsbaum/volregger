# volregger

`volregger` provides retrospective fMRI motion correction with a rigid-first default and optional slice-aware refinement.

## Core API

- `volreg()` main entry point.
- `estimate_motion()` estimation only.
- `apply_motion()` resampling only.
- `motion_qc()`, `framewise_disp()`, `dvars()` for QC.
- `repair_spin_history()` optional post-correction intensity repair.

## Engines

- `engine = "rigid_robust"`: rigid-body estimation with robust templateing, warm-start traversal, and restart seeds.
- `engine = "rigid_spline"`: rigid_robust + spline-smoothed pose traces + slice/packet-time varying application.

## Example

```r
library(volregger)

x <- array(rnorm(32 * 32 * 16 * 50), dim = c(32, 32, 16, 50))
fit <- volreg(
  x,
  engine = "rigid_robust",
  interp_final = "bspline",
  qc = TRUE
)

summary(fit)
plot(fit)
```

## QC interpretation

- `fd`: framewise displacement from motion traces.
- `dvars`, `robust_dvars`: temporal intensity-change metrics.
- `fit_failure`: frame-level fitting warning when optimization does not improve objective.
- `censor_suggest`: conservative downstream censor hint.
- `packet_correction_mag` (rigid_spline): magnitude of spline-vs-rigid adjustment.

## Benchmark harness

Run:

```bash
Rscript inst/benchmarks/benchmark_volregger.R
```

This emits runtime, fixture size, object size, and session information for tracking baseline performance.
