# volregger

> Rigid-first fMRI motion correction in R, with optional acquisition-aware refinement.

`volregger` is an R package for **retrospective motion correction** of 4D fMRI
runs. It estimates six rigid-body parameters per volume (three translations,
three rotations) using a robust, coarse-to-fine optimizer implemented in C++
via Rcpp, resamples the data to a common reference, and reports standard QC
metrics — all from a single `volreg()` call.

A second engine adds slice/packet-timing aware refinement on top of the rigid
estimates, intended for interleaved or multiband acquisitions where a single
pose per volume is too coarse.

## Highlights

- **One-call pipeline**: `volreg()` runs estimation, resampling, and QC end to end.
- **Rigid-first defaults**: robust Huber loss, multi-resolution pyramid, warm-start traversal, optional restart seeds.
- **Acquisition-aware option**: spline-smoothed pose traces with slice/packet-time application via `engine = "rigid_spline"`.
- **Built-in QC**: framewise displacement, DVARS / robust DVARS, fit-failure flags, and a conservative censor suggestion.
- **Composable**: estimate and apply are also exposed separately so you can swap interpolation, share motion across runs, or diagnose intermediates.
- **C++ core**: hot loops are implemented in C++ (Rcpp) and parallelizable via `nthreads`.

## Installation

`volregger` is not on CRAN. Install from a local checkout or directly from a
git remote:

```r
# from a local clone
install.packages("path/to/volregger", repos = NULL, type = "source")

# or with remotes
# remotes::install_github("bbuchsbaum/volregger")
```

Requires R >= 4.1.0 and a working C++ toolchain (Rcpp). Optional suggests:
`RNiftyReg`, `testthat`, `knitr`, `rmarkdown`.

## Quick start

```r
library(volregger)

# 4D fMRI run: x, y, z, time
x <- array(rnorm(32 * 32 * 16 * 50), dim = c(32, 32, 16, 50))

fit <- volreg(
  x,
  engine       = "rigid_robust",
  interp_final = "bspline",
  qc           = TRUE
)

summary(fit)   # snapshot of motion severity and QC
plot(fit)      # six pose traces + framewise displacement
head(coef(fit))         # nt x 6 matrix: tx, ty, tz, rx, ry, rz
fit$corrected           # resampled 4D array
```

For more detail (synthetic phantoms with known motion, tuned controls, the
two-step estimate/apply workflow, and the spline engine) see the
`vignette("volregger")` walkthrough.

## Core API

| Function | Purpose |
|---|---|
| `volreg()` | Main entry point: estimate + resample + QC. |
| `estimate_motion()` | Estimation only; returns pose, mask, cost, flags. |
| `apply_motion()` | Resampling only; apply pose parameters to a 4D array. |
| `motion_qc()` | Compute QC metrics from motion + (optionally) corrected data. |
| `framewise_disp()` | Framewise displacement from a pose matrix. |
| `dvars()` | DVARS (and robust DVARS) from a 4D array. |
| `repair_spin_history()` | Optional intensity repair for through-plane motion. |
| `volreg_control()` | Default control list; override fields and pass to `volreg()`. |

S3 methods on `volregger_fit`: `print()`, `summary()`, `plot()`, `coef()`,
`as.data.frame()`.

## Engines

- **`rigid_robust`** *(default)* — Rigid-body estimation with robust
  templating, coarse-to-fine pyramid, warm-start traversal, and optional
  restart seeds. Fast and appropriate for most data.
- **`rigid_spline`** — `rigid_robust` + spline-smoothed pose traces applied
  with slice/packet-timing awareness. Use this for interleaved or multiband
  acquisitions when you need sub-TR motion estimates. Pass `slice_times` and
  optionally `mb_groups` to `volreg()`.

```r
fit_spline <- volreg(
  x,
  engine      = "rigid_spline",
  slice_times = seq(0, 1, length.out = dim(x)[3])
)
```

## QC interpretation

`fit$qc` (and `as.data.frame(fit)`) include:

- **`fd`** — framewise displacement derived from the motion traces.
- **`dvars`**, **`robust_dvars`** — temporal intensity-change metrics.
- **`fit_failure`** — frame-level flag when the optimizer did not improve the
  objective; useful for diagnosing problem frames.
- **`censor_suggest`** — conservative downstream censor hint (boolean per frame).
- **`packet_correction_mag`** *(rigid_spline only)* — magnitude of the
  spline-vs-rigid pose adjustment per frame.

A common workflow is to threshold `fd` (e.g. > 0.5 mm) and union with
`censor_suggest` to produce a scrubbing mask for downstream analysis.

## Tuning

`volreg_control()` returns the default optimizer settings. Override any field
and pass back via `control = `:

```r
fit <- volreg(
  x,
  control = modifyList(
    volreg_control(),
    list(
      max_iter   = c(60L, 40L, 30L),  # iterations per pyramid level
      n_restarts = 2L,                # extra seeds for hard frames
      nthreads   = 4L                 # parallel resampling
    )
  )
)
```

Frequently tuned fields: `downsample`, `max_iter`, `n_samples`, `huber_k`,
`n_restarts`, `capture_boost`, `use_ic_stencil`, `whiten_modes`, `nthreads`.
See `?volreg_control` for the full list.

## Benchmark harness

A reproducible micro-benchmark is shipped with the package:

```bash
Rscript inst/benchmarks/benchmark_volregger.R
```

It emits runtime, fixture size, object size, and `sessionInfo()` for tracking
baseline performance across changes.

## Documentation

- Vignette: `vignette("volregger")` — narrative walkthrough on a synthetic
  phantom with known injected motion.
- Reference: `?volreg`, `?estimate_motion`, `?apply_motion`,
  `?volreg_control`, `?motion_qc`.

## License

MIT. See `LICENSE` and `DESCRIPTION`.

<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The palette family is provided via `params$family` (default 'violet'). The pkgdown site uses `template: { package: albersdown }`.
<!-- albersdown:theme-note:end -->
