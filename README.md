# volregger

<!-- badges: start -->
[![R-CMD-check](https://github.com/bbuchsbaum/volregger/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/volregger/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bbuchsbaum/volregger/actions/workflows/pkgdown.yaml/badge.svg)](https://bbuchsbaum.github.io/volregger/)
<!-- badges: end -->

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

Requires R >= 4.1.0 and a working C++ toolchain (Rcpp). Optional suggests
worth installing:

- `neuroim2` — needed for NIfTI (`.nii`, `.nii.gz`) input/output, including
  the CLI.
- `knitr`, `rmarkdown` — needed to build the vignette.

`testthat`, `pkgdown`, `pkgload`, and `RNiftyReg` are used during
development and aren't needed at runtime.

## Quick start

```r
library(volregger)

# 4D fMRI run: x, y, z, time
x <- array(rnorm(32 * 32 * 16 * 50), dim = c(32, 32, 16, 50))

# Bare volreg() applies the `fast_fmri` profile by default.
fit <- volreg(x)

summary(fit)            # snapshot of motion severity and QC
plot(fit)               # six pose traces + framewise displacement
head(coef(fit))         # nt x 6 matrix: tx, ty, tz, rx, ry, rz
fit$corrected           # resampled 4D array
```

A toy fixture is exported for cheap experiments and examples:

```r
fit <- volreg(volregger_example(), profile = "fast_native", qc = FALSE)
```

For more detail (synthetic phantoms with known motion, tuned controls, the
two-step estimate/apply workflow, and the spline engine) see the
`vignette("volregger")` walkthrough.

### Transform-only workflow

If you only need per-volume transforms, skip resampling and keep the motion
estimates:

```r
fit <- volregger::volreg(x, resample = FALSE)

pose6 <- coef(fit)                              # nt x 6: tx ty tz rx ry rz
mats  <- volregger::motion_matrices(fit$mats)   # 4 x 4 x nt
```

`motion_matrices()` returns homogeneous transforms mapping moving world
coordinates into the fixed/reference world coordinate frame.
`estimate_motion()` exposes the same values directly as `est$par` and
`est$mats`.

## Core API

| Function | Purpose |
|---|---|
| `volreg()` | Main entry point: estimate + resample + QC. |
| `estimate_motion()` | Estimation only; returns pose, mask, cost, flags. |
| `apply_motion()` | Resampling only; apply pose parameters to a 4D array. |
| `motion_qc()` | Compute QC metrics from motion + (optionally) corrected data. |
| `framewise_disp()`, `dvars()` | Standard motion / intensity QC vectors. |
| `repair_spin_history()` | Optional intensity repair for through-plane motion. |
| `volreg_control()` | Build a control list; every field is documented. |
| `volreg_validate_control()` | Type-check and range-check a control list. |
| `volreg_profile()` | Named optimizer presets (e.g. `fast_fmri`, `slice_spline`). |
| `volreg_report()` / `write_volreg_report()` | Structured QC report bundle. |
| `volregger_cli()`, `install_cli()` | Programmatic CLI entry / installer. |
| `volregger_example()` | Tiny synthetic 4D fixture for examples and tests. |

S3 methods on `volregger_fit`: `print()`, `summary()`, `plot()`, `coef()`,
`as.data.frame()`. Same on `volregger_report` minus `coef()`.

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

## Profiles

Most users never need to touch `volreg_control()` directly — `volreg_profile()`
exposes named presets that bundle a coherent set of optimizer flags. Pass the
profile name (or a `volreg_profile` object) via `profile = `:

```r
fit <- volreg(x, profile = "fast_fmri")          # default for bare volreg(x)
fit <- volreg(x, profile = "ic_nuisance_adaptive_tuned")  # tuned ablation
fit <- volreg(x, profile = "slice_spline",       # multiband-aware engine
              slice_times = seq(0, 1, length.out = dim(x)[3]))
```

Profiles span a spectrum from `dense_baseline` (no IC stencil, no whitening)
through `fast_native` and `fast_fmri` (production defaults) up to the
`ic_nuisance_*` family used for benchmarking. See `?volreg_profile` for the
full list. If you supply both an explicit `control` and a `profile`, the
profile wins for any field it sets and `volreg()` warns about each conflict
so silent overrides don't masquerade as your settings.

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

`volreg_control()` builds a control list. Every field is a documented
argument with a default, so override what you need by passing it as a named
argument:

```r
ctl <- volreg_control(
  max_iter   = c(60L, 40L, 30L),  # iterations per pyramid level
  n_restarts = 2L,                # extra seeds for hard frames
  nthreads   = 4L                 # parallel resampling
)
fit <- volreg(x, control = ctl, profile = NULL)   # opt out of the default profile
```

`volreg()` and `estimate_motion()` run `volreg_validate_control()` on the
final merged control list before reaching C++, so type / range mistakes
surface as readable R errors:

```r
volreg_validate_control(volreg_control(huber_k = -1), strict = FALSE)
#> [1] "huber_k must be a positive scalar"
```

Frequently tuned fields: `downsample`, `max_iter`, `n_samples`, `huber_k`,
`n_restarts`, `capture_boost`, `use_ic_stencil`, `whiten_modes`, `nthreads`.
See `?volreg_control` for the full list (organized into Pyramid /
Optimization / Capture boost / IC stencil / Polish / Temporal / Low-motion
/ Whitening / Nuisance / Engine sections).

## Command line

`volregger` ships a shell wrapper that exposes the full estimate / apply /
report pipeline. NIfTI input/output (`.nii`, `.nii.gz`) goes through the
optional `neuroim2` package; install it once to enable on-disk fMRI I/O.

Install the wrapper into a directory on your `PATH`:

```r
install.packages("neuroim2")        # for NIfTI I/O
volregger::install_cli("~/.local/bin", overwrite = TRUE)
```

```sh
export PATH="$HOME/.local/bin:$PATH"   # one-time, if needed
volregger --help
volregger help run
```

End-to-end pipeline (estimate + apply + report) from a 4D NIfTI:

```sh
volregger run \
  --input  run.nii.gz \
  --mask   brain.nii.gz \
  --out-dir out/ \
  --prefix  sub01 \
  --threads 4
```

Outputs (under `--out-dir`, named with `--prefix`):

- `<prefix>_corrected.nii.gz` — resampled 4D run
- `<prefix>_fit.rds` — saved `volregger_fit`
- `<prefix>_motion.tsv` — 6-column motion params (tx, ty, tz, rx, ry, rz)
- `<prefix>_summary.csv`, `<prefix>_frames.csv`, `<prefix>_matrices.csv`
- `<prefix>_summary.txt`, `<prefix>_report.rds`, `<prefix>_overview.png`

For matrix-only processing, skip corrected image output but keep the report
bundle, including `<prefix>_matrices.csv`:

```sh
volregger run \
  --input sub-01_task-rest_bold.nii.gz \
  --out-dir out \
  --prefix sub-01_task-rest \
  --no-corrected \
  --no-plot
```

Or estimate first and build the report from the saved fit:

```sh
volregger estimate --input run.nii.gz --out-fit fit.rds --out-motion motion.tsv
volregger report --fit fit.rds --out-dir out --prefix run --no-plot
```

Subcommands:

| Command | Purpose |
|---|---|
| `volregger estimate` | Estimate motion only; write fit RDS (and optional motion TSV). |
| `volregger apply` | Apply a saved fit to a 4D run; write a NIfTI or RDS. |
| `volregger run` | Estimate + apply + report bundle in one call. |
| `volregger report` | Build a report bundle from a saved fit. |

Useful flags (most commands accept these):

- `--json` — emit artifact paths as a single JSON object on stdout.
- `--quiet` — write files but suppress stdout summary.
- `--no-plot` (run, report) — skip the overview PNG.
- `--threads N` — set worker threads.
- `--profile NAME` — `volreg()` preset, e.g. `fast_fmri`.
- `--engine rigid_robust|rigid_spline` — choose engine; with `rigid_spline`
  pass `--slice-times` (and optionally `--mb-groups`).

Exit status: `0` success, `1` domain failure (e.g. missing/invalid input
file), `2` usage or runtime error. Run `volregger help <command>` for a
full per-command listing.

If you don't have the package installed yet, the wrapper also runs from a
source checkout via `pkgload::load_all()`:

```sh
./exec/volregger help
```

## Benchmark harness

A reproducible micro-benchmark is shipped with the package:

```bash
Rscript inst/benchmarks/benchmark_volregger.R
```

It emits runtime, fixture size, object size, and `sessionInfo()` for
tracking baseline performance across changes.

The `inst/benchmarks/benchmark_vs_nifreeze_realdata.R` and friends compare
volregger to [nifreeze](https://github.com/nipreps/nifreeze) on real fMRI
runs. The nifreeze checkout is not vendored in the repo; fetch it on
demand with the included script before running those benchmarks:

```bash
inst/benchmarks/fetch_nifreeze.sh
```

## Documentation

- Vignette: `vignette("volregger")` — narrative walkthrough on a synthetic
  phantom with known injected motion.
- Reference: `?volreg`, `?estimate_motion`, `?apply_motion`,
  `?volreg_control`, `?volreg_validate_control`, `?volreg_profile`,
  `?motion_qc`, `?volreg_report`.

## License

MIT. See `LICENSE` and `DESCRIPTION`.

<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The palette family is provided via `params$family` (default 'violet'). The pkgdown site uses `template: { package: albersdown }`.
<!-- albersdown:theme-note:end -->
