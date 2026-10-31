# AGENTS.md — ProFuse

R package combining ProFound (detection), ProFit (profile modelling) and ProSpect (SED modelling) for physical multi-band galaxy fitting.

## Layout
- `R/profuseSingleImage.R` — `profuseFound2Fit`, `profuseDoFit` (single image).
- `R/profuseMultiImage.R` — multi-image, same band.
- `R/profuseMultiBand.R` — multi-band with ProSpect; `.smooth_parm`, `profuseSegimGlobal`.
- `R/profuseAllStar.R` — star/PSF fitting.
- `R/profuseRegen.R` — rebuild PSF convolvers on loaded objects.
- `man/` — hand-written Rd (no roxygen); update by hand when changing arguments.
- `NAMESPACE` — `exportPattern("^[[:alpha:]]+")`; names starting with `.` are internal.
- `tests/testthat/` — testthat (edition 3) tests; internals reached via `ProFuse:::`.

## Commands
- Install: `R CMD INSTALL .`
- Test: `Rscript -e 'devtools::test()'` or `R CMD check --no-manual .`
- CI: `.github/workflows/main.yml` runs `R CMD check --as-cran` on macOS and Ubuntu (devel/release/oldrel).

## Conventions
- Assignment uses `=`; keep the existing style.
- Dependencies (ProFit, ProFound, ProSpect, etc.) are GitHub-sourced (see `Remotes`); add new ones to both `DESCRIPTION` and CI `extra-packages`.
- Optional deps (e.g. `imager`) must be guarded with `requireNamespace`; tests must use `skip_if_not_installed`.
- Full fits are slow and stochastic: keep tests small, deterministic (`set.seed`) and fast; test helpers and small pure functions in preference to full MCMC runs.
- Bump `Version`/`Date` in `DESCRIPTION` for user-visible changes.
