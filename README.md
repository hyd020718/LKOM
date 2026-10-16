# LKOM

**LKOM** (Longitudinal Knockoff with Masked likelihood ratio) jointly analyzes
genetic effects on the trajectory mean and within-subject (WS) variability of
a longitudinal phenotype. It uses paired model-X knockoffs for variant-level
FDR control under the joint null and assigns selected variants to the states
`10`, `01`, and `11`.

This repository is an initial research release (`0.1.0`).

## Core workflow

```r
library(LKOM)

toy <- lkom_simulate_example(seed = 718)

# A fast deterministic fit for a tutorial. The manuscript analyses used the
# Gibbs configuration returned by lkom_control() without these overrides.
fit <- lkom(
  G = toy$G,
  G_tilde = toy$G_tilde,
  Y = toy$Y,
  q = 0.10,
  control = lkom_control(method = "conditional_bf", preset = "fast")
)

print(fit)
lkom_select(fit, selected_only = TRUE)
lkom_discoveries(fit)
```

For an analysis starting from a user's longitudinal phenotype, covariates,
original genotypes, and precomputed knockoff genotypes, see the
**[end-to-end workflow](docs/end-to-end-workflow.md)**. It covers input
formats, subject and variant alignment, WiSER null-model fitting, construction
of the two LKOM response channels, swap-symmetric marginal screening,
LKOM fitting, knockoff+ selection, and extraction of joint and channel-specific
discoveries.

The default `lkom_control()` reproduces the main analysis settings: a
conditional-Bayes-factor warm start followed by 100 burn-in and 300 retained
Gibbs sweeps, with inverse-Wishart updating of the two-channel residual
covariance and grid updating of the slab parameters.

## Symmetric marginal screening

Marginal screening is a strongly recommended computational step in the LKOM
workflow. It reduces the number of variant pairs entering the masked-
likelihood-ratio model and can substantially reduce runtime in genome-wide
analyses. To preserve knockoff validity, screening must treat the original and
knockoff member of every pair symmetrically. Given marginal p-values from
TrajGWAS or another paired scan:

```r
keep <- lkom_screen_pairs(
  p_original = cbind(p_beta, p_tau),
  p_knockoff = cbind(p_tilde_beta, p_tilde_tau),
  threshold = 1e-3
)

fit <- lkom_fit(toy$G, toy$G_tilde, toy$Y, screen = keep)
```

Pairs excluded by screening receive `W = 0`.

Fitting every pair without screening is intended only for toy examples and
small diagnostic analyses.

## Related samples

For small examples, a dense GRM can be supplied directly. For biobank data,
use a sparse three-column GRM file so that LKOM works component by component.

```r
K <- diag(nrow(toy$G))
structure <- lkom_grm_structure(K, subject_ids = rownames(toy$G))

fit_related <- lkom_fit_related(
  toy$G,
  toy$G_tilde,
  toy$Y,
  structure = structure,
  control = lkom_control(method = "conditional_bf", preset = "fast")
)
```

The same GLS-centering and whitening transformation is applied to each member
of an original-knockoff pair, so paired swaps commute with preprocessing.

## WiSER and Julia environment

The exact Julia environment used during development is bundled under
`inst/julia`. After installing the R package, instantiate it once:

```r
julia_dir <- system.file("julia", package = "LKOM")
system2("julia", c(
  paste0("--project=", shQuote(julia_dir)),
  "-e", shQuote("using Pkg; Pkg.instantiate()")
))
```

The locked environment was generated with Julia 1.11.7 and includes
TrajGWAS 0.4.7 and WiSER 0.2.6. Use `lkom_check_dependencies()` to inspect
optional software available on the current system.

## Scope of this package

The package contains the reusable LKOM statistical method. Cohort extraction,
PLINK quality control, knockoff construction for a specific genotype resource,
SLURM submission, and manuscript-wide simulation aggregation remain separate
reproducibility workflows because they depend on local data and cluster paths.

