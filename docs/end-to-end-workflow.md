# End-to-end LKOM analysis from user data

This guide starts with a longitudinal phenotype, covariates, an original
genotype matrix, and a matching model-X knockoff matrix, and ends with the
variant-level LKOM statistics, knockoff+ discoveries, and association-state
classifications.

Swap-symmetric marginal screening is strongly recommended before
masked-likelihood-ratio inference. This step is computationally important: it
removes pairs with little marginal evidence and can greatly reduce the
dimension and runtime of the subsequent MLR model. Analyses should generally
omit screening only for toy examples or very small diagnostic datasets.

## 1. What LKOM expects

LKOM uses two linked data representations.

1. **Long-format phenotype data**, with one row per visit, are used to fit the
   WiSER null model. The required objects are:

   - `id`: subject identifier for every visit;
   - `y`: longitudinal phenotype for every visit;
   - `X`: fixed-effect covariates for the trajectory mean;
   - `Z`: random-effect covariates, such as visit time for a random slope;
   - `W`: fixed-effect covariates for within-subject variability.

2. **Subject-level genotype data**, with one row per subject, are used by
   LKOM:

   - `G`: original genotype matrix;
   - `G_tilde`: paired model-X knockoff genotype matrix.

`G` and `G_tilde` must have identical dimensions and identically ordered SNP
columns. Their row order must ultimately match `null_fit$subject_id`. LKOM does
not construct model-X knockoffs; users must generate a valid knockoff matrix
before fitting the model.

The Julia helper adds intercepts internally. Passing `X`, `Z`, and `W` without
an intercept column is therefore recommended. If a column consisting entirely
of ones is supplied, `lkom_fit_null()` removes it before calling Julia.

## 2. Install LKOM and initialize Julia

```r
remotes::install_github("hyd020718/LKOM")
library(LKOM)

julia_dir <- system.file("julia", package = "LKOM")
system2(
  "julia",
  c(
    paste0("--project=", shQuote(julia_dir)),
    "-e",
    shQuote("using Pkg; Pkg.instantiate(); Pkg.precompile()")
  )
)
```

The Julia environment normally needs to be instantiated only once for a given
Julia version, package environment, and depot.

## 3. Prepare the longitudinal phenotype and covariates

Suppose `long_data` contains one row per visit and the columns `subject_id`,
`phenotype`, `age`, `sex`, `PC1`, ..., `PC10`, and `visit_time`.

```r
id <- as.character(long_data$subject_id)
y  <- as.numeric(long_data$phenotype)

# Mean-model fixed effects. Remove the model.matrix intercept because the
# Julia helper adds an intercept internally.
X <- model.matrix(
  ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 +
    PC6 + PC7 + PC8 + PC9 + PC10,
  data = long_data
)[, -1, drop = FALSE]

# Example random-slope design. A random intercept is added internally.
Z <- matrix(long_data$visit_time, ncol = 1L)
colnames(Z) <- "visit_time"

# Choose the WSV covariates required by the scientific analysis.
W <- model.matrix(
  ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 +
    PC6 + PC7 + PC8 + PC9 + PC10,
  data = long_data
)[, -1, drop = FALSE]

stopifnot(
  length(id) == length(y),
  nrow(X) == length(y),
  nrow(Z) == length(y),
  nrow(W) == length(y),
  all(is.finite(y)),
  !anyNA(id)
)
```

The covariate formulas above are examples rather than required formulas.
Covariate coding should follow the analysis plan for the phenotype under
study. Categorical variables should be converted to numeric design columns,
for example with `model.matrix()`.

## 4. Fit the WiSER null model

```r
null_fit <- lkom_fit_null(
  X = X,
  Z = Z,
  W = W,
  Y = y,
  id = id,
  julia_bin = "julia",
  runs = 10L,
  maxeval = 4000L,
  work_dir = "lkom_wiser_work",
  keep_work_dir = TRUE,
  verbose = TRUE
)
```

Keeping `work_dir` is useful for inspecting `julia_wiser.log` and for
reproducibility. A successful fit contains one record per subject in the order
given by:

```r
subject_order <- as.character(null_fit$subject_id)
length(subject_order)
```

Construct the two subject-level LKOM response channels:

```r
Y_lkom <- lkom_residual_channels(null_fit, transform = TRUE)
stopifnot(identical(rownames(Y_lkom), subject_order))
```

The two columns of `Y_lkom` represent the trajectory-mean (`beta`) and
within-subject-variability (`tau`) channels. `transform = TRUE` applies the
rank-based inverse-normal transformation used in the main LKOM analyses.

## 5. Align original and knockoff genotypes

Both genotype matrices must be numeric subject-by-SNP matrices. Assign stable,
unique subject and SNP identifiers before alignment.

```r
# Example input objects:
# rownames(G) and rownames(G_tilde): subject identifiers
# colnames(G) and colnames(G_tilde): SNP identifiers

stopifnot(
  identical(dim(G), dim(G_tilde)),
  !is.null(rownames(G)),
  !is.null(rownames(G_tilde)),
  !is.null(colnames(G)),
  !is.null(colnames(G_tilde)),
  !anyDuplicated(rownames(G)),
  !anyDuplicated(rownames(G_tilde)),
  !anyDuplicated(colnames(G)),
  identical(colnames(G), colnames(G_tilde))
)

index_G  <- match(subject_order, rownames(G))
index_ko <- match(subject_order, rownames(G_tilde))
if (anyNA(index_G) || anyNA(index_ko)) {
  stop("Some WiSER subjects are absent from G or G_tilde.")
}

G       <- G[index_G, , drop = FALSE]
G_tilde <- G_tilde[index_ko, , drop = FALSE]

rownames(G) <- rownames(G_tilde) <- subject_order
stopifnot(
  identical(rownames(G), rownames(Y_lkom)),
  identical(rownames(G_tilde), rownames(Y_lkom)),
  identical(colnames(G), colnames(G_tilde)),
  all(is.finite(G)),
  all(is.finite(G_tilde))
)
```

This alignment step is essential. Equal matrix dimensions alone do not prove
that the same subjects or variants are in the same order.

## 6. Perform swap-symmetric marginal screening

Marginal screening is strongly recommended for practical LKOM analyses,
particularly at genome-wide scale. Compute marginal p-values for both the
original and knockoff genotype of every SNP, in both response channels, using
the same null model and the same numerical procedure. Arrange them in the SNP
order of `G`:

```r
p_original <- cbind(
  beta = p_beta[colnames(G)],
  tau  = p_tau[colnames(G)]
)
p_knockoff <- cbind(
  beta = p_tilde_beta[colnames(G)],
  tau  = p_tilde_tau[colnames(G)]
)

screen <- lkom_screen_pairs(
  p_original = p_original,
  p_knockoff = p_knockoff,
  threshold = 1e-3
)

table(screen)
```

The rule retains a pair when at least one of its four p-values is below the
threshold. Because the rule is unchanged when the original and knockoff labels
are exchanged, it preserves the swap symmetry required by the knockoff
procedure. The resulting reduction in the number of pairs entering MLR
inference is the main computational benefit of this step.

Marginal p-value calculation is intentionally external to the core package.
The exact implementation used for the manuscript was based on the TrajGWAS
SPA scan. SPA-GRM is a comparator in the related-sample simulations and is not
part of the LKOM fitting API.

For transparency, `lkom()` can technically fit all pairs by setting
`screen = NULL`. This should be reserved for toy examples, unit tests, or very
small diagnostic datasets; it is not the recommended workflow for substantive
or genome-wide analyses.

## 7. Fit LKOM and apply knockoff+

The `paper` preset uses the manuscript configuration: a bidirectional
conditional-Bayes-factor warm start, 100 burn-in sweeps, and 300 retained
Gibbs sweeps.

```r
control <- lkom_control(
  method = "gibbs",
  preset = "paper",
  seed = 718L
)

fit <- lkom(
  G = G,
  G_tilde = G_tilde,
  Y = Y_lkom,
  q = 0.05,
  offset = 1L,
  control = control,
  channel_names = c("beta", "tau"),
  screen = screen
)

print(fit)
summary(fit)
```

`offset = 1L` applies knockoff+, which is the procedure used for FDR control in
the manuscript. `q` is the target variant-level FDR for the joint null of no
association with either channel.

## 8. Obtain the final results

Return one row per variant:

```r
results_all <- lkom_select(fit)
results_selected <- lkom_select(fit, selected_only = TRUE)

head(results_selected)
write.csv(results_all, "lkom_results_all.csv", row.names = FALSE)
write.csv(results_selected, "lkom_discoveries.csv", row.names = FALSE)
saveRDS(fit, "lkom_fit.rds")
```

Important output columns are:

- `W`: signed LKOM knockoff statistic;
- `screened`: whether the pair passed marginal screening;
- `selected`: whether the pair passed the knockoff+ threshold;
- `state`: maximum-posterior association state;
- `posterior_original`: posterior probability favoring the observed original
  orientation;
- posterior probabilities for states `00`, `10`, `01`, and `11`.

Extract the joint and channel-specific discovery sets:

```r
discoveries <- lkom_discoveries(fit)

discoveries$joint  # selected against the joint null
discoveries$beta   # selected and assigned to state 10 or 11
discoveries$tau    # selected and assigned to state 01 or 11
```

State interpretation:

| State | Interpretation |
|---|---|
| `00` | no detected association with either channel |
| `10` | trajectory-mean association |
| `01` | WS-variability association |
| `11` | association with both channels |

## 9. Related samples

For related subjects, construct a GRM structure and use
`lkom_fit_related()`. The function estimates the whitening transformation and
applies the same transformation to the original genotype, knockoff genotype,
and response channels.

```r
structure <- lkom_grm_structure(
  grm = "sparse_grm.tsv",
  subject_ids = subject_order
)

fit_related <- lkom_fit_related(
  G = G,
  G_tilde = G_tilde,
  Y = Y_lkom,
  structure = structure,
  control = control,
  channel_names = c("beta", "tau"),
  screen = screen
)

related_results <- lkom_select(fit_related, q = 0.05, offset = 1L)
related_discoveries <- lkom_discoveries(fit_related, q = 0.05, offset = 1L)
```

For a dense GRM, pass the matrix itself to `lkom_grm_structure()`. For biobank
data, a sparse three-column GRM is recommended.

## 10. Minimum reproducibility record

For each analysis, retain at least:

- the LKOM package version and R `sessionInfo()`;
- the Julia `Project.toml` and `Manifest.toml` versions;
- sample and SNP order files;
- the marginal screening threshold and retained-pair vector;
- `lkom_control()` settings and random seed;
- the target `q` and knockoff offset;
- the fitted `lkom_fit` object and final variant-level table.

```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
saveRDS(
  list(
    subject_order = subject_order,
    snp_order = colnames(G),
    screen = screen,
    control = control,
    q = 0.05,
    offset = 1L
  ),
  "lkom_analysis_metadata.rds"
)
```

## 11. Common errors

- **`G and G_tilde must have identical dimensions`**: the two matrices do not
  contain the same number of subjects and SNPs.
- **Genotype column-name error**: original and knockoff SNP columns are not in
  exactly the same order.
- **Julia package not found**: instantiate the bundled Julia environment and
  ensure that the same project is used when invoking the bundled helper.
- **No discoveries**: `knockoff+` requires enough positive statistics to
  satisfy its data-dependent threshold; an infinite threshold is valid and
  means no variants were selected.
- **Large memory use**: apply swap-symmetric marginal screening externally and
  fit LKOM only to retained pairs.
