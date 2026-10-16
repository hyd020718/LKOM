#' Simulate a small self-contained LKOM example
#'
#' Generates independent-SNP genotype and knockoff matrices together with two
#' residual-derived response channels. This helper is intended for examples,
#' tutorials, and regression tests, not as a replacement for the manuscript's
#' full simulation pipeline.
#'
#' @param n Number of subjects.
#' @param p Number of variant pairs.
#' @param n_beta,n_tau Number of causal variants in each channel.
#' @param overlap Number of causal variants shared by the two channels.
#' @param effect_scale Signal scale.
#' @param seed Random seed used locally.
#' @return A list containing `G`, `G_tilde`, `Y`, and variant-level truth.
#' @export
lkom_simulate_example <- function(
    n = 120L,
    p = 20L,
    n_beta = 4L,
    n_tau = 4L,
    overlap = 2L,
    effect_scale = 0.45,
    seed = 718L) {
  n <- as.integer(n)
  p <- as.integer(p)
  n_beta <- as.integer(n_beta)
  n_tau <- as.integer(n_tau)
  overlap <- as.integer(overlap)
  if (n < 20L || p < 4L || n_beta < 0L || n_tau < 0L ||
      n_beta > p || n_tau > p || overlap < 0L ||
      overlap > min(n_beta, n_tau) || n_beta + n_tau - overlap > p) {
    stop("Invalid example dimensions or causal-overlap specification.",
         call. = FALSE)
  }
  if (!is.finite(effect_scale) || effect_scale <= 0) {
    stop("effect_scale must be positive and finite.", call. = FALSE)
  }

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))

  maf <- stats::runif(p, 0.08, 0.45)
  G <- vapply(maf, function(f) stats::rbinom(n, 2L, f), numeric(n))
  G_tilde <- vapply(maf, function(f) stats::rbinom(n, 2L, f), numeric(n))
  # Extremely unlikely at these dimensions, but make the example deterministic
  # and usable even if a simulated column has zero variance.
  for (j in seq_len(p)) {
    if (stats::var(G[, j]) == 0) G[seq_len(min(3L, n)), j] <- 0:2
    if (stats::var(G_tilde[, j]) == 0) G_tilde[seq_len(min(3L, n)), j] <- 0:2
  }
  snp <- sprintf("rsToy%03d", seq_len(p))
  colnames(G) <- colnames(G_tilde) <- snp
  subject <- sprintf("ID%04d", seq_len(n))
  rownames(G) <- rownames(G_tilde) <- subject

  shared <- if (overlap) seq_len(overlap) else integer()
  beta_only_n <- n_beta - overlap
  tau_only_n <- n_tau - overlap
  beta_only <- if (beta_only_n) overlap + seq_len(beta_only_n) else integer()
  tau_start <- overlap + beta_only_n
  tau_only <- if (tau_only_n) tau_start + seq_len(tau_only_n) else integer()
  beta_index <- c(shared, beta_only)
  tau_index <- c(shared, tau_only)

  beta_effect <- numeric(p)
  tau_effect <- numeric(p)
  if (length(beta_index)) {
    beta_effect[beta_index] <- effect_scale *
      rep(c(-1, 1), length.out = length(beta_index))
  }
  if (length(tau_index)) {
    tau_effect[tau_index] <- effect_scale *
      rep(c(1, -1), length.out = length(tau_index))
  }
  Gs <- scale(G)
  y_beta <- drop(Gs %*% beta_effect + stats::rnorm(n))
  y_tau <- drop(Gs %*% tau_effect + stats::rnorm(n))
  Y <- cbind(beta = as.numeric(scale(y_beta)), tau = as.numeric(scale(y_tau)))
  rownames(Y) <- subject

  state <- rep("00", p)
  state[beta_index] <- "10"
  state[tau_index] <- "01"
  state[intersect(beta_index, tau_index)] <- "11"
  truth <- data.frame(
    SNP = snp,
    state = state,
    causal_beta = seq_len(p) %in% beta_index,
    causal_tau = seq_len(p) %in% tau_index,
    MAF = maf,
    stringsAsFactors = FALSE
  )
  list(G = G, G_tilde = G_tilde, Y = Y, truth = truth, seed = as.integer(seed))
}

