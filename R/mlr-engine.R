# KOALA longitudinal masked-likelihood-ratio (MLR) statistics.
#
# This file adds three interchangeable inference engines on top of the
# existing KOALA/WiSER null-model and TrajGWAS-SPA preparation functions:
#   1) "gibbs"          : latent orientation/state Gibbs sampler;
#   2) "conditional_bf": coordinate conditional Bayes factors (MAP residual);
#   3) "variational"    : damped mean-field variational approximation.
#
# Required functions supplied by koala_wiser_three_methods.R:
#   .koala_rank_int(), koala_prepare_block(), koala_make_block_groups().

.koala_mlr_logsumexp <- function(x) {
  x <- as.numeric(x)
  xmax <- max(x)
  if (!is.finite(xmax)) return(xmax)
  xmax + log(sum(exp(x - xmax)))
}

.koala_mlr_softmax <- function(log_weight) {
  z <- .koala_mlr_logsumexp(log_weight)
  out <- exp(log_weight - z)
  out / sum(out)
}

.koala_mlr_logdet_spd <- function(x, jitter = 1e-10) {
  x <- (as.matrix(x) + t(as.matrix(x))) / 2
  for (k in 0:8) {
    add <- if (k == 0L) 0 else jitter * 10^(k - 1L)
    ch <- try(chol(x + diag(add, nrow(x))), silent = TRUE)
    if (!inherits(ch, "try-error")) {
      return(list(
        logdet = 2 * sum(log(diag(ch))),
        chol = ch,
        matrix = x + diag(add, nrow(x)),
        jitter = add
      ))
    }
  }
  stop("Matrix remains non-positive-definite after numerical ridge adjustment.")
}

.koala_mlr_rank_int <- function(x) {
  if (exists(".koala_rank_int", mode = "function")) {
    return(.koala_rank_int(x))
  }
  x <- as.numeric(x)
  n <- sum(is.finite(x))
  if (n < 2L || any(!is.finite(x))) stop("INT input must contain at least two finite values.")
  stats::qnorm((rank(x, ties.method = "average") - 0.5) / n)
}

.koala_mlr_scale_matrix <- function(x, name, center = TRUE) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (anyNA(x) || any(!is.finite(x))) stop(name, " contains non-finite values.")
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("center must be TRUE or FALSE.")
  }
  location <- if (center) colMeans(x) else rep(0, ncol(x))
  xc <- sweep(x, 2L, location, "-")
  scale <- sqrt(colSums(xc^2) / nrow(xc))
  if (any(!is.finite(scale)) || any(scale <= sqrt(.Machine$double.eps))) {
    bad <- which(!is.finite(scale) | scale <= sqrt(.Machine$double.eps))
    stop(name, " contains zero-variance columns: ", paste(bad, collapse = ","))
  }
  list(
    x = sweep(xc, 2L, scale, "/"),
    center = location,
    scale = scale
  )
}

# Canonical ordering depends only on the unordered pair {g_j, gtilde_j}.
# Therefore swapping the labels cannot alter the data seen by any MLR engine;
# it only reverses observed_orientation and hence the sign of W.
.koala_mlr_canonical_pairs <- function(G, G_tilde, tol = 0) {
  G <- as.matrix(G)
  G_tilde <- as.matrix(G_tilde)
  storage.mode(G) <- "double"
  storage.mode(G_tilde) <- "double"
  if (!identical(dim(G), dim(G_tilde))) stop("G and G_tilde have different dimensions.")
  if (anyNA(G) || anyNA(G_tilde) ||
      any(!is.finite(G)) || any(!is.finite(G_tilde))) {
    stop("G or G_tilde contains non-finite values.")
  }

  p <- ncol(G)
  A <- G
  B <- G_tilde
  observed_orientation <- rep.int(1, p)
  identical_pair <- logical(p)

  for (j in seq_len(p)) {
    d <- G[, j] - G_tilde[, j]
    idx <- which(abs(d) > tol)
    if (!length(idx)) {
      identical_pair[j] <- TRUE
      next
    }
    # Lexicographically smaller vector is canonical A.
    if (G[idx[1L], j] > G_tilde[idx[1L], j]) {
      A[, j] <- G_tilde[, j]
      B[, j] <- G[, j]
      observed_orientation[j] <- -1
    }
  }

  # Direct internal swap check, before any model fitting.
  observed_orientation_swap <- -observed_orientation
  observed_orientation_swap[identical_pair] <- observed_orientation[identical_pair]
  list(
    A = A,
    B = B,
    observed_orientation = observed_orientation,
    identical_pair = identical_pair,
    swap_check = list(
      canonical_data_identical = TRUE,
      orientation_flipped = all(
        observed_orientation_swap[!identical_pair] ==
          -observed_orientation[!identical_pair]
      )
    )
  )
}

.koala_mlr_state_spec <- function(k, effect_sd, effect_rho = 0) {
  effect_sd <- rep(as.numeric(effect_sd), length.out = k)
  if (any(!is.finite(effect_sd)) || any(effect_sd <= 0)) {
    stop("MLR slab effect_sd values must be finite and positive.")
  }
  if (k == 1L) {
    return(list(
      names = c("0", "1"),
      active = list(integer(0), 1L),
      psi = list(NULL, matrix(effect_sd[1L]^2, 1L, 1L))
    ))
  }
  if (k != 2L) stop("The MLR engine supports one or two response channels.")
  if (!is.finite(effect_rho) || abs(effect_rho) >= 1) {
    stop("effect_rho must lie strictly between -1 and 1.")
  }
  psi11 <- matrix(
    c(
      effect_sd[1L]^2,
      effect_rho * effect_sd[1L] * effect_sd[2L],
      effect_rho * effect_sd[1L] * effect_sd[2L],
      effect_sd[2L]^2
    ),
    2L,
    2L
  )
  list(
    names = c("00", "10", "01", "11"),
    active = list(integer(0), 1L, 2L, 1:2),
    psi = list(
      NULL,
      matrix(effect_sd[1L]^2, 1L, 1L),
      matrix(effect_sd[2L]^2, 1L, 1L),
      psi11
    )
  )
}

.koala_mlr_normalize_pi <- function(pi, state_names) {
  pi <- as.numeric(pi)
  if (length(pi) != length(state_names) ||
      any(!is.finite(pi)) || any(pi <= 0)) {
    stop("The MLR state prior has the wrong length or invalid values.")
  }
  pi <- pi / sum(pi)
  names(pi) <- state_names
  pi
}

# Integrated likelihood for R = x b' + E, E_i ~ N(0, Omega),
# with a Gaussian slab on the active coordinates of b.
.koala_mlr_component <- function(u, q, omega_inv, active, psi, sample = FALSE) {
  k <- length(u)
  if (!length(active)) {
    return(list(log_bf = 0, mean = numeric(k), covariance = NULL,
                draw = numeric(k)))
  }
  psi_info <- .koala_mlr_logdet_spd(psi)
  psi_inv <- chol2inv(psi_info$chol)
  J <- q * omega_inv[active, active, drop = FALSE] + psi_inv
  J_info <- .koala_mlr_logdet_spd(J)
  V <- chol2inv(J_info$chol)
  h <- drop((omega_inv %*% u)[active])
  mu <- drop(V %*% h)
  log_bf <- -0.5 * (psi_info$logdet + J_info$logdet) +
    0.5 * sum(h * mu)
  mean_full <- numeric(k)
  mean_full[active] <- mu
  draw_full <- mean_full
  if (sample) {
    z <- stats::rnorm(length(active))
    draw_full[active] <- mu + drop(t(chol(V)) %*% z)
  }
  list(log_bf = log_bf, mean = mean_full, covariance = V,
       draw = draw_full)
}

.koala_mlr_pair_posterior <- function(
    j, C, XtY, Gram, omega_inv, state_spec, pi_state,
    sample_effect = FALSE) {

  p <- nrow(C) / 2L
  rows <- c(j, p + j)
  C_pair <- C[rows, , drop = FALSE]
  # Scores against the residual excluding the complete j-th pair.
  u_pair <- XtY[rows, , drop = FALSE] -
    Gram[rows, , drop = FALSE] %*% C +
    Gram[rows, rows, drop = FALSE] %*% C_pair

  n_state <- length(state_spec$names)
  components <- vector("list", 2L * n_state)
  log_weight <- numeric(2L * n_state)
  index <- 0L
  for (s in 1:2) {
    q <- Gram[rows[s], rows[s]]
    for (h in seq_len(n_state)) {
      index <- index + 1L
      comp <- .koala_mlr_component(
        u = drop(u_pair[s, ]),
        q = q,
        omega_inv = omega_inv,
        active = state_spec$active[[h]],
        psi = state_spec$psi[[h]],
        sample = sample_effect
      )
      components[[index]] <- comp
      log_weight[index] <- log(0.5) + log(pi_state[h]) + comp$log_bf
    }
  }
  probability <- .koala_mlr_softmax(log_weight)
  orientation_index <- rep(1:2, each = n_state)
  state_index <- rep(seq_len(n_state), times = 2L)
  orientation_probability <- c(
    sum(probability[orientation_index == 1L]),
    sum(probability[orientation_index == 2L])
  )
  state_probability <- vapply(seq_len(n_state), function(h) {
    sum(probability[state_index == h])
  }, numeric(1))
  names(state_probability) <- state_spec$names
  log_evidence <- vapply(1:2, function(s) {
    ids <- which(orientation_index == s)
    .koala_mlr_logsumexp(log_weight[ids] - log(0.5))
  }, numeric(1))

  list(
    rows = rows,
    components = components,
    log_weight = log_weight,
    probability = probability,
    orientation_index = orientation_index,
    state_index = state_index,
    orientation_probability = orientation_probability,
    state_probability = state_probability,
    canonical_log_odds = log_evidence[1L] - log_evidence[2L]
  )
}

.koala_mlr_expected_C <- function(post, k) {
  out <- matrix(0, 2L, k)
  for (r in seq_along(post$probability)) {
    s <- post$orientation_index[r]
    out[s, ] <- out[s, ] +
      post$probability[r] * post$components[[r]]$mean
  }
  out
}

.koala_mlr_final_posteriors <- function(
    C, XtY, Gram, omega_inv, state_spec, pi_state, epsilon) {
  p <- nrow(C) / 2L
  n_state <- length(state_spec$names)
  pA <- numeric(p)
  state_prob <- matrix(0, p, n_state,
                       dimnames = list(NULL, state_spec$names))
  log_odds <- numeric(p)
  for (j in seq_len(p)) {
    post <- .koala_mlr_pair_posterior(
      j, C, XtY, Gram, omega_inv, state_spec, pi_state
    )
    pA[j] <- post$orientation_probability[1L]
    state_prob[j, ] <- post$state_probability
    log_odds[j] <- post$canonical_log_odds
  }
  pA <- pmin(pmax(pA, epsilon), 1 - epsilon)
  # Recompute from pA so Gibbs and deterministic backends share clipping.
  log_odds <- stats::qlogis(pA)
  list(pA = pA, state_probability = state_prob, log_odds = log_odds)
}

.koala_mlr_fit_conditional_bf <- function(
    XtY, Gram, omega_inv, state_spec, pi_state, control) {
  p <- nrow(Gram) / 2L
  k <- ncol(XtY)
  C <- matrix(0, 2L * p, k)
  selected_state <- rep.int(1L, p)
  converged <- FALSE
  delta <- Inf
  iterations <- 0L

  for (iter in seq_len(control$max_iter)) {
    pi_used <- pi_state
    delta <- 0
    state_sum <- numeric(length(state_spec$names))
    for (j in seq_len(p)) {
      post <- .koala_mlr_pair_posterior(
        j, C, XtY, Gram, omega_inv, state_spec, pi_used
      )
      chosen <- which.max(post$log_weight)
      C_new <- matrix(0, 2L, k)
      s <- post$orientation_index[chosen]
      C_new[s, ] <- post$components[[chosen]]$mean
      selected_state[j] <- post$state_index[chosen]
      old <- C[post$rows, , drop = FALSE]
      delta <- max(delta, max(abs(C_new - old)))
      C[post$rows, ] <- C_new
      state_sum <- state_sum + post$state_probability
    }
    if (isTRUE(control$learn_pi)) {
      pi_state <- state_sum + control$pi_prior_strength * control$pi_initial
      pi_state <- pi_state / sum(pi_state)
      names(pi_state) <- state_spec$names
    }
    iterations <- iter
    if (delta <= control$tol) {
      converged <- TRUE
      break
    }
  }
  final <- .koala_mlr_final_posteriors(
    C, XtY, Gram, omega_inv, state_spec, pi_state, control$epsilon
  )
  c(final, list(C = C, selected_state = selected_state,
                pi_state = pi_state, converged = converged,
                iterations = iterations, final_delta = delta))
}

.koala_mlr_fit_variational <- function(
    XtY, Gram, omega_inv, state_spec, pi_state, control) {
  p <- nrow(Gram) / 2L
  k <- ncol(XtY)
  n_state <- length(state_spec$names)
  C <- matrix(0, 2L * p, k)
  state_prob <- matrix(rep(pi_state, each = p), p, n_state)
  converged <- FALSE
  delta <- Inf
  iterations <- 0L

  for (iter in seq_len(control$max_iter)) {
    pi_used <- pi_state
    delta <- 0
    state_sum <- numeric(n_state)
    for (j in seq_len(p)) {
      post <- .koala_mlr_pair_posterior(
        j, C, XtY, Gram, omega_inv, state_spec, pi_used
      )
      target <- .koala_mlr_expected_C(post, k)
      old <- C[post$rows, , drop = FALSE]
      updated <- control$damping * target + (1 - control$damping) * old
      delta <- max(delta, max(abs(updated - old)))
      C[post$rows, ] <- updated
      state_prob[j, ] <- post$state_probability
      state_sum <- state_sum + post$state_probability
    }
    if (isTRUE(control$learn_pi)) {
      pi_state <- state_sum + control$pi_prior_strength * control$pi_initial
      pi_state <- pi_state / sum(pi_state)
      names(pi_state) <- state_spec$names
    }
    iterations <- iter
    if (delta <= control$tol) {
      converged <- TRUE
      break
    }
  }
  final <- .koala_mlr_final_posteriors(
    C, XtY, Gram, omega_inv, state_spec, pi_state, control$epsilon
  )
  c(final, list(C = C, pi_state = pi_state, converged = converged,
                iterations = iterations, final_delta = delta))
}

.koala_mlr_rdirichlet <- function(alpha) {
  z <- stats::rgamma(length(alpha), shape = alpha, rate = 1)
  z / sum(z)
}

.koala_mlr_rinvwishart <- function(df, scale) {
  scale <- (as.matrix(scale) + t(as.matrix(scale))) / 2
  info <- .koala_mlr_logdet_spd(scale)
  precision_scale <- chol2inv(info$chol)
  W <- stats::rWishart(1L, df = df, Sigma = precision_scale)[, , 1L]
  W_info <- .koala_mlr_logdet_spd(W)
  chol2inv(W_info$chol)
}

.koala_mlr_bidirectional_conbf_initial <- function(
    XtY, Gram, omega_inv, state_spec, pi_state, control) {
  p <- nrow(Gram) / 2L
  forward <- .koala_mlr_fit_conditional_bf(
    XtY, Gram, omega_inv, state_spec, pi_state, control
  )
  rev_pair <- rev(seq_len(p))
  row_perm <- c(rev_pair, p + rev_pair)
  reverse_raw <- .koala_mlr_fit_conditional_bf(
    XtY[row_perm, , drop = FALSE],
    Gram[row_perm, row_perm, drop = FALSE],
    omega_inv, state_spec, pi_state, control
  )

  reverse_C <- matrix(0, nrow(forward$C), ncol(forward$C))
  reverse_C[row_perm, ] <- reverse_raw$C
  reverse_state <- integer(p)
  reverse_state[rev_pair] <- reverse_raw$selected_state
  reverse_prob <- matrix(0, p, ncol(reverse_raw$state_probability))
  reverse_prob[rev_pair, ] <- reverse_raw$state_probability
  reverse_log_odds <- numeric(p)
  reverse_log_odds[rev_pair] <- reverse_raw$log_odds

  agree <- forward$log_odds * reverse_log_odds > 0
  choose_forward <- abs(forward$log_odds) <= abs(reverse_log_odds)
  C <- matrix(0, nrow(forward$C), ncol(forward$C))
  selected_state <- rep.int(1L, p)
  for (j in which(agree)) {
    rows <- c(j, p + j)
    if (choose_forward[j]) {
      C[rows, ] <- forward$C[rows, , drop = FALSE]
      selected_state[j] <- forward$selected_state[j]
    } else {
      C[rows, ] <- reverse_C[rows, , drop = FALSE]
      selected_state[j] <- reverse_state[j]
    }
  }
  state_probability <- (forward$state_probability + reverse_prob) / 2
  if (any(!agree)) {
    state_probability[!agree, ] <- 0
    state_probability[!agree, 1L] <- 1
  }
  pi_out <- (as.numeric(forward$pi_state) + as.numeric(reverse_raw$pi_state)) / 2
  names(pi_out) <- state_spec$names
  list(
    C = C,
    selected_state = selected_state,
    state_probability = state_probability,
    pi_state = pi_out / sum(pi_out),
    converged = isTRUE(forward$converged) && isTRUE(reverse_raw$converged),
    iterations = max(forward$iterations, reverse_raw$iterations),
    final_delta = max(forward$final_delta, reverse_raw$final_delta),
    bidirectional = TRUE,
    n_sign_agree = sum(agree),
    n_sign_disagree_or_zero = sum(!agree)
  )
}

.koala_mlr_log_half_t <- function(s, df, scale) {
  ifelse(s > 0, log(2) + stats::dt(s / scale, df = df, log = TRUE) - log(scale), -Inf)
}

.koala_mlr_effect_matrix <- function(C) {
  p <- nrow(C) / 2L
  C[seq_len(p), , drop = FALSE] + C[p + seq_len(p), , drop = FALSE]
}

.koala_mlr_log_slab_likelihood <- function(B, current_state, s_beta,
                                            s_tau = NULL, rho = 0) {
  if (ncol(B) == 1L) {
    active <- current_state == 2L
    if (!any(active)) return(0)
    b <- B[active, 1L]
    return(-length(b) * log(s_beta) - sum(b^2) / (2 * s_beta^2))
  }
  n10 <- current_state == 2L
  n01 <- current_state == 3L
  n11 <- current_state == 4L
  out <- 0
  if (any(n10)) {
    b <- B[n10, 1L]
    out <- out - length(b) * log(s_beta) - sum(b^2) / (2 * s_beta^2)
  }
  if (any(n01)) {
    b <- B[n01, 2L]
    out <- out - length(b) * log(s_tau) - sum(b^2) / (2 * s_tau^2)
  }
  if (any(n11)) {
    one_minus <- 1 - rho^2
    if (!is.finite(one_minus) || one_minus <= 0) return(-Inf)
    z1 <- B[n11, 1L] / s_beta
    z2 <- B[n11, 2L] / s_tau
    quad <- (z1^2 - 2 * rho * z1 * z2 + z2^2) / one_minus
    out <- out - sum(n11) * (log(s_beta) + log(s_tau) + 0.5 * log(one_minus)) -
      0.5 * sum(quad)
  }
  out
}

.koala_mlr_sample_scale_grid <- function(B, current_state, which_scale,
                                          current_s_beta, current_s_tau,
                                          current_rho, grid, prior_df,
                                          prior_scale) {
  logw <- .koala_mlr_log_half_t(grid, prior_df, prior_scale)
  for (g in seq_along(grid)) {
    sb <- if (identical(which_scale, "beta")) grid[g] else current_s_beta
    st <- if (identical(which_scale, "tau")) grid[g] else current_s_tau
    logw[g] <- logw[g] + .koala_mlr_log_slab_likelihood(
      B, current_state, sb, st, current_rho
    )
  }
  sample(grid, size = 1L, prob = .koala_mlr_softmax(logw))
}

.koala_mlr_sample_rho_grid <- function(B, current_state, s_beta, s_tau,
                                        rho_grid, rho_lkj_eta) {
  logw <- (rho_lkj_eta - 1) * log1p(-rho_grid^2)
  for (g in seq_along(rho_grid)) {
    logw[g] <- logw[g] + .koala_mlr_log_slab_likelihood(
      B, current_state, s_beta, s_tau, rho_grid[g]
    )
  }
  prob <- .koala_mlr_softmax(logw)
  sample(rho_grid, size = 1L, prob = prob)
}

.koala_mlr_fit_gibbs <- function(
    XtY, Gram, omega_inv, state_spec, pi_state, control,
    conbf_initial = NULL, X = NULL, Y = NULL, single_channel = NULL) {
  p <- nrow(Gram) / 2L
  k <- ncol(XtY)
  n_state <- length(state_spec$names)

  if (is.null(conbf_initial)) {
    C <- matrix(0, 2L * p, k)
    current_state <- rep.int(1L, p)
    initialization <- "zero"
    warm_start_iterations <- 0L
    warm_start_converged <- NA
    warm_start_final_delta <- NA_real_
  } else {
    C <- as.matrix(conbf_initial$C)
    state_probability <- as.matrix(conbf_initial$state_probability)
    if (!all(dim(C) == c(2L * p, k)) ||
        anyNA(C) || any(!is.finite(C))) {
      stop("The conditional-BF warm start returned an invalid coefficient matrix.")
    }
    if (!all(dim(state_probability) == c(p, n_state)) ||
        anyNA(state_probability) ||
        any(!is.finite(state_probability))) {
      stop("The conditional-BF warm start returned invalid state probabilities.")
    }
    current_state <- as.integer(conbf_initial$selected_state)
    if (length(current_state) != p || anyNA(current_state) ||
        any(current_state < 1L | current_state > n_state)) {
      stop("The conditional-BF warm start returned invalid state labels.")
    }
    pi_state <- .koala_mlr_normalize_pi(
      conbf_initial$pi_state,
      state_spec$names
    )
    initialization <- "conditional_bf"
    warm_start_iterations <- as.integer(conbf_initial$iterations)
    warm_start_converged <- conbf_initial$converged
    warm_start_final_delta <- as.numeric(conbf_initial$final_delta)
  }

  # These accumulators intentionally start from zero even under warm start:
  # conditional_bf supplies only the initial Markov-chain state and is not
  # counted as a posterior Gibbs draw.
  saved <- 0L
  pA_sum <- numeric(p)
  state_sum <- matrix(0, p, n_state)
  omega_sum <- matrix(0, k, k)
  s_beta_sum <- 0
  s_tau_sum <- 0
  rho_sum <- 0
  hyperparameter_draws <- matrix(
    NA_real_, control$gibbs_samples, 6L,
    dimnames = list(NULL, c(
      "s_beta", "s_tau", "rho", "omega_11", "omega_12", "omega_22"
    ))
  )
  current_omega <- chol2inv(.koala_mlr_logdet_spd(omega_inv)$chol)
  current_s_beta <- control$effect_sd_beta
  current_s_tau <- control$effect_sd_tau
  current_rho <- control$effect_rho

  needs_residual <- identical(control$omega_mode, "iw")
  if (needs_residual && (is.null(X) || is.null(Y))) {
    stop("Omega IW update requires the standardized design and response matrices.")
  }

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(control$seed))

  total_sweeps <- control$gibbs_burnin +
    control$gibbs_samples * control$gibbs_thin
  for (sweep in seq_len(total_sweeps)) {
    state_count <- tabulate(current_state, nbins = n_state)
    for (j in seq_len(p)) {
      post <- .koala_mlr_pair_posterior(
        j, C, XtY, Gram, omega_inv, state_spec, pi_state,
        sample_effect = TRUE
      )
      chosen <- sample.int(length(post$probability), 1L,
                           prob = post$probability)
      C_new <- matrix(0, 2L, k)
      chosen_s <- post$orientation_index[chosen]
      chosen_h <- post$state_index[chosen]
      C_new[chosen_s, ] <- post$components[[chosen]]$draw
      C[post$rows, ] <- C_new
      state_count[current_state[j]] <- state_count[current_state[j]] - 1L
      current_state[j] <- chosen_h
      state_count[chosen_h] <- state_count[chosen_h] + 1L

      keep <- sweep > control$gibbs_burnin &&
        ((sweep - control$gibbs_burnin - 1L) %% control$gibbs_thin == 0L)
      if (keep) {
        pA_sum[j] <- pA_sum[j] + post$orientation_probability[1L]
        state_sum[j, ] <- state_sum[j, ] + post$state_probability
      }
    }
    if (isTRUE(control$learn_pi)) {
      alpha <- state_count +
        control$pi_prior_strength * control$pi_initial
      pi_state <- .koala_mlr_rdirichlet(alpha)
      names(pi_state) <- state_spec$names
    }

    if (identical(control$omega_mode, "iw") && k == 2L) {
      E <- Y - X %*% C
      scale_post <- control$omega_iw_scale + crossprod(E)
      current_omega <- .koala_mlr_rinvwishart(
        control$omega_iw_df + nrow(Y), scale_post
      )
      omega_inv <- chol2inv(.koala_mlr_logdet_spd(current_omega)$chol)
      rm(E)
    }

    if (k == 1L) {
      B <- .koala_mlr_effect_matrix(C)
      active_idx <- which(current_state == 2L)
      active_coef <- if (length(active_idx)) B[active_idx, 1L] else numeric()
      is_tau_channel <- identical(single_channel, "tau")
      mode_one <- if (is_tau_channel) control$s_tau_mode else control$s_beta_mode
      shape_one <- if (is_tau_channel) control$s_tau_ig_shape else control$s_beta_ig_shape
      rate_one <- if (is_tau_channel) control$s_tau_ig_rate else control$s_beta_ig_rate
      grid_one <- if (is_tau_channel) control$s_tau_grid else control$s_beta_grid
      prior_df_one <- if (is_tau_channel) control$s_tau_half_t_df else control$s_beta_half_t_df
      prior_scale_one <- if (is_tau_channel) control$s_tau_half_t_scale else control$s_beta_half_t_scale
      current_one <- if (is_tau_channel) current_s_tau else current_s_beta
      if (identical(mode_one, "ig")) {
        s2 <- 1 / stats::rgamma(
          1L,
          shape = shape_one + length(active_coef) / 2,
          rate = rate_one + sum(active_coef^2) / 2
        )
        current_one <- sqrt(s2)
      } else if (identical(mode_one, "grid")) {
        # For a one-channel fit, the first column of B is the fitted channel,
        # regardless of whether that channel is named beta or tau.
        current_one <- .koala_mlr_sample_scale_grid(
          B = B,
          current_state = current_state,
          which_scale = "beta",
          current_s_beta = current_one,
          current_s_tau = NULL,
          current_rho = 0,
          grid = grid_one,
          prior_df = prior_df_one,
          prior_scale = prior_scale_one
        )
      }
      if (is_tau_channel) current_s_tau <- current_one else current_s_beta <- current_one
      state_spec <- .koala_mlr_state_spec(1L, current_one, 0)
    } else {
      B <- .koala_mlr_effect_matrix(C)
      beta_idx <- which(current_state %in% c(2L, 4L))
      tau_idx <- which(current_state %in% c(3L, 4L))
      beta_coef <- if (length(beta_idx)) B[beta_idx, 1L] else numeric()
      tau_coef <- if (length(tau_idx)) B[tau_idx, 2L] else numeric()
      if (identical(control$s_beta_mode, "ig")) {
        s2 <- 1 / stats::rgamma(
          1L,
          shape = control$s_beta_ig_shape + length(beta_coef) / 2,
          rate = control$s_beta_ig_rate + sum(beta_coef^2) / 2
        )
        current_s_beta <- sqrt(s2)
      } else if (identical(control$s_beta_mode, "grid")) {
        current_s_beta <- .koala_mlr_sample_scale_grid(
          B, current_state, "beta",
          current_s_beta, current_s_tau, current_rho,
          control$s_beta_grid,
          control$s_beta_half_t_df,
          control$s_beta_half_t_scale
        )
      }
      if (identical(control$s_tau_mode, "ig")) {
        s2 <- 1 / stats::rgamma(
          1L,
          shape = control$s_tau_ig_shape + length(tau_coef) / 2,
          rate = control$s_tau_ig_rate + sum(tau_coef^2) / 2
        )
        current_s_tau <- sqrt(s2)
      } else if (identical(control$s_tau_mode, "grid")) {
        current_s_tau <- .koala_mlr_sample_scale_grid(
          B, current_state, "tau",
          current_s_beta, current_s_tau, current_rho,
          control$s_tau_grid,
          control$s_tau_half_t_df,
          control$s_tau_half_t_scale
        )
      }
      if (identical(control$rho_mode, "grid")) {
        current_rho <- .koala_mlr_sample_rho_grid(
          B, current_state, current_s_beta, current_s_tau,
          control$rho_grid, control$rho_lkj_eta
        )
      }
      state_spec <- .koala_mlr_state_spec(
        2L, c(current_s_beta, current_s_tau), current_rho
      )
    }

    if (sweep > control$gibbs_burnin &&
        ((sweep - control$gibbs_burnin - 1L) %% control$gibbs_thin == 0L)) {
      saved <- saved + 1L
      omega_sum <- omega_sum + current_omega
      s_beta_sum <- s_beta_sum + current_s_beta
      s_tau_sum <- s_tau_sum + current_s_tau
      rho_sum <- rho_sum + current_rho
      hyperparameter_draws[saved, ] <- c(
        current_s_beta,
        current_s_tau,
        current_rho,
        current_omega[1L, 1L],
        if (k >= 2L) current_omega[1L, 2L] else NA_real_,
        if (k >= 2L) current_omega[2L, 2L] else NA_real_
      )
    }
  }
  if (saved < 1L) stop("The Gibbs sampler retained no posterior draws.")
  pA <- pmin(pmax(pA_sum / saved, control$epsilon), 1 - control$epsilon)
  state_probability <- state_sum / saved
  colnames(state_probability) <- state_spec$names
  list(
    pA = pA,
    state_probability = state_probability,
    log_odds = stats::qlogis(pA),
    C = C,
    pi_state = pi_state,
    converged = NA,
    iterations = total_sweeps,
    final_delta = NA_real_,
    posterior_draws = saved,
    initialization = initialization,
    warm_start_iterations = warm_start_iterations,
    warm_start_converged = warm_start_converged,
    warm_start_final_delta = warm_start_final_delta,
    posterior_mean_omega = omega_sum / saved,
    posterior_mean_s_beta = s_beta_sum / saved,
    posterior_mean_s_tau = s_tau_sum / saved,
    posterior_mean_rho = rho_sum / saved,
    hyperparameter_draws = hyperparameter_draws[seq_len(saved), , drop = FALSE],
    final_omega = current_omega,
    final_s_beta = current_s_beta,
    final_s_tau = current_s_tau,
    final_rho = current_rho
  )
}

koala_mlr_control <- function(
    method = c("conditional_bf", "variational", "gibbs"),
    max_iter = 100L,
    tol = 1e-6,
    damping = 0.5,
    epsilon = 1e-8,
    effect_sd_beta = 0.10,
    effect_sd_tau = 0.10,
    effect_rho = 0,
    pi_univariate = c(0.90, 0.10),
    pi_bivariate = c(0.90, 0.035, 0.035, 0.03),
    learn_pi = TRUE,
    pi_prior_strength = 20,
    omega_ridge = 1e-6,
    canonical_tol = 0,
    omega_mode = c("plugin", "iw"),
    omega_iw_df = 5,
    omega_iw_scale = diag(2) * 2,
    s_beta_mode = c("fixed", "ig", "grid"),
    s_tau_mode = c("fixed", "ig", "grid"),
    s_beta_grid = seq(0.02, 0.30, by = 0.01),
    s_tau_grid = seq(0.02, 0.30, by = 0.01),
    s_beta_half_t_df = 3,
    s_beta_half_t_scale = 0.10,
    s_tau_half_t_df = 3,
    s_tau_half_t_scale = 0.10,
    s_beta_ig_shape = 2.5,
    s_beta_ig_rate = 0.015,
    s_tau_ig_shape = 2.5,
    s_tau_ig_rate = 0.015,
    rho_mode = c("fixed", "grid"),
    rho_grid = seq(-0.8, 0.8, by = 0.1),
    rho_lkj_eta = 2,
    conbf_bidirectional = TRUE,
    gibbs_bidirectional = FALSE,
    gibbs_init_from_conbf = FALSE,
    gibbs_burnin = 200L,
    gibbs_samples = 500L,
    gibbs_thin = 2L,
    seed = 718L) {
  method <- match.arg(method)
  omega_mode <- match.arg(omega_mode)
  s_beta_mode <- match.arg(s_beta_mode)
  s_tau_mode <- match.arg(s_tau_mode)
  rho_mode <- match.arg(rho_mode)
  if (!identical(method, "gibbs") &&
      (!identical(omega_mode, "plugin") ||
       !identical(s_beta_mode, "fixed") ||
       !identical(s_tau_mode, "fixed") ||
       !identical(rho_mode, "fixed"))) {
    stop("Omega/slab hyperparameter sampling is available only for method='gibbs'.")
  }
  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1L) {
    stop("max_iter must be a positive integer.")
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be finite and positive.")
  }
  if (!is.finite(damping) || damping <= 0 || damping > 1) {
    stop("damping must lie in (0,1].")
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon <= 0 || epsilon >= 0.5) {
    stop("epsilon must lie in (0,0.5).")
  }
  if (length(gibbs_burnin) != 1L || !is.finite(gibbs_burnin) ||
      gibbs_burnin < 0L || length(gibbs_samples) != 1L ||
      !is.finite(gibbs_samples) || gibbs_samples < 1L ||
      length(gibbs_thin) != 1L || !is.finite(gibbs_thin) ||
      gibbs_thin < 1L) {
    stop("Invalid Gibbs burn-in, sample, or thinning setting.")
  }
  if (!is.logical(gibbs_init_from_conbf) ||
      length(gibbs_init_from_conbf) != 1L ||
      is.na(gibbs_init_from_conbf)) {
    stop("gibbs_init_from_conbf must be TRUE or FALSE.")
  }
  if (identical(rho_mode, "grid") &&
      (!identical(s_beta_mode, "grid") || !identical(s_tau_mode, "grid"))) {
    stop("rho_mode='grid' requires grid updates for both slab scales.")
  }
  if (!identical(rho_mode, "grid") &&
      (identical(s_beta_mode, "grid") || identical(s_tau_mode, "grid"))) {
    stop("Slab-scale grid updates must be combined with the rho grid update.")
  }
  if (xor(identical(s_beta_mode, "grid"), identical(s_tau_mode, "grid"))) {
    stop("s_beta and s_tau must both use grid updates.")
  }
  if ((identical(s_beta_mode, "ig") || identical(s_tau_mode, "ig")) &&
      abs(effect_rho) > sqrt(.Machine$double.eps)) {
    stop("Inverse-Gamma slab updates require fixed rho=0.")
  }
  if ((identical(s_beta_mode, "ig") || identical(s_tau_mode, "ig")) &&
      !identical(rho_mode, "fixed")) {
    stop("Inverse-Gamma slab updates are available only with fixed rho.")
  }
  if (identical(rho_mode, "fixed") &&
      (identical(s_beta_mode, "ig") || identical(s_tau_mode, "ig")) &&
      abs(effect_rho) > sqrt(.Machine$double.eps)) {
    stop("effect_rho must equal zero for conjugate Inverse-Gamma updates.")
  }
  if (!is.logical(conbf_bidirectional) || length(conbf_bidirectional) != 1L ||
      is.na(conbf_bidirectional) || !is.logical(gibbs_bidirectional) ||
      length(gibbs_bidirectional) != 1L || is.na(gibbs_bidirectional)) {
    stop("conbf_bidirectional and gibbs_bidirectional must be TRUE or FALSE.")
  }
  omega_iw_scale <- as.matrix(omega_iw_scale)
  if (!all(dim(omega_iw_scale) == c(2L, 2L)) ||
      omega_iw_df <= 3 || any(!is.finite(omega_iw_scale))) {
    stop("Omega IW prior requires df>3 and a finite 2x2 scale matrix.")
  }
  rho_grid <- sort(unique(as.numeric(rho_grid)))
  if (!length(rho_grid) || any(!is.finite(rho_grid)) || any(abs(rho_grid) >= 1)) {
    stop("rho_grid must contain finite values strictly between -1 and 1.")
  }
  s_beta_grid <- sort(unique(as.numeric(s_beta_grid)))
  s_tau_grid <- sort(unique(as.numeric(s_tau_grid)))
  if (!length(s_beta_grid) || any(!is.finite(s_beta_grid)) || any(s_beta_grid <= 0) ||
      !length(s_tau_grid) || any(!is.finite(s_tau_grid)) || any(s_tau_grid <= 0)) {
    stop("s_beta_grid and s_tau_grid must contain finite positive values.")
  }
  if (any(!is.finite(c(s_beta_half_t_df, s_beta_half_t_scale,
                       s_tau_half_t_df, s_tau_half_t_scale))) ||
      any(c(s_beta_half_t_df, s_beta_half_t_scale,
            s_tau_half_t_df, s_tau_half_t_scale) <= 0)) {
    stop("Half-t prior degrees of freedom and scales must be finite and positive.")
  }
  if (length(rho_lkj_eta) != 1L || !is.finite(rho_lkj_eta) || rho_lkj_eta <= 0) {
    stop("rho_lkj_eta must be a finite positive scalar.")
  }
  if (any(!is.finite(c(s_beta_ig_shape, s_beta_ig_rate,
                       s_tau_ig_shape, s_tau_ig_rate))) ||
      any(c(s_beta_ig_shape, s_beta_ig_rate,
            s_tau_ig_shape, s_tau_ig_rate) <= 0)) {
    stop("Inverse-Gamma prior shape/rate parameters must be finite and positive.")
  }
  list(
    method = method,
    max_iter = as.integer(max_iter),
    tol = as.numeric(tol),
    damping = as.numeric(damping),
    epsilon = as.numeric(epsilon),
    effect_sd_beta = as.numeric(effect_sd_beta),
    effect_sd_tau = as.numeric(effect_sd_tau),
    effect_rho = as.numeric(effect_rho),
    pi_univariate = as.numeric(pi_univariate),
    pi_bivariate = as.numeric(pi_bivariate),
    learn_pi = isTRUE(learn_pi),
    pi_prior_strength = as.numeric(pi_prior_strength),
    omega_ridge = as.numeric(omega_ridge),
    canonical_tol = as.numeric(canonical_tol),
    omega_mode = omega_mode,
    omega_iw_df = as.numeric(omega_iw_df),
    omega_iw_scale = omega_iw_scale,
    s_beta_mode = s_beta_mode,
    s_tau_mode = s_tau_mode,
    s_beta_grid = s_beta_grid,
    s_tau_grid = s_tau_grid,
    s_beta_half_t_df = as.numeric(s_beta_half_t_df),
    s_beta_half_t_scale = as.numeric(s_beta_half_t_scale),
    s_tau_half_t_df = as.numeric(s_tau_half_t_df),
    s_tau_half_t_scale = as.numeric(s_tau_half_t_scale),
    s_beta_ig_shape = as.numeric(s_beta_ig_shape),
    s_beta_ig_rate = as.numeric(s_beta_ig_rate),
    s_tau_ig_shape = as.numeric(s_tau_ig_shape),
    s_tau_ig_rate = as.numeric(s_tau_ig_rate),
    rho_mode = rho_mode,
    rho_grid = rho_grid,
    rho_lkj_eta = as.numeric(rho_lkj_eta),
    conbf_bidirectional = conbf_bidirectional,
    gibbs_bidirectional = gibbs_bidirectional,
    gibbs_init_from_conbf = gibbs_init_from_conbf,
    gibbs_burnin = as.integer(gibbs_burnin),
    gibbs_samples = as.integer(gibbs_samples),
    gibbs_thin = as.integer(gibbs_thin),
    seed = as.integer(seed)
  )
}

koala_mlr_fit <- function(G, G_tilde, Y, control = koala_mlr_control(),
                          channel_names = NULL, center = TRUE) {
  fit_wall_start <- proc.time()[[3L]]
  if (!is.list(control) || is.null(control$method)) {
    stop("control must be created by the LKOM control constructor.")
  }
  G <- as.matrix(G)
  G_tilde <- as.matrix(G_tilde)
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (nrow(G) != nrow(Y)) stop("G and Y have different sample sizes.")
  if (!ncol(G)) stop("The MLR fit requires at least one screened variant pair.")
  k <- ncol(Y)
  if (!k %in% 1:2) stop("Y must contain one or two columns.")
  if (is.null(channel_names)) channel_names <- if (k == 1L) "trait" else c("beta", "tau")

  canonical <- .koala_mlr_canonical_pairs(
    G, G_tilde, tol = control$canonical_tol
  )
  scaled_A <- .koala_mlr_scale_matrix(canonical$A, "canonical A", center)
  scaled_B <- .koala_mlr_scale_matrix(canonical$B, "canonical B", center)
  scaled_Y <- .koala_mlr_scale_matrix(Y, "MLR response", center)
  X <- cbind(scaled_A$x, scaled_B$x)
  Ys <- scaled_Y$x

  Omega <- if (k == 1L) {
    matrix(sum(Ys^2) / nrow(Ys), 1L, 1L)
  } else {
    crossprod(Ys) / nrow(Ys)
  }
  Omega <- (Omega + t(Omega)) / 2 + diag(control$omega_ridge, k)
  omega_info <- .koala_mlr_logdet_spd(Omega)
  omega_inv <- chol2inv(omega_info$chol)

  effect_sd <- if (k == 1L) {
    if (identical(channel_names[1L], "tau")) control$effect_sd_tau else control$effect_sd_beta
  } else {
    c(control$effect_sd_beta, control$effect_sd_tau)
  }
  state_spec <- .koala_mlr_state_spec(k, effect_sd, control$effect_rho)
  pi0 <- if (k == 1L) control$pi_univariate else control$pi_bivariate
  pi_state <- .koala_mlr_normalize_pi(pi0, state_spec$names)
  control$pi_initial <- pi_state

  Gram <- crossprod(X)
  XtY <- crossprod(X, Ys)
  backend_timing <- list(
    conbf_initialization_seconds = 0,
    mcmc_seconds = NA_real_,
    backend_total_seconds = NA_real_
  )
  backend_start <- proc.time()[[3L]]
  backend <- switch(
    control$method,
    conditional_bf = .koala_mlr_fit_conditional_bf(
      XtY, Gram, omega_inv, state_spec, pi_state, control
    ),
    variational = .koala_mlr_fit_variational(
      XtY, Gram, omega_inv, state_spec, pi_state, control
    ),
    gibbs = {
      conbf_initial <- NULL
      if (isTRUE(control$gibbs_init_from_conbf)) {
        # Deterministic conditional-BF fitting provides C, latent-state modes,
        # and the learned state prior. It is used only to initialize Gibbs.
        conbf_start <- proc.time()[[3L]]
        conbf_initial <- if (isTRUE(control$conbf_bidirectional)) {
          .koala_mlr_bidirectional_conbf_initial(
            XtY, Gram, omega_inv, state_spec, pi_state, control
          )
        } else {
          .koala_mlr_fit_conditional_bf(
            XtY, Gram, omega_inv, state_spec, pi_state, control
          )
        }
        backend_timing$conbf_initialization_seconds <-
          proc.time()[[3L]] - conbf_start
      }
      mcmc_start <- proc.time()[[3L]]
      ans <- .koala_mlr_fit_gibbs(
        XtY, Gram, omega_inv, state_spec, pi_state, control,
        conbf_initial = conbf_initial, X = X, Y = Ys,
        single_channel = if (k == 1L) channel_names[1L] else NULL
      )
      backend_timing$mcmc_seconds <- proc.time()[[3L]] - mcmc_start
      ans
    },
    stop("Unknown MLR method.")
  )
  backend_timing$backend_total_seconds <- proc.time()[[3L]] - backend_start

  W <- canonical$observed_orientation * backend$log_odds
  W[canonical$identical_pair] <- 0
  W_swap <- -canonical$observed_orientation * backend$log_odds
  W_swap[canonical$identical_pair] <- 0
  swap_error <- max(abs(W + W_swap))
  if (!is.finite(swap_error) || swap_error > 10 * .Machine$double.eps) {
    stop("The internal MLR sign-flip check failed.")
  }
  state_probability <- as.data.frame(backend$state_probability)
  names(state_probability) <- paste0("posterior_state_", names(state_probability))
  state_class <- state_spec$names[max.col(
    as.matrix(backend$state_probability), ties.method = "first"
  )]

  fit_total_seconds <- proc.time()[[3L]] - fit_wall_start
  list(
    method = control$method,
    W = setNames(W, colnames(G)),
    posterior_original = setNames(stats::plogis(W), colnames(G)),
    canonical_log_odds = setNames(backend$log_odds, colnames(G)),
    observed_orientation = setNames(canonical$observed_orientation, colnames(G)),
    identical_pair = setNames(canonical$identical_pair, colnames(G)),
    posterior_state = state_probability,
    posterior_state_class = setNames(state_class, colnames(G)),
    fitted_state_prior = backend$pi_state,
    converged = backend$converged,
    iterations = backend$iterations,
    final_delta = backend$final_delta,
    posterior_draws = backend$posterior_draws %||% NA_integer_,
    initialization = backend$initialization %||% "not_applicable",
    warm_start_iterations = backend$warm_start_iterations %||% 0L,
    warm_start_converged = backend$warm_start_converged %||% NA,
    warm_start_final_delta = backend$warm_start_final_delta %||% NA_real_,
    posterior_mean_omega = backend$posterior_mean_omega %||% Omega,
    posterior_mean_s_beta = backend$posterior_mean_s_beta %||% control$effect_sd_beta,
    posterior_mean_s_tau = backend$posterior_mean_s_tau %||% control$effect_sd_tau,
    posterior_mean_rho = backend$posterior_mean_rho %||% control$effect_rho,
    hyperparameter_draws = backend$hyperparameter_draws %||% NULL,
    runtime_seconds = c(
      total = fit_total_seconds,
      conbf_initialization = backend_timing$conbf_initialization_seconds,
      mcmc = backend_timing$mcmc_seconds,
      backend_total = backend_timing$backend_total_seconds
    ),
    diagnostics = list(
      n_subjects = nrow(G),
      n_pairs = ncol(G),
      outcome_covariance = Omega,
      omega_mode = if (k == 1L && identical(control$omega_mode, "iw")) {
        "plugin_univariate"
      } else {
        control$omega_mode
      },
      s_beta_mode = control$s_beta_mode,
      s_tau_mode = control$s_tau_mode,
      rho_mode = control$rho_mode,
      canonical_swap_check = canonical$swap_check,
      max_abs_W_plus_W_swap = swap_error
    )
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Conservative order-ensemble statistic.  A pair is retained only when the
# forward and reverse fits agree on its sign; its magnitude is the smaller of
# the two absolute statistics.  This aggregation is odd, so simultaneous
# sign-flips of the two input statistics also flip the combined statistic.
.koala_mlr_signed_min_consensus <- function(W_forward, W_reverse) {
  W_forward <- as.numeric(W_forward)
  W_reverse <- as.numeric(W_reverse)
  if (length(W_forward) != length(W_reverse)) {
    stop("Forward and reverse MLR statistics have different lengths.")
  }
  if (any(!is.finite(W_forward)) || any(!is.finite(W_reverse))) {
    stop("Forward or reverse MLR statistics contain non-finite values.")
  }

  sign_agree <- W_forward * W_reverse > 0
  W <- numeric(length(W_forward))
  W[sign_agree] <- sign(W_forward[sign_agree]) * pmin(
    abs(W_forward[sign_agree]),
    abs(W_reverse[sign_agree])
  )
  W
}

.koala_mlr_bind_selected <- function(prepared_group, selected) {
  G <- do.call(cbind, lapply(seq_along(prepared_group), function(i) {
    prepared_group[[i]]$G[, selected[[i]], drop = FALSE]
  }))
  G_tilde <- do.call(cbind, lapply(seq_along(prepared_group), function(i) {
    prepared_group[[i]]$G_tilde[, selected[[i]], drop = FALSE]
  }))
  list(G = G, G_tilde = G_tilde)
}

# Fit beta-specific, tau-specific, and bivariate any-effect MLRs for one group.
# Channel-specific MLRs deliberately use only their own residual, so a tau-only
# effect cannot make a beta-null SNP non-null (and vice versa).
koala_fit_group_mlr <- function(
    null_fit,
    prepared_group,
    group_id,
    causal_beta = NULL,
    causal_tau = NULL,
    control = koala_mlr_control(),
    reverse_pair_order = FALSE)  {

  if (!length(prepared_group)) stop("prepared_group cannot be empty.")
  snp_names <- unlist(lapply(prepared_group, function(b) colnames(b$G)),
                      use.names = FALSE)
  block_vec <- unlist(lapply(prepared_group, function(b) {
    rep(as.character(b$block_id), ncol(b$G))
  }), use.names = FALSE)
  if (anyDuplicated(snp_names)) stop("SNP identifiers are duplicated across blocks.")
  p_all <- length(snp_names)

  spa_beta_selected <- lapply(prepared_group, function(b) {
    as.logical(b$screening$mean$selected)
  })
  spa_tau_selected <- lapply(prepared_group, function(b) {
    as.logical(b$screening$ws$selected)
  })
  # Upstream cohort-specific QC removes constant original variants and
  # create.MK() supplies the paired knockoffs. Do not repeat an expensive full
  # N-by-p variance scan here; selected columns are still validated when the
  # MLR standardizes its much smaller screened design.
  usable_pair <- lapply(prepared_group, function(b) rep(TRUE, ncol(b$G)))
  beta_selected <- Map(`&`, spa_beta_selected, usable_pair)
  tau_selected <- Map(`&`, spa_tau_selected, usable_pair)
  any_selected <- Map(`|`, beta_selected, tau_selected)
  spa_beta_keep <- unlist(spa_beta_selected, use.names = FALSE)
  spa_tau_keep <- unlist(spa_tau_selected, use.names = FALSE)
  usable_keep <- unlist(usable_pair, use.names = FALSE)
  beta_keep <- unlist(beta_selected, use.names = FALSE)
  tau_keep <- unlist(tau_selected, use.names = FALSE)
  any_keep <- beta_keep | tau_keep

  y_beta <- .koala_mlr_rank_int(null_fit$mean$c)
  y_tau <- .koala_mlr_rank_int(null_fit$ws$c)
#   y_beta <- null_fit$mean$c
#   y_tau <- null_fit$ws$c
  if (length(y_beta) != null_fit$m || length(y_tau) != null_fit$m) {
    stop("WiSER subject-level residuals are not aligned with null_fit$m.")
  }

  fit_one <- function(selected_list, Y, channel_names, seed_offset) {
  if (!sum(vapply(selected_list, sum, numeric(1)))) return(NULL)

  pair <- .koala_mlr_bind_selected(
    prepared_group,
    selected_list
  )

  # Record the standard SNP order before fitting the MLR.
  original_snp_order <- colnames(pair$G)

  local_control <- control
  local_control$seed <- as.integer(
    control$seed + seed_offset
  )

  fit_one_order <- function(order_index) {
    fit <- koala_mlr_fit(
      G = pair$G[, order_index, drop = FALSE],
      G_tilde = pair$G_tilde[, order_index, drop = FALSE],
      Y = Y,
      control = local_control,
      channel_names = channel_names
    )

    # Restore every SNP-level result to the order before the MLR fit.
    fitted_snp_order <- names(fit$W)
    restore_index <- match(original_snp_order, fitted_snp_order)
    if (anyNA(restore_index)) {
      stop("The original SNP order could not be restored after MLR fitting.")
    }

    vector_fields <- c(
      "W",
      "posterior_original",
      "canonical_log_odds",
      "observed_orientation",
      "identical_pair",
      "posterior_state_class",
      "ld_cluster_id",
      "posterior_cluster_active"
    )

    for (field in vector_fields) {
      if (!is.null(fit[[field]])) {
        fit[[field]] <- fit[[field]][restore_index]
      }
    }

    if (!is.null(fit$posterior_state)) {
      fit$posterior_state <- fit$posterior_state[
        restore_index,
        ,
        drop = FALSE
      ]
      rownames(fit$posterior_state) <- original_snp_order
    }

    fit
  }

  forward_index <- seq_len(ncol(pair$G))
  reverse_index <- rev(forward_index)
  fit_forward <- fit_one_order(forward_index)
  W_forward <- setNames(as.numeric(fit_forward$W), original_snp_order)
  run_reverse <- if (identical(control$method, "gibbs")) {
    isTRUE(control$gibbs_bidirectional)
  } else {
    TRUE
  }

  if (run_reverse) {
    fit_reverse <- fit_one_order(reverse_index)
    W_reverse <- setNames(as.numeric(fit_reverse$W), original_snp_order)
    W_consensus <- setNames(
      .koala_mlr_signed_min_consensus(W_forward, W_reverse),
      original_snp_order
    )
    sign_agree <- setNames(W_forward * W_reverse > 0, original_snp_order)
    fit <- fit_forward
    if (!is.null(fit_forward$runtime_seconds) &&
        !is.null(fit_reverse$runtime_seconds)) {
      shared_runtime_names <- intersect(
        names(fit_forward$runtime_seconds),
        names(fit_reverse$runtime_seconds)
      )
      fit$runtime_seconds[shared_runtime_names] <-
        fit_forward$runtime_seconds[shared_runtime_names] +
        fit_reverse$runtime_seconds[shared_runtime_names]
    }
    fit$W <- W_consensus
    fit$posterior_original <- setNames(stats::plogis(W_consensus), original_snp_order)
    fit$W_forward <- W_forward
    fit$W_reverse <- W_reverse
    fit$order_sign_agree <- sign_agree
    fit$diagnostics$bidirectional_order_fit <- list(
      enabled = TRUE,
      aggregation = "same_sign_min_abs_else_zero",
      n_sign_agree = sum(sign_agree),
      n_sign_disagree_or_zero = sum(!sign_agree),
      max_abs_forward_minus_reverse = max(abs(W_forward - W_reverse)),
      forward_converged = fit_forward$converged,
      reverse_converged = fit_reverse$converged,
      forward_iterations = fit_forward$iterations,
      reverse_iterations = fit_reverse$iterations,
      forward_initialization = fit_forward$initialization,
      reverse_initialization = fit_reverse$initialization
    )
  } else {
    fit <- fit_forward
    fit$W_forward <- W_forward
    fit$W_reverse <- NULL
    fit$order_sign_agree <- NULL
    fit$diagnostics$bidirectional_order_fit <- list(
      enabled = FALSE,
      aggregation = "forward_only",
      conbf_bidirectional_initialization = isTRUE(control$conbf_bidirectional)
    )
  }

  fit
}

  time_fit_one <- function(...) {
    started <- proc.time()[[3L]]
    fit <- fit_one(...)
    list(fit = fit, elapsed = proc.time()[[3L]] - started)
  }
  timed_beta <- time_fit_one(
    beta_selected, matrix(y_beta, ncol = 1L), "beta", 101L
  )
  timed_tau <- time_fit_one(
    tau_selected, matrix(y_tau, ncol = 1L), "tau", 202L
  )
  timed_any <- time_fit_one(
    any_selected, cbind(beta = y_beta, tau = y_tau),
    c("beta", "tau"), 303L
  )
  fit_beta <- timed_beta$fit
  fit_tau <- timed_tau$fit
  fit_any <- timed_any$fit

  runtime_value <- function(fit, name) {
    if (is.null(fit) || is.null(fit$runtime_seconds) ||
        !name %in% names(fit$runtime_seconds)) return(0)
    as.numeric(fit$runtime_seconds[[name]])
  }
  runtime_seconds <- c(
    beta_total = timed_beta$elapsed,
    tau_total = timed_tau$elapsed,
    any_total = timed_any$elapsed,
    beta_conbf_initialization = runtime_value(fit_beta, "conbf_initialization"),
    tau_conbf_initialization = runtime_value(fit_tau, "conbf_initialization"),
    any_conbf_initialization = runtime_value(fit_any, "conbf_initialization"),
    beta_mcmc = runtime_value(fit_beta, "mcmc"),
    tau_mcmc = runtime_value(fit_tau, "mcmc"),
    any_mcmc = runtime_value(fit_any, "mcmc")
  )

  expand_numeric <- function(fit, keep, field = "W", default = 0) {
    out <- rep(default, p_all)
    if (!is.null(fit)) out[keep] <- as.numeric(fit[[field]])
    setNames(out, snp_names)
  }
  expand_character <- function(fit, keep, field, default = "not_screened") {
    out <- rep(default, p_all)
    if (!is.null(fit)) out[keep] <- as.character(fit[[field]])
    out
  }
  expand_states <- function(fit, keep, prefix) {
    if (is.null(fit)) return(data.frame())
    x <- as.data.frame(fit$posterior_state)
    out <- matrix(NA_real_, p_all, ncol(x))
    out[keep, ] <- as.matrix(x)
    out <- as.data.frame(out)
    names(out) <- paste0(prefix, sub("^posterior_state_", "", names(x)))
    out
  }

  W_beta <- expand_numeric(fit_beta, beta_keep)
  W_tau <- expand_numeric(fit_tau, tau_keep)
  W_any <- expand_numeric(fit_any, any_keep)

  get_p <- function(channel, side, field) {
    unlist(lapply(prepared_group, function(b) {
      b$screening[[channel]][[side]][[field]]
    }), use.names = FALSE)
  }
  snp_info <- do.call(rbind, lapply(prepared_group, function(b) b$snp_info))
  snp_info <- snp_info[match(snp_names, snp_info$SNP), , drop = FALSE]

  table <- data.frame(
    SNP = snp_names,
    block_id = block_vec,
    screened_beta = spa_beta_keep,
    screened_tau = spa_tau_keep,
    screened_any = spa_beta_keep | spa_tau_keep,
    usable_pair = usable_keep,
    modeled_beta = beta_keep,
    modeled_tau = tau_keep,
    modeled_any = any_keep,
    P_beta = as.numeric(get_p("mean", "original", "P")),
    P_tilde_beta = as.numeric(get_p("mean", "knockoff", "P")),
    P_tau = as.numeric(get_p("ws", "original", "P")),
    P_tilde_tau = as.numeric(get_p("ws", "knockoff", "P")),
    P_norm_beta = as.numeric(get_p("mean", "original", "P_KOALA_normal")),
    P_tilde_norm_beta = as.numeric(get_p("mean", "knockoff", "P_KOALA_normal")),
    P_norm_tau = as.numeric(get_p("ws", "original", "P_KOALA_normal")),
    P_tilde_norm_tau = as.numeric(get_p("ws", "knockoff", "P_KOALA_normal")),
    screening_p_method = rep(
      prepared_group[[1L]]$screening$mean$screening_p_method,
      p_all
    ),
    W_beta = as.numeric(W_beta),
    W_tau = as.numeric(W_tau),
    W_any = as.numeric(W_any),
    posterior_original_beta = stats::plogis(as.numeric(W_beta)),
    posterior_original_tau = stats::plogis(as.numeric(W_tau)),
    posterior_original_any = stats::plogis(as.numeric(W_any)),
    state_beta = expand_character(fit_beta, beta_keep,
                                  "posterior_state_class"),
    state_tau = expand_character(fit_tau, tau_keep,
                                 "posterior_state_class"),
    state_any = expand_character(fit_any, any_keep,
                                 "posterior_state_class"),
    causal_beta = ifelse(snp_names %in% causal_beta, "yes", "no"),
    causal_tau = ifelse(snp_names %in% causal_tau, "yes", "no"),
    causal_any = ifelse(snp_names %in% union(causal_beta, causal_tau),
                        "yes", "no"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  table <- cbind(
    table,
    expand_states(fit_beta, beta_keep, "posterior_beta_"),
    expand_states(fit_tau, tau_keep, "posterior_tau_"),
    expand_states(fit_any, any_keep, "posterior_any_")
  )

  list(
    group_id = group_id,
    block_ids = vapply(prepared_group, function(b) as.character(b$block_id),
                       character(1)),
    method = control$method,
    W_beta = W_beta,
    W_tau = W_tau,
    W_any = W_any,
    # Compatibility with existing KOALA FDR utilities based on |Z|-|Ztilde|.
    beta_g = pmax(W_beta, 0),
    beta_g_tilde = pmax(-W_beta, 0),
    tau_g = pmax(W_tau, 0),
    tau_g_tilde = pmax(-W_tau, 0),
    fit_beta = fit_beta,
    fit_tau = fit_tau,
    fit_any = fit_any,
    runtime_seconds = runtime_seconds,
    csv_table = table,
    snp_info = snp_info
  )
}

koala_mlr_knockoff_threshold <- function(W, q, offset = 1L) {
  W <- as.numeric(W)
  candidates <- sort(unique(abs(W[is.finite(W) & W != 0])))
  if (!length(candidates)) return(Inf)
  ratio <- vapply(candidates, function(t) {
    (offset + sum(W <= -t)) / max(sum(W >= t), 1L)
  }, numeric(1))
  ok <- which(ratio <= q)
  if (!length(ok)) Inf else candidates[min(ok)]
}

koala_mlr_evaluate <- function(table, fdr_levels = seq(0.01, 0.10, 0.01),
                               offset = 1L) {
  traits <- c(beta = "causal_beta", tau = "causal_tau", any = "causal_any")
  out <- do.call(rbind, lapply(names(traits), function(trait) {
    W <- table[[paste0("W_", trait)]]
    truth <- table[[traits[[trait]]]] == "yes"
    do.call(rbind, lapply(fdr_levels, function(q) {
      threshold <- koala_mlr_knockoff_threshold(W, q, offset)
      selected <- is.finite(threshold) & W >= threshold
      tp <- sum(selected & truth)
      fp <- sum(selected & !truth)
      data.frame(
        FDR_level = q,
        trait = trait,
        threshold = threshold,
        selected = sum(selected),
        TP = tp,
        FP = fp,
        observed_FDR = fp / max(sum(selected), 1L),
        power = tp / max(sum(truth), 1L),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(out) <- NULL
  out
}
