# Public LKOM interface -----------------------------------------------------

#' Configure LKOM inference
#'
#' Creates the control list used by [lkom_fit()]. The `paper` preset matches
#' the main simulation configuration: Gibbs sampling, a conditional-Bayes-
#' factor warm start, inverse-Wishart residual-covariance updates, and joint
#' grid updates of the slab scales and correlation.
#'
#' @param method Inference engine: `"gibbs"`, `"conditional_bf"`, or
#'   `"variational"`.
#' @param preset Either `"paper"` or a lighter `"fast"` preset.
#' @param seed Integer random seed used locally by the Gibbs sampler.
#' @param ... Named arguments passed to the internal MLR control constructor.
#' @return An LKOM control list.
#' @export
lkom_control <- function(
    method = c("gibbs", "conditional_bf", "variational"),
    preset = c("paper", "fast"),
    seed = 718L,
    ...) {
  method <- match.arg(method)
  preset <- match.arg(preset)
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(names(dots) == ""))) {
    stop("All arguments in ... must be named.", call. = FALSE)
  }

  defaults <- list(method = method, seed = as.integer(seed))
  if (identical(method, "gibbs")) {
    defaults <- c(defaults, if (identical(preset, "paper")) {
      list(
        gibbs_init_from_conbf = TRUE,
        conbf_bidirectional = TRUE,
        gibbs_bidirectional = FALSE,
        gibbs_burnin = 100L,
        gibbs_samples = 300L,
        gibbs_thin = 1L,
        pi_bivariate = c(0.90, 0.04, 0.04, 0.02),
        omega_mode = "iw",
        s_beta_mode = "grid",
        s_tau_mode = "grid",
        rho_mode = "grid"
      )
    } else {
      list(
        gibbs_init_from_conbf = TRUE,
        conbf_bidirectional = TRUE,
        gibbs_bidirectional = FALSE,
        gibbs_burnin = 20L,
        gibbs_samples = 50L,
        gibbs_thin = 1L,
        omega_mode = "plugin",
        s_beta_mode = "fixed",
        s_tau_mode = "fixed",
        rho_mode = "fixed"
      )
    })
  } else {
    defaults <- c(defaults, list(
      omega_mode = "plugin",
      s_beta_mode = "fixed",
      s_tau_mode = "fixed",
      rho_mode = "fixed"
    ))
  }
  args <- utils::modifyList(defaults, dots, keep.null = TRUE)
  do.call(koala_mlr_control, args)
}

.lkom_validate_design <- function(G, G_tilde, Y) {
  G <- as.matrix(G)
  G_tilde <- as.matrix(G_tilde)
  Y <- as.matrix(Y)
  storage.mode(G) <- "double"
  storage.mode(G_tilde) <- "double"
  storage.mode(Y) <- "double"
  if (!identical(dim(G), dim(G_tilde))) {
    stop("G and G_tilde must have identical dimensions.", call. = FALSE)
  }
  if (nrow(G) != nrow(Y)) {
    stop("G, G_tilde, and Y must contain the same subjects.", call. = FALSE)
  }
  if (!ncol(G) || !nrow(G)) stop("The genotype matrices are empty.", call. = FALSE)
  if (!ncol(Y) %in% 1:2) stop("Y must have one or two columns.", call. = FALSE)
  if (any(!is.finite(G)) || any(!is.finite(G_tilde)) || any(!is.finite(Y))) {
    stop("G, G_tilde, and Y must contain only finite values.", call. = FALSE)
  }
  ids <- colnames(G)
  if (is.null(ids)) ids <- sprintf("SNP%05d", seq_len(ncol(G)))
  if (anyNA(ids) || any(ids == "") || anyDuplicated(ids)) {
    stop("Genotype column names must be non-missing and unique.", call. = FALSE)
  }
  if (!is.null(colnames(G_tilde)) && !identical(ids, colnames(G_tilde))) {
    stop("G and G_tilde column names must be identically ordered.", call. = FALSE)
  }
  colnames(G) <- colnames(G_tilde) <- ids
  list(G = G, G_tilde = G_tilde, Y = Y, snp_ids = ids)
}

#' Fit the LKOM masked-likelihood-ratio model
#'
#' Fits a one- or two-channel masked model to paired original and knockoff
#' genotypes. Any preliminary screen must be symmetric in the two members of
#' each original-knockoff pair; [lkom_screen_pairs()] constructs such a rule
#' from externally computed marginal p-values.
#'
#' @param G Numeric subject-by-variant matrix of original genotypes.
#' @param G_tilde Numeric matrix of paired model-X knockoffs.
#' @param Y Numeric response-channel matrix with one or two columns.
#' @param control Output from [lkom_control()].
#' @param channel_names Names of the response channels.
#' @param screen Logical vector identifying pairs retained by a symmetric
#'   preliminary screen. By default all pairs are fitted.
#' @param center Whether the engine should center matrices before scaling.
#'   Use `FALSE` only for already GLS-centered and whitened related-sample data.
#' @return An object of class `lkom_fit`.
#' @export
lkom_fit <- function(
    G,
    G_tilde,
    Y,
    control = lkom_control(),
    channel_names = NULL,
    screen = NULL,
    center = TRUE) {
  dat <- .lkom_validate_design(G, G_tilde, Y)
  p <- ncol(dat$G)
  if (is.null(channel_names)) {
    channel_names <- if (ncol(dat$Y) == 2L) c("beta", "tau") else "trait"
  }
  if (length(channel_names) != ncol(dat$Y) || anyNA(channel_names)) {
    stop("channel_names must name every column of Y.", call. = FALSE)
  }
  colnames(dat$Y) <- channel_names
  if (is.null(screen)) screen <- rep(TRUE, p)
  if (!is.logical(screen) || length(screen) != p || anyNA(screen)) {
    stop("screen must be a complete logical vector with one value per pair.",
         call. = FALSE)
  }

  if (!any(screen)) {
    out <- list(
      method = control$method,
      W = stats::setNames(numeric(p), dat$snp_ids),
      posterior_original = stats::setNames(rep(0.5, p), dat$snp_ids),
      canonical_log_odds = stats::setNames(numeric(p), dat$snp_ids),
      observed_orientation = stats::setNames(rep(NA_integer_, p), dat$snp_ids),
      identical_pair = stats::setNames(rep(FALSE, p), dat$snp_ids),
      posterior_state = data.frame(),
      posterior_state_class = stats::setNames(rep("not_screened", p), dat$snp_ids),
      screened = stats::setNames(screen, dat$snp_ids),
      channel_names = channel_names,
      control = control,
      diagnostics = list(n_subjects = nrow(dat$G), n_pairs = p,
                         n_pairs_screened = 0L),
      call = match.call()
    )
    class(out) <- "lkom_fit"
    return(out)
  }

  engine <- koala_mlr_fit(
    G = dat$G[, screen, drop = FALSE],
    G_tilde = dat$G_tilde[, screen, drop = FALSE],
    Y = dat$Y,
    control = control,
    channel_names = channel_names,
    center = center
  )
  selected_ids <- dat$snp_ids[screen]
  expand_num <- function(x, default = 0) {
    z <- stats::setNames(rep(default, p), dat$snp_ids)
    z[selected_ids] <- as.numeric(x)
    z
  }
  expand_chr <- function(x, default = "not_screened") {
    z <- stats::setNames(rep(default, p), dat$snp_ids)
    z[selected_ids] <- as.character(x)
    z
  }
  state_small <- as.data.frame(engine$posterior_state, check.names = FALSE)
  state_full <- matrix(
    NA_real_, p, ncol(state_small),
    dimnames = list(dat$snp_ids, names(state_small))
  )
  state_full[screen, ] <- as.matrix(state_small)

  engine$W <- expand_num(engine$W)
  engine$posterior_original <- expand_num(engine$posterior_original, 0.5)
  engine$canonical_log_odds <- expand_num(engine$canonical_log_odds)
  engine$observed_orientation <- expand_num(engine$observed_orientation, NA_real_)
  engine$identical_pair <- as.logical(expand_num(engine$identical_pair))
  names(engine$identical_pair) <- dat$snp_ids
  engine$posterior_state <- as.data.frame(state_full, check.names = FALSE)
  engine$posterior_state_class <- expand_chr(engine$posterior_state_class)
  engine$screened <- stats::setNames(screen, dat$snp_ids)
  engine$channel_names <- channel_names
  engine$control <- control
  engine$diagnostics$n_pairs_total <- p
  engine$diagnostics$n_pairs_screened <- sum(screen)
  engine$call <- match.call()
  class(engine) <- "lkom_fit"
  engine
}

#' Fit LKOM and apply knockoff+ selection
#'
#' @inheritParams lkom_fit
#' @param q Target false discovery rate.
#' @param offset Knockoff threshold offset; `1` gives knockoff+.
#' @return An `lkom_fit` object with threshold and selection fields.
#' @export
lkom <- function(G, G_tilde, Y, q = 0.05, offset = 1L,
                 control = lkom_control(), channel_names = NULL,
                 screen = NULL, center = TRUE) {
  fit <- lkom_fit(
    G, G_tilde, Y, control = control, channel_names = channel_names,
    screen = screen, center = center
  )
  fit$q <- q
  fit$offset <- offset
  fit$threshold <- lkom_threshold(fit$W, q = q, offset = offset)
  fit$selected <- stats::setNames(
    is.finite(fit$threshold) & fit$W >= fit$threshold,
    names(fit$W)
  )
  fit
}

#' Knockoff or knockoff+ data-dependent threshold
#'
#' @param W Numeric knockoff statistics.
#' @param q Target FDR in `(0,1)`.
#' @param offset `1` for knockoff+ and `0` for the original knockoff filter.
#' @return The selected threshold, or `Inf` when no threshold qualifies.
#' @export
lkom_threshold <- function(W, q = 0.05, offset = 1L) {
  if (length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1) {
    stop("q must be a finite scalar in (0,1).", call. = FALSE)
  }
  if (length(offset) != 1L || !offset %in% c(0L, 1L)) {
    stop("offset must be 0 or 1.", call. = FALSE)
  }
  koala_mlr_knockoff_threshold(W, q = q, offset = as.integer(offset))
}

#' Construct a swap-symmetric marginal screen
#'
#' @param p_original Vector or variant-by-channel matrix of original p-values.
#' @param p_knockoff Matching vector or matrix of knockoff p-values.
#' @param threshold Marginal screening threshold.
#' @return A logical vector, one value per original-knockoff pair.
#' @export
lkom_screen_pairs <- function(p_original, p_knockoff, threshold = 1e-3) {
  po <- as.matrix(p_original)
  pk <- as.matrix(p_knockoff)
  if (!identical(dim(po), dim(pk)) || !nrow(po)) {
    stop("Original and knockoff p-values must have identical dimensions.",
         call. = FALSE)
  }
  if (any(!is.finite(po)) || any(!is.finite(pk)) ||
      any(po < 0 | po > 1) || any(pk < 0 | pk > 1)) {
    stop("All p-values must be finite and lie in [0,1].", call. = FALSE)
  }
  if (length(threshold) != 1L || !is.finite(threshold) ||
      threshold <= 0 || threshold > 1) {
    stop("threshold must be in (0,1].", call. = FALSE)
  }
  apply(cbind(po, pk), 1L, min) < threshold
}

#' Obtain posterior association-state classifications
#'
#' @param object An `lkom_fit` object or a posterior-state probability matrix.
#' @return A named character vector of maximum-posterior states.
#' @export
lkom_state_class <- function(object) {
  if (inherits(object, "lkom_fit")) return(object$posterior_state_class)
  prob <- as.matrix(object)
  if (!nrow(prob) || !ncol(prob) || any(!is.finite(prob))) {
    stop("Posterior-state probabilities must form a finite nonempty matrix.",
         call. = FALSE)
  }
  states <- sub("^posterior_state_", "", colnames(prob))
  if (is.null(states) || any(states == "")) states <- as.character(seq_len(ncol(prob)))
  stats::setNames(states[max.col(prob, ties.method = "first")], rownames(prob))
}

#' Tabulate LKOM selections
#'
#' @param object An `lkom_fit` object.
#' @param q Target FDR if the object was created by [lkom_fit()].
#' @param offset Knockoff threshold offset.
#' @param selected_only Return only selected variants.
#' @return A variant-level data frame.
#' @export
lkom_select <- function(object, q = object$q %||% 0.05,
                        offset = object$offset %||% 1L,
                        selected_only = FALSE) {
  if (!inherits(object, "lkom_fit")) stop("object must be an lkom_fit.", call. = FALSE)
  threshold <- object$threshold %||% lkom_threshold(object$W, q, offset)
  selected <- is.finite(threshold) & object$W >= threshold
  ans <- data.frame(
    SNP = names(object$W),
    W = as.numeric(object$W),
    screened = as.logical(object$screened),
    selected = as.logical(selected),
    state = as.character(object$posterior_state_class),
    posterior_original = as.numeric(object$posterior_original),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(object$posterior_state)) {
    ans <- cbind(ans, object$posterior_state)
  }
  attr(ans, "threshold") <- threshold
  attr(ans, "q") <- q
  if (isTRUE(selected_only)) ans <- ans[ans$selected, , drop = FALSE]
  rownames(ans) <- NULL
  ans
}

#' Return joint and channel-specific discovery sets
#'
#' A selected variant contributes to the mean channel only when classified as
#' state `10` or `11`, and to the WS-variability channel only when classified
#' as `01` or `11`.
#'
#' @inheritParams lkom_select
#' @return A list containing `joint`, `beta`, and `tau` SNP identifiers.
#' @export
lkom_discoveries <- function(object, q = object$q %||% 0.05,
                             offset = object$offset %||% 1L) {
  tab <- lkom_select(object, q = q, offset = offset)
  joint <- tab$SNP[tab$selected]
  state <- tab$state
  list(
    joint = joint,
    beta = tab$SNP[tab$selected & state %in% c("10", "11", "1")],
    tau = tab$SNP[tab$selected & state %in% c("01", "11")]
  )
}

#' @export
print.lkom_fit <- function(x, ...) {
  threshold <- x$threshold %||% NA_real_
  selected <- if (!is.null(x$selected)) sum(x$selected) else NA_integer_
  cat("LKOM fit\n")
  cat("  subjects:", x$diagnostics$n_subjects, "\n")
  cat("  variant pairs:", x$diagnostics$n_pairs_total %||% length(x$W), "\n")
  cat("  screened pairs:", sum(x$screened), "\n")
  cat("  engine:", x$method, "\n")
  if (is.finite(threshold)) cat("  knockoff threshold:", signif(threshold, 5), "\n")
  if (!is.na(selected)) cat("  discoveries:", selected, "\n")
  invisible(x)
}

#' @export
summary.lkom_fit <- function(object, ...) {
  tab <- lkom_select(object)
  list(
    method = object$method,
    n_subjects = object$diagnostics$n_subjects,
    n_pairs = length(object$W),
    n_screened = sum(object$screened),
    threshold = attr(tab, "threshold"),
    n_selected = sum(tab$selected),
    state_counts = table(tab$state[tab$selected])
  )
}

