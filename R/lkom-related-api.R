# Related-sample and WiSER interfaces --------------------------------------

.read_fam <- function(plink_prefix) {
  path <- paste0(plink_prefix, ".fam")
  if (!file.exists(path)) stop("Missing PLINK FAM file: ", path, call. = FALSE)
  x <- utils::read.table(path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 6L) stop("PLINK FAM file must contain six columns.", call. = FALSE)
  names(x)[1:6] <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
  x
}

#' Construct an efficient GRM eigensystem
#'
#' `grm` may be a dense matrix for small analyses or a three-column sparse-GRM
#' file (`ID1`, `ID2`, `Value`) for biobank analyses. Sparse input is split
#' into connected components, avoiding a dense whole-cohort eigendecomposition.
#'
#' @param grm Dense GRM matrix or path to a sparse GRM file.
#' @param subject_ids Subject identifiers in analysis order.
#' @param psd_tolerance Numerical tolerance for negative eigenvalues.
#' @param psd_action Either clip small/negative eigenvalues or stop.
#' @param ... Additional arguments used by the sparse-GRM constructor.
#' @return An object of class `lkom_grm_structure`.
#' @export
lkom_grm_structure <- function(
    grm,
    subject_ids = NULL,
    psd_tolerance = 1e-7,
    psd_action = c("clip", "error"),
    ...) {
  psd_action <- match.arg(psd_action)
  if (is.character(grm) && length(grm) == 1L) {
    if (is.null(subject_ids)) {
      stop("subject_ids are required for a sparse GRM file.", call. = FALSE)
    }
    out <- koala_related_grm_structure(
      sparse_grm_file = grm,
      subject_ids = subject_ids,
      psd_tolerance = psd_tolerance,
      psd_action = psd_action,
      ...
    )
  } else {
    K <- as.matrix(grm)
    storage.mode(K) <- "double"
    if (!nrow(K) || nrow(K) != ncol(K) || any(!is.finite(K))) {
      stop("A dense GRM must be a finite nonempty square matrix.", call. = FALSE)
    }
    if (is.null(subject_ids)) {
      subject_ids <- rownames(K)
      if (is.null(subject_ids)) subject_ids <- sprintf("ID%05d", seq_len(nrow(K)))
    }
    subject_ids <- as.character(subject_ids)
    if (length(subject_ids) != nrow(K) || anyNA(subject_ids) || anyDuplicated(subject_ids)) {
      stop("subject_ids must uniquely identify every GRM row.", call. = FALSE)
    }
    if (!is.null(rownames(K))) {
      idx <- match(subject_ids, rownames(K))
      if (anyNA(idx)) stop("subject_ids cannot be matched to GRM row names.", call. = FALSE)
      K <- K[idx, idx, drop = FALSE]
    }
    K <- (K + t(K)) / 2
    eg <- eigen(K, symmetric = TRUE)
    if (min(eg$values) < -psd_tolerance && identical(psd_action, "error")) {
      stop("The GRM is not positive semidefinite; minimum eigenvalue = ",
           signif(min(eg$values), 6), ".", call. = FALSE)
    }
    repair <- data.frame(
      component = 1L,
      size = nrow(K),
      min_eigenvalue = min(eg$values),
      n_clipped = sum(eg$values < 0),
      stringsAsFactors = FALSE
    )
    eg$values[eg$values < 0] <- 0
    out <- list(
      subject_ids = subject_ids,
      n_families = 1L,
      n_components = 1L,
      component_sizes = nrow(K),
      psd_repair = repair,
      families = list(list(
        index = seq_len(nrow(K)),
        ids = subject_ids,
        values = eg$values,
        vectors = eg$vectors
      ))
    )
  }
  class(out) <- c("lkom_grm_structure", "list")
  out
}

#' Estimate the related-sample whitening parameter
#'
#' @param Y One- or two-channel subject-level response matrix.
#' @param structure Output from [lkom_grm_structure()].
#' @param eta_upper Upper bound for the GRM variance proportion.
#' @param optimize_tol Tolerance used in profile optimization.
#' @return A whitening fit.
#' @export
lkom_fit_whitening <- function(Y, structure, eta_upper = 0.999,
                               optimize_tol = 1e-6) {
  if (!inherits(structure, "lkom_grm_structure") &&
      !is.list(structure$families)) {
    stop("structure must be produced by lkom_grm_structure().", call. = FALSE)
  }
  out <- koala_related_fit_whitening(
    Y = Y, structure = structure, eta_upper = eta_upper,
    optimize_tol = optimize_tol
  )
  class(out) <- c("lkom_whitening", "list")
  out
}

#' Apply GLS centering and GRM whitening
#'
#' @param x Subject-aligned vector or matrix.
#' @param whitening Output from [lkom_fit_whitening()].
#' @param gls_center Remove a GLS intercept before applying `V^{-1/2}`.
#' @return A whitened matrix.
#' @export
lkom_whiten <- function(x, whitening, gls_center = TRUE) {
  koala_related_whiten(x, whitening, gls_center = gls_center)
}

#' Fit LKOM after related-sample whitening
#'
#' The same estimated linear transformation is applied to each member of every
#' original-knockoff pair, preserving paired swaps.
#'
#' @inheritParams lkom_fit
#' @param structure A GRM structure from [lkom_grm_structure()].
#' @param eta_upper Upper bound for profile estimation of `eta`.
#' @param optimize_tol Profile-optimization tolerance.
#' @return An `lkom_fit` with its whitening fit attached.
#' @export
lkom_fit_related <- function(
    G,
    G_tilde,
    Y,
    structure,
    control = lkom_control(),
    channel_names = NULL,
    screen = NULL,
    eta_upper = 0.999,
    optimize_tol = 1e-6) {
  subject_ids <- as.character(structure$subject_ids)
  if (length(subject_ids) != nrow(as.matrix(G))) {
    stop("The GRM structure and genotype matrices have different sample sizes.",
         call. = FALSE)
  }
  align_rows <- function(x, label) {
    x <- as.matrix(x)
    if (is.null(rownames(x))) return(x)
    idx <- match(subject_ids, rownames(x))
    if (anyNA(idx)) {
      stop(label, " row names do not contain all GRM subjects.", call. = FALSE)
    }
    x[idx, , drop = FALSE]
  }
  G <- align_rows(G, "G")
  G_tilde <- align_rows(G_tilde, "G_tilde")
  Y <- align_rows(Y, "Y")
  dat <- .lkom_validate_design(G, G_tilde, Y)
  whitening <- lkom_fit_whitening(
    dat$Y, structure, eta_upper = eta_upper, optimize_tol = optimize_tol
  )
  Gw <- lkom_whiten(dat$G, whitening, gls_center = TRUE)
  Gtw <- lkom_whiten(dat$G_tilde, whitening, gls_center = TRUE)
  Yw <- lkom_whiten(dat$Y, whitening, gls_center = TRUE)
  colnames(Gw) <- colnames(Gtw) <- dat$snp_ids
  fit <- lkom_fit(
    Gw, Gtw, Yw, control = control, channel_names = channel_names,
    screen = screen, center = FALSE
  )
  fit$whitening <- whitening
  fit$related_samples <- TRUE
  fit
}

#' Extract subject-level response channels from a WiSER null fit
#'
#' @param null_fit Output from [lkom_fit_null()] or the underlying WiSER helper.
#' @param transform Apply the rank-based inverse-normal transformation.
#' @return A subject-by-two response matrix with `beta` and `tau` columns.
#' @export
lkom_residual_channels <- function(null_fit, transform = TRUE) {
  if (is.null(null_fit$mean$c) || is.null(null_fit$ws$c)) {
    stop("null_fit does not contain both residual-derived channels.", call. = FALSE)
  }
  if (!is.logical(transform) || length(transform) != 1L || is.na(transform)) {
    stop("transform must be TRUE or FALSE.", call. = FALSE)
  }
  fun <- if (transform) .koala_rank_int else as.numeric
  out <- cbind(beta = fun(null_fit$mean$c), tau = fun(null_fit$ws$c))
  if (!is.null(null_fit$subject_id)) rownames(out) <- as.character(null_fit$subject_id)
  out
}

#' Fit the WiSER null model through the bundled Julia helper
#'
#' @param X,Z,W,Y,id Longitudinal-model inputs accepted by the original LKOM
#'   null-model implementation.
#' @param julia_bin Julia executable.
#' @param julia_project Julia environment directory. The bundled locked
#'   environment is used by default.
#' @param julia_script Path to the bundled Julia helper.
#' @param ... Additional arguments passed to the null-model fitter.
#' @return A WiSER null-fit object containing the two subject-level channels.
#' @export
lkom_fit_null <- function(
    X, Z, W, Y, id,
    julia_bin = "julia",
    julia_project = system.file("julia", package = "LKOM"),
    julia_script = file.path(julia_project, "koala_wiser_null.jl"),
    ...) {
  if (!nzchar(julia_project) || !dir.exists(julia_project)) {
    stop("The bundled Julia project could not be located.", call. = FALSE)
  }
  old_project <- Sys.getenv("JULIA_PROJECT", unset = NA_character_)
  Sys.setenv(JULIA_PROJECT = julia_project)
  on.exit({
    if (is.na(old_project)) Sys.unsetenv("JULIA_PROJECT")
    else Sys.setenv(JULIA_PROJECT = old_project)
  }, add = TRUE)
  fit_wiser_null(
    X = X, Z = Z, W = W, Y = Y, id = id,
    julia_bin = julia_bin, julia_script = julia_script, ...
  )
}

#' Check optional LKOM software dependencies
#'
#' @return A data frame describing available commands, R packages, and bundled
#'   Julia metadata.
#' @export
lkom_check_dependencies <- function() {
  r_packages <- c("data.table", "glmnet", "GRAB", "KnockoffScreen")
  r_rows <- lapply(r_packages, function(pkg) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    data.frame(
      component = paste0("R package: ", pkg),
      available = ok,
      version = if (ok) as.character(utils::packageVersion(pkg)) else NA_character_,
      path = if (ok) system.file(package = pkg) else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  commands <- c("julia", "plink2", "gcta64")
  command_rows <- lapply(commands, function(cmd) {
    path <- Sys.which(cmd)
    data.frame(
      component = paste0("command: ", cmd),
      available = nzchar(path),
      version = NA_character_,
      path = if (nzchar(path)) unname(path) else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  project <- system.file("julia", "Project.toml", package = "LKOM")
  manifest <- system.file("julia", "Manifest.toml", package = "LKOM")
  bundled <- data.frame(
    component = c("bundled Julia Project.toml", "bundled Julia Manifest.toml"),
    available = c(file.exists(project), file.exists(manifest)),
    version = c(NA_character_, "Julia 1.11.7"),
    path = c(project, manifest),
    stringsAsFactors = FALSE
  )
  do.call(rbind, c(r_rows, command_rows, list(bundled)))
}

