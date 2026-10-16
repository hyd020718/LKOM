# Utilities for KOALA simulations with related subjects.

.koala_related_rank_int <- function(x, c = 0.5) {
  x <- as.numeric(x)
  if (length(x) < 2L || any(!is.finite(x))) {
    stop("INT input must contain at least two finite values.")
  }
  stats::qnorm((rank(x, ties.method = "average") - c) /
                 (length(x) + 1 - 2 * c))
}

koala_related_subject_ids <- function(plink_prefix) {
  fam <- .read_fam(plink_prefix)
  ids <- trimws(as.character(fam$IID))
  if (anyNA(ids) || any(ids == "") || anyDuplicated(ids)) {
    stop("The source FAM contains missing, empty, or duplicated IID values.")
  }
  ids
}

.koala_related_read_table <- function(path, required_names = NULL) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  x <- utils::read.table(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = ""
  )
  if (!is.null(required_names) && !all(required_names %in% names(x))) {
    if (ncol(x) < length(required_names)) {
      stop("File has fewer than ", length(required_names), " columns: ", path)
    }
    names(x)[seq_along(required_names)] <- required_names
  }
  x
}

# Build eigensystems for arbitrary connected components of a sparse GRM.  This
# supports a mixture of unrelated singletons, pairs, and larger real pedigrees
# without materializing or diagonalizing a dense 10,000 x 10,000 matrix.
koala_related_grm_structure <- function(
    sparse_grm_file,
    subject_ids,
    FamMode = NULL,
    family_ids = NULL,
    psd_tolerance = 1e-7,
    psd_action = c("clip", "error")) {

  subject_ids <- as.character(subject_ids)
  psd_action <- match.arg(psd_action)
  if (anyDuplicated(subject_ids)) stop("subject_ids contains duplicates.")

  sparse <- .koala_related_read_table(
    sparse_grm_file, c("ID1", "ID2", "Value")
  )
  sparse <- sparse[, c("ID1", "ID2", "Value"), drop = FALSE]
  sparse$ID1 <- trimws(as.character(sparse$ID1))
  sparse$ID2 <- trimws(as.character(sparse$ID2))
  sparse$Value <- suppressWarnings(as.numeric(sparse$Value))
  if (anyNA(sparse$ID1) || anyNA(sparse$ID2) ||
      any(!is.finite(sparse$Value))) {
    stop("SparseGRM_list.txt contains invalid IDs or GRM values.")
  }
  if (length(setdiff(union(sparse$ID1, sparse$ID2), subject_ids))) {
    stop("Sparse GRM contains IDs that are absent from the source FAM.")
  }

  # Union-find over nonzero off-diagonal entries recovers the exact block
  # structure of the sparse GRM, including singleton components.
  n <- length(subject_ids)
  id_index <- stats::setNames(seq_len(n), subject_ids)
  parent <- seq_len(n)
  rank <- integer(n)
  find_root <- function(x) {
    y <- x
    while (parent[[y]] != y) y <- parent[[y]]
    while (parent[[x]] != x) {
      z <- parent[[x]]
      parent[[x]] <<- y
      x <- z
    }
    y
  }
  union_root <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra == rb) return(invisible(NULL))
    if (rank[[ra]] < rank[[rb]]) {
      parent[[ra]] <<- rb
    } else {
      parent[[rb]] <<- ra
      if (rank[[ra]] == rank[[rb]]) rank[[ra]] <<- rank[[ra]] + 1L
    }
    invisible(NULL)
  }
  off <- sparse$ID1 != sparse$ID2 & abs(sparse$Value) > psd_tolerance
  if (any(off)) {
    for (z in which(off)) {
      union_root(unname(id_index[[sparse$ID1[[z]]]]),
                 unname(id_index[[sparse$ID2[[z]]]]))
    }
  }
  roots <- vapply(seq_len(n), find_root, integer(1L))
  family_indices <- split(seq_len(n), roots, drop = TRUE)

  eigensystems <- vector("list", length(family_indices))
  names(eigensystems) <- paste0("component_", seq_along(family_indices))
  repair_log <- vector("list", length(family_indices))
  for (f in seq_along(family_indices)) {
    idx <- family_indices[[f]]
    ids <- subject_ids[idx]
    K <- diag(1, length(idx))
    dimnames(K) <- list(ids, ids)
    rows <- sparse$ID1 %in% ids & sparse$ID2 %in% ids
    sf <- sparse[rows, , drop = FALSE]
    if (nrow(sf)) {
      i <- match(sf$ID1, ids)
      j <- match(sf$ID2, ids)
      for (z in seq_len(nrow(sf))) {
        K[i[z], j[z]] <- sf$Value[z]
        K[j[z], i[z]] <- sf$Value[z]
      }
    }
    K <- (K + t(K)) / 2
    eg <- eigen(K, symmetric = TRUE)
    if (min(eg$values) < -psd_tolerance) {
      if (psd_action == "error") {
        stop("GRM component ", f, " is non-PSD (minimum eigenvalue=",
             signif(min(eg$values), 6), ").")
      }
    }
    repair_log[[f]] <- data.frame(
      component = f, size = length(idx), min_eigenvalue = min(eg$values),
      n_clipped = sum(eg$values < 0), stringsAsFactors = FALSE
    )
    eg$values[eg$values < 0] <- 0
    eigensystems[[f]] <- list(
      index = idx,
      ids = ids,
      values = eg$values,
      vectors = eg$vectors
    )
  }

  list(
    subject_ids = subject_ids,
    FamMode = FamMode,
    n_families = length(family_indices),
    n_components = length(family_indices),
    component_sizes = vapply(family_indices, length, integer(1L)),
    psd_repair = do.call(rbind, repair_log),
    families = eigensystems
  )
}

.koala_related_vinv <- function(x, structure, eta) {
  x <- as.matrix(x)
  out <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  for (fam in structure$families) {
    U <- fam$vectors
    d <- 1 - eta + eta * fam$values
    out[fam$index, ] <- U %*%
      sweep(crossprod(U, x[fam$index, , drop = FALSE]), 1L, d, "/")
  }
  out
}

.koala_related_logdet_v <- function(structure, eta) {
  sum(vapply(structure$families, function(fam) {
    sum(log(1 - eta + eta * fam$values))
  }, numeric(1L)))
}

koala_related_fit_whitening <- function(
    Y,
    structure,
    eta_upper = 0.999,
    optimize_tol = 1e-6) {

  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  n <- nrow(Y)
  k <- ncol(Y)
  if (n != length(structure$subject_ids) || !k %in% 1:2) {
    stop("Whitening requires an aligned n x 1 or n x 2 response.")
  }
  if (any(!is.finite(Y))) stop("Whitening response contains non-finite values.")
  eta_upper <- as.numeric(eta_upper)
  if (!is.finite(eta_upper) || eta_upper <= 0 || eta_upper >= 1) {
    stop("eta_upper must be a finite value in (0,1).")
  }

  one <- matrix(1, n, 1L)
  objective <- function(eta, details = FALSE) {
    VinvY <- .koala_related_vinv(Y, structure, eta)
    Vinv1 <- .koala_related_vinv(one, structure, eta)
    c_eta <- sum(Vinv1)
    if (!is.finite(c_eta) || c_eta <= 0) return(Inf)
    mu <- colSums(VinvY) / c_eta
    R <- sweep(Y, 2L, mu, "-")
    VinvR <- .koala_related_vinv(R, structure, eta)
    S <- crossprod(R, VinvR)
    Omega <- (S + t(S)) / (2 * (n - 1))
    ch <- try(chol(Omega), silent = TRUE)
    if (inherits(ch, "try-error")) return(Inf)
    value <- (k / 2) * (
      .koala_related_logdet_v(structure, eta) + log(c_eta)
    ) +
      (n - 1) * sum(log(diag(ch)))
    if (details) {
      list(
        value = value, mu = mu, Omega = Omega,
        c_eta = c_eta, vinv_one = drop(Vinv1)
      )
    } else {
      value
    }
  }

  interior <- stats::optimize(
    objective,
    interval = c(0, eta_upper),
    tol = optimize_tol
  )
  candidates <- c(0, interior$minimum, eta_upper)
  values <- vapply(candidates, objective, numeric(1L))
  eta_hat <- candidates[which.min(values)]
  detail <- objective(eta_hat, details = TRUE)
  if (!is.list(detail) || !is.finite(detail$value)) {
    stop("Profile REML failed to produce a finite whitening estimate.")
  }

  list(
    eta = eta_hat,
    objective = detail$value,
    mu = detail$mu,
    Omega = detail$Omega,
    c_eta = detail$c_eta,
    vinv_one = detail$vinv_one,
    eta_upper = eta_upper,
    optimize_tol = optimize_tol,
    structure = structure
  )
}

koala_related_whiten <- function(x, whitening, gls_center = TRUE) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  structure <- whitening$structure
  if (nrow(x) != length(structure$subject_ids)) {
    stop("The object to whiten is not aligned to the GRM subjects.")
  }
  if (any(!is.finite(x))) stop("The object to whiten contains non-finite values.")
  if (!is.logical(gls_center) || length(gls_center) != 1L || is.na(gls_center)) {
    stop("gls_center must be TRUE or FALSE.")
  }
  if (!ncol(x)) return(x)

  # GLS intercept removal must precede V^{-1/2}.  Ordinary centering after
  # whitening is not equivalent when V^{-1/2} 1 is not a constant vector,
  # as occurs in multi-generation pedigrees.
  if (gls_center) {
    vinv_one <- whitening$vinv_one
    if (is.null(vinv_one)) {
      vinv_one <- drop(.koala_related_vinv(
        matrix(1, nrow(x), 1L), structure, whitening$eta
      ))
    }
    c_eta <- sum(vinv_one)
    if (length(vinv_one) != nrow(x) || any(!is.finite(vinv_one)) ||
        !is.finite(c_eta) || c_eta <= 0) {
      stop("Invalid GLS intercept weights in the whitening fit.")
    }
    mu <- drop(crossprod(vinv_one, x)) / c_eta
    x <- sweep(x, 2L, mu, "-")
  }

  out <- matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
  eta <- whitening$eta
  for (fam in structure$families) {
    U <- fam$vectors
    d_inv_sqrt <- 1 / sqrt(1 - eta + eta * fam$values)
    out[fam$index, ] <- U %*%
      sweep(crossprod(U, x[fam$index, , drop = FALSE]), 1L, d_inv_sqrt, "*")
  }
  out
}

koala_related_response_matrix <- function(null_fit, apply_int = TRUE) {
  if (is.null(null_fit$mean$c) || is.null(null_fit$ws$c)) {
    stop("The WiSER null fit does not contain both subject-level residual channels.")
  }
  if (!is.logical(apply_int) || length(apply_int) != 1L || is.na(apply_int)) {
    stop("apply_int must be TRUE or FALSE.")
  }
  transform_channel <- if (apply_int) {
    .koala_related_rank_int
  } else {
    function(x) {
      x <- as.numeric(x)
      if (length(x) < 2L || any(!is.finite(x))) {
        stop("Residual channels must contain at least two finite values.")
      }
      x
    }
  }
  Y <- cbind(
    beta = transform_channel(null_fit$mean$c),
    tau = transform_channel(null_fit$ws$c)
  )
  rownames(Y) <- as.character(null_fit$subject_id)
  Y
}

