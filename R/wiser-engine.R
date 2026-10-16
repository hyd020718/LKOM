# Minimal R interface to the bundled Julia WiSER null-model helper.

.koala_as_matrix <- function(x, n, name) {
  if (is.null(x)) return(matrix(numeric(), nrow = n, ncol = 0L))
  if (is.vector(x) && !is.list(x)) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (nrow(x) != n) {
    stop(sprintf("%s must have %d rows; found %d.", name, n, nrow(x)))
  }
  if (any(!is.finite(x))) stop(name, " contains non-finite values.")
  x
}

.koala_strip_intercept <- function(x, name, tol = 1e-12, verbose = TRUE) {
  if (!ncol(x)) return(x)
  intercept <- vapply(seq_len(ncol(x)), function(j) {
    max(abs(x[, j] - 1)) <= tol
  }, logical(1L))
  if (sum(intercept) > 1L) stop(name, " contains multiple intercept columns.")
  if (any(intercept)) {
    if (isTRUE(verbose)) {
      message(name, ": removing the supplied intercept; Julia adds one internally.")
    }
    x <- x[, !intercept, drop = FALSE]
  }
  x
}

.koala_read_vec <- function(path, n = NULL) {
  if (!file.exists(path)) stop("Missing Julia output: ", path)
  x <- scan(path, what = numeric(), quiet = TRUE)
  if (!is.null(n) && length(x) != n) {
    stop("Unexpected length in ", basename(path), ": expected ", n,
         ", found ", length(x), ".")
  }
  x
}

.koala_read_mat <- function(path, nr, nc) {
  x <- .koala_read_vec(path)
  if (length(x) != nr * nc) {
    stop("Unexpected element count in ", basename(path), ": expected ",
         nr * nc, ", found ", length(x), ".")
  }
  matrix(x, nrow = nr, ncol = nc, byrow = TRUE)
}

.koala_rank_int <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || any(!is.finite(x))) {
    stop("INT input must be a nonempty finite numeric vector.")
  }
  stats::qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
}

fit_wiser_null <- function(
    X, Z, W, Y, id,
    julia_bin = "julia",
    julia_script = "koala_wiser_null.jl",
    runs = 10L,
    maxeval = 4000L,
    work_dir = NULL,
    keep_work_dir = FALSE,
    verbose = TRUE) {
  Y <- as.numeric(Y)
  n <- length(Y)
  if (length(id) != n) stop("id and Y must have the same length.")
  if (any(!is.finite(Y))) stop("Y contains non-finite values.")
  if (anyNA(id)) stop("id cannot contain missing values.")

  X <- .koala_strip_intercept(.koala_as_matrix(X, n, "X"), "X", verbose = verbose)
  Z <- .koala_strip_intercept(.koala_as_matrix(Z, n, "Z"), "Z", verbose = verbose)
  W <- .koala_strip_intercept(.koala_as_matrix(W, n, "W"), "W", verbose = verbose)
  if (!file.exists(julia_script)) stop("Missing Julia helper: ", julia_script)
  julia_script <- normalizePath(julia_script, winslash = "/", mustWork = TRUE)

  id_chr <- as.character(id)
  subject_id <- unique(id_chr)
  m <- length(subject_id)
  sid_index <- match(id_chr, subject_id)
  internal_id <- paste0("LKOM_", sid_index)
  ord <- order(sid_index, seq_len(n))

  dat <- data.frame(
    id = internal_id[ord],
    y = Y[ord],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  add_design <- function(dat, design, prefix) {
    if (!ncol(design)) return(dat)
    design <- design[ord, , drop = FALSE]
    colnames(design) <- paste0(prefix, seq_len(ncol(design)))
    cbind(dat, design)
  }
  dat <- add_design(dat, X, "X_")
  dat <- add_design(dat, Z, "Z_")
  dat <- add_design(dat, W, "W_")

  auto_work_dir <- is.null(work_dir)
  if (auto_work_dir) work_dir <- tempfile("lkom_wiser_null_")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  input_csv <- file.path(work_dir, "longitudinal_data.csv")
  out_dir <- file.path(work_dir, "null_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(dat, input_csv, row.names = FALSE, quote = TRUE)
  log_file <- file.path(work_dir, "julia_wiser.log")

  if (isTRUE(verbose)) {
    message(sprintf(
      "Fitting WiSER null model: %d observations from %d subjects.", n, m
    ))
  }
  args <- c(
    shQuote(julia_script),
    shQuote(normalizePath(input_csv, winslash = "/", mustWork = TRUE)),
    shQuote(normalizePath(out_dir, winslash = "/", mustWork = TRUE)),
    as.character(as.integer(runs)),
    as.character(as.integer(maxeval))
  )
  status <- system2(
    command = julia_bin,
    args = args,
    stdout = log_file,
    stderr = log_file
  )
  if (!identical(status, 0L)) {
    log_text <- if (file.exists(log_file)) {
      paste(readLines(log_file, warn = FALSE), collapse = "\n")
    } else {
      ""
    }
    stop("Julia/WiSER null-model fitting failed. Log:\n", log_text)
  }

  meta_file <- file.path(out_dir, "meta.csv")
  if (!file.exists(meta_file)) stop("Julia did not create meta.csv.")
  meta <- utils::read.csv(meta_file, check.names = FALSE)
  required <- c("m", "q_beta", "r_tau", "q_w", "q_z")
  if (!all(required %in% names(meta)) || !nrow(meta)) {
    stop("Julia meta.csv is incomplete.")
  }
  m_julia <- as.integer(meta$m[1L])
  q_beta <- as.integer(meta$q_beta[1L])
  r_tau <- as.integer(meta$r_tau[1L])
  q_w <- as.integer(meta$q_w[1L])
  q_z <- as.integer(meta$q_z[1L])
  if (m_julia != m) {
    stop("R and Julia subject counts differ: ", m, " versus ", m_julia, ".")
  }

  julia_ids <- scan(
    file.path(out_dir, "subject_ids.tsv"), what = character(), quiet = TRUE
  )
  expected_ids <- paste0("LKOM_", seq_len(m))
  idx <- match(julia_ids, expected_ids)
  if (length(idx) != m || anyNA(idx) || anyDuplicated(idx)) {
    stop("Julia subject order cannot be mapped back to the R subject order.")
  }
  inv_idx <- match(seq_len(m), idx)

  beta_c <- .koala_read_vec(file.path(out_dir, "beta_c.tsv"), m)[inv_idx]
  beta_d <- .koala_read_vec(file.path(out_dir, "beta_d.tsv"), m)[inv_idx]
  beta_K <- .koala_read_mat(file.path(out_dir, "beta_K.tsv"), m, q_beta)[
    inv_idx, , drop = FALSE
  ]
  beta_M <- .koala_read_mat(file.path(out_dir, "beta_M.tsv"), q_beta, q_beta)
  tau_c <- .koala_read_vec(file.path(out_dir, "tau_c.tsv"), m)[inv_idx]
  tau_d <- .koala_read_vec(file.path(out_dir, "tau_d.tsv"), m)[inv_idx]
  tau_K <- .koala_read_mat(file.path(out_dir, "tau_K.tsv"), m, r_tau)[
    inv_idx, , drop = FALSE
  ]
  tau_M <- .koala_read_mat(file.path(out_dir, "tau_M.tsv"), r_tau, r_tau)

  beta_hat <- .koala_read_vec(file.path(out_dir, "beta_hat.tsv"), q_beta)
  tau_hat <- .koala_read_vec(file.path(out_dir, "tau_hat.tsv"), q_w)
  Sigma_hat <- .koala_read_mat(file.path(out_dir, "Sigma_hat.tsv"), q_z, q_z)
  eta_labels_file <- file.path(out_dir, "eta_labels.csv")
  eta_labels <- if (file.exists(eta_labels_file)) {
    utils::read.csv(eta_labels_file, stringsAsFactors = FALSE)$parameter
  } else {
    character()
  }
  spa_state_file <- file.path(out_dir, "trajgwas_spa_state.jls")
  if (!file.exists(spa_state_file)) {
    stop("The Julia output does not contain trajgwas_spa_state.jls.")
  }

  cleanup_holder <- NULL
  if (!keep_work_dir && auto_work_dir) {
    cleanup_holder <- new.env(parent = emptyenv())
    cleanup_holder$path <- work_dir
    reg.finalizer(
      cleanup_holder,
      function(e) {
        if (!is.null(e$path) && dir.exists(e$path)) {
          unlink(e$path, recursive = TRUE, force = TRUE)
        }
      },
      onexit = TRUE
    )
  }
  out <- list(
    subject_id = subject_id,
    n = n,
    m = m,
    mean = list(c = beta_c, d = beta_d, K = beta_K, M = beta_M),
    ws = list(c = tau_c, d = tau_d, K = tau_K, M = tau_M),
    null_estimates = list(
      beta = beta_hat,
      tau = tau_hat,
      Sigma_gamma = Sigma_hat,
      eta_labels = eta_labels
    ),
    spa = list(
      state_file = spa_state_file,
      julia_bin = julia_bin,
      julia_script = julia_script,
      julia_subject_index = idx,
      julia_subject_id = julia_ids,
      cleanup_holder = cleanup_holder
    ),
    work_dir = if (keep_work_dir || !auto_work_dir) work_dir else NULL,
    julia_log = if (file.exists(log_file)) {
      readLines(log_file, warn = FALSE)
    } else {
      character()
    }
  )
  class(out) <- c("lkom_wiser_null", "list")
  out
}

