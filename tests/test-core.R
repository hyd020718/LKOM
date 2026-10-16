library(LKOM)

# Known knockoff+ threshold.
W_known <- c(a = 5, b = 4, c = 3, d = -1)
stopifnot(identical(lkom_threshold(W_known, q = 0.5, offset = 1L), 3))
stopifnot(is.infinite(lkom_threshold(c(1, -1), q = 0.1, offset = 1L)))

# The screening rule is invariant to exchanging original and knockoff p-values.
p_org <- cbind(beta = c(0.2, 1e-4, 0.4), tau = c(0.3, 0.5, 0.8))
p_ko <- cbind(beta = c(0.7, 0.9, 0.6), tau = c(5e-4, 0.4, 0.9))
screen_1 <- lkom_screen_pairs(p_org, p_ko, threshold = 1e-3)
screen_2 <- lkom_screen_pairs(p_ko, p_org, threshold = 1e-3)
stopifnot(identical(screen_1, screen_2), identical(screen_1, c(TRUE, TRUE, FALSE)))

# Maximum-posterior state classification.
state_prob <- rbind(
  rs1 = c(`00` = 0.1, `10` = 0.7, `01` = 0.1, `11` = 0.1),
  rs2 = c(`00` = 0.1, `10` = 0.1, `01` = 0.2, `11` = 0.6)
)
stopifnot(identical(unname(lkom_state_class(state_prob)), c("10", "11")))

# Full statistic flip-sign property under paired column exchanges.
toy <- lkom_simulate_example(
  n = 60, p = 8, n_beta = 2, n_tau = 2, overlap = 1, seed = 21
)
ctl_det <- lkom_control(
  method = "conditional_bf", preset = "fast", max_iter = 30L,
  learn_pi = TRUE
)
fit <- lkom_fit(toy$G, toy$G_tilde, toy$Y, control = ctl_det)
swap <- c(1L, 4L, 7L)
G_swap <- toy$G
Gt_swap <- toy$G_tilde
G_swap[, swap] <- toy$G_tilde[, swap]
Gt_swap[, swap] <- toy$G[, swap]
fit_swap <- lkom_fit(G_swap, Gt_swap, toy$Y, control = ctl_det)
expected <- fit$W
expected[swap] <- -expected[swap]
stopifnot(max(abs(fit_swap$W - expected)) < 1e-10)
stopifnot(max(abs(as.matrix(fit_swap$posterior_state) -
                      as.matrix(fit$posterior_state)), na.rm = TRUE) < 1e-10)

# A fixed Gibbs seed is reproducible and does not alter the caller's RNG state.
ctl_gibbs <- lkom_control(
  method = "gibbs", preset = "fast", seed = 2026L,
  gibbs_burnin = 2L, gibbs_samples = 5L, gibbs_thin = 1L,
  gibbs_init_from_conbf = FALSE
)
set.seed(991)
rng_before <- .Random.seed
fit_g1 <- lkom_fit(toy$G[, 1:6], toy$G_tilde[, 1:6], toy$Y,
                   control = ctl_gibbs)
stopifnot(identical(rng_before, .Random.seed))
fit_g2 <- lkom_fit(toy$G[, 1:6], toy$G_tilde[, 1:6], toy$Y,
                   control = ctl_gibbs)
stopifnot(identical(fit_g1$W, fit_g2$W))
stopifnot(identical(fit_g1$posterior_state_class,
                    fit_g2$posterior_state_class))

# Whitening centers an identity-GRM example and commutes with paired swaps.
K <- diag(nrow(toy$G))
rownames(K) <- colnames(K) <- rownames(toy$G)
grm <- lkom_grm_structure(K, subject_ids = rownames(toy$G))
white <- lkom_fit_whitening(toy$Y, grm)
Yw <- lkom_whiten(toy$Y, white)
stopifnot(max(abs(colMeans(Yw))) < 1e-10)
Gw <- lkom_whiten(toy$G, white)
Gtw <- lkom_whiten(toy$G_tilde, white)
stopifnot(isTRUE(all.equal(Gtw[, swap], lkom_whiten(G_swap, white)[, swap])))
stopifnot(isTRUE(all.equal(Gw[, swap], lkom_whiten(Gt_swap, white)[, swap])))

