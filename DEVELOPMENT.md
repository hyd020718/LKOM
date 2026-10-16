# LKOM 0.1.0 development status

## Completed

- Installable R source package with a documented public API.
- Independent-sample dual-channel MLR using the latest shared core supplied in
  the independent-sample simulation archive.
- Knockoff and knockoff+ thresholds.
- Swap-symmetric screening interface.
- Four-state posterior classification and state-aware beta/tau discoveries.
- Sparse-component and dense-GRM whitening for related samples.
- Bundled Julia 1.11.7 environment with TrajGWAS 0.4.7 and WiSER 0.2.6.
- Small deterministic example generator.
- Base-R tests for flip-sign invariance, thresholds, state classification,
  reproducible Gibbs seeds, preservation of the caller RNG state, and
  whitening/swap commutation.
- Local package installation, examples, tests, and `R CMD check`: Status OK.

## Validation environment

The package was checked under R 4.3.3 because that was the local validation
runtime available during packaging. The supplied analysis environment reports
R 4.5.3 on Rocky Linux 9.3. No compiled R code is used, so this difference is
not expected to affect installation; the target-server smoke test remains
recommended.

