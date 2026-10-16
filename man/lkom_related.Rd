\name{lkom_related}
\alias{lkom_grm_structure}
\alias{lkom_fit_whitening}
\alias{lkom_whiten}
\alias{lkom_fit_related}
\title{Related-sample whitening for LKOM}
\usage{
lkom_grm_structure(grm, subject_ids = NULL, psd_tolerance = 1e-7,
  psd_action = c("clip", "error"), ...)

lkom_fit_whitening(Y, structure, eta_upper = 0.999,
  optimize_tol = 1e-6)

lkom_whiten(x, whitening, gls_center = TRUE)

lkom_fit_related(G, G_tilde, Y, structure, control = lkom_control(),
  channel_names = NULL, screen = NULL, eta_upper = 0.999,
  optimize_tol = 1e-6)
}
\arguments{
  \item{grm}{Dense GRM or path to a sparse three-column GRM file.}
  \item{subject_ids}{Identifiers in analysis order.}
  \item{psd_tolerance}{Tolerance for negative eigenvalues.}
  \item{psd_action}{Clip negative eigenvalues or stop.}
  \item{Y}{One- or two-channel subject-level response matrix.}
  \item{structure}{GRM structure returned by \code{lkom_grm_structure()}.}
  \item{eta_upper,optimize_tol}{Bounds and tolerance for profile estimation.}
  \item{x}{A subject-aligned object to whiten.}
  \item{whitening}{A fitted whitening object.}
  \item{gls_center}{Whether to remove a GLS intercept before whitening.}
  \item{G,G_tilde}{Original and knockoff genotype matrices.}
  \item{control,channel_names,screen}{Arguments passed to \code{lkom_fit()}.}
  \item{...}{Additional sparse-GRM construction arguments.}
}
\description{
Construct componentwise GRM eigensystems, estimate the variance proportion,
and apply the same linear transformation to original and knockoff genotypes.
}
\value{
A GRM structure, whitening fit, whitened matrix, or related-sample LKOM fit.
}

