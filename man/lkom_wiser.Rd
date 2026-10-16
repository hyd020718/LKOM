\name{lkom_wiser}
\alias{lkom_fit_null}
\alias{lkom_residual_channels}
\title{WiSER null-model interface and LKOM response channels}
\usage{
lkom_fit_null(X, Z, W, Y, id, julia_bin = "julia",
  julia_project = system.file("julia", package = "LKOM"),
  julia_script = file.path(julia_project, "koala_wiser_null.jl"), ...)

lkom_residual_channels(null_fit, transform = TRUE)
}
\arguments{
  \item{X,Z,W,Y,id}{Longitudinal null-model inputs.}
  \item{julia_bin}{Julia executable.}
  \item{julia_project}{Directory containing the locked Julia environment.}
  \item{julia_script}{Path to the LKOM Julia helper.}
  \item{null_fit}{Fitted WiSER null object.}
  \item{transform}{Apply the rank-based inverse-normal transformation.}
  \item{...}{Additional null-model fitting arguments.}
}
\description{
Fit the longitudinal WiSER null model through the bundled Julia helper and
extract aligned trajectory-mean and within-subject-variability response
channels.
}
\value{A null-fit object or a subject-by-two response matrix.}

