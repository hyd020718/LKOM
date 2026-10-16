\name{lkom_utilities}
\alias{lkom_simulate_example}
\alias{lkom_check_dependencies}
\title{LKOM examples and dependency checks}
\usage{
lkom_simulate_example(n = 120L, p = 20L, n_beta = 4L, n_tau = 4L,
  overlap = 2L, effect_scale = 0.45, seed = 718L)

lkom_check_dependencies()
}
\arguments{
  \item{n,p}{Numbers of subjects and variants.}
  \item{n_beta,n_tau,overlap}{Causal counts and their overlap.}
  \item{effect_scale}{Toy signal scale.}
  \item{seed}{Local simulation seed.}
}
\description{
Create a small self-contained example for tutorials and tests, or inspect
optional software used by extended LKOM workflows.
}
\value{A toy data list or a dependency-status data frame.}

