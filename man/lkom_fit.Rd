\name{lkom_fit}
\alias{lkom_fit}
\alias{lkom}
\alias{lkom_control}
\alias{print.lkom_fit}
\alias{summary.lkom_fit}
\title{Fit the LKOM masked-likelihood-ratio model}
\usage{
lkom_control(method = c("gibbs", "conditional_bf", "variational"),
  preset = c("paper", "fast"), seed = 718L, ...)

lkom_fit(G, G_tilde, Y, control = lkom_control(), channel_names = NULL,
  screen = NULL, center = TRUE)

lkom(G, G_tilde, Y, q = 0.05, offset = 1L,
  control = lkom_control(), channel_names = NULL, screen = NULL,
  center = TRUE)

\method{print}{lkom_fit}(x, ...)

\method{summary}{lkom_fit}(object, ...)
}
\arguments{
  \item{G}{Subject-by-variant matrix of original genotypes.}
  \item{G_tilde}{Matching matrix of model-X knockoff genotypes.}
  \item{Y}{Subject-by-channel response matrix with one or two columns.}
  \item{control}{Control list returned by \code{lkom_control()}.}
  \item{channel_names}{Names for the one or two response channels.}
  \item{screen}{Logical vector produced by a swap-symmetric screen.}
  \item{center}{Whether to center before scaling. Set to false only after GLS centering and whitening.}
  \item{q}{Target false discovery rate.}
  \item{offset}{Zero for knockoff or one for knockoff+.}
  \item{method}{Inference engine.}
  \item{preset}{Manuscript or fast demonstration settings.}
  \item{seed}{Local random seed used by Gibbs sampling.}
  \item{x,object}{An object of class \code{lkom_fit}.}
  \item{...}{Additional control arguments or unused print arguments.}
}
\value{
\code{lkom_control()} returns a control list. \code{lkom_fit()} and
\code{lkom()} return objects of class \code{lkom_fit}; the latter also stores
the knockoff threshold and selected variants.
}
\description{
Fit LKOM to paired original and knockoff genotypes. The two-channel fit assigns
posterior probabilities to states 00, 10, 01, and 11.
}
\examples{
toy <- lkom_simulate_example(n = 60, p = 8, n_beta = 2, n_tau = 2,
                             overlap = 1)
ctl <- lkom_control(method = "conditional_bf", preset = "fast")
fit <- lkom(toy$G, toy$G_tilde, toy$Y, q = 0.2, control = ctl)
print(fit)
}

