\name{lkom_selection}
\alias{lkom_threshold}
\alias{lkom_screen_pairs}
\alias{lkom_state_class}
\alias{lkom_select}
\alias{lkom_discoveries}
\title{Screen, select, and classify LKOM variants}
\usage{
lkom_threshold(W, q = 0.05, offset = 1L)

lkom_screen_pairs(p_original, p_knockoff, threshold = 1e-3)

lkom_state_class(object)

lkom_select(object, q = object$q \%||\% 0.05,
  offset = object$offset \%||\% 1L, selected_only = FALSE)

lkom_discoveries(object, q = object$q \%||\% 0.05,
  offset = object$offset \%||\% 1L)
}
\arguments{
  \item{W}{Numeric knockoff statistics.}
  \item{q}{Target false discovery rate.}
  \item{offset}{Zero for knockoff or one for knockoff+.}
  \item{p_original,p_knockoff}{Matching variant-by-channel p-value objects.}
  \item{threshold}{Marginal screening threshold.}
  \item{object}{An LKOM fit or a matrix of posterior state probabilities.}
  \item{selected_only}{Whether to return only discoveries.}
}
\description{
Utilities for symmetric marginal screening, knockoff thresholding, and
post-selection state classification. Mean-channel discoveries must have state
10 or 11, while variability-channel discoveries must have state 01 or 11.
}
\value{
The functions return, respectively, a numeric threshold, a logical screen, a
state vector, a variant-level data frame, or a list of discovery sets.
}

