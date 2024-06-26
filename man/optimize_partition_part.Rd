% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimize_partition_part.R
\name{optimize_partition_part}
\alias{optimize_partition_part}
\title{Title}
\usage{
optimize_partition_part(
  counts,
  levels,
  ess,
  regular,
  min_score_improv = 0,
  P = as.list(1:nrow(counts) - 1),
  ldag_consistent = T,
  verbose = FALSE
)
}
\arguments{
\item{counts}{}

\item{levels}{}

\item{ess}{}

\item{P}{}

\item{verbose}{}

\item{nlev}{}

\item{lkappa}{}
}
\description{
Title
}
\examples{
levels <- list(0:1, 0:1)
counts <- cbind(c(10, 100, 100, 100), rep(10, 4))
optimize_partition_part(counts, levels, ess = 1, verbose = T)
}
