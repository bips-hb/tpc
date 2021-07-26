\name{dat_sim}

\alias{dat_sim, true_sim}

\docType{data}

\title{Simulated data with a partial ordering}

\description{A simple graph and corresponding dataset used in the examples
illustrating \code{\link{tpc}}.}

\usage{data("dat_sim")
Rgraphviz::plot(true_sim) }

\format{
  A data frame with 1000 observations on 9 numerical variables
  \code{A1}, \code{B1}, \code{C1}, \code{A2}, \code{B2}, \code{C2},
  \code{A3}, \code{B3}, \code{C3}, simulated by
  drawing from a multivariate distribution according to the DAG \code{true_sim}.
}

\examples{
data("dat_sim")
head(dat_sim)
Rgraphviz::plot(true_sim)
}
