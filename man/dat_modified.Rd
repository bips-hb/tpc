\name{dat_modified}

\alias{dat_cohort_dis}
\alias{dat_miss}

\docType{data}

\title{Modified simulated cohort data}

\description{Modified versions of \code{\link{dat_cohort}}. \code{dat_miss}
contains missing values. In \code{dat_cohort_dis}, numerical
variables have been discretised according to the terciles.}

\usage{
data("dat_miss")
data("dat_cohort_dis")
}

\references{
Andrews RM, Foraita R, Witte J (2021). A practical guide to causal discovery
with cohort data. Working paper.
}

\examples{
data("dat_miss")
head(dat_miss)

data("dat_cohort_dis")
head(dat_cohort_dis)
}
