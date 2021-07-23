\name{tskeleton}

\alias{tskeleton}

\title{Estimate the Skeleton of a DAG while Accounting for a Partial Ordering}

\description{Like \code{pcalg::\link[pcalg]{skelton}}, but takes into account a
user-specified partial node ordering. The conditional independence between
\code{x} and \code{y} given \code{S} is not tested if any variable in \code{S}
lies in the future of both \code{x} and \code{y}.}

\usage{
tskeleton(suffStat, indepTest, alpha, labels, p,
          method = c("stable", "original"), m.max = Inf,
          fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
          tiers = NULL, verbose = FALSE)
}

\arguments{
  \item{suffStat}{A list of sufficient statistics, containing all necessary
  elements for the conditional independence decisions in the function
  \code{indepTest}.}

  \item{indepTest}{Predefined \code{\link[base]{function}} for testing
  conditional independence. It is internally called as
  \code{indepTest(x,y,S,suffStat)}, and tests conditional independence of
  \code{x} and \code{y} given \code{S}. Here, \code{x} and \code{y} are
  variables, and \code{S} is a (possibly empty) vector of variables (all
  variables are denoted by their (integer) column positions in the adjacency
  matrix). \code{suffStat} is a list, see the argument above. The return value
  of \code{indepTest} is the p-value of the test for conditional independence.}

  \item{alpha}{significance level (number in \emph{(0,1)} for the individual
  conditional independence tests.}

  \item{labels}{(optional) character vector of variable (or "node") names.
  Typically preferred to specifying \code{p}.}

  \item{p}{(optional) number of variables (or nodes). May be specified if
  \code{labels} are not, in which case \code{labels} is set to \code{1:p}.}

  \item{method}{Character string specifying method; the default, "stable"
  provides an \emph{order-independent} skeleton, see 'Details' below.}

  \item{m.max}{Maximal size of the conditioning sets that are considered in the
  conditional independence tests.}

  \item{fixedGaps}{logical \emph{symmetric} matrix of dimension p*p. If entry
  \code{[i,j]} is true, the edge \emph{i-j} is removed before starting the
  algorithm. Therefore, this edge is guaranteed to be \emph{absent} in the
  resulting graph.}

  \item{fixedEdges}{a logical \emph{symmetric} matrix of dimension p*p. If entry
  \code{[i,j]} is true, the edge \emph{i-j} is never considered for removal.
  Therefore, this edge is guaranteed to be \emph{present} in the resulting
  graph.}

  \item{NAdelete}{logical needed for the case \code{indepTest(*)} returns
  \code{NA}. If it is true, the corresponding edge is deleted, otherwise not.}

  \item{tiers}{Numeric vector specifying the tier / time point for each
  variable. Must be of length 'p', if specified, or have the same length as
  'labels',if specified. A smaller number corresponds to an earlier tier / time
  point. Conditional independence testing is restricted such that if x is in
  tier t(x) and y is in t(y), only those variables are allowed in the
  conditioning set whose tier is not larger than max(t(x), t(y)).}

  \item{verbose}{if \code{TRUE}, detailed output is provided.}
}

\details{See \code{pcalg::\link[pcalg]{skeleton}} for further information on the
skeleton algorithm.}

\return{An object of class "pcAlgo" (see \code{pcalg::\link[pcalg]{pcAlgo}})
containing an estimate of the skeleton of the underlying DAG, the conditioning
sets (sepset) that led to edge removals and several other parameters.}

\author{Original code by Markus Kalisch, Martin Maechler, Alain Hauser and
Diego Colombo. Modifications by Janine Witte.}


\examples{
# load numeric cohort data
data(dat_cohort_num)


n <- nrow(gmG8$x)
V <- colnames(gmG8$x) # labels aka node names

## estimate skeleton
skel.fit <- skeleton(suffStat = list(C = cor(dat_), n = n),
                     indepTest = mixCItest, alpha = 0.01, labels = V,
                     verbose = TRUE)
if (require(Rgraphviz)) {
  ## show estimated Skeleton
  par(mfrow=c(1,2))
  plot(skel.fit, main = "Estimated Skeleton")
  plot(gmG8$g, main = "True DAG")
}
}
