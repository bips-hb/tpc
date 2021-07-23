\name{ida_invalid}

\alias{ida_invalid}

\title{IDA adjustment sets with invalid graph}

\description{Determines adjustment sets compatible with a partially directed
graph that is not a valid CPDAG or MPDAG.}

\usage{ida_invalid(x.pos, y.pos, graphEst, method = NULL, verbose = TRUE,
plot = TRUE)}

\arguments{
\item{x.pos}{Position of exposure variable in graphEst@nodes.}

\item{y.pos}{Position of outcome varialbe in graphEst@nodes.}

\item{graphEst}{Estimated invalid partially directed graph. Usually obtained by
pc.fit@graph, where pc.fit is the output of pcalg::pc.}

\item{method}{Character string specifying the method.\\
"local": returns all compatible parent sets of x.\\
"optimal": returns all compatible optimal adjustment sets relative to (x,y).\\
Adjustment for the optimal adjustment set is often more efficient, but
determining the optimal adjustment sets is more sensitive to faulty graph
estimates than determining the parent sets of x.}

\item{verbose}{Logical. If TRUE, details are printed to the console.}

\item{plot}{Logical. If TRUE, a plot is are produced.}
}

\details{The IDA algorithm (Maathuis et al. 2019, implemented as
\code{pcalg::\link[pcalg]{ida}}) requires a valid CPDAG or MPDAG as an input.
\code{ida_invalid} is an alternative for the case that the graph returned by PC
or another causal discovery algorithm is invalid. It determines all local or
optimal adjustment sets compatible with the invalid graph by trying out all
possible orientations of undirected edges in relevant parts of the graph while
ignoring the rest of the graph.}

\value{A list of adjustment sets compatible with graphEst. 'NULL' stands for
the empty set. '0' means that the effect is zero (this occurs if y is among the
parents of x in local adjustment, or if there is no path from x to y in optimal
adjustment).}

\references{Maathuis MH, Kalisch M, & Buehlmann P (2009). Estimating
high-dimensional intervention effects from observational data. The Annals of
Statistics, 37(6A), 3133-3164.}

\examples{}

