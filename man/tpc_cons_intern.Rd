\name{tpc_cons_intern}

\alias{tpc_cons_intern}

\title{Utility for Conservative and Majority Rule in tpc}

\description{Like \code{pcalg::\link[pcalg]{pc.cons.intern}}, but takes into
account the user-specified partial node/variable ordering.}

\usage{
tpc_cons_intern(sk, suffStat, indepTest, alpha, version.unf = c(NA, NA),
                maj.rule = FALSE, forbEdge s= NULL, tiers = NULL,
                context.all = NULL, context.tier = NULL, verbose = FALSE))
}

\arguments{
\item{sk}{A skeleton object as returned from
\code{pcalg::\link[pcalg]{skeleton}}.}

\item{suffStat}{Sufficient statistic: List containing all relevant elements for
the conditional independence decisions.}

\item{indepTest}{Pre-defined \code{\link[base]{function}} for testing
conditional independence. The function is internally called as
\code{indepTest(x,y,S,suffStat)}, and tests conditional independence of \code{x}
and \code{y} given \code{S}. Here, \code{x} and \code{y} are variables, and
\code{S} is a (possibly empty) vector of variables (all variables are denoted by
their (integer) column positions in the adjacency matrix). The return value of
\code{indepTest} is the p-value of the test for conditional independence.}

\item{alpha}{significance level for the individual conditional independence
tests.}

\item{version.unf}{Vector of length two. If \code{version.unf[2]==1}, the
inititial separating set found by the PC/FCI algorithm is added to the set of
separating sets; if \code{version.unf[2]==2}, it is not added. In the latter
case, if the set of separating sets is empty, the triple is marked as
unambiguous if \code{version.unf[1]==1}, and as ambiguous if
\code{version.unf[1]==2}.}

\item{maj.rule}{Logical indicating if the triples are checked for ambiguity
using the majority rule idea, which is less strict than the standard
conservative method.}

\item{forbEdges}{A logical matrix of dimension p*p. If [i,j] is TRUE, then the
directed edge i->j is forbidden. If both [i,j] and [j,i] are TRUE, then any type
of edge between i and j is forbidden.}

\item{tiers}{Numeric vector specifying the tier / time point for each variable.
A smaller number corresponds to an earlier tier / time point.}

\item{context.all}{Numeric or character vector. Specifies the positions or names
of global context variables. Global context variables have no incoming edges,
i.e. no parents, and are themselves parents of all non-context variables in the
graph.}

\item{context.tier}{Numeric or character vector. Specifies the positions or
names of tier-specific context variables. Tier-specific context variables have
no incoming edges, i.e. no parents, and are themselves parents of all
non-context variables in the same tier.}

\item{verbose}{Logical asking for detailed output.}
}

\details{
See \code{pcalg::\link[pcalg]{pc.cons.intern}} for further information on the
majority and conservative approaches to learning v-structures.

Specifying a tier for each variable using the \code{tier} argument has the
following effects:

1) Only those triples x-y-z are considered as potential v-structures that
satisfy t(y)=max(t(x),t(z)). This allows for three constellations: either y is
in the same tier as x and both are later than z, or y is in the same tier as z
and both are later than x, or all three are in the same tier. Triples where y is
earlier than one or both of x and z need not be considered, as y being a
collider would be against the partial ordering. Triples where y is later than
both x and z will be oriented later in the pc algorithm and are left out here to
minimise the number of conditional independence tests.

2) Conditional independence testing is restricted such that if x is in tier t(x)
and y is in t(y), only those variables are allowed in the conditioning set whose
tier is not larger than max(t(x), t(y)).

Context variables specified via \code{context.all} or \code{context.tier} are
not considered as candidate colliders or candidate parents of colliders.
}

\value{
\item{unfTripl}{numeric vector of triples coded as numbers (via
\code{pcalg::\link[pcalg:pc_cons_intern]{triple2numb}})
that were marked as ambiguous.}

\item{sk}{The updated skeleton-object (separating sets might have been
updated).}
}

\author{Original code by Markus Kalisch and Diego Colombo. Modifications by
Janine Witte.}
