\name{tpc}

\alias{tpc}

\title{PC Algorithm Accounting for a Partial Node Ordering}

\description{Like \code{pcalg::\link[pcalg]{pc}}, but takes into account a
user-specified partial ordering of the nodes/variables. This has two effects:
1) The conditional independence between \code{x} and \code{y} given \code{S} is
ot tested if any variable in \code{S} lies in the future of both \code{x} and
\code{y}; 2) edges cannot be oriented from a higher-order to a lower-order node.
In addition, the user may specify individual forbidden edges and context
variables.}

\usage{
tpc(suffStat, indepTest, alpha, labels, p,
    forbEdges = NULL, m.max = Inf,
    conservative = FALSE, maj.rule = TRUE,
    tiers = NULL, context.all = NULL, context.tier = NULL,
    verbose = FALSE)
}

\arguments{
\item{suffStat}{A \code{\link[base]{list}} of sufficient statistics, containing
all necessary elements for the conditional independence decisions in the
function \code{indepTest}.}

\item{indepTest}{A \code{\link[base]{function}} for testing conditional
independence. It is internally called as \code{indepTest(x,y,S,suffStat)}, and
tests conditional independence of \code{x} and \code{y} given \code{S}. Here,
\code{x} and \code{y} are variables, and \code{S} is a (possibly empty) vector
of variables (all variables are denoted by their (integer) column positions in
the adjacency matrix). \code{suffStat} is a list, see the argument above. The
return value of \code{indepTest} is the p-value of the test for conditional
independence.}

\item{alpha}{significance level (number in \emph{(0,1)} for the individual
conditional independence tests.}

\item{labels}{(optional) character vector of variable (or "node") names.
Typically preferred to specifying \code{p}.}

\item{p}{(optional) number of variables (or nodes). May be specified if
\code{labels} are not, in which case \code{labels} is set to \code{1:p}.}

\item{forbEdges}{A logical matrix of dimension p*p. If [i,j] is TRUE, then the
directed edge i->j is forbidden. If both [i,j] and [j,i] are TRUE, then any type
of edge between i and j is forbidden.}

\item{m.max}{Maximal size of the conditioning sets that are considered in the
conditional independence tests.}

\item{conservative}{Logical indicating if conservative PC should be used.
Defaults to FALSE. See \code{\link[pcalg]{pc}} for details.}

\item{maj.rule}{Logical indicating if the majority rule should be used. Defaults
to TRUE. See \code{\link[pcalg]{pc}} for details.}

\item{tiers}{Numeric vector specifying the tier / time point for each variable.
Must be of length 'p', if specified, or have the same length as 'labels', if
specified. A smaller number corresponds to an earlier tier / time point.}

\item{context.all}{Numeric or character vector. Specifies the positions or names
of global context variables. Global context variables have no incoming edges,
i.e. no parents, and are themselves parents of all non-context variables in the
graph.}

\item{context.tier}{Numeric or character vector. Specifies the positions or
names of tier-specific context variables. Tier-specific context variables have
no incoming edges, i.e. no parents, and are themselves parents of all
non-context variables in the same tier.}

\item{verbose}{if \code{TRUE}, detailed output is provided.}
}

\details{See \code{pcalg::\link[pcalg]{pc}} for further information on the PC
algorithm.

Specifying a tier for each variable using the \code{tier} argument has the
following effects: 1) In the skeleton phase and v-structure learing phases,
conditional independence testing is restricted such that if x is in tier t(x)
and y is in t(y), only those variables are allowed in the conditioning set whose
tier is not larger than max(t(x), t(y)). 2) Following the v-structure phase, all
edges that were found between two tiers are directed into the direction of the
higher-order tier. If context variables are specified using \code{context.all}
and/or \code{context.tier}, the corresponding orientations are added in this
step.}

\item{return}{An object of \code{\link[base]{class}} "\code{pcAlgo}" (see
\code{\link[pcalg]{pcAlgo}}) containing an estimate of the equivalence class of
the underlying DAG.}

\item{author}{Original code by Markus Kalisch, Martin Maechler, and Diego
Colombo. Modifications by Janine Witte.}

\references{
%% ~put references to the literature/web site here ~
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
