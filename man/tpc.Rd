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

\value{An object of \code{\link[base]{class}} "\code{pcAlgo}" (see
\code{\link[pcalg]{pcAlgo}}) containing an estimate of the equivalence class of
the underlying DAG.}

\author{Original code by Markus Kalisch, Martin Maechler, and Diego
Colombo. Modifications by Janine Witte.}

\examples{
## load simulated cohort data
data(dat_sim)
n <- nrow(dat_sim)
lab <- colnames(dat_sim)

## estimate skeleton without taking background
## information into account
tpc.fit <- tpc(suffStat = list(C = cor(dat_sim), n = n),
               indepTest = gaussCItest, alpha = 0.01,
               labels = lab)
pc.fit <- pcalg::pc(suffStat = list(C = cor(dat_sim), n = n),
                   indepTest = gaussCItest, alpha = 0.01,
                   labels = lab,
                   maj.rule = TRUE, solve.conf = TRUE)
identical(pc.fit@graph, tpc.fit@graph) # TRUE

## estimate skeleton with temporal ordering as
## background information
tiers <- rep(c(1,2,3), times=c(3,3,3))
tpc.fit2 <- tpc(suffStat = list(C = cor(dat_sim), n = n),
                indepTest = gaussCItest, alpha = 0.01,
                labels = lab, tiers = tiers)

if (require(Rgraphviz)) {
  ## compare estimated CPDAGs
  data("true_sim")
  par(mfrow = c(1,3))
  plot(true_sim, main = "True DAG")
  plot(tpc.fit, main = "PC estimate")
  plot(tpc.fit2, main = "tPC estimate")
}

## require that there is no edge between A1 and A1, and
## that any edge between A2 and B2 or A2 and C2 is
## directed away from A2
forb <- matrix(FALSE, nrow=9, ncol=9)
rownames(forb) <- colnames(forb) <- lab
forb["A1","A3"] <- forb["A3","A1"] <- TRUE
forb["B2","A2"] <- TRUE
forb["C2","A2"] <- TRUE

tpc.fit3 <- tpc(suffStat = list(C = cor(dat_sim), n = n),
                indepTest = gaussCItest, alpha = 0.01,
                labels = lab,
                forbEdges = forb, tiers = tiers)

if (require(Rgraphviz)) {
  ## compare estimated CPDAGs
  data("true_sim")
  par(mfrow = c(1,3))
  plot(tpc.fit2, main = "old tPC estimate")
  plot(tpc.fit3, main = "new tPC estimate")
}

## force edge from A1 to all other nodes measured at time 1
## into the graph (note that the edge from A1 to A2 is then
## forbidden)
tpc.fit4 <- tpc(suffStat = list(C = cor(dat_sim), n = n),
                indepTest = gaussCItest, alpha = 0.01,
                labels = lab,
                tiers = tiers, context.tier = "A1")

if (require(Rgraphviz)) {
  ## compare estimated CPDAGs
  data("true_sim")
  par(mfrow = c(1,1))
  plot(tpc.fit4, main = "alternative tPC estimate")
}

## force edge from A1 to all other nodes into the graph
tpc.fit5 <- tpc(suffStat = list(C = cor(dat_sim), n = n),
                indepTest = gaussCItest, alpha = 0.01,
                labels = lab,
                tiers = tiers, context.all = "A1")

if (require(Rgraphviz)) {
  ## compare estimated CPDAGs
  data("true_sim")
  par(mfrow = c(1,1))
  plot(tpc.fit5, main = "alternative tPC estimate")
}


}
