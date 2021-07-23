\name{MeeksRules}
\alias{MeeksRules}

\title{Last Step of tPC Algorithm: Apply Meek's rules}

\description{This is a modified version of
\code{pcalg::\link[pcalg:udag2pdag]{udag2pdagRelaxed}}. It applies Meek's rules
to the partially oriented graph obtained after orienting edges between time
points / tiers.
}

\usage{
MeeksRules(gInput, verbose = FALSE, unfVect = NULL, solve.confl = FALSE,
           rules = rep(TRUE, 4))
}

\arguments{
  \item{gInput}{"pcAlgo"-object containing skeleton and conditional indepedence
  information.}

  \item{verbose}{FALSE: No output; TRUE: Details}

  \item{unfVect}{vector containing numbers that encode ambiguous triples (as
  returned by \code{\link{tpc.cons.intern}}). This is needed in the conservative
  and majority rule PC algorithms.}

  \item{solve.confl}{if \code{TRUE}, the orientation rules work with lists for
  candidate sets and allow bi-directed edges to resolve conflicting edge
  orientations. Note that therefore the resulting object is order-independent
  but might not be a PDAG because bi-directed edges can be present.}

  \item{rules}{vector of length 4 containing \code{TRUE} or \code{FALSE} for
  each rule. \code{TRUE} in position i means that rule i (Ri) will be applied.
  By default, all rules are used.}
}

\details{If \code{unfVect = NULL} (no ambiguous triples), the four orientation
rules are applied to each eligible structure until no more edges can be
oriented. Otherwise, unfVect contains the numbers of all ambiguous triples in
the graph as determined by \code{\link{tpc.cons.intern}}. Then the orientation
rules take this information into account. For example, if a -> b - c and <a,b,c>
is an unambigous triple and a non-v-structure, then rule 1 implies b -> c. On
the other hand, if a -> b - c but <a,b,c> is an ambiguous triple, then the edge
b - c is not oriented.

If \code{solve.confl = FALSE}, earlier edge orientations are overwritten by
later ones.

If \code{solv.confl = TRUE}, both the v-structures and the orientation rules
work with lists for the candidate edges and allow bi-directed edges if there are
conflicting orientations. For example, two v-structures a -> b <- c and
b -> c <- d then yield a -> b <-> c <-d. This option can be used to get an
order-independent version of the PC algorithm (see Colombo and Maathuis (2014)).

We denote bi-directed edges, for example between two variables i and j, in the
adjacency matrix M of the graph as M[i,j]=2 and M[j,i]=2. Such edges should be
interpreted as indications of conflicts in the algorithm, for example due to
errors in the conditional independence tests or violations of the faithfulness
assumption.
}

\value{An object of class \code{\link[pcalg]{pcAlgo-class}}
}
\references{C. Meek (1995). Causal inference and causal explanation with
background knowledge. In: Proceedings of the Eleventh Conference on Uncertainty
in Artificial Intelligence (UAI-95), pp. 403-411. Morgan Kaufmann Publishers.

D. Colombo and M.H. Maathuis (2014). Order-independent constraint-based causal
structure learning. Journal of Machine Learning Research 15:3741-3782.}

\author{Original code by Markus Kalisch, modifications by Janine Witte.}

\examples{}
