#' Estimate the Skeleton of a DAG while Accounting for a Partial Ordering
#'
#' Like \code{pcalg::\link[pcalg]{skeleton}}, but takes a user-specified partial node
#' ordering into account. The conditional independence
#' between \code{x} and \code{y} given \code{S} is not tested if any variable in
#' \code{S} lies in the future of both \code{x} and \code{y}.
#'
#' @param suffStat A list of sufficient statistics, containing all necessary elements for
#' the conditional independence decisions in the function \code{indepTest}.
#' @param indepTest Predefined \code{\link[base]{function}} for testing conditional
#' independence. It is internally called as \code{indepTest(x,y,S,suffStat)}, and tests
#' conditional independence of \code{x} and \code{y} given \code{S}. Here, \code{x} and
#' \code{y} are variables, and \code{S} is a (possibly empty) vector of variables (all
#' variables are denoted by their (integer) column positions in the adjacency matrix).
#' \code{suffStat} is a list, see the argument above. The return value of \code{indepTest}
#' is the p-value of the test for conditional independence.
#' @param alpha Significance level (number in \emph{(0,1)} for the individual conditional
#' independence tests.
#' @param labels (optional) character vector of variable (or "node") names.
#' Typically preferred to specifying \code{p}.
#' @param p (optional) number of variables (or nodes). May be specified if \code{labels}
#' are not, in which case \code{labels} is set to \code{1:p}.
#' @param method Character string specifying method; the default, "stable" provides an
#' \emph{order-independent} skeleton, see 'Details' below.
#' @param m.max Maximal size of the conditioning sets that are considered in the
#' conditional independence tests.
#' @param fixedGaps logical \emph{symmetric} matrix of dimension \code{p*p}. If entry
#' \code{[i,j]} is true, the edge \emph{i-j} is removed before starting the
#' algorithm. Therefore, this edge is guaranteed to be \emph{absent} in the
#' resulting graph.
#' @param fixedEdges a logical \emph{symmetric} matrix of dimension \code{p*p}. If entry
#' \code{[i,j]} is true, the edge \emph{i-j} is never considered for removal.
#' Therefore, this edge is guaranteed to be \emph{present} in the resulting graph.
#' @param NAdelete logical needed for the case \code{indepTest(*)} returns \code{NA}.
#' If it is true, the corresponding edge is deleted, otherwise not.
#' @param tiers Numeric vector specifying the tier / time point for each variable.
#' Must be of length 'p', if specified, or have the same length as 'labels', if specified.
#' A smaller number corresponds to an earlier tier / time point. Conditional independence
#' testing is restricted such that if \code{x} is in tier \code{t(x)} and \code{y} is
#' in \code{t(y)}, only those variables are allowed in the conditioning set whose tier is
#' not larger than \code{t(x)}.
#' @param verbose if \code{TRUE}, detailed output is provided.
#'
#' @details See \code{pcalg::\link[pcalg]{skeleton}} for further information on the
#' skeleton algorithm.
#'
#' @return An object of class "pcAlgo" (see \code{pcalg::\link[pcalg]{pcAlgo}})
#' containing an estimate of the skeleton of the underlying DAG, the conditioning
#' sets (sepset) that led to edge removals and several other parameters.
#'
#' @author
#' Original code by Markus Kalisch, Martin Maechler, Alain Hauser and Diego Colombo.
#' Modifications by Janine Witte.
#'
#' @importFrom methods as new
#'
#' @export
#'
#' @examples
#' # load simulated cohort data
#' data("dat_sim")
#' n <- nrow(dat_sim)
#' lab <- colnames(dat_sim)
#' # estimate skeleton without taking background information into account
#' tskel.fit <- tskeleton(suffStat = list(C = cor(dat_sim), n = n),
#'                        indepTest = gaussCItest, alpha = 0.01, labels = lab)
#' skel.fit <- pcalg::skeleton(suffStat = list(C = cor(dat_sim), n = n),
#'                             indepTest = gaussCItest, alpha = 0.01, labels = lab)
#'                             identical(skel.fit@graph, tskel.fit@graph) # TRUE
#'
#'# estimate skeleton with temporal ordering as background information
#' tiers <- rep(c(1,2,3), times=c(3,3,3))
#' tskel.fit2 <- tskeleton(suffStat = list(C = cor(dat_sim), n = n),
#'                        indepTest = gaussCItest, alpha = 0.01, labels = lab, tiers = tiers)
#'
#' # in this case, the skeletons estimated with and without
#' # background knowledge are identical, but fewer conditional
#' # independence tests were performed when background
#' # knowledge was taken into account
#' identical(tskel.fit@graph, tskel.fit2@graph) # TRUE
#' tskel.fit@n.edgetests
#' tskel.fit2@n.edgetests
#'
#'

tskeleton <- function (suffStat, indepTest, alpha, labels, p,
                       method = c("stable", "original"), m.max = Inf,
                       fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                       tiers = NULL, verbose = FALSE) {

     cl <- match.call()
     if (!missing(p))
        stopifnot(is.numeric(p), length(p <- as.integer(p)) ==
                     1, p >= 2)
     if (missing(labels)) {
        if (missing(p))
           stop("need to specify 'labels' or 'p'")
        labels <- as.character(seq_len(p))
     }   else {
        stopifnot(is.character(labels))
        if (missing(p))
           p <- length(labels)
        else if (p != length(labels))
           stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
     }
     seq_p <- seq_len(p)
     method <- match.arg(method)
     if (is.null(fixedGaps)) {
        G <- matrix(TRUE, nrow = p, ncol = p)
     } else if (!identical(dim(fixedGaps), c(p, p)))
        stop("Dimensions of the dataset and fixedGaps do not agree.")
     else if (!identical(fixedGaps, t(fixedGaps)))
        stop("fixedGaps must be symmetric")
     else G <- !fixedGaps
     diag(G) <- FALSE
     #################################################
     ## if no tiers are specified, everything is tier 0
     if (is.null(tiers)) {
        tiers <- rep(0, p)
     } else {
        ## check if 'tiers' are correctly specified
        if (!is.numeric(tiers)) {stop("'tiers' must be a numeric vector")}
        if (length(tiers) != p) {stop("length of 'tiers' does not match 'p' or length of 'labels'")}
     }
     #################################################
     if (any(is.null(fixedEdges))) {
        fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
     }
     else if (!identical(dim(fixedEdges), c(p, p)))
        stop("Dimensions of the dataset and fixedEdges do not agree.")
     else if (!identical(fixedEdges, t(fixedEdges)))
        stop("fixedEdges must be symmetric")

    pval <- NULL
    # seq_p is just the vector 1:p
    # sepset is a list of p lists with p elements each,
    # so each element represents an edge, and each edge is represented twice
    sepset <- lapply(seq_p, function(.) vector("list", p))
    # pMax is a matrix with one p-value per edge, at the beginning all p-values
    # are -Inf
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    # ord is the size of the conditioning set
    ord <- 0L
    # n.edgetests is for recording how many cond. ind. tests have been conducted
    # in total
    n.edgetests <- numeric(1)
    # G is a (pxp)-matrix (each entry represents an edge and each edge is
    # represented twice); at the beginning, all elements are TRUE except for the
    # diagonal
    while (!done && any(G) && ord <= m.max) {
      # done is FALSE if for every remaining edge, the number of neighbours is
      # smaller than the new ord
       n.edgetests[ord1 <- ord + 1L] <- 0
       done <- TRUE
       # ind is a two-column matrix, each row represents an edge (indices of
       # both endpoints)
       ind <- which(G, arr.ind = TRUE)
       # the next command just reorders ind
       ind <- ind[order(ind[, 1]), ]
       # how many edges are remaining?
       remEdges <- nrow(ind)
       if (verbose)
          cat("Order=", ord, "; remaining edges:", remEdges,
              "\n", sep = "")
       if (method == "stable") {
         # G is split into p vectors, each vector respresenting the neighbours
         # of one node
          G.l <- split(G, gl(p, p))
       }
       for (i in 1:remEdges) {
         # every edge is visited twice, so that each of the endpoints gets to
         # be the node whose neighbours are considered
          if (verbose && (verbose >= 2 || i%%100 == 0))
             cat("|i=", i, "|iMax=", remEdges, "\n")
          # endpoint 1 of current edge
          x <- ind[i, 1]
          # endpoint 2 of current edge
          y <- ind[i, 2]
          if (G[y, x] && !fixedEdges[y, x]) {
             # only edges are considered that are still in the current skeleton
             # else go to next remaining edge
             # in nbrsBool, the neighbours of the current node are TRUE
             nbrsBool <- if (method == "stable")
                G.l[[x]] #
             else G[, x]
             #################################################
             # this excludes neighbours in a later tier than x from the
             # conditioning set
             nbrsBool[tiers > tiers[x]] <- FALSE
             #################################################
             nbrsBool[y] <- FALSE
             # nbrs contains the indices of all eligible neighbours
             nbrs <- seq_p[nbrsBool]
             length_nbrs <- length(nbrs)
             # next steps only possible if there are enough neighbours to form
             # conditioning sets of cardinality length_nbrs
             if (length_nbrs >= ord) { #else go to next remaining edge
                if (length_nbrs > ord)
                  # done is reset to FALSE if for any node with remaining edges,
                  # the number of neighbours is at least as large as the order
                  # that will come next
                   done <- FALSE
                S <- seq_len(ord)
                repeat { # the repeat loop goes over all subsets of the
                  # neighbours with length ord
                   n.edgetests[ord1] <- n.edgetests[ord1] +
                      1
                   pval <- indepTest(x, y, nbrs[S], suffStat)
                   if (verbose)
                      cat("x=", x, " y=", y, " S=", nbrs[S],
                          ": pval =", pval, "\n")
                   if (is.na(pval))
                     pval <- as.numeric(NAdelete)
                   # pMax is the maximum p-value of all the tests conditioning
                   # on different subsets of the neighbours
                   # what is pMax for?
                   if (pMax[x, y] < pval)
                      pMax[x, y] <- pval
                   if (pval >= alpha) {
                      G[x, y] <- G[y, x] <- FALSE
                      sepset[[x]][[y]] <- nbrs[S]
                      break # exit repeat loop (?)
                   }
                   else {
                      # chose ord elements from the neighbours as the new S
                      nextSet <- getNextSet(length_nbrs, ord,
                                            S)
                      if (nextSet$wasLast)
                         break
                      S <- nextSet$nextSet
                   }
                } # end repeat
             } # end if
          } # end if
       } # end for
       ord <- ord + 1L
    }
    for (i in 1:(p - 1)) {
       for (j in 2:p) pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j, i])
    }
   Gobject <- if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
   } else {
      colnames(G) <- rownames(G) <- labels
      as(G, "graphNEL")
   }
   new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
       max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
       sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}
