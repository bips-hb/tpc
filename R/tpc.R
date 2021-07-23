#' Use PC while Accounting for a Partial Ordering
#'
#' Like \code{pcalg::\link[pcalg]{pc}}, but takes into account a user-specified
#' partial ordering. The conditional independence between \code{x} and \code{y}
#' given \code{S} is only tested if all variables in \code{S} precede at least
#' one of \code{x} and \code{y} in the partial ordering, or are in the same tier
#' as \code{x} and \code{y}. Edges cannot be oriented from a higher-order to a
#' lower-order node. In addition, context variables can be specified
#' (see below for details).
#'
#' @param suffStat A \code{\link[base]{list}} of sufficient statistics,
#' containing all necessary elements for the conditional independence decisions
#' in the function \code{indepTest}.
#' @param indepTest  A \code{\link[base]{function}} for testing conditional
#' independence. It is internally called as \code{indepTest(x,y,S,suffStat)},
#' and tests conditional independence of \code{x} and \code{y} given \code{S}.
#' Here, \code{x} and \code{y} are variables, and \code{S} is a (possibly empty)
#' vector of variables (all variables are denoted by their (integer) column
#' positions in the adjacency matrix). \code{suffStat} is a list, see the
#' argument above. The return value of \code{indepTest} is the p-value of the
#' test for conditional independence.
#' @param alpha significance level (number in \emph{(0,1)} for the individual
#' conditional independence tests.
#' @param labels (optional) character vector of variable (or "node") names.
#' Typically preferred to specifying \code{p}.
#' @param p (optional) number of variables (or nodes). May be specified if
#' \code{labels} are not, in which case \code{labels} is set to \code{1:p}.
#' @param fixedGaps logical \emph{symmetric} matrix of dimension p*p. If entry
#' \code{[i,j]} is true, the edge \emph{i-j} is removed before starting the
#' algorithm. Therefore, this edge is guaranteed to be \emph{absent} in the
#' resulting graph.
#' @param fixedEdges A logical \emph{symmetric} matrix of dimension p*p.
#' If entry \code{[i,j]} is true, the edge \emph{i-j} is never considered for
#' removal. Therefore, this edge is guaranteed to be \emph{present} in the
#' resulting graph.
#' @param m.max Maximal size of the conditioning sets that are considered in the
#' conditional independence tests.
#' @param skel.method  Character string specifying method; the default,
#' "\code{stable}" provides an \emph{order-independent} skeleton,
#' see \code{\link{skeleton}}.
#' @param conservative  Logical indicating if the conservative PC is used. In
#' this case, only option \code{u2pd = "relaxed"} is supported. Note that
#' therefore the resulting object might not be extendable to a DAG. See details
#' for more information.
#' @param maj.rule  Logical indicating that the triples shall be checked for
#' ambiguity using a majority rule idea, which is less strict than the
#' conservative PC algorithm. For more information, see details.
#' @param NAdelete If indepTest returns \code{NA} and this option is \code{TRUE}, the corresponding edge is deleted. If this option is \code{FALSE}, the edge is not deleted.
#' @param verbose  if \code{TRUE}, detailed output is provided.
#' @param tiers Numeric vector specifying the tier / time point for each variable.
#' Must be of length 'p', if specified, or have the same length as 'labels',
#' if specified. A smaller number corresponds to an earlier tier / time point.
#' Conditional independence testing is restricted such that if x is in tier t(x)
#' and y is in t(y), only those variables are allowed in the conditioning set
#' whose tier is not larger than max(t(x), t(y)).
#' @param context.all Numeric or character vector. Specifies the positions or
#' names of global context variables. Global context variables have no incoming
#' edges, i.e. no parents, and are themselves parents of all non-context
#' variables in the graph.
#' @param context.tier  Numeric or character vector. Specifies the positions or
#' names of tier-specific context variables. Tier-specific context variables
#' have no incoming edges, i.e. no parents, and are themselves parents of all
#' non-context variables in the same tier.
#'
#' @details See \code{pcalg::\link[pcalg]{pc}} for further information on the PC
#' algorithm.
#'
#' Specifying a tier for each variable using the \code{tier} argument has the
#' following effects: 1) In the skeleton phase and v-structure learing phases,
#' conditional independence testing is restricted such that if x is in tier t(x)
#' and y is in t(y), only those variables are allowed in the conditioning set
#' whose tier is not larger than max(t(x), t(y)). 2) Following the v-structure
#' phase, all edges that were found between two tiers are directed into the
#' direction of the higher-order tier. If context variables are specifiec using
#' \code{context.all} and/or \code{context.tier}, the corresponding orientations
#' are added in this step.
#'
#' A note regarding Meek's rules:\cr
#' Note that only rules 1-3 of Meek's rules are implemented, rule 4 is missing.
#' This is only appropriate when PC is used without background knowledge,
#' but here we have a partial ordering, so it would be more appropriate to use
#' all four rules and we should implement that should we ever plan to make this
#' a real package accessible for people outside our group.
#'
#' @return An object of \code{\link[base]{class}} "\code{pcAlgo}"
#' (see \code{\link[pcalg]{pcAlgo}}) containing an estimate of the equivalence
#' class of the underlying DAG.
#'
#' @author Original code by Markus Kalisch, Martin Maechler, and Diego Colombo.
#' Modifications by Janine Witte.
#'
#' @export
#'
#' @examples
#' ## Load Gaussian data from pcalg package
#' data(gmG)
#' n <- nrow(gmG8$x)
#' V <- colnames(gmG8$x)
#'
#' ## estimate CPDAG
#' tpc.fit <- tpc(suffStat = list(C = cor(gmG8$x), n = n),
#'                indepTest = gaussCItest,
#'                alpha = 0.01, labels = V,
#'                tiers = c(1,1,2,2,3,3,3,3),
#'                context.all = "Author",
#'                context.tier = "V5")
#'
#' ## estimate CPDAG using data with missing values
#' daten <- mice::windspeed[,1]
#' for(i in 2:ncol(windspeed)) daten <- c(daten, windspeed[,i])
#' daten[sample(1:length(daten), 260)] <- NA
#' daten <- matrix(daten, ncol = 6)
#'
#' ## Using gaussCItestMI
#' imp <- mice(daten)
#' tpc.fit2 <- tpc(suffStat = imp,
#'                indepTest = gaussCItestMI,
#'                 alpha = 0.01, labels = colnames(imp$data),
#'                 tiers = c(1,1,2,2,2,2))
#'
#' ## Using mixMItest
#' suffStat <- mice::complete(imp, action = "all")
#' tpc.fit3 <- tpc(suffStat = suffStat,
#'                indepTest = mixMItest,
#'                 alpha = 0.01, labels = colnames(imp$data),
#'                 tiers = c(1,1,2,2,2,2))
#'

tpc <- function (suffStat, indepTest, alpha, labels, p, forbEdges = NULL,
                 m.max = Inf,
                 skel.method = c("stable", "stable.parallel"),
                 conservative = FALSE,
                 maj.rule = TRUE, ### set this to TRUE
                 verbose = FALSE,
                 tiers = NULL, ### new argument
                 context.all = NULL, ### new argument
                 context.tier = NULL, ### new argument
                 numCores = NULL ### new argument
){

  ### forbEdges: A logical matrix of dimension p*p. If [i,j] is TRUE, then the
  ###            directed edge i->j is forbidden. If both [i,j] and [j,i] are
  ###            TRUE, then any type of edge between i and j is forbidden.
  ###
  ### tiers:  Numeric vector specifying the tier / time point for each variable.
  ###         Must be of length 'p', if specified, or have the same length as
  ###         'labels',if specified. A smaller number corresponds to an earlier
  ###         tier / time point. Conditional independence testing is restricted
  ###         such that if x is in tier t(x) and y is in t(y), only those
  ###         variables are allowed in the conditioning set whose tier is not
  ###         larger than max(t(x), t(y)).
  ###
  ### context.all: Numeric or character vector. Specifies the positions or names of
  ###         global context variables. Global context variables have no incoming edges, i.e. no
  ###         parents, and are themselves parents of all non-context variables in the
  ###         graph.
  ###
  ### context.tier: Numeric or character vector. Specifies the positions or names of
  ###         tier-specific context variables. Tier-specific context variables have
  ###         no incoming edges, i.e. no parents, and are themselves parents of all
  ###         non-context variables in the same tier.



  cl <- match.call()
  if (!missing(p)) {
    stopifnot(is.numeric(p), length(p <- as.integer(p)) ==
                1, p >= 2)
  }
  if (missing(labels)) {
    if (missing(p))
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else {
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    } else if (p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else message("No need to specify 'p', when 'labels' is given")
  }

  if (is.null(tiers)) {
    ## if no tiers are specified, everything is tier 1
    tiers <- rep(1, p)
  } else {
    ## check if 'tiers' are correctly specified
    if (!is.numeric(tiers)) {stop("'tiers' must be a numeric vector")}
    if (length(tiers) != p) {stop("length of 'tiers' does not match 'p' or length of 'labels'")}
  }

  if (!is.null(context.all)) {
    if (is.character(context.all)) {
      if (!all(context.all %in% labels)) {stop("'context.all' includes variable names not in 'labels'")}
      context.all <- which(labels %in% context.all)
    }

    if (is.numeric(context.all)) {
      if (!all(context.all %in% (1:p))) {stop("'context.all' contains elements that are smaller than 1 or larger than 'p'")}
      if (!all(tiers[context.all]==min(tiers))) {stop("'context.all' variables must be in the first tier")}
    } else {
      stop("'context.all' must be an integer vector or character vector")
    }
  }

  if (!is.null(context.tier)) {
    if (is.character(context.tier)) {
      if (!all(context.tier %in% labels)) {stop("'context.tier' includes variable names not in 'labels'")}
      context.tier <- which(labels %in% context.tier)
    }
    if (is.numeric(context.tier)) {
      if (!all(context.tier %in% 1:p)) {stop("'context.tier' contains elements that are smaller than 1 or larger than 'p'")}
    } else {
      stop("'context.tier' must be a numeric or character vector")
    }
  }

  if ( !is.null(context.tier) & !is.null(context.all) ) {
    if (length(intersect(context.tier, context.all)) > 0) {stop(paste("The following variables are in both 'context.tier' and 'context.all': ", paste(intersect(context.tier, context.all), collapse=",")))}
  }


  skel.method <- match.arg(skel.method)
  if (conservative && maj.rule) {
    stop("Choose either conservative PC or majority rule PC!")
  }
  if ((!conservative) && (!maj.rule)) {
    stop("Choose one of conservative PC and majority rule PC!")
  }

  ## generate fixedEdges and fixedGaps according to context.all and context.tier
  fixedEdges <- matrix(FALSE, p, p)
  fixedGaps <- matrix(FALSE, p, p)

  ## context.all
  if (!is.null(context.all)) {
    for (i in context.all) {
      fixedEdges[i, ] <- TRUE
      fixedEdges[ ,i] <- TRUE
      for (j in context.all) {
        fixedEdges[i,j] <- FALSE
        fixedEdges[j,i] <- FALSE
        fixedGaps[i,j] <- TRUE
        fixedGaps[j,i] <- TRUE
      }
    }
  }
  ## context.tier
  if (!is.null(context.tier)) {
    for (i in context.tier) {
      for (j in c(context.tier, context.all)) {
        fixedGaps[i,j] <- TRUE
        fixedGaps[j,i] <- TRUE
      }
      k <- tiers[i]
      fixedEdges[i,tiers==k] <- TRUE
      fixedEdges[tiers==k,i] <- TRUE
      fixedGaps[i,tiers!=k] <- TRUE
      fixedGaps[tiers!=k,i] <- TRUE
    }
  }



  if (!is.null(forbEdges)) {

    ## check if forbEdges contradicts context.tier or context.all
    checkMatrix <- fixedEdges
    for (i in context.all) {
      checkMatrix[ ,i] <- FALSE
    }
    for (i in context.tier) {
      k <- tiers[i]
      checkMatrix[tiers==k,i] <- FALSE
    }

    if ( sum(checkMatrix*forbEdges) > 0 ) {
      ConflictList <- which((checkMatrix*forbEdges) > 0, arr.ind=TRUE)
      ConflictList[ ,1] <- labels[ConflictList[ ,1]]
      ConflictList[ ,2] <- labels[as.numeric(ConflictList[ ,2])]
      colnames(ConflictList) <- NULL
      cat("Note: 'forbEdges' overrules 'context.tier' and 'context.all'.
          Edges between the following pairs of nodes are forbidden by
          'forbEdges' even though 'context.tier' and/or 'context.all'
          suggest they should be present:\n")
      print(ConflictList)
    }

    ## matrix of forbidden adjacencies (no type of edge allowed between a and b):
    forbAdj <- forbEdges * t(forbEdges)

    ## modify fixedEdges and fixedGaps according to forbEdges
    fixedGaps <- ( fixedGaps + forbAdj ) > 0
    fixedEdges <- ( fixedEdges - forbAdj ) > 0

    ## generate list of forbidden arrows (where the other direction is allowed)
    forbArrows <- forbEdges - forbEdges * t(forbEdges)
    forbArrowsL <- which(forbArrows > 0, arr.ind=TRUE)
    forbArrowList <- lapply(seq_len(nrow(forbArrowsL)), function(i) forbArrowsL[i,])
  } else {
    forbArrowList <- list()
  }


  if (skel.method=="stable.parallel") {
    if (is.null(numCores)) {stop("Please specify 'numCores'.")}
    skel <- tskeleton.parallel(suffStat, indepTest, alpha, labels = labels,
                               method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                               m.max = m.max, verbose = verbose, tiers = tiers,
                               numCores = numCores)
  } else {
    skel <- tskeleton(suffStat, indepTest, alpha, labels = labels,
                      method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                      m.max = m.max, verbose = verbose, tiers = tiers)
  }

  skel@call <- cl
  if (graph::numEdges(skel@graph) == 0) {
    return(skel)
  }
  ## step II, orientation of v-structures:
  skelII <- tpc.cons.intern(skel, suffStat, indepTest, alpha, version.unf = c(2, 1), maj.rule = maj.rule,
                            verbose = verbose, tiers=tiers, context.all=context.all, context.tier=context.tier,
                            forbEdges = forbEdges)

  ## step III, orientation of edges between tiers:
  gIII <- as(skelII$sk@graph, "matrix")
  for (t in unique(tiers)) {
    gIII[tiers>t, tiers==t] <- 0
  }
  ## context variables:
  for (i in context.all) {
    gIII[ ,i] <- 0
  }
  for (i in context.tier) {
    k <- tiers[i]
    gIII[tiers==k,i] <- 0
  }
  ## edges in forbArrowList:
  for ( r in forbArrowList ) {
    gIII[r[1],r[2]] <- 0
  }

  # step IV, Meek's rules
  skelIII <- skelII$sk
  skelIII@graph <- as(gIII, "graphNEL")
  MeekRules(skelIII, verbose = verbose, unfVect = skelII$unfTripl, solve.confl = TRUE)
}
