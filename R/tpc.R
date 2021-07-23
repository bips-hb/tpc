

tpc <- function (suffStat, indepTest, alpha, labels, p,
                 forbEdges = NULL, m.max = Inf,
                 conservative = FALSE, maj.rule = TRUE,
                 tiers = NULL, context.all = NULL, context.tier = NULL,
                 verbose = FALSE){

  cl <- match.call()
  if (!missing(p)) {
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
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
    if (length(intersect(context.tier, context.all)) > 0) {
      stop(paste("The following variables are in both 'context.tier' and 'context.all': ",
                 paste(intersect(context.tier, context.all), collapse=",")))
    }
  }

  skel.method <- "stable"
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

    ## matrix of forbidden adjacencies (no type of edge allowed between a and
    ## b):
    forbAdj <- forbEdges * t(forbEdges)

    ## modify fixedEdges and fixedGaps according to forbEdges
    fixedGaps <- ( fixedGaps + forbAdj ) > 0
    fixedEdges <- ( fixedEdges - forbAdj ) > 0

    ## generate list of forbidden arrows (where the other direction is allowed)
    forbArrows <- forbEdges - forbEdges * t(forbEdges)
    forbArrowsL <- which(forbArrows > 0, arr.ind=TRUE)
    forbArrowList <- lapply(seq_len(nrow(forbArrowsL)),
                            function(i) forbArrowsL[i,])
  } else {
    forbArrowList <- list()
  }


  skel <- tskeleton(suffStat, indepTest, alpha, labels = labels,
                    method = skel.method, m.max = m.max,
                    fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                    tiers = tiers, verbose = verbose)

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
  MeekRules(skelIII, verbose = verbose, unfVect = skelII$unfTripl,
            solve.confl = TRUE)
}
