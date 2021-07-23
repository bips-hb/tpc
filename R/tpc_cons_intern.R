
tpc.cons.intern <- function(
  sk, suffStat, indepTest, alpha, version.unf = c(NA, NA),
  maj.rule = FALSE, forbEdges=NULL, tiers=NULL,
  context.all=NULL, context.tier=NULL, verbose = FALSE) {

  orientConflictCollider <- function(pdag, x, y, z) {
    if (pdag[x, y] == 1) {
      pdag[y, x] <- 0
    }
    else {
      pdag[x, y] <- pdag[y, x] <- 2
    }
    if (pdag[z, y] == 1) {
      pdag[y, z] <- 0
    }
    else {
      pdag[z, y] <- pdag[y, z] <- 2
    }
    pdag
  }

  g <- as(sk@graph, "matrix")
  stopifnot(all(g == t(g)))
  p <- as.numeric(dim(g)[1])
  if (is.null(tiers)) { tiers <- rep(1, p) }
  # initiate list of ambiguous triples
  unfTripl <- vers <- rep(NA, min(p * p, 1e+05))
  # the counter counts the number of ambiguous triples
  counter <- 0
  if (sum(g) > 0) {
    # indices of edges (each edge is listed twice)
    ind <- which(g == 1, arr.ind = TRUE)
    # context variables can never be colliders
    ind <- ind[!(ind[ ,2] %in% union(context.all, context.tier)), ,drop=FALSE]
    # context varialbes do not need to be considered as potential collider parents
    ind <- ind[!(ind[ ,1] %in% union(context.all, context.tier)), ,drop=FALSE]

    tripleMatrix <- matrix( , nrow=0, ncol=4)
    # go through all edges until tripleMatrix contains all unshielded triples
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1] # candidate 'left-hand' parent of collider
      b <- ind[i, 2] # candidate collider
      allC <- setdiff(which(g[b, ] == 1), a) # candidate 'right-hand' parents of collider
      newC <- allC[g[a, allC] == 0] # only those where candidate parents are not neighbours
      newC <- setdiff(newC, union(context.all, context.tier)) # exclude context variables
      tmpMatrix <- cbind(rep(a, length(newC)), rep(b, length(newC)), newC, rep(0, length(newC))) # added the zeros
      # add these triples to tripleMatrix
      tripleMatrix <- rbind(tripleMatrix, tmpMatrix)
      colnames(tripleMatrix) <- c("", "", "", "") # added a forth ""
    }
    if ((m <- nrow(tripleMatrix))>0) {
      # delete duplicates by only keeping those triples where candidate 'left-hand' parent has smaller node ID than candidate 'right-hand' parent
      deleteDupl <- logical(m)
      for (i in seq_len(m)) if (tripleMatrix[i, 1] > tripleMatrix[i, 3])
        deleteDupl[i] <- TRUE
      if (any(deleteDupl))
        tripleMatrix <- tripleMatrix[!deleteDupl, , drop = FALSE]

      # delete triples where forbEdges[a,b]==1 (a->b forbidden) or forbEdges[c,b]==TRUE (b<-c forbidden)
      if (!is.null(forbEdges)) {
        m2 <- nrow(tripleMatrix)
          deleteForb <- logical(m2)
          for (i in seq_len(m2)) {
            if ( forbEdges[tripleMatrix[i,1], tripleMatrix[i,2]] ) {
              deleteForb[i] <- TRUE
            }
            if ( forbEdges[tripleMatrix[i,3], tripleMatrix[i,2]] ) {
              deleteForb[i] <- TRUE
            }
          }
          if ( any(deleteForb) ) {
            tripleMatrix <- tripleMatrix[!deleteForb, , drop=FALSE]
          }
      }


      # now go through all triples in earnest
      # mark the ones ot be oriented with '1' in column 4
      # add the ambiguous ones to unfTripl
      for (i in seq_len(nrow(tripleMatrix))) {
        # add extra place holders if necessary
        if (counter + 1L == length(unfTripl)) {
          n.xtra <- min(p * p, 1e+05)
          new.len <- counter + 1L + n.xtra
          length(unfTripl) <- new.len
        }

        a <- tripleMatrix[i, 1]
        b <- tripleMatrix[i, 2] # candidate collider
        c <- tripleMatrix[i, 3]
        ############################# new code
        # go to next triple if the candidate collider is a context variable, or against the time structure
        if (b %in% union(context.all, context.tier)) {next}
        if (tiers[b]!=max(tiers[a],tiers[c])) {next}
        ############################# end new code

        # find all neighbours of candidate parents for the new independence tests
        nbrsA <- which(g[, a] != 0)
        nbrsC <- which(g[, c] != 0)

        ############################# new code
        # neighbours in the future of both candidate parents are excluded
        nbrsA <- intersect(nbrsA, which(tiers<=max(tiers[a],tiers[c])))
        nbrsC <- intersect(nbrsC, which(tiers<=max(tiers[a],tiers[c])))
        ############################# end new code
        if (verbose) {
          cat("\nTriple:", a, b, c, "and sepset by skelet:",
              unique(sk@sepset[[a]][[c]], sk@sepset[[c]][[a]]),
              "\n")
        }
        # checkTriple goes through all subsets of nbrsA and nbrsC
        r.abc <- tcheckTriple(a, b, c, nbrsA, nbrsC, sk@sepset[[a]][[c]],
                             sk@sepset[[c]][[a]], suffStat = suffStat, indepTest = indepTest,
                             alpha = alpha, version.unf = version.unf, maj.rule = maj.rule,
                             verbose = verbose)
        if (r.abc$decision == 3) { # decision 3 is ambiguous
          # add ambiguous triple to unfTripl
          counter <- counter + 1
          unfTripl[counter] <- triple2numb(p, a, b, c)
          vers[counter] <- r.abc$version
        }
        if ((version.unf[1] == 2) && (r.abc$version == 2) && (r.abc$decision != 3)) {
          # this triple is added to unfTripl because a and c are not independent given any subset of their neighbours
          counter <- counter + 1
          unfTripl[counter] <- triple2numb(p, a, b, c)
          vers[counter] <- r.abc$version
        }
        ############################# new code
        # mark this triple as 'to be oriented' if it meets the assumptions
        if ( (r.abc$decision == 1) &&
             ((version.unf[1] == 1) | (version.unf[1]==2) & r.abc$version==1) ){
          tripleMatrix[i,4] <- 1
        }
        ############################# end new code
      }
    }

    ############################# new code
    tripleMatrix <- tripleMatrix[tripleMatrix[ ,4]==1, , drop=FALSE]
    if (nrow(tripleMatrix) > 0) {
      pdag <- g
      for (i in 1:nrow(tripleMatrix)) {
        x <- tripleMatrix[i,1]
        y <- tripleMatrix[i,2]
        z <- tripleMatrix[i,3]
        pdag <- orientConflictCollider(pdag, x, y, z)
      }
      sk@graph <- as(pdag, "graphNEL")
    }
    ############################# end new code
  }
  length(unfTripl) <- counter
  return(list(unfTripl=unfTripl, sk=sk))
}
