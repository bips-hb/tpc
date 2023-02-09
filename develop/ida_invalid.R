#' IDA Adjustment Sets with Invalid Graph
#'
#' Determines adjustment sets compatible with a partially directed graph that is not a
#' valid CPDAG or MPDAG.
#'
#' @param x.pos Position of exposure variable in `graphEst@nodes`.
#' @param y.pos Position of outcome variable in `graphEst@nodes`.
#' @param graphEst Estimated invalid partially directed graph. Usually obtained by
#' `pc.fit@graph`, where \code{pc.fit} is the output of [pcalg::pc()].
#' @param method Character string specifying the method.
#' \code{"local"}: returns all compatible parent sets of x.
#' \code{"optimal"}: returns all compatible optimal adjustment sets relative to \code{(x,y)}.
#' @param verbose Logical. If TRUE, details are printed to the console.
#' @param plot Logical. If TRUE, a plot is are produced.
#'
#' @details The IDA algorithm (Maathuis et al. 2019, implemented as
#' \code{pcalg::\link[pcalg]{ida}}) requires a valid CPDAG or MPDAG as an input.
#' \code{ida_invalid} is an alternative for the case that the graph returned by PC
#' or another causal discovery algorithm is invalid. It determines all local or
#' optimal adjustment sets compatible with the invalid graph by trying out all
#' possible orientations of undirected edges in relevant parts of the graph while
#' ignoring the rest of the graph.
#'
#' Adjustment for the optimal adjustment set is often more efficient, but determining the
#' optimal adjustment sets is more sensitive to faulty graph estimates than determining
#' the parent sets of x.
#'
#' Note that \code{ida_invalid} does not have any theoretical guarantees.
#'
#' @return A list of adjustment sets compatible with \code{graphEst}. \code{NULL} stands for
#' the empty set. \code{0} means that the effect is zero (this occurs if \code{y} is among the
#' parents of \code{x} in local adjustment, or if there is no path from \code{x} to \code{y}
#' in optimal adjustment).
#'
#' @references Maathuis MH, Kalisch M, & Buehlmann P (2009). Estimating high-dimensional
#' intervention effects from observational data. The Annals of Statistics, 37(6A), 3133-3164.
#'
#' @import igraph
#'
#' @export
#'
#' @examples
#' # load simulated cohort data
#' data(dat_sim)
#' n <- nrow(dat_sim)
#' lab <- colnames(dat_sim)
#'
#' # estimate skeleton with temporal ordering as background information
#' tiers <- rep(c(1,2,3), times = c(3,3,3))
#' tpc.fit <- tpc(suffStat = list(C = cor(dat_sim), n = n),
#'                indepTest = gaussCItest, alpha = 0.01, labels = lab,
#'                tiers = tiers)
#'
#' if(require("Rgraphviz", character.only = TRUE, quietly = TRUE) & packageVersion("igraph") < "1.4.0"){
#' # compare estimated CPDAGs
#'  data("true_sim")
#'  par(mfrow = c(1,1))
#'  plot(tpc.fit, main = "tPC estimate")
#'  }
#'
#' # check if the output is a valid CPDAG or MPDAG
#' isValidGraph(as(tpc.fit, "amat"), type="cpdag")
#' isValidGraph(as(tpc.fit, "amat"), type="pdag")
#'
#' # determine possibly valid local adjustment sets for the effect of A1 on B1
#' ida_invalid(x.pos = 1, y.pos = 2, graphEst = tpc.fit@graph, method = "local", verbose = TRUE)
#' pcalg::ida(x.pos = 1, y.pos = 2, mcov = cor(dat_sim), graphEst = tpc.fit@graph,
#'            method = "local", type = "pdag")
#'
#' # determine possibly valid optimal adjustment sets for the effect of A1 on A3
#' ida_invalid(x.pos = 1, y.pos = 7, graphEst = tpc.fit@graph, method = "optimal")


ida_invalid <- function(x.pos, y.pos, graphEst, method = NULL,
                       verbose = TRUE, plot = TRUE) {
  pc.obj <- prep.graph.ida(graphEst)

  if ( !identical(method,"local") & !identical(method,"optimal") ) {
    stop("'method' must be one of 'local', 'optimal'")
  }

  stopifnot(x.pos == as.integer(x.pos), y.pos == as.integer(y.pos),
            length(x.pos) == 1, length(y.pos) == 1, x.pos != y.pos)

   if (verbose) {
    x <- pc.obj$lab[x.pos]
    y <- pc.obj$lab[y.pos]
   }

  if ( method == "local" ) {
    ### possible parents of X
    possPa <-
      sort(unlist(adjacent_vertices(pc.obj$ig, x.pos, mode="in")))

    ### definite parents of X
    Pa <- sort(unlist(adjacent_vertices(pc.obj$dg, x.pos, mode="in")))

    ### ambiguous parents of X
    amPa <- setdiff(possPa, Pa)

    if (verbose){
      cat("The exposure x =", x, "has", length(Pa), "definite parents and",
          length(amPa), "ambiguous parents.\n")
      }
    if (length(Pa) > 0){
        Pa_string <- paste(pc.obj$lab[Pa], collapse=", ")
        Pa_string <- paste("(", Pa_string, ").\n")
        if(verbose) cat("The definite parents are", Pa_string)
    }
    if (length(amPa) > 0){
        amPa_string <- paste(pc.obj$lab[amPa], collapse=", ")
        amPa_string <- paste("(", amPa_string, ").\n")
        if(verbose) cat("The ambiguous parents are", amPa_string)
    }

    ### plot subgraph showing all possible parents
    if ( plot ) {
      sub <- induced_subgraph(pc.obj$ig, vids=c(possPa, x.pos))
      col <- rep("white", length(possPa))
      col[names(V(sub)) %in% pc.obj$lab[amPa]] <- "yellow"
      col[names(V(sub)) %in% pc.obj$lab[x.pos]] <- "green"
      col[names(V(sub)) %in% pc.obj$lab[y.pos]] <- "red"

      plot(sub, vertex.color=col)
       legend('topleft', legend=c("exposure","parents","ambiguous"),
                       pt.bg=c("green", "white", "yellow"), pch=21)
    }

    ### check for directed cycles among the definite parents
    if ( length(Pa) > 1 ) {
      am_dir_Pa <- pc.obj$am_dir[c(Pa,x.pos),c(Pa,x.pos)]
      ok_Pa <- noCycles_pcalg(am_dir_Pa)
      if (!ok_Pa) {stop("There is a directed cycle in the relevant part of the graph (x, definite parents). Remove the cycle and try again.")}
    }

    ### If the outcome is among the definite parents, stop here. The effect is
    ### then zero.

    if ( y.pos %in% Pa ) {
      if (verbose) {
        cat("Since the outcome = ", y, "is among the definite parents, the causal effect of",
            x, "on", y, "is zero.\n")
      }
      return( list(0) )
    }

    ### If there are no ambiguous (yellow) parents, stop here. The unambiguous
    ### local adjustment set is then 'Pa'.

    if(length(amPa) == 0){
      if ( length(Pa) == 0 ) { Pa <- NULL }
      if (verbose) {
        if ( is.null(Pa) ) {
          cat("The only local adjustment set compatible with graphEst is the empty set.\n")
        } else {
          cat("They form the only local adjustment set compatible with graphEst.\n")
        }
      }
      return( list(unname(Pa)) )
    }

    ### create list of undirected edges:
    un_list <- which(pc.obj$am_un_tri[c(possPa,x.pos), c(possPa,x.pos)]==1,
                     arr.ind=TRUE)
    un_list <- apply(un_list, 2, function(i){pc.obj$lab[c(possPa,x.pos)][i]})
    if ( is.vector(un_list) ) {
      un_list <- t(un_list)
    }

    ### go through the list and try all combinations of orientations
    decimals <- 1:2^nrow(un_list)
    m <- lapply(decimals,function(x){
      as.logical(intToBits(x))[1:nrow(un_list)]
    })

    adj_sets <- lapply(m, function(k){
      am2 <- pc.obj$am
      for (i in 1:nrow(un_list)) {
        am2[un_list[i,1], un_list[i,2]] <- k[i]
        am2[un_list[i,2], un_list[i,1]] <- 1 - k[i]
      }
      # 0 = edge from left to right in un_list
      # 1 = edge from right to left in un_list
      # check if this newly created adjacency matrix contains directed cycles
      am2_dir <- am2 - am2 * t(am2)
      ok <- noCycles_pcalg(am2_dir)
      if (!ok) {return(NA)}

      adj_set <- which(am2_dir[x, ]==1)
      if ( length(adj_set)==0 ) {adj_set <- NULL}
      return(unname(adj_set))
    })

    ### delete elements corresponding to graphs with directed cycles
    cycl <- sapply(adj_sets, anyNA)
    adj_sets[cycl] <- NULL

    ### remove duplicates
    adj_sets <- unique(adj_sets)

    ### remove sets including the outcome (and add 0 to the list of adj. sets)
    if ( y.pos %in% unlist(adj_sets) ) {
      del <- sapply(adj_sets, function(i) {y.pos %in% i})
      adj_sets[del] <- NULL
      adj_sets <- c(adj_sets, 0)
    }

    if ( verbose ) {
      if ( length(adj_sets)==1 ) {
        cat("There is 1 local adjustment set compatible with graphEst.\n")
      } else {
        cat("There are", length(adj_sets),
          "local adjustment sets compatible with graphEst.\n")
      }
    }

    return(adj_sets)

  } else if ( method == "optimal" ) {

    ### possible mediators
    possMedPaths <- all_simple_paths(pc.obj$ig, from=x.pos, to=y.pos,
                                     mode="out")
    possMedNodes <- sort(setdiff(unlist(possMedPaths), c(x.pos,y.pos)))

    ### definite mediators
    MedPaths <- all_simple_paths(pc.obj$dg, from=x, to=y, mode="out")
    MedNodes <- sort(setdiff(unlist(MedPaths), c(x.pos,y.pos)))

    ### ambiguous mediators
    amMedNodes <- setdiff(possMedNodes, MedNodes)

    ### ambiguous parents of possible mediators (undirected edges)
    possPaMed <- adjacent_vertices(pc.obj$ug, union(possMedNodes, y.pos),
                                   mode="in")
    possPaMed <- sort(setdiff(unlist(possPaMed), c(x.pos,y.pos,possMedNodes)))

    ### all relevant nodes
    rel <- sort(unique(c(possMedNodes, possPaMed, x.pos, y.pos)))

### 24.08.22
    if (verbose){
      cat("There are", length(MedNodes),
          "definite mediators and", length(amMedNodes),
          "ambiguous mediators of the effect of exposure", x, "on outcome",
          y, ".\n")}

    if (length(MedNodes) > 0){
        MedNodes_string <- paste(pc.obj$lab[MedNodes], collapse=", ")
        MedNodes_string <- paste("(", MedNodes_string, ").\n")
        if(verbose) cat("The definite mediators are", MedNodes_string)
    }
    if (length(amMedNodes) > 0){
        amMedNodes_string <- paste(pc.obj$lab[amMedNodes], collapse=", ")
        amMedNodes_string <- paste("(", amMedNodes_string, ").\n")
        if(verbose) cat("The ambiguous mediators are", amMedNodes_string)
    }
    if (length(possPaMed) > 0) {
        possPaMed_string <- paste(pc.obj$lab[possPaMed], collapse=", ")
        possPaMed_string <- paste("(", possPaMed_string, ").\n")
        if(verbose){
          cat("Further, the definite and ambiguous mediators have",
              length(possPaMed), "ambiguous parents:", possPaMed_string)
        }
    }

    if ( plot ) {
      ### definite parents of mediators
      PaMed <- adjacent_vertices(pc.obj$dg, union(MedNodes, y.pos), mode="in")
      PaMed <- sort(setdiff(unlist(PaMed), c(x.pos,possMedNodes)))

      rel2 <- c(rel, PaMed)

      ### subgraph showing all possible mediators and their (ambiguous) parents
      sub <- induced_subgraph(pc.obj$ig, vids=rel2)
      col <- rep("gray", length(rel2))
      col[names(V(sub)) %in% pc.obj$lab[amMedNodes]] <- "yellow"
      col[names(V(sub)) %in% pc.obj$lab[possPaMed]] <- "yellow"
      col[names(V(sub)) %in% pc.obj$lab[PaMed]] <- "white"
      col[names(V(sub)) %in% pc.obj$lab[c(x.pos,y.pos)]] <- "green"

      plot(sub, vertex.color=col)
          legend('topleft',
                legend=c("exp./outc.", "mediators", "parents med./outc.", "ambiguous"),
                 pt.bg=c("green", "gray", "white", "yellow"), pch=21)
    }

    ### check if there are cycles among the definite mediators
    am_dir_Med <- pc.obj$am_dir[c(MedNodes,x.pos,y.pos),c(MedNodes,x.pos,y.pos)]
    ok_Med <- noCycles_pcalg(am_dir_Med)
    if (!ok_Med) {stop("There is a directed cycle in the relevant part of the graph (x, y, definite mediators). Remove the cycle and try again.")}

    ### create list of undirected edges:
    un_list <- which(pc.obj$am_un_tri[rel, rel]==1, arr.ind=TRUE)

    if ( nrow(un_list)==0 ) {
      allPa <- adjacent_vertices(pc.obj$ig, union(MedNodes,y.pos), mode="in")
      allPa <- sort(unique(unlist(allPa)))
      if (x.pos %in% allPa) {
        O_set <- setdiff(allPa, c(MedNodes, x.pos))
      } else {
        O_set <- 0
      }

      if (verbose) {
        cat("There is 1 optimal adjustment set compatible with graphEst.\n")
      }

      return( list(unname(O_set)) )
    }

    un_list <- apply(un_list, 2, function(i){pc.obj$lab[rel][i]})
    if ( is.vector(un_list) ) {
      un_list <- t(un_list)
    }

    ### go through the list and try all combinations of orientations
    decimals <- 1:2^nrow(un_list)
    m <- lapply(decimals,function(i){
      as.logical(intToBits(i))[1:nrow(un_list)]
    })

    adj_sets <- lapply(m, function(k){
      am <- pc.obj$am
      for (i in 1:nrow(un_list)) {
        am[un_list[i,1], un_list[i,2]] <- k[i]
        am[un_list[i,2], un_list[i,1]] <- 1 - k[i]
      }
      # 0 = edge from left to right in un_list
      # 1 = edge from right to left in un_list
      # check if this newly created adjacency matrix contains directed cycles
      am_dir <- am - am * t(am)
      ok <- noCycles_pcalg(am_dir)
      if (!ok) {return(NA)}

      ig2 <- graph_from_adjacency_matrix(t(am))
      # mediators = causal nodes (forbidden) including x and y
      MedPaths2 <- all_simple_paths(ig2, from=x.pos, to=y.pos, mode="out")
      MedNodes2 <- unique(unlist(MedPaths2))
      if ( length(MedNodes2)==0 ) {return(0)}
      # remove x
      MedNodes2 <- sort(setdiff(MedNodes2, x.pos))
      # parents of mediators
      Pa <- sort(
        setdiff(unique(unlist(adjacent_vertices(ig2, MedNodes2, mode="in"))),
                c(MedNodes2,x.pos)))
      if ( length(Pa)==0 ) {Pa <- NULL}
      return(unname(Pa))
    })

    ### delete elements corresponding to graphs with directed cycles
    cycl <- sapply(adj_sets, anyNA)
    adj_sets[cycl] <- NULL

    ### remove duplicates
    adj_sets <- unique(adj_sets)

    if (verbose){
      if (length(adj_sets) == 1) {
        cat("There is 1 optimal adjustment set compatible with graphEst.\n")
      } else {
        cat("There are", length(adj_sets),
          "optimal adjustment sets compatible with graphEst.\n")
      }
    }
    return(unname(adj_sets))
  }
}



prep.graph.ida <- function(graphEst){
  ## graphEst   Estimated invalid graph. Usually obtained by pc.fit@graph, where
  ##            pc.fit is the output of pcalg::pc.

  ### check if input is a graphNEL object
  if (!inherits(graphEst, "graphNEL")) {
    stop("The input is not a graphNEL object")
  }

  ### adjacency matrix
  am <- wgtMatrix(graphEst)
  am[am==2] <- 1
  am[am > 0 & am < 1] <- 1

  ### check if input is a valid pdag
  if ( isValidGraph(am, type="pdag") ) {
    message("The input is a valid pdag! Consider using pcalg::ida instead\n")
  }

  ### convert into an igraph object
  ig <- graph_from_adjacency_matrix(t(am))

  ### directed subgraph
  am_dir <- am - am * t(am)
  dg <- graph_from_adjacency_matrix(t(am_dir))

  ### undirected subgraph
  am_un <- am * t(am)
  ug <- graph_from_adjacency_matrix(t(am_un))

  ### triangular adjacency matrix of undirected subgraph
  am_un_tri <- am_un
  for (i in 1:nrow(am_un_tri)) {
    am_un_tri[i, 1:i] <- 0
  }

  output <- list(am = am, am_dir = am_dir, dg = dg, ig = ig, ug = ug,
                 am_un_tri = am_un_tri, lab = graphEst@nodes)
  return(output)
}






# this is a not exported function from the pcalg package, version 2.7-6.
noCycles_pcalg <- function (mat)
{
  ok <- TRUE
  for (i in 1:length(mat[1, ])) {
    pa.i <- which(mat[i, ] != 0 & mat[, i] == 0)
    if (length(pa.i) != 0) {
      for (j in 1:length(pa.i)) {
        pos.anc <- setdiff(possAn(m = mat, x = pa.i[j],
                                  ds = FALSE, possible = TRUE), pa.i[j])
        if (i %in% pos.anc)
          ok <- FALSE
      }
    }
  }
  ok
}
