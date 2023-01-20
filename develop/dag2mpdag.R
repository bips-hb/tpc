#################################################################################
### filename:           dag2mpdag.R                                           ###
###                                                                           ###
### project:            DFG Causal Discovery                                  ###
###                                                                           ###
### function:           Determines the maximally oriented partially directed  ###
###                     acyclic graph (MPDAG) given a DAG and temporal back-  ###
###                     ground information / a partial node ordering.         ###
###                                                                           ###
### author:             Maria Geers,  Janine Witte                            ###
###                                                                           ###
### R version:          4.0.1                                                 ###
#################################################################################

#' Determine the MPDAG of a DAG Given a Partial Node Ordering
#'
#' Takes as input a DAG and a partial node ordering and returns the maximally
#' oriented partially directed acyclic graph (MPDAG).
#'
#' @param dag A \code{graphNEL} or \code{bn} object.
#' @param tiers Numeric vector specifying the tier / time point for each variable.
#' The length must equal the number of nodes in cpdag.
#' @param context.tier Numeric vector specifying the positions of variables that
#' are parents of all other variables in their tier (except for the other
#' exotier and the exoall variables).
#' @param context.all Numeric vector specifying the positions of variables that
#' are parents of all other variables (except for the other exoall and the
#' exotier variables).
#'
#' @details The Markov equivalence class of a DAG can be uniquely represented by
#' a completed partially directed acyclic graph (CPDAG, see Andersen et al., 1997).
#' The subset of a DAG's Markov equivalence class that takes a known partial
#' node ordering into account can be uniquely represented by a maximally oriented
#' partially directed acyclic graph (MPDAG, see Perkovic et al., 2017).
#' This function takes as input a DAG and a partial node ordering, and outputs
#' the true MPDAG. The true MPDAG is the target graph of the tpc algorithm in
#' the sense that if tpc was given oracle independence information and a correct
#' partial node ordering, then the output would be the true MPDAG.
#'
#' @return A \code{graphAM} object.
#'
#' @references
#' E. Perkovic, M. Kalisch and M.H. Maathuis (2017). Interpreting and using CPDAGs
#' with background knowledge. \emph{Proceedings of the Thirty-Third Conference on
#' Uncertainty in Artificial Intelligence (UAI-17)}, page ID 120.
#'
#' S.A. Andersson, D. Madigan, and M.D. Perlman (1997). A characterization of
#' Markov equivalence classes for acyclic digraphs.
#' \emph{The Annals of Statistics} 25(2):505-541.
#'
#' @author Janine Witte, Maria Geers
# #' @export
#'
#' @examples
#' dag2mpdag(trueDAG,
#' tiers = c(1,1,2,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6),
#' exotier=c("age_t0", "age_t1", "age_t2"),
#' exoall=c("sex", "country"))

dag2mpdag <- function (dag, tiers=NULL, context.tier=NULL, context.all=NULL) {
### dag           A 'graphNEL' or 'bn' object.
### tiers         Numeric vector specifying the tier / time point for each variable.
###               The length must equal the number of nodes in cpdag.
### context.all   Integer or character vector. Specifies the positions or names of
###               context variables. Context variables have no incoming edges, i.e. no
###               parents, and are themselves parents of all non-context variables in the
###               graph.
### context.tier  Integer or character vector. Specifies the positions or names of
###               tier-specific context variables. Tier-specific context variables have
###               no incoming edges, i.e. no parents, and are themselves parents of all
###               non-context variables in the same tier.

  if (class(dag)=="bn") {
    dag <- as.graphNEL(dag)
  }

  # check that no node is in both context.all and context.tier
  if ( length(intersect(context.all, context.tier)) >0 ) {stop(paste("The following variables are in both 'context.tier' and 'context.all': ", paste(intersect(context.tier, context.all), collapse=",")))}

  # check if all nodes in context.tier and context.all are orphans
  M <- t(as(dag, "matrix"))
  npa <- apply(M[c(context.all, context.tier), , drop=FALSE], 1, sum)
  if (sum(npa)>0) {stop("Only nodes without parents in 'dag' can be specified as context variables")}

  if (is.null(tiers)) {
    tiers <- rep(1, length(dag@nodes))
  } else if ( length(tiers)!= length(dag@nodes) ) {
    stop("length of tiers must equal number of nodes in 'dag'")
  }

  # determine CPDAG
  cpdag <- dag2cpdag(dag)

  # determine undirected edges in cpdag
  W <- t(as(cpdag, "matrix"))
  Wu <- ( W + t(W) )==2
  B <- which(Wu, arr.ind = TRUE)

  # only keep those undirected edges where we know the orientation from the partial
  # ordering or from the information in context.all and context.tier
  dir1 <- apply(B, 1, function(i) { tiers[i[2]] > tiers[i[1]] })
  dir2 <- apply(B, 1, function(i) { (tiers[i[1]]==tiers[i[2]]) & (i[1] %in% context.tier) & !(i[2] %in% union(context.tier,context.all)) })
  dir3 <- apply(B, 1, function(i) { (tiers[i[1]]==tiers[i[2]]) & (i[1] %in% context.all) & !(i[2] %in% union(context.tier,context.all)) })
  A <- B[(dir1+dir2+dir3)>0, , drop=FALSE]

  # labels
  lbl <- colnames(W)
  xl <- lbl[A[ ,1]]
  yl <- lbl[A[ ,2]]

  V <- addBgKnowledge(gInput=W, xl, yl, verbose=TRUE)
  if (is.null(V)) {stop("The partial node ordering is not compatible with the DAG.")}
  getGraph(t(V))
}
