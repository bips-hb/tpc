#' Tiered PC Algorithm
#'
#' Constraint-based causal discovery using the PC algorithm while accounting for a
#' partial node ordering, e.g. a partial temporal ordering when the data were
#' collected in different waves of a cohort study.
#'
#' @docType package
#'
#' @keywords internal
#' @aliases tpc-package
"_PACKAGE"

utils::globalVariables(c("remEdges"))

## usethis namespace: start
#' @importFrom graph numEdges
#' @importFrom graphics legend
#' @importFrom igraph adjacent_vertices all_simple_paths graph_from_adjacency_matrix induced_subgraph
#' @importFrom methods as new
#' @importFrom parallel clusterEvalQ clusterExport makeCluster parLapply stopCluster
#' @import pcalg
#' @import Rfast
#' @importFrom utils combn
#  use_import_from("parallel", 'makeCluster')
## usethis namespace: end
NULL
