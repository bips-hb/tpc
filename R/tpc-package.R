#' Temporal PC Algorithm
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
#' @import pcalg
#' @importFrom graph numEdges
#' @importFrom methods as
#' @importFrom methods new
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom utils combn
#  use_import_from("parallel", 'makeCluster')
## usethis namespace: end
NULL
