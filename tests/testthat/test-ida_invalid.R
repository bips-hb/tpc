test_that("IDA_invalid produces lists", {

  library(graph)
  library(pcalg)
  library(tpc)

  dag.amat <- matrix(rep(0,25), nrow = 5)

  while(isValidGraph(dag.amat, type = "cpdag")){
    possible_p <- c(11, 20, 30)
    possible_neighb_size <- c(3:10)

    p <-  sample(possible_p, 1)
    neigh <- sample(possible_neighb_size, 1)
    prob <- round(neigh / (p-1), 3)

    g <- pcalg:::randomDAG(p, prob)

    dag.amat <- t(as(g,"matrix"))
    dag.amat[dag.amat != 0] <- 1
  }


  knoten <- sample(p, 2, replace = FALSE)
  out.o <- suppressMessages(ida_invalid(x.pos = knoten[1], y.pos = knoten[2], graphEst = g,
                                        method = "optimal", verbose = FALSE, plot = FALSE))
  out.l <- suppressMessages(ida_invalid(x.pos = knoten[1], y.pos = knoten[2], graphEst = g,
                                        method = "local", verbose = FALSE, plot = FALSE))


  ## the output should be a list, even if no adjustment is required
  expect_true(all(is.list(out.o), is.list(out.l)))

  ## include 1 -> 2 -> 5 <- 3 <- 1: in this case, node 1 should be selected into the
  ## adjustment set in both the optimal and the local method
  g <- suppressWarnings(addEdge("1","2",g,1))
  g <- suppressWarnings(addEdge("1","3",g,1))
  g <- suppressWarnings(addEdge("2","5",g,1))
  g <- suppressWarnings(addEdge("3","5",g,1))
  out.o2 <- suppressMessages(ida_invalid(x.pos = 2, y.pos = 5, graphEst = g,
                                        method = "optimal", verbose = FALSE, plot = FALSE))
  out.l2 <- suppressMessages(ida_invalid(x.pos = 2, y.pos = 5, graphEst = g,
                                        method = "local", verbose = FALSE, plot = FALSE))

  expect_true(all(any(c(1,3) %in% unlist(out.o2)), 1 %in% unlist(out.l2)))

})
