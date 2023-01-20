test_that("pc and tpc produce equal results", {

  suppressPackageStartupMessages(library(pcalg))
  library(tpc)

  possible_p <- c(11, 20, 30)
  possible_neighb_size <- c(3:10)

  p <-  sample(possible_p, 1)
  neigh <- sample(possible_neighb_size, 1)
  prob <- round(neigh / (p-1), 3)

  g <- pcalg:::randomDAG(p, prob)
  n <- 200
  daten <- pcalg::rmvDAG(n, g, errDist = "normal")

  tpc.fit <- tpc(suffStat = list(C = cor(daten), n = n),
                 indepTest = gaussCItest, alpha = 0.05, labels = g@nodes)
  pc.fit <- pcalg::pc(suffStat = list(C = cor(daten), n = n),
                      indepTest = gaussCItest, alpha = 0.05, labels = g@nodes,
                      maj.rule = TRUE, solve.conf = TRUE)

  expect_equal(tpc.fit@graph, pc.fit@graph)
})
