test_that("Meek rules work", {

  library(MASS)
  library(tpc)
  library(pcalg)

  R <- cov(swiss)
  colnames(R) <- rownames(R) <- rep("",6)

  daten <- data.frame(MASS::mvrnorm(n = 100,
                                    mu = rep(0,6),
                                    Sigma = R))
  names(daten) <- paste0("X", 1:6)

  suffStat <- list(C = cor(daten), n = nrow(daten))
  alpha <- runif(1,0.01,0.1)

  sk1.fit <- pcalg::skeleton(suffStat = suffStat,
                             labels = names(daten),
                             indepTest = gaussCItest, alpha = alpha)

  tskel1.fit <- tskeleton(suffStat = suffStat,
                          labels = names(daten),
                          indepTest = gaussCItest, alpha = alpha)

  s1 <- MeekRules(sk1.fit)
  t1 <- MeekRules(tskel1.fit)

  expect_true(all.equal(s1@graph@edgeL, t1@graph@edgeL))
})
