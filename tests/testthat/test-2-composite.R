context("composite")
require(testthat)
# generate data
set.seed(142857)
n <- 2000
theta <- rnorm(n)
x1 <- runif(n)
theta <- theta + x1 * 0.2
# function to generate data
gen <- function(paramTab, theta, key) {
  # get parameter values from items
  a <- paramTab$slope
  b <- paramTab$difficulty
  c <- paramTab$guessing
  t(sapply(theta, function(thetai) {
    rbinom(length(a),1, c + (1 - c) / (1 + exp(-1.7 * a * (thetai - b))))
  }))
}
# paramTab from NAEP
dichotParamTab <- structure(list(ItemID = structure(1:13,
                                              .Label = c("m017401", "m017701", "m017901", "m018201", "m018401", "m018501", "m018601", "m020001", "m020501", "m046301", "m046501", "m051501", "n202831"), class = "factor"),
                      slope = c(0.25, 1, 1.15, 0.52, 1.11, 1.64, 0.78, 0.72, 0.72, 0.89, 0.92, 1.2, 0.75),
                      difficulty = c(-5.16, -1.01, -0.93, -1.21, -1.03, 0.34, 0.9, -0.49, -0.62, -1.07, -0.23, 1.22, -2.58),
                      guessing = c(0.19, 0.16, 0.15, 0.03, 0.24, 0.26, 0.12, 0, 0, 0.28, 0.33, 0.2, 0.25),
                      D=rep(1.7,13),
                      ScorePoints = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
                      MODEL = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "3pl", class = "factor")), class = "data.frame", row.names = c(NA, -13L))
dichotParamTab$test <- "composite"
dichotParamTab$subtest <- c("subA","subA","subA","subA","subA","subB","subB","subB","subB","subB","subB","subB","subA")
dichotParamTab$omitted <- rep(8,13)

# generate and then format the data
mat <- gen(dichotParamTab, theta, key=dichotParamTab$ItemID)
colnames(mat) <- c("m017401", "m017701", "m017901", "m018201", "m018401", "m018501", 
"m018601", "m020001", "m020501", "m046301", "m046501", "m051501", 
"n202831")
mat <- data.frame(mat)
rownames(mat) <- paste0("pseudo-student",1:nrow(mat))
mat$origwt <- 1
nperstratum <- 10
nstrata <- length(theta)/nperstratum
mat$repgrp1 <- rep(1:nstrata, each=nperstratum)
mat$jkunit <- rep(rep(1:2, each=nperstratum/2), nstrata)
stuItems <- mat[,1:13]
stuItems$oppID <- factor(rownames(mat), levels=rownames(mat))
mat$oppID <- factor(rownames(mat), levels=rownames(mat))
stuItems <- reshape(data=stuItems, varying=c(dichotParamTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=dichotParamTab$ItemID, timevar="key")
rownames(stuItems) <- NULL
stuDat <- mat[, c('origwt', 'repgrp1', 'jkunit')]
stuDat$oppID <- rownames(stuDat)
############### test functions ###############
testDat <- data.frame(test=c("composite", "composite"),
           subtest=c("subA","subB"),
           location=c(277.1563,278),
           scale=c(37.7297, 38),
           subtestWeight=c(0.3,0.7))
mat$x1 <- stuDat$x1 <- x1
stuDat$origwt <- mat$origwt <- runif(nrow(stuDat)) * 4 * abs(stuDat$x1 + 3)

# tests:

test_that("composite", {
  suppressWarnings(mml1 <- mml(composite ~ 1, stuItems=stuItems, stuDat=stuDat, 
                               dichotParamTab=dichotParamTab, Q=34, idVar="oppID", 
                               testScale=testDat, composite=TRUE, 
                               strataVar="repgrp1", PSUVar="jkunit"))
  expect_is(mml1, "mmlCompositeMeans")
  mml1s <- summary(mml1, nint=1000)
  expect_is(mml1s, "summary.mmlCompositeMeans")
})

context("composite with regressor")
test_that("composite with regressors", {
  # non-composite
  suppressWarnings(mml2 <- mml(subA ~ x1,stuItems=stuItems, stuDat=stuDat, 
                               dichotParamTab=dichotParamTab, Q=34, idVar="oppID", 
                               testScale=testDat, weightVar="origwt", composite=FALSE))
  # partial Taylor agrees with AM
  mml2T <- summary(mml2, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit")
  expect_equal(mml2T$coef[,1], # check estimates
                 c(`(Intercept)` = 275.8414, x1 = 13.0923, `Population SD` =  38.0835),
                 tol=2e2*sqrt(.Machine$double.eps))
  expect_equal(mml2T$coef[,2], # check estimates
                 c(`(Intercept)` = 2.817651, x1 = 4.907722, `Population SD` = 2.020495),
                 tol=2e2*sqrt(.Machine$double.eps))
  # non-composite but named composite; treated as "overall"
  testDat <- data.frame(test=c("composite", "composite", "composite"),
             subtest=c("subA","subB", NA),
             location=c(277.1563,278, 278),
             scale=c(37.7297, 38, 38),
             subtestWeight=c(0.3,0.7, NA))

  suppressWarnings(mml3 <- mml(composite~x1,stuItems=stuItems, stuDat=stuDat, 
                               dichotParamTab=dichotParamTab, Q=34, idVar="oppID", 
                               testScale=testDat, weightVar="origwt", composite=FALSE, 
                               strataVar="repgrp1", PSUVar="jkunit", minNode=-5, maxNode=5))
  mml3T <- summary(mml3, varType="Taylor")

  suppressWarnings(mml3c <- mml(composite~x1,stuItems=stuItems, stuDat=stuDat, 
                                dichotParamTab=dichotParamTab, Q=34, idVar="oppID", 
                                testScale=testDat, weightVar="origwt", composite=TRUE, 
                                strataVar="repgrp1", PSUVar="jkunit", minNode=-5, maxNode=5))
  mml3cT <- summary(mml3c, varType="Taylor")

  # AM results:
  expect_equal(mml3T$coef[,1],
                 c(`(Intercept)` = 277.626052117303, x1 = 9.98571505425182, `Population SD` = 37.2516385710825),
                 tolerance=1e-5)

  expect_equal(mml3T$coef[,2], # check estimates
                 c(`(Intercept)` = 2.17194822389352, x1 = 3.88465057751056, `Population SD` = 1.07353915893749),
                 tolerance=1e-4)

  expect_equal(mml3T$coef[,1], mml3cT$coef[,1], tol=0.02)

})

# PV generation
context("Composite PV generation")
test_that("Composite PV generation", {
  set.seed(2)
  stuDat$x2 <- factor(sample(1:3, nrow(stuDat), replace=TRUE), 1:3, LETTERS[1:3])
  suppressWarnings(mml2 <- mml(composite ~ x2, stuItems=stuItems, stuDat=stuDat, 
                               dichotParamTab=dichotParamTab, Q=34, idVar="oppID", 
                               testScale=testDat, composite=TRUE, strataVar="repgrp1",  PSUVar="jkunit"))
  expect_is(mml2, "mmlCompositeMeans")
  mml2tt <- summary(mml2, varType="Taylor") # new Taylor
  expect_is(mml2tt, "summary.mmlCompositeMeans")

  pvsA <- drawPVs(mml2tt, npv=20L, newStuDat=stuDat, newStuItems=stuItems, verbose=FALSE)
  expect_is(pvsA, "data.frame")
  # make sure each subtest + composite has the correct number of columns/PVs
  expect_equal(sum(grepl("subA_", colnames(pvsA))), 20)
  expect_equal(sum(grepl("subB_", colnames(pvsA))), 20)
  expect_equal(sum(grepl("composite_", colnames(pvsA))), 20)

})
