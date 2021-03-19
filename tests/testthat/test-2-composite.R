skip_if_not(capabilities()["long.double"])
context("composite works")
require(testthat)
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
dichotParamTab$subtest <- c(2,2,2,2,2,3,3,3,3,3,3,3,2)
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
           subtest=c(2,3),
           location=c(277.1563,278),
           scale=c(37.7297, 38),
           subtestWeight=c(0.3,0.7))
mat$x1 <- stuDat$x1 <- x1
stuDat$origwt <- mat$origwt <- runif(nrow(stuDat)) * 4 * abs(stuDat$x1 + 3)

# write "mat" here to compare to AM
# save sav in wide version
if(FALSE) {
  # do this to test against AM
  saveSav(stuItems, mat, "testComp", "~/", idVar="oppID")
  paramTab2AMdct(dichotParamTab=dichotParamTab,
                 polyParamTab=NULL,
                 testScale=testDat,
                 weightVar="origwt",
                 stratumVar="repgrp1", 
                 PSUVar="jkunit",
                 saveFilePath="~/",
                 testName="testComp")
}



# tests:
stuItemsSplit <- split(stuItems, stuItems$oppID)
mml1 <- mml(composite ~ 1, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, composite=TRUE)
expect_is(mml1, "mmlCompositeMeans")
mml1s <- summary(mml1, gradientHessian=TRUE)
expect_is(mml1s, "summary.mmlCompositeMeans")

context("composite with regressor")

# composite
mml1 <- mml(composite ~ x1, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, weightVar="origwt", strataVar="repgrp1", PSUVar="jkunit", minNode=-5, maxNode=5)
mml1s <- summary(mml1, gradientHessian=TRUE, varType="consistent")
mml1Robust <- summary(mml1, varType="robust")
mml1Taylor <- summary(mml1, varType="Taylor", gradientHessian=TRUE)
# regression tests, no external ref
expect_equal(mml1s$coef[,2],
             c(`(Intercept)` = 0.8821387089, x1 = 1.4734883920, `Population SD` = NA),
             tolerance=sqrt(.Machine$double.eps) * 2000)

expect_equal(mml1Robust$coef[,2],
             c(`(Intercept)` = 0.3405808015, x1 = 0.5490783818, `Population SD` = NA),
             tolerance=sqrt(.Machine$double.eps)*200)

# AM results
# Parameter Name  Estimate  Standard Error  t Statistic p > |t| 
# Constant     277.670       2.576      -0.030       0.976  
# X1       9.707       4.419       2.197       0.029  
# Root MSE      37.554  (--)  
expect_equal(mml1Taylor$coef[,1],
             c(`(Intercept)` = 277.670, x1 = 9.707, `Population SD` = 37.554),
             tolerance=0.001)

# AM and Dire do not use the same equations
expect_equal(mml1Taylor$coef[1:2,2],
             c(`(Intercept)` = 2.561, x1 = 4.393),
             tolerance=0.005)

# non-composite
mml2 <- mml(2~x1,stuItems=stuItemsSplit, stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, weightVar="origwt", composite=FALSE)
mml2T <- summary(mml2, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit")
expect_equal(mml2T$coef[,1], # check estimates
               c(`(Intercept)` = 275.774229447017, x1 = 13.0517930716158, `Population SD` = 37.9227262346112),
               tol=20*sqrt(.Machine$double.eps))
expect_equal(mml2T$coef[,2], # check estimates
               c(`(Intercept)` = 2.64675830633024, x1 = 4.56193953861948, `Population SD` = 1.97288870616125),
               tol=2000*sqrt(.Machine$double.eps))
# non-composite but named composite; treated as "overall"
testDat <- data.frame(test=c("composite", "composite", "composite"),
           subtest=c(2,3, NA),
           location=c(277.1563,278, 200),
           scale=c(37.7297, 38, 90),
           subtestWeight=c(0.3,0.7, NA))

mml3 <- mml(composite~x1,stuItems=stuItemsSplit, stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, weightVar="origwt", composite=FALSE, strataVar="repgrp1", PSUVar="jkunit", minNode=-5, maxNode=5)

mml3T <- summary(mml3, varType="Taylor", gradientHessian=TRUE)

# Iterations: 21 
# Log Likelihood: -87057.6
#
# Adjusted Wald Test
# F(1,200) = 6.6087
# p(F > f) = 0.0108748  
# Dependent Variable: composite, overall  
# Parameter Name  Estimate  Standard Error  t Statistic p > |t| 
# Constant  199.11432422  5.14372445  38.71014595 0.00000000  
# X1  23.65037709 9.19983795  2.57073844  0.01089412  
# Root MSE  88.22756305 2.54241308  --- --- 
# AM Statistical Software Beta Version 0.06.04. (c) The American Institutes for Research and Jon Cohen  


# AM results:
expect_equal(mml3T$coef[,1],
               c(`(Intercept)` = 199.11432422, x1 = 23.65037709, `Population SD` = 88.22756305),
               tolerance=1e-5)

expect_equal(mml3T$coef[,2], # check estimates
               c(`(Intercept)` = 5.14372445, x1 = 9.19983795, `Population SD` = 2.54241308),
               tolerance=1e-4)
