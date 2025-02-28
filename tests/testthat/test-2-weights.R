# calculating these accurately requires long doubles
require(testthat)
context("simple univariate") 
### generate data
set.seed(142857)
n <- 2000
theta <- rnorm(n)
x1 <- runif(n)
theta <- theta + x1 * 0.2
# function to generate data
gen <- function(dParamsTab, theta, key) {
  # get parameter values from items
  a <- dParamsTab$slope
  b <- dParamsTab$difficulty
  c <- dParamsTab$guessing
  D <- dParamsTab$D
  t(sapply(theta, function(thetai) {
    rbinom(length(a),1, c + (1 - c) / (1 + exp(-D * a * (thetai - b))))
  }))
}
# paramTab from NAEP
dichotParamTab <- structure(list(ItemID = structure(1:13, .Label = c("m017401", "m017701", "m017901", "m018201", "m018401", "m018501", "m018601", "m020001", "m020501", "m046301", "m046501", "m051501", "n202831"), class = "factor"),
                                 test = c("composite", "composite", "composite", "composite", "composite", "composite", "composite", "composite", "composite", "composite", "composite", "composite", "composite"),
                                 subtest = c("num", "num", "num", "num", "num", "num", "num", "num", "alg", "alg", "alg", "alg", "alg"),
                                 slope = c(0.25, 1, 1.15, 0.52, 1.11, 1.64, 0.78, 0.72, 0.72, 0.89, 0.92, 1.2, 0.75),
                                 difficulty = c(-5.16, -1.01, -0.93, -1.21, -1.03, 0.34, 0.9, -0.49, -0.62, -1.07, -0.23, 1.22, -2.58),
                                 guessing = c(0.19, 0.16, 0.15, 0.03, 0.24, 0.26, 0.12, 0, 0, 0.28, 0.33, 0.2, 0.25),
                                 D = rep(1.7,13),
                                 missingCode = rep(8,13),
                                 missingValue = rep("c",13)),
                                 row.names = c(NA, -13L), class = "data.frame")

# generate and then format the data
mat <- gen(dichotParamTab, theta, key=dichotParamTab$ItemID)
colnames(mat) <- c("m017401", "m017701", "m017901", "m018201", "m018401",
                   "m018501", "m018601", "m020001", "m020501", "m046301",
                   "m046501", "m051501", "n202831")
testDat <- data.frame(location=c(277.1563),
                      scale=c(37.7297))

testDat <- data.frame(test=c("composite", "composite","composite") ,
                      subtest=c("num", "alg",NA),
                      location=c(277.1563, 280.2948,277.1563),
                      scale=c(37.7297, 36.3887, 37.7297),
                      subtestWeight=c(0.3,0.7,1))

mat <- data.frame(mat)
rownames(mat) <- paste0("pseudo-student",1:nrow(mat))
mat$origwt <- 1
nperstratum <- 10
nstrata <- length(theta)/nperstratum
mat$repgrp1 <- rep(1:nstrata, each=nperstratum)
mat$jkunit <- rep(rep(1:2, each=nperstratum/2), nstrata)
stuItems <- mat[,1:13]
stuItems$oppID <- factor(rownames(mat), levels=rownames(mat))
stuItems <- reshape(data=stuItems, varying=c(dichotParamTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=dichotParamTab$ItemID, timevar="key")
rownames(stuItems) <- NULL
stuDat0 <- mat
stuDat0$oppID <- rownames(stuDat0)
############### test functions ###############
mat$x1 <- stuDat0$x1 <- x1
stuDat0$origwt <- mat$origwt <- 1

stuDatW <- stuDat0
stuDat2 <- stuDat0
stuDat3 <- stuDat0

stuDat2$oppID <- paste0(stuDat2$oppID, "r1")
stuDat3$oppID <- paste0(stuDat3$oppID, "r3")

stuDat0 <- rbind(stuDat0, stuDat2[stuDat2$x1 < 0.3, ])
stuDat0 <- rbind(stuDat0, stuDat3[stuDat3$x1 < 0.5, ])

stuDatW$origwt <- ifelse(stuDatW$x1 < 0.3, 3,
	                ifelse(stuDatW$x1 < 0.5, 2, 1))

# tests:
test_that("weighted vs unwweighted", {
  stuDatU <- stuDat0
  suppressWarnings(mml1qnuw <- mml(composite ~ x1, stuDat=stuDatU, dichotParamTab=dichotParamTab,
                  Q=34, minNode=-4, maxNode=4, idVar="oppID", testScale=testDat,
                  composite=FALSE, strataVar="repgrp1", PSUVar="jkunit", weightVar="origwt",
                  optimizer="QN"))
  suppressWarnings(mml1qnw <- mml(composite ~ x1, stuDat=stuDatW, dichotParamTab=dichotParamTab,
                 Q=34, minNode=-4, maxNode=4, idVar="oppID", testScale=testDat,
                 composite=FALSE, strataVar="repgrp1", PSUVar="jkunit", weightVar="origwt",
                 optimizer="QN"))
  expect_equal(coef(mml1qnuw), coef(mml1qnw), tolerance=1e-4)

  suppressWarnings(mml1emuw <- mml(composite ~ x1, stuDat=stuDatU, dichotParamTab=dichotParamTab,
                  Q=34, minNode=-4, maxNode=4, idVar="oppID", testScale=testDat,
                  composite=FALSE, strataVar="repgrp1", PSUVar="jkunit", weightVar="origwt",
                  optimizer="EM"))
  suppressWarnings(mml1emw <- mml(composite ~ x1, stuDat=stuDatW, dichotParamTab=dichotParamTab,
                 Q=34, minNode=-4, maxNode=4, idVar="oppID", testScale=testDat,
                 composite=FALSE, strataVar="repgrp1", PSUVar="jkunit", weightVar="origwt",
                 optimizer="EM"))
  expect_equal(coef(mml1emuw), coef(mml1emw), tolerance=1e-4)
  expect_equal(coef(mml1emw), coef(mml1qnw), tolerance=1e-4)
})
