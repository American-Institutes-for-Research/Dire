# calculating these accurately requires long doubles
skip_if_not(capabilities()["long.double"])
require(testthat)
context("simple univariate") 
### generate data
set.seed(142857)
n <- 2000
theta <- rnorm(n)
x1 <- runif(n)
theta <- theta + x1 * 0.2
mean(theta)
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
stuDat <- mat[, c('origwt', 'repgrp1', 'jkunit')]
stuDat$oppID <- rownames(stuDat)
############### test functions ###############
mat$x1 <- stuDat$x1 <- x1
stuDat$origwt <- mat$origwt <- runif(nrow(stuDat)) * 4 * abs(stuDat$x1 + 3)

# tests:
test_that("simple univariate", {
  mml1 <<- mml(composite ~ 1, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, composite=FALSE, strataVar="repgrp1", PSUVar="jkunit")
  expect_is(mml1, "mmlMeans")
  mml1s <<- summary(mml1, gradientHessian=TRUE)
  mml1Robust <<- summary(mml1, varType="robust")
  mml1Cluster <<- summary(mml1, varType="cluster", clusterVar="repgrp1")
  mml1Taylor <<- summary(mml1, varType="Taylor")
})


context("mml errors") 
test_that("mml errors", {
  expect_error(mml1 <- mml(composite ~ 1, stuItems=stuItems[1:100,], stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat), "pseudo-student101") 
  expect_error(mml1 <- mml(composite ~ 1, stuItems=stuItems, stuDat=stuDat[-11,], dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat), "pseudo-student11")
  stuItems2 <- stuItems
  stuItems2$score[stuItems2$key == "m046501"] <- sample(0:2, sum(stuItems2$key == "m046501"), replace=TRUE)
  expect_error(mml1 <- mml(composite ~ 1, stuItems=stuItems2, stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat), "score points inconsistent")
})

context("mml summary and variance estimation") 
test_that("mml summary and variance estimation", {
  expect_equal(mml1$coefficients, c(`(Intercept)` = 0.116721240655054, "Population SD" = 0.972350274042811))
  # gradientHessian is a quirk of AM, set it to TRUE to compare the apples to apples
  expect_is(mml1s, "summary.mmlMeans")
  # values from AM with convergence set to 1E-12
  mml1s_coef_REF <- structure(c(281.560163073, 36.6864807567,
                                0.974354739545, 0.93732233951),
                              .Dim = c(2L,2L), .Dimnames = list(c("(Intercept)", "Population SD"),
                                                                c("Estimate", "StdErr")))
  # compare means
  expect_equal(mml1s$coefficients[,1], mml1s_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*20)
  # compare var estimates
  expect_equal(mml1s$coefficients[,2], mml1s_coef_REF[,2], tolerance=sqrt(.Machine$double.eps)*20)

  mml1R_SE_REF <- c(`(Intercept)` = 0.974716304084688, `Population SD` = 0.913234113680424)
  expect_equal(mml1Robust$coef[,"StdErr"], mml1R_SE_REF, tolerance=sqrt(.Machine$double.eps)*20)

  mml1SPrintRef <- c(
    "Call:",
    "mml(formula = composite ~ 1, stuItems = stuItems, stuDat = stuDat, ",
    "    idVar = \"oppID\", dichotParamTab = dichotParamTab, testScale = testDat, ",
    "    Q = 34, composite = FALSE, strataVar = \"repgrp1\", PSUVar = \"jkunit\")",
    "Summary Call:",
    "summary.mmlMeans(object = mml1, varType = \"robust\")",
    "",
    "Summary:",
    "            Estimate  StdErr t.value",
    "(Intercept)  281.560   0.975     289",
    "",
    "Residual Variance Estimate:",
    "              Estimate StdErr", 
    "Population SD       37   0.91",
    "",
    "Convergence = converged", 
    "LogLike = -12438.48", 
    "Observations = 2000", 
    "Weighted observations = 2000", 
    "location = 277.1563 scale = 37.7297")

  withr::with_options(list(digits=2), # accuracy is not the point here, the output is.
                      mml1SPrint <- capture.output(mml1Robust))
  mml1SPrint <- mml1SPrint[!grepl("^Iterations", mml1SPrint)]
  expect_equal(mml1SPrint, mml1SPrintRef)

  # REF values from AM
  mml1T_SE_REF <- c(`(Intercept)` = 0.967524019094, `Population SD` = 0.873102164238)
  expect_equal(mml1Taylor$coef[,"StdErr"], mml1T_SE_REF, tolerance=0.02)

  mml1C_SE_REF <- c(`(Intercept)` = 0.958703029704291, `Population SD` = 0.881032943695999)
  expect_equal(mml1Cluster$coef[,"StdErr"], mml1C_SE_REF, tolerance=sqrt(.Machine$double.eps)*20)

  mml2 <<- mml(composite ~x1, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, composite=FALSE)
  expect_is(mml2, "mmlMeans")
  expect_equal(mml2$coefficients, c(`(Intercept)` = 0.00592827512675005, x1 = 0.224670345051881, "Population SD" = 0.970126544720125), tolerance=sqrt(.Machine$double.eps)*200)
  mml2s <<- summary(mml2, gradientHessian=TRUE)
  # results form AM
  mml2s_coef_REF <- structure(c(277.379973218, 8.47675551939, 36.602582391, 
                                1.91227846027, 3.42771590088, 0.938701562756),
                              .Dim = c(3L, 2L), .Dimnames = list(c("(Intercept)", "x1", "Population SD"),
                                                                 c("Estimate", "StdErr")))
  expect_equal(mml2s$coefficients[,1], mml2s_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*10)
  expect_equal(mml2s$coefficients[,2], mml2s_coef_REF[,2], tolerance=sqrt(.Machine$double.eps)*10)
  mml2R_SE_REF <- c(`(Intercept)` =  1.93496054639559, x1 = 3.31985602943675, `Population SD` = 0.910098075346815)
  mml2Robust <<- summary(mml2, varType="robust")
  expect_equal(mml2Robust$coef[,"StdErr"], mml2R_SE_REF, tolerance=(.Machine$double.eps)^0.25)

  # gradientHessian=TRUE is what AM uses
  mml2Taylor <<- summary(mml2, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit", gradientHessian=TRUE)
  mml2T_SE_REF <- c(`(Intercept)` = 1.91342587986, x1 = 3.35043639323, `Population SD` = 0.872891128653)
  expect_equal(mml2Taylor$coef[,"StdErr"], mml2T_SE_REF, tolerance=0.001)
  mml2Cluster <<- summary(mml2, varType="cluster", clusterVar="repgrp1")
  mml2C_SE_REF <- c(`(Intercept)` = 1.88505689357269, x1 = 3.33275071678814, `Population SD` = 0.87139149522642)
  expect_equal(mml2Cluster$coef[,"StdErr"], mml2C_SE_REF, tolerance=(.Machine$double.eps)^0.25)
})

context("weighted case")
test_that("weighted case", {
  mml2W <- mml(composite~x1,stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", weightVar="origwt", testScale=testDat, composite=FALSE)
  mml2s_coef_REF <- structure(c(276.78360779, 9.90868559817, 36.9630033481, 
                                2.1549072924, 3.85424777114, 1.06630035114),
                              .Dim = c(3L, 2L), .Dimnames = list(c("(Intercept)", "x1", "Population SD"),
                                                                 c("Estimate", "StdErr")))
  mml2WTaylor <- summary(mml2W, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit", gradientHessian=TRUE)
  expect_equal(mml2WTaylor$coefficients[,1], mml2s_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*20)
  # compare var estimates
  expect_equal(mml2WTaylor$coefficients[,2], mml2s_coef_REF[,2], tolerance=2*(.Machine$double.eps)^0.25)
})

context("factor")
test_that("factor", {
  set.seed(142857)
  b <- table(stuItems$oppID, stuItems$score)
  stuDat2 <- stuDat
  stuDat2$multiLevel <- factor(sample(c(1:5,NA), nrow(stuDat2), prob=c(0.5,0.3,0.10,0.03,0.02,0.05), replace=TRUE), 1:5, LETTERS[1:5])
  stuDat2$mlM <- mat$mlM <- 0 + stuDat2$multiLevel %in% NA
  stuDat2$mlE <- mat$mlE <- 0 + (stuDat2$multiLevel %in% "E" | stuDat2$oppID %in% rownames(b[b[,2] == 13,])[1:80])
  stuDat2$mlA <- mat$mlA <- 0 + (stuDat2$multiLevel %in% "A" & !stuDat2$mlE)
  stuDat2$mlB <- mat$mlB <- 0 + (stuDat2$multiLevel %in% "B" & !stuDat2$mlE)
  stuDat2$mlC <- mat$mlC <- 0 + (stuDat2$multiLevel %in% "C" & !stuDat2$mlE)
  stuDat2$mlD <- mat$mlD <- 0 + (stuDat2$multiLevel %in% "D" & !stuDat2$mlE)
  stuDat2$mlM <- mat$mlM <- 0 + (stuDat2$multiLevel %in% NA & !stuDat2$mlE)
  stuItemsFactor <- stuItems
  naitems <- c("m017401", "m017701", "m017901", "m018201", "m018401", "m018501", "m018601", "m020001", "m020501", "m046301", "m046501")
  stuItemsFactor$score[stuItems$oppID %in% rownames(b[b[,2] == 13,])[1:20] & stuItems$key %in% naitems] <- NA
  matFactor <- mat
  for(nai in 1:length(naitems)) {
    matFactor[rownames(matFactor) %in% rownames(b[b[,2] == 13,])[1:20],naitems[nai]] <- NA
  }
  
  mmlml <- mml(composite~mlA + mlB + mlC + mlD + mlE, stuItems=stuItemsFactor, stuDat=stuDat2, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", weight="origwt", testScale=testDat, composite=FALSE)
  mmlmlsTaylor <- summary(mmlml, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit", gradientHessian=TRUE)
  mmlmls_coef_REF <- structure(c(287.342509, -8.764645, -7.781229, -18.667426, -1.636193,  55.627185, 34.659638,
                                 4.839807, 5.110929,  5.236816,  5.593414,  7.638541,  6.295088, 1.017022),
                               .Dim = c(7L, 2L), .Dimnames = list(c("(Intercept)", "mlA", "mlB", "mlC", "mlD", "mlE", "Population SD"),
                                                                  c("Estimate", "StdErr")))
  
  expect_equal(mmlmlsTaylor$coefficients[,1], mmlmls_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*200)
  expect_equal(mmlmlsTaylor$coefficients[,2], mmlmls_coef_REF[,2], tolerance=2*(.Machine$double.eps)^0.25)
})

context("missing data elements")
test_that("missing data elements", {

  mat1 <- mat
  mat1$m017401[1:1000] <- NA
  stuItems2 <- mat1[,1:13]
  stuItems2$oppID <- factor(rownames(mat1), levels=rownames(mat1))
  stuItems2 <- reshape(data=stuItems2, varying=c(dichotParamTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=dichotParamTab$ItemID, timevar="key")

  mml1B <- mml(composite ~ 1, stuItems=stuItems2, stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, composite=FALSE)
  mml1Bs <- summary(mml1B, gradientHessian=TRUE)

  # AM results with NA
  mml1Bs_coef_REF <- structure(c(281.55602373, 36.6399666578,
                                 0.973407969332, 0.940574435195),
                               .Dim = c(2L,2L), .Dimnames = list(c("(Intercept)", "Population SD"),
                                                            c("Estimate", "StdErr")))

  expect_equal(mml1Bs$coefficients[,2], mml1Bs_coef_REF[,2], tolerance=20*sqrt(.Machine$double.eps))

  # AM Taylor results with NA
  mml1BsT_coef_REF <- structure(c(281.55602373, 36.6399666578,
                                  0.96547714892, 0.877505088603),
                                .Dim = c(2L,2L), .Dimnames = list(c("(Intercept)", "Population SD"),
                                                                  c("Estimate", "StdErr")))

  mml1BTaylor <- summary(mml1B, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit", gradientHessian=TRUE)
  expect_equal(mml1BTaylor$coefficients[,1], mml1BsT_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*20)
  expect_equal(mml1BTaylor$coefficients[,2], mml1BsT_coef_REF[,2], tolerance=10*(.Machine$double.eps)^0.25)

  # missing data coded 8
  mat8 <- mat
  set.seed(142857) # make sure sample is the same
  randomRow <- sample(1:2000,1000,replace=TRUE) 
  randomCol <- sample(1:13,1000,replace=TRUE)
  # randomly set 1000 values to the missing value of 8
  for(i in 1:1000) {
    mat8[randomRow[i], randomCol[i]] <- 8
  }
  stuItems8 <- mat8[,1:13]
  stuItems8$oppID <- factor(rownames(mat8), levels=rownames(mat8))
  stuItems8 <- reshape(data=stuItems8, varying=c(dichotParamTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=dichotParamTab$ItemID, timevar="key")
  # it seems AM uses a missing value of c, so use that, note: weighted
  mml1A <- mml(composite~1, stuItems=stuItems8, stuDat=stuDat, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", weightVar="origwt", testScale=testDat, composite=FALSE)
  mml1As <- summary(mml1A, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit", gradientHessian=TRUE)
  mml1AsT_coef_REF <- structure(c(276.603996367644, 33.27212524293007,
                                  0.9749103126644552, 0.9559105270104056),
                                .Dim = c(2L,2L), .Dimnames = list(c("(Intercept)", "Population SD"),
                                                                  c("Estimate", "StdErr")))

  expect_equal(mml1As$coefficients[,1], mml1AsT_coef_REF[,1], tolerance=sqrt(.Machine$double.eps)*20)
  expect_equal(mml1As$coefficients[,2], mml1AsT_coef_REF[,2], tolerance=10*(.Machine$double.eps)^0.25)

  #missing data rows
  mat2 <- mat[1:1000,]
  mat2[,1:13] <- NA
  mat2 <- rbind(mat, mat2[1:1000,])
  stuItems3 <- mat2[,1:13]
  stuItems3$oppID <- factor(rownames(mat2), levels=rownames(mat2))
  stuItems3 <- reshape(data=stuItems3, varying=c(dichotParamTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=dichotParamTab$ItemID, timevar="key")

  rownames(mat) <- paste0("pseudo-student", 1:nrow(mat))
  stuDat3 <- rbind(stuDat, stuDat[1:1000, ])
  stuDat3$oppID <- rownames(mat2)
  mml1C <- mml(composite ~ 1, stuItems=stuItems3, stuDat=stuDat3, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, composite=FALSE)

  # test fit
  expect_equal(coef(mml1), coef(mml1C), tolerance=sqrt(.Machine$double.eps)*20)
  # test standard errors
  mml1Cs <- summary(mml1C, gradientHessian=TRUE)
  expect_equal(mml1s$coefficients, mml1Cs$coefficients, tolerance=(.Machine$double.eps)^0.25)
  mml1CRobust <- summary(mml1C, varType="robust")
  expect_equal(mml1CRobust$coefficients, mml1Robust$coefficients, tolerance=sqrt(.Machine$double.eps)*200)
  mml1CTaylor <- summary(mml1C, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit")
  expect_equal(mml1CTaylor$coefficients, mml1Taylor$coefficients, tolerance=sqrt(.Machine$double.eps)*20)
  mml1CCluster <- summary(mml1C, varType="cluster", clusterVar="repgrp1")
  expect_equal(mml1CCluster$coefficients, mml1Cluster$coefficients, tolerance=sqrt(.Machine$double.eps)*20)

  #missing data rows, with regressor

  mml2C <<- mml(composite~x1,stuItems=stuItems3, stuDat=stuDat3, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, composite=FALSE)

  # test fit
  expect_equal(coef(mml2), coef(mml2C), tolerance=sqrt(.Machine$double.eps)*20)
  # test standard errors
  mml2Cs <- summary(mml2C, gradientHessian=TRUE)
  expect_equal(mml2s$coefficients, mml2Cs$coefficients, tolerance=(.Machine$double.eps)^0.25)
  mml2CRobust <- summary(mml2C, varType="robust")
  expect_equal(mml2CRobust$coefficients, mml2Robust$coefficients, tolerance=sqrt(.Machine$double.eps)*200)
  mml2CTaylor <- summary(mml2C, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit", gradientHessian=TRUE)
  expect_equal(mml2CTaylor$coefficients, mml2Taylor$coefficients, tolerance=(.Machine$double.eps)^0.25)
  mml2CCluster <- summary(mml2C, varType="cluster", clusterVar="repgrp1")
  expect_equal(mml2CCluster$coefficients, mml2Cluster$coefficients, tolerance=sqrt(.Machine$double.eps)*200)
})

context("unsorted rows")
test_that("unsorted rows", {
  skip_on_cran()
  # this jumble should be undone. Most data will enter like this rather than pre-sorted
  jumble <- function(df) {
    df[sample(1:nrow(df),nrow(df)),]
  }
  #missing data rows
  mat2 <- mat[1:1000,]
  mat2[,1:13] <- NA
  mat2 <- rbind(mat, mat2[1:1000, ])
  stuItems3 <- mat2[ , 1:13]
  stuItems3$oppID <- factor(rownames(mat2), levels=rownames(mat2))
  stuItems3 <- reshape(data=stuItems3, varying=c(dichotParamTab$ItemID), idvar=c("oppID"), direction="long", v.names="score", times=dichotParamTab$ItemID, timevar="key")
  stuItems3$oppID <-  as.numeric(gsub("[a-zA-Z]|[-]", "", as.character(stuItems3$oppID)))
  stuItems3 <- jumble(stuItems3)

  rownames(mat) <- paste0("pseudo-student", 1:nrow(mat))
  stuDat3 <- rbind(stuDat, stuDat[1:1000,])
  stuDat3$oppID <- rownames(mat2)
  stuDat3$oppID <-  as.numeric(gsub("[a-zA-Z]|[-]", "", as.character(stuDat3$oppID)))
  stuDat3 <- jumble(stuDat3)
  dichotParamTab <- jumble(dichotParamTab)
  testDat <- jumble(testDat)
  mml2Cjumble <- mml(composite~x1,stuItems=stuItems3, stuDat=stuDat3, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, composite=FALSE)
  mml2Cjumble$scale <- testDat$scale[2]
  mml2Cjumble$location <- testDat$location[2]

  expect_equal(coef(mml2C), coef(mml2Cjumble), tolerance=sqrt(.Machine$double.eps)*20)
})

context("Draw PVs")
test_that("Draw PVs", {

  stuDatA <- (stuDat[order(as.character(stuDat$oppID)),])[1:1800,]
  stuItemsA <- subset(stuItems, oppID %in% stuDatA$oppID)
  mmlz <- mml(num ~ x1, stuItems=stuItemsA, stuDat=stuDatA, dichotParamTab=dichotParamTab, Q=34, idVar="oppID", testScale=testDat, composite=FALSE)
  set.seed(3)
  pv1 <- drawPVs.mmlMeans(mmlz, npv=1L, stochasticBeta=FALSE, newStuDat=NULL, newStuItems=NULL)
  set.seed(3)
  pv2 <- drawPVs.mmlMeans(mmlz, npv=1L, stochasticBeta=FALSE, newStuDat=stuDat, newStuItems=stuItemsA)
  # PVs (for first PV, when new data comes alphabetically at the end) are the same if new data is used or not
  expect_equal(pv1$data[1:1800,2], pv2$data[1:1800,2])
  expect_true(is.numeric(pv2$data[1801:2000, 2]))
  expect_true(!any(is.na(pv2$data[1801:2000, 2])))

  set.seed(3)
  pv3 <- drawPVs.mmlMeans(summary(mmlz), npv=3L, stochasticBeta=TRUE, newStuDat=stuDat, newStuItems=stuItemsA, returnPosterior=FALSE)
  expect_is(pv3, "DirePV")
  expect_is(pv3$data, "data.frame")
  expect_equal(dim(pv3$data), c(2000, 4))
  # TODO: add tests of dims


  set.seed(3)
  pv4 <- drawPVs.mmlMeans(mmlz, npv=3L, stochasticBeta=TRUE, newStuDat=stuDat, newStuItems=stuItemsA, returnPosterior=TRUE)
  expect_identical(names(pv4), c("posterior", "X", "rr1"))
  expect_equal(dim(pv4$posterior), c(2000, 1+3*2))
  expect_equal(dim(pv4$X), c(2000, 2))
  expect_equal(dim(pv4$rr1), c(34,2000))
})