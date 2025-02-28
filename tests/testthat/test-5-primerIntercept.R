context("NAEPPrimer composite")
skip_if_not_installed("EdSurvey")
require(EdSurvey)

dichotParamTab <- structure(list(ItemID = structure(1:11,
                                                    .Label = c("m109801", "m020001", "m111001", "m046301", "m046501", "m051501", "m111601", "m111301", "m111201", "m110801", "m110101"), class = "factor"),
                                 test = rep("composite",11),
                                 subtest = c(rep("num",6),rep("alg",5)),
                                 slope = c(0.96, 0.69, 0.83, 0.99, 1.03, 0.97, 1.45, 0.59, 0.34, 0.18, 1.20),
                                 difficulty = c(-0.37, -0.55, 0.85 , -0.97, -0.14, 1.21, 0.53, -1.84, -0.46, 2.43, 0.70),
                                 guessing = c(0.16,0.00,0.17,0.31,0.37,0.18, 0.28, 0.15, 0.09, 0.05, 0.18),
                                 D = rep(1.7, 11),
                                 MODEL = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "3pl", class = "factor")), class = "data.frame", row.names = c(NA, -11L))

polyParamTab <- structure(list(ItemID = structure(c(1:2), .Label = c("m0757cl", "m066501"), class = "factor"),
                               test = rep("composite",2),
                               subtest = c(rep("alg",2)),
                               slope = c( 0.43, 0.52 ), # could also be called "a"
                               itemLocation = c( -1.21, -0.96), # could also not exist, added to d1 ... dn
                               d0 = c(0,0),
                               d1 = c(2.38, -0.56), 
                               d2 = c(-0.57, 0.56),
                               d3 = c(-1.18, NA),
                               D = c(1.7, 1.7),
                               scorePoints = c(3L, 2L) # number of score points, read d1 to d(n-1)
), class = "data.frame", row.names=c(1L,2L))

## Test Data information
testDat <- data.frame(
  subtest=c("num", "alg"),
  location=c(277.1563, 280.2948),
  scale=c(37.7297, 36.3887),
  subtestWeight=c(0.3,0.7),
  test=c("composite","composite"))

## Read-in NAEP Primer and Subset to items in Parameter Tabs 
sdf <- readNAEP(system.file("extdata/data", "M36NT2PM.dat", package = "NAEPprimer"))
items <- c(as.character(dichotParamTab$ItemID), as.character(polyParamTab$ItemID))
edf <- EdSurvey::getData(data=sdf, varnames=c(items, "dsex", "origwt", "repgrp1", "jkunit"))

edf$stu_id <- paste0("S",1:nrow(edf))

# apply super basic scoring to Primer 
for(i in 1:length(items)) {
  coli <- items[i]
  # save the original
  edf[,paste0(coli,"raw")] <- edf[,coli]
  if( coli %in% dichotParamTab$ItemID) {
    edf[,coli] <- ifelse(grepl("[ABCDE]", edf[,paste0(coli,"raw")]), 0, NA)
    edf[,coli] <- ifelse(grepl("Incorrect", edf[,paste0(coli,"raw")]), 0, edf[,coli])
    edf[,coli] <- ifelse(grepl("*", edf[,paste0(coli,"raw")], fixed=TRUE), 1, edf[,coli])
  } else {
    edf[,coli] <- ifelse(grepl("Incorrect", edf[,paste0(coli,"raw")], fixed=TRUE), 0, NA)
    edf[,coli] <- ifelse(grepl("None correct", edf[,paste0(coli,"raw")], fixed=TRUE), 0, edf[,coli])
    edf[,coli] <- ifelse(grepl("Incorrect", edf[,paste0(coli,"raw")], fixed=TRUE), 0, edf[,coli])
    edf[,coli] <- ifelse(grepl("One correct", edf[,paste0(coli,"raw")], fixed=TRUE), 1, edf[,coli])
    edf[,coli] <- ifelse(grepl("Partial", edf[,paste0(coli,"raw")], fixed=TRUE), 1, edf[,coli])
    edf[,coli] <- ifelse(grepl("Two correct", edf[,paste0(coli,"raw")], fixed=TRUE), 2, edf[,coli])
    edf[,coli] <- ifelse(grepl("Correct", edf[,paste0(coli,"raw")], fixed=TRUE), 2, edf[,coli])
    edf[,coli] <- ifelse(grepl("Three correct", edf[,paste0(coli,"raw")], fixed=TRUE), 3, edf[,coli])
  }
  # print results
  edf[,paste0(coli,"raw")] <- NULL # delete original
}

# Separate items and student data from eachother 
stuDat <- edf[,c("stu_id", "dsex", "origwt", "repgrp1", "jkunit", items)]


context("NAEP primer, composite as weighted sum of subscales")
test_that("NAEP primer, composite as weighted sum of subscales", {
  suppressWarnings(mmlA <- mml(alg ~ 1, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
                               Q=34, idVar="stu_id", testScale=testDat, composite=FALSE, weightVar="origwt",
                               minNode = -5, maxNode = 5,
                               strataVar="repgrp1", PSUVar="jkunit"))
  mmlATaylor <- summary(mmlA, varType="Taylor")
  
  expect_equal(unname(mmlATaylor$coef[,1]), c(278.334734151, 32.173725366), tolerance=(.Machine$double.eps)^0.25)
  expect_equal(unname(mmlATaylor$coef[1,2]), c(1.035753180), tolerance=5*(.Machine$double.eps)^0.25)
  expect_equal(mmlA$LogLik, -16394.5, tol=1E-5)
  
  # numeracy score
  suppressWarnings(mmlN <- mml(num ~ 1, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
                               Q=34, idVar="stu_id", testScale=testDat, composite=FALSE, weightVar="origwt",
                               minNode = -5, maxNode = 5, strataVar="repgrp1", PSUVar="jkunit"))
  
  # composite
  suppressWarnings(mmlC <- mml(composite ~ 1, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
              Q=34, idVar="stu_id", testScale=testDat, composite=TRUE, weightVar="origwt",
              minNode = -5, maxNode = 5, 
              strataVar="repgrp1", PSUVar="jkunit"))
  # check that composite coefficients really are composite
  # only [1] because the SD is not formed in this way
  expect_equal(0.3* coef(mmlN)[1] + 0.7*coef(mmlA)[1], coef(mmlC)[1])
  
  expect_equal(unname(coef(mmlC)[1]), 278.484489655905, tolerance=1e-6)
  expect_equal(unname(coef(mmlC)[2]), 33.7446804726567, tolerance=1e-4)
  set.seed(555)
  mmlCTaylor <- summary(mmlC)
  expect_equal(mmlCTaylor$coef[1,2], c(0.946295009278728), 
               tolerance=20*(.Machine$double.eps)^0.25)
})

