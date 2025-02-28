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
dichotParamTab$missingCode <- 8
dichotParamTab$missingValue <- 0

## GPCM
polyParamTab <- structure(list(ItemID = structure(c(1:2), .Label = c("m0757cl", "m066501"), class = "factor"),
                               test = rep("composite",2),
                               subtest = c(rep("alg",2)),
                               slope = c( 0.43, 0.52 ), # could also be called "a"
                               itemLocation = c( -1.21, -0.96), # could also not exist, added to d1 ... dn
                               d0 = rep(0,2),
                               d1 = c(2.38, -0.56), 
                               d2 = c(-0.57, 0.56),
                               d3 = c(-1.18, NA),
                               D = c(1.7, 1.7),
                               scorePoints = c(3L, 2L) # number of score points, read d1 to d(n-1)
), class = "data.frame", row.names=c(1L,2L))
polyParamTab$missingCode <- 8
polyParamTab$missingValue <- 0
# read-in NAEP Primer data 
sdf <- readNAEP(system.file("extdata/data", "M36NT2PM.dat", package = "NAEPprimer"))
items <- c(as.character(dichotParamTab$ItemID), as.character(polyParamTab$ItemID))
edf <- EdSurvey::getData(data=sdf, varnames=c(items, "dsex", "origwt", "repgrp1", "jkunit", "sdracem"), dropOmittedLevels = FALSE)
edf$sid <- paste0("S",1:nrow(edf))
# apply super basic scoring
for(i in seq_along(items)) {
  coli <- items[i]
  # save the original
  edf[,paste0(coli,"raw")] <- edf[,coli]
  if( coli %in% dichotParamTab$ItemID) {
    edf[,coli] <- ifelse(grepl("[ABCDE]", edf[,paste0(coli,"raw")]), 0, NA)
    # for collapsed items
    edf[,coli] <- ifelse(grepl("Incorrect", edf[,paste0(coli,"raw")], fixed=TRUE), 0, edf[,coli])
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

# dummy code
indepVars <- c("dsex","sdracem")
dummyVars <- c() 
for (indep in indepVars) {
  dummyMatrix <- model.matrix(~get(indep) - 1, data = edf)
  edf <- cbind(edf, dummyMatrix)
  levs <- levels(edf[[indep]])
  colnames(edf) <- sub('^get\\(indep\\)', indep, colnames(edf)) # change to variableLEVEL
  edf <- edf[,!colnames(edf)%in% c(indep)] # drop original column and first dummy
  # record dummy names for subset
  newNames <- paste0(indep,levs)
  dummyVars <- c(dummyVars, newNames)
}

# AM doesn't like certain overlap of variable names so make a correction here 
colnames(edf)[colnames(edf) == "sdracemAmer Ind/Alaska Natv"] <- "sdracemInd"
dummyVars[dummyVars == "sdracemAmer Ind/Alaska Natv"] <- "sdracemInd"

# make stuItems (long) and stuDat 
stuItems <- stuItems1 <- edf[,c("sid", items)]
stuItems <- reshape(data=stuItems, varying=c(items), idvar=c("sid"), direction="long", v.names="score", times=items, timevar="key")
stuDat <- edf[,c("sid", "origwt", "repgrp1", "jkunit", dummyVars, items)]

# make testDat 
testDat <- data.frame(test=c("composite", "composite", "composite") ,
                      subtest=c("num", "alg", NA),
                      location=c(277.1563, 280.2948, sum(c(0.3,0.7) * c(277.1563, 280.2948))),
                      scale=c(37.7297, 36.3887, sum(c(0.3,0.7) * c(37.7297, 36.3887))),
                      subtestWeight=c(0.3,0.7, NA))


context("NAEPPrimer Regression Non-composite Algebra")
test_that("NAEPPrimer Regression Non-composite Algebra", {
  scoreTest <<- getFromNamespace("scoreTest", "Dire")
  mmlu <<- getFromNamespace("mmlu", "Dire")
  scored_test <- scoreTest(stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
                                 Q=34, idVar="sid", testScale=testDat, weightVar="origwt",
                                 minNode = -4, maxNode = 4, 
                                 strataVar="repgrp1", PSUVar="jkunit")
  mmlA2 <- mmlu(alg ~ dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                  sdracemInd + sdracemOther, scored_test, optimizer="EM")
  mmlA2Taylor <- summary(mmlA2, varType="Taylor")
  # EM gives different results, so tolerance is quite large
  expect_equal(unname(mmlA2Taylor$coefficients[,1]), 
               c(284.0244475898, 4.6425240064, -26.4952368029,-20.3300619109,
                 13.4723776156, -24.0877598983, -2.9570604020, 29.9108170108), 
               tolerance=0.01)
  expect_equal(unname(mmlA2Taylor$coef[,2]),
               c(1.5620600420, 1.6651582440,3.1195434525, 2.9454209419,
                 4.7778038300, 7.1463851501, 10.8693889849, 1.1684294466),
               tolerance=0.01)
  mmlA3 <- mmlu(alg ~ dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                  sdracemInd + sdracemOther, scored_test, optimizer="QN")
  mmlA3Taylor <- summary(mmlA3, varType="Taylor")
  expect_equal(unname(mmlA3Taylor$coefficients[,1]), 
               c(284.0244475898, 4.6425240064, -26.4952368029,-20.3300619109,
                 13.4723776156, -24.0877598983, -2.9570604020, 29.9108170108), 
               tolerance=2e4*sqrt(.Machine$double.eps))
  expect_equal(unname(mmlA3Taylor$coef[,2]),
               c(1.5620600420, 1.6651582440,3.1195434525, 2.9454209419,
                 4.7778038300, 7.1463851501, 10.8693889849, 1.1684294466),
               tolerance=5*(.Machine$double.eps)^0.25)
  
  # algebra
  suppressWarnings(mmlA <- mml(alg ~ dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                sdracemInd + sdracemOther, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
              Q=34, idVar="sid", testScale=testDat, composite=FALSE, weightVar="origwt",
              minNode = -4, maxNode = 4, 
              strataVar="repgrp1", PSUVar="jkunit"))
  mmlATaylor <- summary(mmlA, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit")
  
  # similarly, the PV results should be approximately equal
  set.seed(2)
  pvA <- drawPVs(mmlATaylor)
  lmA <- summary(lm(alg_dire1 ~  dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                      sdracemInd + sdracemOther, data=pvA, weights=pvA$origwt))
  expect_equal(coef(mmlATaylor)[1:6, 1], lmA$coef[1:6, 1], 0.1)
  
  # check the posterior means agree with the MML estimate
  pvApos <- mmlATaylor$posteriorEsts
  stuDatApos <- merge(stuDat, pvApos, by.x="sid", by.y="id", all.x=TRUE, all.y=FALSE)
  lmApos <- summary(lm(mu ~  dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                         sdracemInd + sdracemOther, data=stuDatApos, weights=stuDatApos$origwt))
  # this should be exact
  expect_equal(mmlATaylor$rawCoef[1:6], lmApos$coef[1:6, 1], tol=2e3*sqrt(.Machine$double.eps))
  
  expect_equal(unname(mmlATaylor$coefficients[,1]), 
               c(284.0244475898, 4.6425240064, -26.4952368029,-20.3300619109,
                 13.4723776156, -24.0877598983, -2.9570604020, 29.9108170108), 
               tolerance=2e3*sqrt(.Machine$double.eps))
  expect_equal(unname(mmlATaylor$coef[,2]),
               c(1.5620600420, 1.6651582440,3.1195434525, 2.9454209419,
                 4.7778038300, 7.1463851501, 10.8693889849, 1.1684294466),
               tolerance=5*(.Machine$double.eps)^0.25)
  
})

context("NAEPPrimer Regression Non-composite Numeric")
test_that("NAEPPrimer Regression Non-composite Numeric", {
  # Numeric 
  mmlN <- mml(num ~ dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + sdracemInd + sdracemOther, 
              stuDat=stuDat, 
              dichotParamTab=dichotParamTab, 
              polyParamTab=polyParamTab,
              Q=34, 
              idVar="sid", 
              testScale=testDat, 
              composite=FALSE, 
              weightVar="origwt",
              minNode = -4, 
              maxNode = 4,  
              strataVar="repgrp1", 
              PSUVar="jkunit")
  mmlNTaylor <- summary(mmlN, varType="Taylor", strataVar="repgrp1", PSUVar="jkunit")
  
  # compare estimate and SE 
  expect_equal(unname(mmlNTaylor$coefficients[,1]), 
               c(289.634387565144, -2.35572533192638, -28.7203722454655, -24.8520550271598, 
                 5.52373326613286, -17.5779673562621, -5.39532958283866, 34.6819277835647), 
               tolerance=200*sqrt(.Machine$double.eps))
  expect_equal(unname(mmlNTaylor$coef[,2]),
               c(1.56628867157975, 1.90197738244694, 2.85274599070158, 3.38928411761147, 
                 5.93325627970376, 8.21153414630517, 15.1023665060971, 1.29272963212231),
               tolerance=3*(.Machine$double.eps)^0.25)
  
})

context("NAEPPrimer Regression Composite, reduced RI")
test_that("NAEPPrimer Regression Composite", {
  # Composite 
  # Reduced information
  set.seed(555)
  mmlC_RI <- mml(composite ~ dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                `sdracemInd` + sdracemOther, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
                 Q=34, idVar="sid", testScale=testDat, composite=TRUE, weightVar="origwt",
                 minNode = -4, maxNode = 4, 
                 strataVar="repgrp1", PSUVar="jkunit",
                 retainedInformation = 0.8, verbose=FALSE)
  
  mmlC_RI_Taylor <- summary(mmlC_RI)
  
  expect_equal(unname(mmlC_RI_Taylor$coefficients[-7,1]), 
               c(285.328160709463, 735.844284546021, 705.874589968883,
                 1434.28529703307, 250.074669776892, 31.4135429802756), 
               tolerance=2*sqrt(.Machine$double.eps))
  
})
