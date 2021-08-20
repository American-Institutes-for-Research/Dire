\dontrun{
require(EdSurvey)

# 1) make param tab for dichotomous items
dichotParamTab <- data.frame(ItemID = c("m109801", "m020001", "m111001",
                                        "m046301", "m046501", "m051501",
                                        "m111601", "m111301", "m111201",
                                        "m110801", "m110101"),
                             test = rep("composite",11),
                             subtest = c(rep("num",6),rep("alg",5)),
                             slope = c(0.96, 0.69, 0.83,
                                       0.99, 1.03, 0.97,
                                       1.45, 0.59, 0.34,
                                       0.18, 1.20),
                             difficulty = c(-0.37, -0.55,  0.85,
                                            -0.97, -0.14,  1.21,
                                             0.53, -1.84, -0.46,
                                             2.43,  0.70),
                             guessing = c(0.16, 0.00, 0.17,
                                          0.31, 0.37, 0.18,
                                          0.28, 0.15, 0.09,
                                          0.05, 0.18),
                             D = rep(1.7, 11),
                             MODEL = rep("3pl", 11))

# param tab for GPCM items
polyParamTab <- data.frame(ItemID = factor(c("m0757cl", "m066501")),
                           test = rep("composite",2),
                           subtest = rep("alg",2),
                           slope = c(0.43, 0.52), # could also be called "a"
                           itemLocation = c(-1.21, -0.96), # added to d1 ... dn
                           d1 = c(2.38, -0.56), 
                           d2 = c(-0.57, 0.56),
                           d3 = c(-1.18, NA),
                           D = c(1.7, 1.7),
                           scorePoints = c(4L, 3L)) # number of score points, read d1 to d(n-1)
# read-in NAEP Primer data 
sdf <- readNAEP(system.file("extdata/data", "M36NT2PM.dat", package = "NAEPprimer"))
# read in these items
items <- c(as.character(dichotParamTab$ItemID), as.character(polyParamTab$ItemID))
# dsex, student sex
# origwt, full sample weights
# repgrp1, stratum indicator
# jkunit, PSU indicator
edf <- getData(data=sdf, varnames=c(items, "dsex", "origwt", "repgrp1", "jkunit", "sdracem"),
               omittedLevels = FALSE, returnJKreplicates=FALSE)
# make up a student ID
edf$sid <- paste0("S",1:nrow(edf))
# apply simplified scoring
for(i in 1:length(items)) {
  coli <- items[i]
  # save the original
  rawcol <- paste0(coli,"raw")
  edf[,rawcol] <- edf[,coli]
  if( coli %in% dichotParamTab$ItemID) {
    edf[,coli] <- ifelse(grepl("[ABCDE]", edf[,rawcol]), 0, NA)
    edf[,coli] <- ifelse(grepl("*", edf[,rawcol]), 1, edf[,coli])
  } else {
    # scale for m066501
    edf[,coli] <- ifelse(grepl("Incorrect", edf[,rawcol]), 0, NA)
    edf[,coli] <- ifelse(grepl("Partial", edf[,rawcol]), 1, edf[,coli])
    edf[,coli] <- ifelse(grepl("Correct", edf[,rawcol]), 2, edf[,coli])
    # scale for m0757cl
    edf[,coli] <- ifelse(grepl("None correct", edf[,rawcol]), 0, edf[,coli])
    edf[,coli] <- ifelse(grepl("One correct", edf[,rawcol]), 1, edf[,coli])
    edf[,coli] <- ifelse(grepl("Two correct", edf[,rawcol]), 2, edf[,coli])
    edf[,coli] <- ifelse(grepl("Three correct", edf[,rawcol]), 3, edf[,coli])
  }
  edf[,rawcol] <- NULL # delete original
}

# stuItems has one row per student/item combination
stuItems <- edf[,c("sid", items)]
stuItems <- reshape(data=stuItems, varying=c(items), idvar=c("sid"),
                    direction="long", v.names="score", times=items, timevar="key")
# stuDat is one row per student an contains student-level information
stuDat <- edf[,c("sid", "origwt", "repgrp1", "jkunit", "dsex", "sdracem")]

# testDat shows scaling and weights for subtests, an overall score can be treated as a subtest
testDat <- data.frame(test=c("composite", "composite") ,
                      subtest=c("num", "alg"),
                      location=c(277.1563, 280.2948),
                      scale=c(37.7297, 36.3887),
                      subtestWeight=c(0.3,0.7))

# estimate a regression for Algebra subscale
mmlA <- mml(alg ~ dsex,
            stuItems=stuItems, stuDat=stuDat,
            dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
            testScale=testDat,
            idVar="sid", weightVar="origwt", # these are column names on stuDat
            strataVar="repgrp1", PSUVar="jkunit")
# summary, with Taylor standard errors
mmlAs <- summary.mmlMeans(mmlA, varType="Taylor")


# estimate a regression for Numeracy subscale
mmlN <- mml(num ~ dsex,
            stuItems=stuItems, stuDat=stuDat,
            dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
            testScale=testDat,
            idVar="sid", weightVar="origwt", # these are column names on stuDat
            strataVar="repgrp1", PSUVar="jkunit")
# summary, with Taylor standard errors
mmlNs <- summary.mmlMeans(mmlN, varType="Taylor")
mmlNs

set.seed(2)
head(pvd <- drawPVs.mmlMeans(mmlA))
head(pvs <- drawPVs.mmlMeans(summary.mmlMeans(mmlA, varType="Taylor"), stochasticBeta=T))
sd(apply(pvd[,-1], 2, mean))
sd(apply(pvs[,-1], 2, mean)) # should be larger, accounts for uncertainty in beta


# composite regression 
mmlC <- mml(composite ~ dsex ,
            stuItems=stuItems, stuDat=stuDat,
            dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
            testScale=testDat,
            idVar="sid", weightVar="origwt", # these are column names on stuDat
            strataVar="repgrp1", PSUVar="jkunit")
# summary, with Taylor standard errors
summary(mmlC, varType="Taylor")

set.seed(2)
head(pvd <- drawPVs.mmlCompositeMeans(mmlC))
mmlCs <- summary.mmlCompositeMeans(mmlC, varType="Taylor")
head(pvs <- drawPVs.mmlCompositeMeans(mmlCs, stochasticBeta=TRUE))


}