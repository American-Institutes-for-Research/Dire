skip_on_cran()
require(EdSurvey)
# paramTab$a/slope <-
# paramTab$g/asymptote/guessing <- 
# paramTab$d/difficulty <- 
# param$D <- 1.7
dichotParamTab <- structure(list(ItemID = structure(1:11,
                                                    .Label = c("m109801", "m020001", "m111001", "m046301", "m046501", "m051501", "m111601", "m111301", "m111201", "m110801", "m110101"), class = "factor"),
                                 test = rep("composite",11),
                                 subtest = c(rep("num",6),rep("alg",5)),
                                 slope = c(0.96, 0.69, 0.83, 0.99, 1.03, 0.97, 1.45, 0.59, 0.34, 0.18, 1.20),
                                 difficulty = c(-0.37, -0.55, 0.85 , -0.97, -0.14, 1.21, 0.53, -1.84, -0.46, 2.43, 0.70),
                                 guessing = c(0.16,0.00,0.17,0.31,0.37,0.18, 0.28, 0.15, 0.09, 0.05, 0.18),
                                 D = rep(1.7, 11),
                                 MODEL = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "3pl", class = "factor")), class = "data.frame", row.names = c(NA, -11L))

## GPCM
# a / slope
# d1,2,...,(param-1)
# b/itemLocation, then map to d1, d2
# D / D
polyParamTab <- structure(list(ItemID = structure(c(1:2), .Label = c("m0757cl", "m066501"), class = "factor"),
                               test = rep("composite",2),
                               subtest = c(rep("alg",2)),
                               slope = c( 0.43, 0.52 ), # could also be called "a"
                               itemLocation = c( -1.21, -0.96), # could also not exist, added to d1 ... dn
                               d1 = c(2.38, -0.56), 
                               d2 = c(-0.57, 0.56),
                               d3 = c(-1.18, NA),
                               D = c(1.7, 1.7),
                               scorePoints = c(4L, 3L) # number of score points, read d1 to d(n-1)
), class = "data.frame", row.names=c(1L,2L))
# read-in NAEP Primer data 
sdf <- readNAEP(system.file("extdata/data", "M36NT2PM.dat", package = "NAEPprimer"))
items <- c(as.character(dichotParamTab$ItemID), as.character(polyParamTab$ItemID))
edf <- getData(data=sdf, varnames=c(items, "dsex", "origwt", "repgrp1", "jkunit", "sdracem"), omittedLevels = FALSE)
edf$sid <- paste0("S",1:nrow(edf))
# apply super basic scoring
for(i in 1:length(items)) {
  coli <- items[i]
  # save the original
  edf[,paste0(coli,"raw")] <- edf[,coli]
  if( coli %in% dichotParamTab$ItemID) {
    edf[,coli] <- ifelse(grepl("[ABCDE]", edf[,paste0(coli,"raw")]), 0, NA)
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

# dummify
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
stuDat <- edf[,c("sid", "origwt", "repgrp1", "jkunit", dummyVars)]

# make testDat 
testDat <- data.frame(test=c("composite", "composite", "composite") ,
                      subtest=c("num", "alg", NA),
                      location=c(277.1563, 280.2948, sum(c(0.3,0.7) * c(277.1563, 280.2948))),
                      scale=c(37.7297, 36.3887, sum(c(0.3,0.7) * c(37.7297, 36.3887))),
                      subtestWeight=c(0.3,0.7, NA))

context("NAEPPrimer Regression Non-composite Algebra")
# algebra
mmlA <- mml(alg ~ dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
              sdracemInd + sdracemOther, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
            Q=34, idVar="sid", testScale=testDat, composite=FALSE, weightVar="origwt",
            minNode = -4, maxNode = 4, 
            strataVar="repgrp1", PSUVar="jkunit")
mmlATaylor <- summary(mmlA, varType="Taylor", gradientHessian=TRUE, strataVar="repgrp1", PSUVar="jkunit")

# similarly, the PV results should be approximately equal
set.seed(2)
pvA <- drawPVs(mmlA, construct="alg")
stuDatA <- merge(stuDat, pvA[["data"]], by.x="sid", by.y="id", all.x=TRUE, all.y=FALSE)
lmA <- summary(lm(alg_dire1 ~  dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                  sdracemInd + sdracemOther, data=stuDatA, weights=stuDatA$origwt))
expect_equal(coef(mmlATaylor)[1:6, 1], lmA$coef[1:6, 1], 0.1)

# check stochastic beta works
pvAs <- drawPVs(mmlATaylor, construct="alg", stochasticBeta=TRUE)
stuDatAs <- merge(stuDat, pvAs[["data"]], by.x="sid", by.y="id", all.x=TRUE, all.y=FALSE)
lmAs <- summary(lm(alg_dire1 ~  dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                   sdracemInd + sdracemOther, data=stuDatAs, weights=stuDatA$origwt))
expect_equal(coef(mmlATaylor)[1:6, 1], lmAs$coef[1:6, 1], 0.1)

# check the posterior means agree with the MML estimate
pvApos <- drawPVs(mmlA, construct="alg", returnPosterior =TRUE)
pvApos <- pvApos[["posterior"]]
pvApos$mus <- pvApos$mu * mmlA$scale + mmlA$location
stuDatApos <- merge(stuDat, pvApos, by.x="sid", by.y="id", all.x=TRUE, all.y=FALSE)
lmApos <- summary(lm(mus ~  dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                     sdracemInd + sdracemOther, data=stuDatApos, weights=stuDatApos$origwt))
# this should be exact
expect_equal(coef(mmlATaylor)[1:6, 1], lmApos$coef[1:6, 1])

# MML Regression	
# Run completed on Wednesday, December 09, 2020. 09:54:09 AM
# 
# Selection:         	ALL
# Observations:      	8826
# Strata Variable:   	REPGRP1
# Cluster Variable:  	JKUNIT
# Weight Variable:   	ORIGWT
# 
# Iterations: 37 
# Log Likelihood: -14495.6
# 
# Adjusted Wald Test
# F(6,57) = 22.4731
# p(F > f) = 2.17826e-013	
# Dependent Variable: composite, alg	
# Parameter Name	Estimate	Standard Error	t Statistic	p > |t|	
# Constant	284.0244475898	1.5620600420	181.8268440084 0.0000000000	
# DSEXFEMA	4.6425240064	  1.6651582440	2.7880377274	0.0070480262	
# SDRACEMB	-26.4952368029	3.1195434525	-8.4933058976	0.0000000000	
# SDRACEMH	-20.3300619109	2.9454209419	-6.9022602582	0.0000000033	
# SDRACEMA	 13.4723776156	4.7778038300	2.8197845904	0.0064602479	
# SDRACEMI	-24.0877598983	7.1463851501	-3.3706215649	0.0013014462	
# SDRACEMO	 -2.9570604020	10.8693889849	-0.2720539679	0.7864836263	
# Root MSE	 29.9108170108	1.1684294466	---	---	
#   AM Statistical Software Beta Version 0.06.04. (c) The American Institutes for Research and Jon Cohen	
# compare estimate and weight 
expect_equal(unname(mmlATaylor$coefficients[,1]), 
             c(284.0244475898, 4.6425240064, -26.4952368029,-20.3300619109,
               13.4723776156, -24.0877598983, -2.9570604020, 29.9108170108), 
             tolerance=50*sqrt(.Machine$double.eps))
expect_equal(unname(mmlATaylor$coef[,2]),
             c(1.5620600420, 1.6651582440,3.1195434525, 2.9454209419,
               4.7778038300, 7.1463851501, 10.8693889849, 1.1684294466),
             tolerance=5*(.Machine$double.eps)^0.25)



context("NAEPPrimer Regression Non-composite Numeric")
# Numeric 
mmlN <- mml(num ~ dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + sdracemInd + sdracemOther, 
            stuItems=stuItems, 
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
mmlNTaylor <- summary(mmlN, varType="Taylor", gradientHessian=TRUE, strataVar="repgrp1", PSUVar="jkunit")


# MML Regression	
# Run completed on Wednesday, December 09, 2020. 10:00:59 AM
# 
# Selection:         	ALL
# Observations:      	6378
# Strata Variable:   	REPGRP1
# Cluster Variable:  	JKUNIT
# Weight Variable:   	ORIGWT
# 
# Iterations: 12 
# Log Likelihood: -8647.33
# 
# Adjusted Wald Test
# F(6,57) = 25.0037
# p(F > f) = 2.63123e-014	
# Dependent Variable: composite, num	
# Parameter Name	Estimate	Standard Error	t Statistic	p > |t|	
# Constant	295.4743460911	1.5484819999	190.8154864667	0.0000000000	
# DSEXFEMA	-2.6901100945	1.7664688656	-1.5228743324	0.1328644564	
# SDRACEMB	-26.8599727157	2.5390553373	-10.5787267893	0.0000000000	
# SDRACEMH	-22.3144730714	2.8487217696	-7.8331528582	0.0000000001	
# SDRACEMA	4.9481736285	5.9500552409	0.8316180990	0.4088183655	
# SDRACEMI	-17.7347803395	7.2513021373	-2.4457373316	0.0173294853	
# SDRACEMO	-4.9037488612	13.9236962237	-0.3521872915	0.7259035127	
# Root MSE	29.6919958757	1.3943320344	---	---	
#   AM Statistical Software Beta Version 0.06.04. (c) The American Institutes for Research and Jon Cohen	

# compare estimate and SE 
expect_equal(unname(mmlNTaylor$coefficients[,1]), 
             c(295.4743460911, -2.6901100945, -26.8599727157,-22.3144730714,
               4.9481736285, -17.7347803395, -4.9037488612, 29.6919958757), 
             tolerance=200*sqrt(.Machine$double.eps))
expect_equal(unname(mmlNTaylor$coef[,2]),
             c(1.5484819999, 1.7664688656,2.5390553373, 2.8487217696,
               5.9500552409, 7.2513021373, 13.9236962237, 1.3943320344),
             tolerance=3*(.Machine$double.eps)^0.25)


context("NAEPPrimer Regression Composite")
# Composite 
mmlC <- mml(composite ~ dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
              `sdracemInd` + sdracemOther, stuItems=stuItems, stuDat=stuDat, dichotParamTab=dichotParamTab, polyParamTab=polyParamTab,
            Q=34, idVar="sid", testScale=testDat, composite=TRUE, weightVar="origwt",
            minNode = -4, maxNode = 4, 
            strataVar="repgrp1", PSUVar="jkunit")
mmlCTaylor <- summary(mmlC, varType="Taylor", gradientHessian=TRUE)

### similarly, the PV results should be approximately equal
# draw the PVs
pvC <- drawPVs(mmlC)
stuDatC <- merge(stuDat, pvC[["data"]], by.x="sid", by.y="id", all.x=TRUE, all.y=FALSE)
# first for a subscale
lmA2 <- summary(lm(alg_dire1 ~  dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                   sdracemInd + sdracemOther, data=stuDatC, weights=stuDatC$origwt))
# this is the above, just Algebra mml estimate
expect_equal(coef(mmlATaylor)[1:6, 1], lmA2$coef[1:6, 1], 0.05)
# now use full composite
lmC2 <- summary(lm(composite_dire1 ~  dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                   sdracemInd + sdracemOther, data=stuDatC, weights=stuDatC$origwt))
# this is the above, just Algebra mml estimate
expect_equal(coef(mmlCTaylor)[1:6, 1], lmC2$coef[1:6, 1], 0.05)

# now draw with stochastic beta
pvCs <- drawPVs(mmlCTaylor, stochasticBeta=TRUE)
stuDatCs <- merge(stuDat, pvCs[["data"]], by.x="sid", by.y="id", all.x=TRUE, all.y=FALSE)
# first for a subscale
lmA2s <- summary(lm(alg_dire1 ~  dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                    sdracemInd + sdracemOther, data=stuDatCs, weights=stuDatC$origwt))
expect_equal(coef(mmlATaylor)[1:6, 1], lmA2s$coef[1:6, 1], 0.05)
# now composite
lmC2s <- summary(lm(composite_dire1 ~  dsexFemale + sdracemBlack + sdracemHispanic + `sdracemAsian/Pacific Island` + 
                    sdracemInd + sdracemOther, data=stuDatCs, weights=stuDatC$origwt))
expect_equal(coef(mmlCTaylor)[1:6, 1], lmC2s$coef[1:6, 1], 0.05)


# MML Composite Regression	
# Run completed on Wednesday, December 09, 2020. 12:54:46 PM
# 
# Selection:         	ALL
# Observations:      	16915
# Strata Variable:   	REPGRP1
# Cluster Variable:  	JKUNIT
# Weight Variable:   	ORIGWT
# 
# Adjusted Wald Test
# F(6,57) = 22.1717
# p(F > f) = 2.83107e-013	
# Composite Results	
# Parameter Name	Estimate	Standard Error	t Statistic	p > |t|	
#   Constant	287.4594172767	1.5217659456	5.3268139977	0.0000014991	
# DSEXFEMA	2.4427336887	1.6550653899	1.4759137033	0.1450151434	
# SDRACEMB	-26.6046575811	2.8827206928	-9.2290098195	0.0000000000	
# SDRACEMH	-20.9253852827	2.8492971588	-7.3440515737	0.0000000006	
# SDRACEMA	10.9151163178	5.0000511649	2.1830009250	0.0328517195	
# SDRACEMI	-22.1818659549	7.0094097999	-3.1645839790	0.0024164736	
# SDRACEMO	-3.5410669629	11.4857615243	-0.3083005820	0.7588921216	
# Root MSE	29.1509253186	(--)	
# AM Statistical Software Beta Version 0.06.04. (c) The American Institutes for Research and Jon Cohen	

# Dependent Variable: composite, alg	
# Parameter Name	Estimate	Standard Error	
# Constant	     0.102	     0.043	
# DSEXFEMA	     0.128	     0.046	
# SDRACEMB	    -0.728	     0.086	
# SDRACEMH	    -0.559	     0.081	
# SDRACEMA	     0.370	     0.131	
# SDRACEMI	    -0.662	     0.196	
# SDRACEMO	    -0.081	     0.299	
# Root MSE	     0.822	     0.032	
# Dependent Variable: composite, num	
# Parameter Name	Estimate	Standard Error	
# Constant	     0.486	     0.041	
# DSEXFEMA	    -0.071	     0.047	
# SDRACEMB	    -0.712	     0.067	
# SDRACEMH	    -0.591	     0.076	
# SDRACEMA	     0.131	     0.158	
# SDRACEMI	    -0.470	     0.192	
# SDRACEMO	    -0.130	     0.369	
# Root MSE	     0.787	     0.037	

# DIRE to AM
expect_equal(unname(mmlCTaylor$coefficients[-8,1]), 
             c(287.4594172767, 2.4427336887, -26.6046575811, -20.9253852827, 
               10.9151163178, -22.1818659549, -3.5410669629), 
             tolerance=50*sqrt(.Machine$double.eps)) 

# AM gets wrong correlations, reported for posterity
# expect_equal(unname(mmlCTaylor$coefficients[1:7,2]), 
#              c(1.5217659456, 1.6550653899, 2.8827206928, 
#                2.8492971588, 5.0000511649, 7.0094097999,11.4857615243), 
#              tolerance=5*(.Machine$double.eps)^0.25) 
