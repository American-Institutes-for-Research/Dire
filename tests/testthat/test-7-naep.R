
# Notes on Running this Test: 
# This test is meant to compare Direct Estimation and AM results on real 
# NAEP data. It can only be run in environments with access to the Confidential-Internal 
# drive, or I-drive such as SPP. Paths are constructed under the assumption
# the test is being run on SPP and may need to be changed if using another 
# environment with access to real NAEP data. 

# Data used for test: 
# 2017 NAEP Math G4 data (M48NT1AT.dat)
# 2014 NAEP Civic G8 Data (C45NT2AT.dat)

# SPP locations: 
# The AM results are saved in SPP: 
# - I:\NCES\NCES_Dev\EdSurvey Project\Direct Estimation\QC\NAEP Direct Estimation QC.xlsx.
# Data files: 
# - I:\NCES\NCES_Dev\EdSurvey Project\Direct Estimation\QC\QC Data Files


if(FALSE) {
  # Set-Up
  context("Real NAEP Test")
  require(profvis)
  require(EdSurvey)
  library(tictoc)
  context("2017 NAEP Math G4 data")
  # 2017 NAEP Math G4 data (M48NT1AT.dat) 
  ### get NAEP Primer for tests 
  sdf <- readNAEP("I:/NCES/NCES_Dev/EdSurvey Project/Direct Estimation/QC/QC Data Files/2017 Math Y48MAT/Data/M48NT1AT.dat")
  searchSDF("race",sdf, levels = TRUE)
  # call mml.sdf for composite with dummy coded variable
  # saving AM files to working directory 
  # Note here that I have to specify a dctPath IRTparams couldn't find this test
  tic("intercept composite")
  mmlMcomp <- mml.sdf(composite~1, 
                      sdf, 
                      weightVar='origwt', 
                      verbose=T,
                      dctPath = "I:/NCES/NCES_Dev/EdSurvey Project/Direct Estimation/QC/QC Data Files/2017 Math Y48MAT/AM/M48NT1AT.dct",
                      composite = TRUE
                      
  )
  toc()
  
  tic("sex and race composite")
  mmlMcompSR <- mml.sdf(composite~dsex + sdracem, 
                        sdf, 
                        weightVar='origwt', 
                        verbose=T,
                        dctPath = "I:/NCES/NCES_Dev/EdSurvey Project/Direct Estimation/QC/QC Data Files/2017 Math Y48MAT/AM/M48NT1AT.dct",
                        composite = TRUE
                        
  )
  toc()
  
  
  
  mmlCsummary <- summary(mmlC)
  saveRDS(mmlC,"I:/NCES/NCES_Dev/EdSurvey Project/Direct Estimation/QC/DE Results/mml/mmlMcomp.rds")
  saveRDS(mmlCompSR,"I:/NCES/NCES_Dev/EdSurvey Project/Direct Estimation/QC/DE Results/mml/mmlMcompSR.rds")
  mmlC <- readRDS("I:/NCES/NCES_Dev/EdSurvey Project/Direct Estimation/QC/DE Results/mml/mmlMcomp.rds")
  mmlCompSR <- readRDS("I:/NCES/NCES_Dev/EdSurvey Project/Direct Estimation/QC/DE Results/mml/mmlMcompSR.rds")
  
  
  context("2014 NAEP Civics G8")
  # 2014 NAEP Math G4 data (C45NT2AT.dat) 
  ### get NAEP Primer for tests 
  sdf <- readNAEP("I:/NCES/NCES_Dev/EdSurvey Project/Direct Estimation/QC/QC Data Files/2014 Social Studies Y45HGC/Data/C45NT2AT.dat")
  searchSDF("race",sdf, levels = TRUE)
  # run intercept model
  tic("intercept")
  mmlC <- mml.sdf(civics~1, 
                  sdf, 
                  weightVar='origwt', 
                  verbose=T,
                  dctPath = "I:/NCES/NCES_Dev/EdSurvey Project/Direct Estimation/QC/QC Data Files/2014 Social Studies Y45HGC/AM/C45NT2AT.dct",
                  composite = FALSE
                  
  )
  toc()
  
  tic("race and sex civics model")
  mmlcSR <- mml.sdf(civics~dsex + sdracem, 
                    sdf, 
                    weightVar='origwt', 
                    verbose=T,
                    dctPath = "I:/NCES/NCES_Dev/EdSurvey Project/Direct Estimation/QC/QC Data Files/2014 Social Studies Y45HGC/AM/C45NT2AT.dct",
                    composite = FALSE
                    
  )
  toc()
  
}