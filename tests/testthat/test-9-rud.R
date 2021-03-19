if(FALSE) {
  library(foreach)
  library(iterators)
  library(doParallel)
  library(numDeriv)
  library(minqa)
  registerDoParallel(cores = 4)
  set.seed(142857)
  
  # Returns the expected value (probability of correct response)
  pl <- function(theta, b) {
    1 / (1 + exp(-1.7* (theta - b)))
  }
  
  
  ######################### Test data #########################
  library(EdSurvey)
  library(dplyr)
  
  ### functions ###
  getScoreCard <- function(filename, itemColsClean, adjustedData) {
    # given fr2 file, list of question items, 
    # produce a score card of question id, answer label, points awarded
    
    # read NAEP machine readable file.
    # this has layout information on it
    t <- try(mrcFile <- readLines(filename), silent=TRUE)
    
    # Split mrcFile by "\n"
    mrcFile <- strsplit(mrcFile , "\n", fixed=T)
    mrcFile <- unlist(mrcFile)
    
    # read in the variables from the file, this is based on the information ETS
    # shared with us
    variableName <- trimws(substring(mrcFile, 1, 8)) # name of the variable
    variableName <- unlist(lapply(variableName, tolower))
    
    # indices/lines of question item variables
    indices <-  match(itemColsClean, variableName)
    itemLines <- mrcFile[indices]
    
    # initialize empty dataframe to build
    scoreCard <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(scoreCard) <- c('key', 'answer', 'score')
    
    # for mulitple choice answer questions
    # temporary scoring where all else set to "0", until figure out how to deal with "8" and "."
    scoreDictMult <- list(resCat=c("Multiple", "Not Reached", "Missing", "Omitted", "Illegible", "Non-Rateable", "Off Task"),
                          point=c(".", ".", ".", ".", ".", ".", "."))
    #scoreDictMult <- list(resCat=c("Multiple", "Not Reached", "Missing", "Omitted", "Illegible", "Non-Rateable", "Off Task"),
    #                      point=c("0", "0", "0", "0", "0", "0", "0"))
    # scoreDictMult <- list(resCat=c("Multiple", "Not Reached", "Missing", "Omitted", "Illegible", "Non-Rateable", "Off Task"), 
    #                       point=c("8", ".", ".", "8", "0", "0", "0"))
    
    # for constructed answer questions
    # temporary scoring where all else set to "0", until figure out how to deal with "8" and "."
    scoreDictConst <- list(resCat=c("Multiple", "Not Reached", "Missing", "Omitted", "Illegible", "Non-Rateable", "Off Task"), 
                           point=c(".", ".", ".", ".", ".", ".", "."))
    # scoreDictConst <- list(resCat=c("Multiple", "Not Reached", "Missing", "Omitted", "Illegible", "Non-Rateable", "Off Task"), 
    #                        point=c("0", "0", "0", "0", "0", "0", "0"))
    # scoreDictConst <- list(resCat=c("Multiple", "Not Reached", "Missing", "Omitted", "Illegible", "Non-Rateable", "Off Task"), 
    #                        point=c("8", ".", ".", "0", "0", "0", "0"))
    
    # items with scores that have been adjusted
    adjustedItems <- adjustedData$NAEP.ID
    
    # get item key, answers, and points per item
    for (line in itemLines) {
      # initialize portion of scorecard for item 
      subScoreCard <- data.frame(matrix(ncol = 3, nrow = 9))
      colnames(subScoreCard) <- c('key', 'answer', 'score')
      
      # get key (question id)
      itemId <- tolower(substring(line, 1,7))
      
      # get points awarded by splitting string (like 0010)
      points <- unlist(strsplit(substring(line, 71,90), ' '))[1]
      
      # switch to new score points if item scoring has been adjusted
      if (itemId %in% adjustedItems) {
        pointsNew <- trimws(gsub(';','', adjustedData[adjustedData$NAEP.ID==itemId, 'To']))
        # check that number of points match old numbers
        if (nchar(points) != nchar(pointsNew)) {
          stop(paste0('Number of adjusted scores for item ', itemId, ' is different from original.'))
        } else {
          points <- pointsNew
        }
      }
      
      # split and process points
      if (nchar(points) != 0) {
        points <- unlist(strsplit(points,''))
        points <- replace(NA * c(1:9), c(1:length(points)), points)
      } else {
        points <- NA * c(1:9)
      }
      
      # get answers (like A, B *, C)
      firstLabel <- getLabel(line, 93,118)
      secondLabel <- getLabel(line, 121, 146)
      thirdLabel <- getLabel(line, 149, 174)
      fourthLabel <- getLabel(line, 177, 202)
      fifthLabel <- getLabel(line, 205, 230)
      sixthLabel <- getLabel(line, 233, 258)
      seventhLabel <- getLabel(line, 261,286)
      eighthLabel <- getLabel(line, 289,314)
      ninthLabel <- getLabel(line, 317,342)
      
      labels <- c(firstLabel, secondLabel, thirdLabel, fourthLabel, fifthLabel, sixthLabel, seventhLabel, eighthLabel, ninthLabel)
      
      # set points for other types of answers (like omitted, illegible, etc.)
      if (firstLabel %in% c('A', 'A*', 'A *')) {
        # this is a multiple choice question
        for (l in 1:length(labels)) {
          # find the corresponding point score in dict
          newScore <- scoreDictMult$point[match(labels[l], scoreDictMult$resCat)]
          if (! is.na(newScore)) {
            points[l] <- newScore
          }
        }
      } else {
        # this is a constructed answer question
        for (l in 1:length(labels)) {
          newScore <- scoreDictMult$point[match(labels[l], scoreDictMult$resCat)]
          if (! is.na(newScore)) {
            points[l] <- newScore
          }
        }
      }
      
      # build sub scorecard and bind 
      subScoreCard$answer <- labels
      subScoreCard$score <- points
      subScoreCard$key <- itemId
      scoreCard <- rbind(scoreCard, subScoreCard)
    }
    return(scoreCard)
  }
  
  getLabel <- function(line, first, last) {
    #get answer labels (A, B *, C, Correct, Omitted, etc.)
    part <- substring(line, first, last) #get the part of string needed with first and last index
    return (trimws(gsub("^(.+?)(?=\\d).*", "\\1", part, perl = TRUE))) #grab label (everything before digits) and trim white space
  }
  
  
  ### readins ###
  # data
  sdf <- readNAEP('R:/data/NAEP/Math/2003/NAEP 2003 Math G8/data/M34NT2AT.dat')
  data <- getData(sdf, varnames=colnames(sdf), addAttributes = TRUE,
                  omittedLevels = FALSE, defaultConditions = FALSE)
  data <- data[1:2,]
  
  # file format for point scoring
  fr2 <- 'R:/data/NAEP/Math/2003/NAEP 2003 Math G8/select/parms/M34NT2AT.fr2'
  
  # item parameters
  params <- read.csv('//dc2fs/DC4Work/ESSIN Task 14/NAEP R Program/Sunjoo/Scoring/master.csv', stringsAsFactors = FALSE)
  
  # item adjustments
  adjust <- read.csv('//dc2fs/DC4Work/ESSIN Task 14/NAEP R Program/Sunjoo/Scoring/adjustments.csv', stringsAsFactors = FALSE)
  
  
  ### paramTab input ###
  # Get test info
  subject <- sdf$subject
  glevel <- sdf$gradeLevel
  year <- sdf$year
  
  # Switch subject wording to match param table data (may need to add more here for other subjects)
  dict <- list(sdf=c('Mathematics','Reading'), param=c('math','read'))
  idx <- grep(subject, dict$sdf)
  subject <- dict$param[idx]
  
  level <- strsplit(glevel, " ")
  levType <- tolower(level[[1]][1])
  lev <- as.integer(level[[1]][2])
  
  # From the master params file, filter only needed params for the subject, grade, year
  paramTab <- params %>% filter(Subject==subject, Level.Type==levType, Level==lev, Year==year)
  paramTab <- paramTab[, c('NAEP.ID','aj','bj','cj','dj1','dj2','dj3','dj4','dj5')]
  paramTab$NAEP.ID <- tolower(paramTab$NAEP.ID)
  
  # Rename, add columns, make numeric
  colnames(paramTab) <- c('ItemID','P0','P1','P2','P3','P4','P5','P6','P7')
  paramTab$P8 <- NA
  paramTab$P9 <- NA
  paramTab$P10 <- NA
  paramTab[paste0('P', 0:10)] <- sapply(paramTab[paste0('P', 0:10)], as.numeric)
  
  # If P2 is missing, P3 becomes P2 (etc.)
  paramTab[is.na(paramTab$P2), paste0('P',2:9)] <- paramTab[is.na(paramTab$P2), paste0('P',3:10)]
  
  # for making any adjustments to items (deleted items, change in scoring)
  adjustData <- adjust[adjust$Level==lev & adjust$Level.Type==levType & adjust$Subject==subject & adjust$Year==year,] #filter from master-adjustment file
  deletedItems <- tolower(adjustData[adjustData$Adjustment=='Deleted','NAEP.ID']) #deleted item keys
  adjustedData <- adjustData[adjustData$Adjustment=='Collapsed',c('NAEP.ID','From','To')] #adjust item key and adjustments
  adjustedData$NAEP.ID <- tolower(adjustedData$NAEP.ID)
  
  # get question item IDs 
  itemCols <- paramTab$ItemID
  itemColsClean <- itemCols[!itemCols %in% deletedItems] #take out deleted items
  
  # take out deleted items in paramTab
  paramTab <- paramTab[paramTab$ItemID %in% itemColsClean, ]
  
  # get score card (giving item id, answer, points for answer)
  sCard <- getScoreCard(fr2, itemColsClean, adjustedData)
  sCard$score <- as.numeric(sCard$score)
  sCard <- sCard[complete.cases(sCard),] # take out rows with NA ScorePoints
  colnames(sCard) <- c('key', 'answer', 'ScorePoints')
  
  # extra configuration for partial scores (make labeling uniform)
  partialKeys <- sCard[sCard$answer=='Partial','key']
  partialNums <- table(partialKeys) 
  twoPartials <- names(partialNums[partialNums==2]) #items with 2 partial scores
  for (item in twoPartials) {
    sCard[sCard$key==item & sCard$answer=='Partial', 'answer'] <- c('Partial1','Partial2')
  }
  threePartials <-names(partialNums[partialNums==3]) #items with 3 partial scores
  for (item in threePartials) {
    sCard[sCard$key==item & sCard$answer=='Partial', 'answer'] <- c('Partial1','Partial2','Partial3')
  }
  
  # get max vals per item
  keyMaxVals <- aggregate(ScorePoints ~ key, sCard, max)
  colnames(keyMaxVals) <- c('ItemID', 'ScorePoints')
  
  # merge max vals into paramTab
  paramTab <- merge(paramTab, keyMaxVals, by='ItemID')
  
  # add MODEL type (though this is not used anywhere in the function)
  paramTab$MODEL <- ifelse(paramTab$ScorePoints==1, '3pl', 'pcl')
  
  ### stuItems input ###
  # Subset student data to only question items (take out admin questions, etc)
  itemDF <- data[itemColsClean]
  
  # Find questions answered per student (not NA for question)
  rowColNotNA <- which(!is.na(itemDF), arr.ind=TRUE)
  preStuItems <- as.data.frame(rowColNotNA)
  colnames(preStuItems) <- c('oppID', 'key') #student's row number in data is its oppID
  
  # Change index to question names, sort by oppID
  preStuItems$key <- itemColsClean[preStuItems$key]
  preStuItems <- preStuItems[order(preStuItems$oppID), ]
  
  # Get student's actual answers for each question answered for point-scoring
  ans <- apply(itemDF, 1, function(x) x[!is.na(x)])
  ans <- unlist(ans)
  preStuItems$answer <- ans
  row.names(preStuItems) <- NULL
  
  # merge student data with scorecard to get point scores for answers
  stuItems <- merge(preStuItems, sCard, by=c('key', 'answer'), all=TRUE)
  
  # cases where score card answers and student answers don't match exactly are likely wrong answers
  stuItems[is.na(stuItems$ScorePoints), 'ScorePoints'] <- 0
  
  # drop answer column -- no longer needed, change colname to match
  stuItems <- stuItems[, c('key', 'oppID', 'ScorePoints')]
  colnames(stuItems) <- c('key', 'oppID', 'score')

  
  ### stuDat input ###
  stuDat <- data[unique(stuItems$oppID), c('dsex', 'srace', 'iep', 'origwt','repgrp1','jkunit')] 
  stuDat$oppID <- row.names(stuDat) #again, row number is student's id
  
  
  ### final cleanup ###
  # keep only complete cases in stuDat (drop row if null value somewhere)
  stuDat <- stuDat[complete.cases(stuDat), ]
  
  # take out students with incomplete cases from stuItems
  stuItems <- stuItems[stuItems$oppID %in% stuDat$oppID, ]
  
  # take out items from paramTab if only answered by deleted students
  paramTab <- paramTab[paramTab$ItemID %in% stuItems$key, ]
  
  
  ### subset to numbers subtest
  library(readxl)
  num <- read_xls('AM_DCT_math8_number.xls')
  num_id <- tolower(num$ID)
  stuItems <- stuItems %>% filter(key %in% num_id)
  paramTab <- paramTab %>% filter(ItemID %in% num_id)
  
  ### take out omitted answer level
  stuItems <- stuItems %>% filter(key!='m066401')
  paramTab <- paramTab %>% filter(ItemID!='m066401')
  
  ### only 1 student
  stuItems <- stuItems %>% filter(oppID == 1)
  paramTab <- paramTab %>% filter(ItemID %in% stuItems$key)
  stuDat <- stuDat[1,]
  
  
  ############### test functions ###############

  resultRegIntercept <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType = 'regression')
  resultWeightRegIntercept <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType = 'regression', weightVar = 'origwt')
  resultWeightPlusRegIntercept <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType = 'regression', weightVar = 'origwt', clusterVar = 'jkunit', stratavar = 'repgrp1')
  
  #resultPopMean <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType = 'regression', ind.vars=c('dsex', 'srace', 'iep'))  
  #resultReg <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType = 'regression', ind.vars=c('dsex', 'srace', 'iep'))
  #system.time(resultReg <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType = 'regression', ind.vars = c('dsex', 'srace')))
  #resultWeightReg <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType = 'regression', ind.vars = c('dsex', 'srace', 'iep'), weightVar = 'origwt')
  #resultWeightPlusReg <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType = 'regression', ind.vars = c('dsex', 'srace', 'iep'), weightVar = 'origwt', clusterVar = 'jkunit', stratavar = 'repgrp1')
  
  #naep primer has items 
  
  summary(resultWeightReg)
  
  #Examples of variance methods
  # library(EdSurvey)
  # 
  # # All this is just to borrow jackknife replicate weights 
  # 
  # cntl <- readPISA(paste0(edsurveyHome, "PISA/2012"), countries = "USA", verbose=FALSE)
  # 
  # EdSurvey:::getPSUVar(cntl)# "var_unit"
  # EdSurvey:::getStratumVar(cntl)# "wvarstrr"
  # 
  # data <- getData(cntl,c("schoolid","pv1math","st29q03","sc14q02","st04q01",
  #                        "escs","w_fschwt","w_fstuwt", "var_unit", "wvarstrr"), 
  #                 omittedLevels = FALSE,   addAttributes = FALSE)
  # 
  # stuDat <- cbind(stuDat, data) # add weights to stuDat
  # 
  # rep_weight_list <- names(data)[grep("w_fstr",names(data))]
  # 
  # ## this example shows consistent variance estimation
  # # 120 s
  # system.time(rwg1 <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType='regression', ind.vars=c('catVar', 'numVar'), varType="consistent", weightVar='weightVar'))
  # rwg1$SE
  # 
  # ## this example shows robust variance estimation
  # # 120 s
  # system.time(rwg2 <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType='regression', ind.vars=c('catVar', 'numVar'), varType="robust", weightVar='weightVar'))
  # data.frame(consistent=rwg1$SE,
  #            robust=rwg2$SE)
  # 
  # ## this example shows cluster variance estimation
  # # 120 s
  # system.time(rwg3 <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType='regression', ind.vars=c('catVar', 'numVar'), varType="cluster", clusterVar=c("wvarstrr", "var_unit") , weightVar='weightVar'))
  # data.frame(consistent=rwg1$SE,
  #            robust=rwg2$SE,
  #            cluster=rwg3$SE)
  # 
  # ## this example shows taylor series without  weight
  # # 100 s
  # system.time(rwg4 <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType='regression', ind.vars=c('catVar', 'numVar'), varType="Taylor", stratavar="wvarstrr", PSUVar="var_unit", singletonFix="use mean"))
  # data.frame(consistent=rwg1$SE,
  #            robust=rwg2$SE,
  #            cluster=rwg3$SE,
  #            taylor=rwg4$SE)
  # 
  # ## this example shows replicate variance estimation. 
  # # 2500s (41m)
  # system.time(rwg5 <- mml(split(stuItems, stuItems$oppID), stuDat, paramTab, regType='regression', ind.vars=c('catVar', 'numVar'), varType="replicate", weightVar='weightVar', rep_weight=rep_weight_list, jkSumMultiplier=cntl$jkSumMultiplier))
  # data.frame(consistent=rwg1$SE,
  #            robust=rwg2$SE,
  #            cluster=rwg3$SE,
  #            taylor=rwg4$SE,
  #            jk=rwg5$SE)
  # 
  
}



