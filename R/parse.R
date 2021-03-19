#' Export arguments from an mml call to generate an AM dct file
#' @param dichotParamTab see \code{\link{mml}}
#' @param polyParamTab see \code{\link{mml}}
#' @param testScale see \code{\link{mml}}
#' @param weightVar see \code{\link{mml}}
#' @param stratumVar see \code{\link{mml}}
#' @param PSUVar see \code{\link{mml}}
#' @param testName the name of the file to save, dct is appended
#' @param saveFilePath the directory to output the dct file
#' 
#' Only works with \code{.sav} file from \code{\link{saveSav}}.
#'
#' @seealso \link{saveSav}, \link{mml}
#' 
#' @importFrom stats reshape
#' @export
paramTab2AMdct <- function(dichotParamTab,
                           polyParamTab,
                           testScale,
                           weightVar="origwt",
                           stratumVar="repgrp1",
                           PSUVar="jkunit",
                           testName="IRTparams",
                           saveFilePath=getwd()) {
  # make testScale have expected columns
  if (nrow(testScale) == 0) {
    stop('Please provide testScale.')
  }
  if ( (is.null(dichotParamTab) || nrow(dichotParamTab)==0) &
       (is.null(polyParamTab) || nrow(polyParamTab)==0)) {
    stop('Please provide one or both paramTabs.')
  }
  # fix NULL to zero length to simplify future tests
  if(is.null(dichotParamTab)) {
    dichotParamTab <- polyParamTab[NULL,]
  }
  # fix NULL to zero length to simplify future tests
  if(is.null(polyParamTab)) {
    polyParamTab <- dichotParamTab[NULL,]
  }
  if(!"test" %in% colnames(testScale)) {
    testScale$test <- "test"
  }
  if(!"subtest" %in% colnames(testScale)) {
    testScale$subtest <- paste0("subtest",1:nrow(testScale))
  }
  if (nrow(dichotParamTab) > 0) {
    if(!"missingCode" %in% colnames(dichotParamTab)) {
      dichotParamTab$missingCode <- 99
    }
    if(!"D" %in% colnames(dichotParamTab)) {
      dichotParamTab$D <- 1
    }
    if(!"test" %in% colnames(dichotParamTab)) {
      dichotParamTab$test <- "test"
    }
    if(!"subtest" %in% colnames(dichotParamTab)) {
      dichotParamTab$subtest <- paste0("subtest1")
    }
    dichotParamTab$testi <- ""
    dichotParamTab$subtesti <- ""
  }
  if (nrow(polyParamTab) > 0) {
    if(!"missingCode" %in% colnames(polyParamTab)) {
      polyParamTab$missingCode <- 99
    }
    if(!"D" %in% colnames(polyParamTab)) {
      polyParamTab$D <- 1
    }
    if(!"test" %in% colnames(polyParamTab)) {
      polyParamTab$test <- "test"
    }
    if(!"subtest" %in% colnames(polyParamTab)) {
      polyParamTab$subtest <- paste0("subtest1")
    }
    polyParamTab$testi <- ""
    polyParamTab$subtesti <- ""
  }
  
  if(!"subtestWeight" %in% colnames(testScale)) {
    testScale$subtestWeight <- 1
  }

  # create file path if it doesn't exist
  saveFilePath <- normalizePath(saveFilePath, winslash="/")
  fi <- file.info(saveFilePath)
  if(is.na(fi$isdir)) {
    dir.create(saveFilePath, recursive=TRUE)
  }
  sink(paste0(saveFilePath,'/',testName,'.dct'))
  # write top block, using testScale
  utest <- unique(testScale$test)[1]
  for(testi in 1:length(utest)) {
    test <- utest[testi]
    tdi <- testScale[testScale$test == test,]
    cat(paste0("test = \"", test, "\" id=",testi," scale=1 location=0\n") )
    if (nrow(dichotParamTab) > 0) {
      dichotParamTab$testi <- ifelse(dichotParamTab$test == test, testi, dichotParamTab$testi)
    }
    if (nrow(polyParamTab) > 0) {
      polyParamTab$testi <- ifelse(polyParamTab$test == test, testi, polyParamTab$testi)
    }
    usub <- tdi$subtest[!is.na(tdi$subtest)]
    for(subtesti in 1:length(usub)) {
      subtest <- usub[subtesti]
      cat(paste0("  subtest = \"", subtest, "\" id=",subtesti," scale=", tdi$scale[subtesti], " location=", tdi$location[subtesti]," weight=", tdi$subtestWeight[subtesti],"\n") )
      if (nrow(dichotParamTab) > 0) {
        dichotParamTab$subtesti <- ifelse(dichotParamTab$test == test & dichotParamTab$subtest == subtest, subtesti, dichotParamTab$subtesti)
      }
      if (nrow(polyParamTab) > 0) {
        polyParamTab$subtesti <- ifelse(polyParamTab$test == test & polyParamTab$subtest == subtest, subtesti, polyParamTab$subtesti)
      }
    }
  }
  # write variables
  cat("\nVariables\n\n")
  
  checkNA <- function(data, colnames) {
    for(colname in colnames)
      if(any(is.na(data[,colname]))) {
        sink()
        stop("NA value in column ", dQuote(colname), ".")
      }
    }
  
  if (nrow(dichotParamTab) > 0) {
  dpt <- dichotParamTab
  checkNA(dpt, c("slope", "difficulty", "guessing", "D"))
  
  for(i in 1:nrow(dichotParamTab)) {
    cat(paste0(dpt$ItemID[i], " IRM=3pl ipa=",dpt$slope[i],
               " ipb=", dpt$difficulty[i], " ipc=", dpt$guessing[i],
               " scale=", dpt$D[i], " ontest=", dpt$testi[i],
               " onsubtest=", dpt$subtesti[i]))
    if(is.na(dpt$missingCode[i])) {
      cat("\n")
    } else {
      cat(paste0(" omitted=", dpt$missingCode[i],"\n"))
    }
   }
  }
  
  if (nrow(polyParamTab) > 0) {
    ppt <- polyParamTab
    checkNA(ppt, c("slope", "itemLocation", "D"))
    for(i in 1:nrow(polyParamTab)) {
      cat(paste0(ppt$ItemID[i], " ontest=", ppt$testi[i],
                 " onsubtest=", ppt$subtesti[i],
                 " IRM=PCL ipa=",ppt$slope[i],
                 " ipb=", ppt$itemLocation[i]))
      j <- 1
      dvar <- paste0("d",j)
      while(dvar %in% colnames(ppt) && !is.na(ppt[i,dvar])) {
        cat(paste0(" ipd",j,"=", ppt[i,dvar], " "))
        j <- j + 1
        dvar <- paste0("d",j)
      }
      cat(paste0(" scale=", ppt$D[i]))
      if(is.na(dpt$missingCode[i])) {
        cat("\n")
      } else {
        cat(paste0(" omitted=", dpt$missingCode[i],"\n"))
      }
    }
  }

  cat("\n")
  cat(paste0(weightVar," type=R design=W\n"))
  cat(paste0(stratumVar," type=R design=S\n"))
  cat(paste0(PSUVar," type=R design=C\n"))
  sink()
}

#' Save a scored file as sav to use in AM
#' @param stuItems see \code{\link{mml}}
#' @param stuDat a data frame 
#' @param testName a character; the name of the save file path, \code{.sav} is appended to the file
#' @param saveFilePath a character; the directory to save in
#' @param idVar a character of the column name of the student identifier on \code{stuItems} and \code{stuDat}
#' @param formula a formula; when present performs all variable transformations R typically performs on a
#'                formula before exporting. This is useful because AM does not dummy code variables for you.
#' @param na.action see \code{na.action}
#'
#' Saves a \code{.sav} formatted file that can be read in to AM. This speeds testing Dire against known software.
#'
#' @author Sun-joo Lee, Paul Bailey
#' @importFrom haven write_sav
#' @export
saveSav <- function(stuItems, stuDat, testName, saveFilePath, idVar="sid", formula=NULL, na.action=options('na.action')$na.action) {
  # test inputs
  if(!inherits(stuItems, "data.frame")) stop(paste0("The argument ", dQuote("stuItems"), " must be a ", dQuote("data.frame"), " object"))
  if(!inherits(stuDat, "data.frame")) stop(paste0("The argument ", dQuote("stuDat"), " must be a ", dQuote("data.frame"), " object"))
  if(!inherits(testName, "character")) stop(paste0("The argument ", dQuote("testName"), " must be a length 1 character."))
  if(length(saveFilePath) > 1) stop(paste0("The argument ", dQuote("saveFilePath"), " must be a length 1 character."))
  if(!inherits(saveFilePath, "character")) stop(paste0("The argument ", dQuote("saveFilePath"), " must be a length 1 character."))
  if(length(idVar) > 1) stop(paste0("The argument ", dQuote("idVar"), " must be a length 1 character."))
  if(!inherits(idVar, "character")) stop(paste0("The argument ", dQuote("idVar"), " must be a length 1 character."))
  if(!idVar %in% colnames(stuDat)) stop(paste0("The argument ", dQuote("idVar"), " must be the name of a column in ", dQuote("stuDat")))
  if(!idVar %in% colnames(stuItems)) stop(paste0("The argument ", dQuote("idVar"), " must be the name of a column in ", dQuote("stuItems")))
  if(!inherits(formula, "formula")) stop(paste0("The argument ", dQuote("formula"), " must be a ", dQuote("formula"), " object"))
  if(length(testName) > 1) stop(paste0("The argument ", dQuote("testName"), " must be a length 1 character."))
  if(!inherits(get(na.action), "function")) stop(paste0("The argument ", dQuote("na.action"), " must be the name of a valid na.action function. See ?na.action"))
  # if the user supplies a formula, generate all variables for that
  if(!is.null(formula)) {
    # na.action before we start, save to reset later
    na0 <-  options('na.action')$na.action
    # update to pass NA values
    options(na.action='na.pass')
    # get all the column we will create from the formula
    modelMatrix <- model.matrix(formula, data=stuDat)
    # add these values
    stuDat <- cbind(stuDat, modelMatrix)
    # reset na.action
    options(na.action=na0)
  }
  # perform na.action
  stuDat <- do.call(na.action, list(object=stuDat))
  stuItemsWideAgain <- reshape(data=stuItems, idvar=idVar, timevar='key', direction='wide')
  # remove items from stuDat
  stuDat <- stuDat[ , c(idVar, colnames(stuDat)[!colnames(stuDat) %in% c(colnames(stuItemsWideAgain), unique(stuItems$key))])]
  # merge on items
  scored_df <- merge(stuDat, stuItemsWideAgain, by=idVar, all.x=TRUE) # merge back by student id
  scored <- scored_df[ , colnames(scored_df)[!colnames(scored_df) %in% 'sid']]  
  colnames(scored) <- gsub('score.', '', colnames(scored))
  # create file path if it doesn't exist
  if(!dir.exists(saveFilePath)) {
    dir.create(saveFilePath)
    print(paste0(saveFilePath, ' has been created.'))
  }

  names(scored) <- gsub("-|\\.|\\/|'|\\[|\\]", "", names(scored))
  names(scored) <- gsub(" ",  "_", names(scored))
  names(scored) <- gsub("[(]|[)]",  "", names(scored))
browser()
  write_sav(scored, paste0(saveFilePath,'/', testName, '.sav'))
}

#' Format AM dct File for Use with Dire
#' @description
#' Takes an \code{AM dct} file and formats it for use with the \code{mml} method
#' as \code{paramTab}.
#' 
#' @param dct a file location from which to read the \code{dct} file
#' @param mml a logical for if the paramTab is being used in \code{mml.sdf}
#' 
#' @return a \code{data.frame} in a format suitable for use with \code{mml} as
#' a \code{paramTab}.
#' 
#' @author Sun-Joo Lee 
#' @export
parseNAEPdct <- function(dct, mml=TRUE) {
  # build 
  dichotParamTab <- data.frame()
  polyParamTab <- data.frame()
  testDat <- data.frame()
  
  # correct factors
  subjects<-c('Mathematics', 'Reading', 'Science', 'Geography', 'History','Writing', 'Civics', 'Economics', 'Arts')
  subtests<-c('algebra', 'geometry', 'measurement', 'number', 'data','information', 'literary', 'earth', 'life', 'physical','environment', 'space', 'spatial', 'cultures', 'democracy','technology', 'world', 'task', 'ltt', 'ltt_bridge','international', 'market', 'national', 'vocabulary', 'univariate', 'one', 'overall')
  
  dct_obj <- readLines(dct)
  for (line in dct_obj) {
    line <- trimws(tolower(line))
    line <- gsub(' = ', '=', line)
    line <- gsub(' =', '=', line)
    line <- gsub('= ', '=', line)
    if(grepl('irm=3pl', line)){
      # parse dichot 
      line_split <- unlist(strsplit(line, " "))
      naepid <- line_split[1]
      ipa <- as.numeric(unlist(strsplit(line_split[grepl('ipa', line_split)], "="))[2])
      ipb <- as.numeric(unlist(strsplit(line_split[grepl('ipb', line_split)], "="))[2])
      ipc <- as.numeric(unlist(strsplit(line_split[grepl('ipc', line_split)], "="))[2])
      subtest <- as.numeric(unlist(strsplit(line_split[grepl('onsubtest', line_split)], "="))[2])
      
      params_row <- as.data.frame(
        list(
          ItemID = tolower(naepid),
          subtest = subtest,
          slope = ipa,
          difficulty = ipb,
          guessing = ipc
        )
      )
      dichotParamTab <- rbind(dichotParamTab, params_row)
      
    } else if (grepl('irm=pcl', line)){
      # parse poly
      line_split <- unlist(strsplit(line, " "))
      naepid <- line_split[1]
      ipa <- as.numeric(unlist(strsplit(line_split[grepl('ipa', line_split)], "="))[2])
      ipb <- as.numeric(unlist(strsplit(line_split[grepl('ipb', line_split)], "="))[2])
      
      ds <- c()
      for (dnum in 1:5) {
        d <- as.numeric(unlist(strsplit(line_split[grepl(paste0('ipd',dnum), line_split)], "="))[2])
        if (length(d) == 0) {
          d <- NA
        }
        ds <- c(ds, d)
      }
      subtest <- as.numeric(unlist(strsplit(line_split[grepl('onsubtest', line_split)], "="))[2])
      
      params_row <- as.data.frame(
        list(
          ItemID = naepid,
          subtest = subtest,
          slope = ipa,
          itemLocation = ipb,
          d1 = ds[1],
          d2 = ds[2],
          d3 = ds[3],
          d4 = ds[4],
          d5 = ds[5]
        )
      )
      polyParamTab <- rbind(polyParamTab, params_row)
      
    } else if (grepl('^subtest', line)) {
      # parse testDat
      line_split <- tolower(unlist(strsplit(line, " ")))
      line_split <- gsub('[[:punct:]]+$', '',line_split)
      line_split <- gsub('^[[:punct:]]+', '',line_split)
      subtest_name <- intersect(line_split, subtests)
      
      ### This part will need to be refined to capture all variations of subtest names ###
      if (length(subtest_name)!=1) {
        subtest_name <- unlist(strsplit(line_split[1], '='))[2]
        subtest_name <- tolower(gsub('[[:punct:]]+', '', subtest_name))
        if (subtest_name == 'numbers') {
          subtest_name <- 'number'
        }
      }
      if ((subtest_name == 'space') & (grepl('earth',line))) {
        subtest_name <- 'earth'
      }
      
      if ((subtest_name == 'one') | (subtest_name == 'overall')) {
        subtest_name <- 'univariate'
      } 
      ### Up to here ###
      
      subtest_num <- as.numeric(unlist(strsplit(line_split[grepl('id', line_split)],"="))[2])
      subtest_scale <- as.numeric(unlist(strsplit(line_split[grepl('scale', line_split)],"="))[2])
      subtest_loc <- as.numeric(unlist(strsplit(line_split[grepl('location', line_split)],"="))[2])
      subtest_wgt <- as.numeric(unlist(strsplit(line_split[grepl('weight', line_split)],"="))[2])
      if (identical(subtest_wgt, numeric(0))) {
        subtest_wgt <- 0
      }
      
      subtest_row <- as.data.frame(
        list(
          subtest = subtest_name,
          subtestNumber = subtest_num,
          scale = subtest_scale,
          location = subtest_loc,
          subtestWeight = subtest_wgt
        )
      )
      testDat <- rbind(testDat, subtest_row)
      
    } else if (grepl('^test', line)) {
      line_split <- trimws(gsub('[[:punct:] ]+',' ',(unlist(strsplit(unlist(strsplit(line, " ")), '=')))))
      nums <- grep('[[:digit:]]+', line_split, value = T)
      theyear <- as.numeric(nums[1])
      thelevel <- as.numeric(nums[2])
      if (('age' %in% line_split) & ('grade' %in% line_split)) {
        warning('Unclear whether grade or age of students.')
        theleveltype <- ''
      } else if ('age' %in% line_split) {
        theleveltype <- 'age'
      } else if ('grade' %in% line_split) {
        theleveltype <- 'grade'
      } else {
        theleveltype <- ''
      }
      if (('national' %in% line_split) & ('state' %in% line_split)) {
        warning('Unclear whether National or State test.')
        theregion <- ''
      } else if ('national' %in% line_split) {
        theregion <- 'national'
      } else if ('state' %in% line_split) {
        theregion <- 'state'
      } else {
        theregion <- ''
      }
      thesubject <- intersect(line_split, tolower(subjects))
      if (length(thesubject) == 1) {
        thesubject <- thesubject
      } else {
        for (part in line_split) {
          subject <- grep(part, tolower(subjects), value = T)
          if (length(subject) == 1) {
            thesubject <- subject
            break
          } else {
                thesubject <- ''
              }
            }
          }
        }
      } 
    

  # change subtest from number to label
  dichotParamTab$subtest <- unlist(lapply(dichotParamTab$subtest, function(x) testDat[testDat$subtestNumber==x, 'subtest']))
  polyParamTab$subtest <- unlist(lapply(polyParamTab$subtest, function(x) testDat[testDat$subtestNumber==x, 'subtest']))
  testDat$subtestNumber <- NULL 
  
  # for feeding into mml.sdf
  if (mml) {
    if (nrow(dichotParamTab) > 0) {
      dichotParamTab$D <- 1.7
      dichotParamTab$missingValue <- 'c'
      dichotParamTab$missingCode <- 8
    }
    if (nrow(polyParamTab) > 0) {
      polyParamTab$D <- 1.7
      polyParamTab$missingValue <- 'c'
      polyParamTab$missingCode <- 8
    }
    
  } else {
    # for making NAEPirtparams
    if (nrow(dichotParamTab) > 0) {
      dichotParamTab$level <- thelevel
      dichotParamTab$levelType <- theleveltype
      dichotParamTab$assessmentCode <- paste0(toupper(substr(theregion, 1, 1)), substr(theregion, 2, nchar(theregion)))
      dichotParamTab$subject <- paste0(toupper(substr(thesubject, 1, 1)), substr(thesubject, 2, nchar(thesubject)))
      dichotParamTab$year <- theyear
    }
    
    if (nrow(polyParamTab) > 0) {
      polyParamTab$level <- thelevel
      polyParamTab$levelType <- theleveltype
      polyParamTab$assessmentCode <- paste0(toupper(substr(theregion, 1, 1)), substr(theregion, 2, nchar(theregion)))
      polyParamTab$subject <- paste0(toupper(substr(thesubject, 1, 1)), substr(thesubject, 2, nchar(thesubject)))
      polyParamTab$year <- theyear
    }
    
    if (nrow(testDat) > 0) {
      testDat$level <- thelevel
      testDat$assessmentCode <- paste0(toupper(substr(theregion, 1, 1)), substr(theregion, 2, nchar(theregion)))
      testDat$subject <- paste0(toupper(substr(thesubject, 1, 1)), substr(thesubject, 2, nchar(thesubject)))
      testDat$year <- theyear
    }
  }
  
  # return list
  ret <- list()
  ret$dichotParamTab <- dichotParamTab
  ret$polyParamTab <- polyParamTab
  ret$testDat <- testDat
  
  return(ret)
}

# helper
getPart <- function(le, element) {
  res <- lapply(le, function(e) {
    if(tolower(e[1]) == tolower(element)) {
      return(e[2])
    }
    return(NULL)
  })
  res <- unlist(res)
  if(is.null(res)) {
    return(NA)
  }
  suppressWarnings(rr <- as.numeric(res))
  if(is.na(rr)) {
    return(res)
  }
  return(rr)
}

# helper
csv2paramList <- function(pModel = c('GPCM', 'GRM', 'PCM'), dat, subsetItems = NULL){
  # If converting modelParameters tab from config, only grab
  # relevant variables and rename ItemID to ITEM_ID
  if (all(c('MODEL', 'ItemID', 'ScorePoints', paste0('P', 0:10)) %in% names(dat))) {
    dat <- dat[, c('ItemID', 'MODEL', 'ScorePoints', paste0('P', 0:10))]
    names(dat)[1] <- 'ITEM_ID'
  }
  
  params2use <- with(dat, dat[order(dat$ITEM_ID), ])

  if(!is.null(subsetItems)){
    dat <- dat[which(dat$ITEM_ID %in% subsetItems), ]
    params2use <- dat[match(subsetItems, dat$ITEM_ID), ] #sort according to item order input
  }
  
  pModel <- toupper(pModel)
  
  names(params2use) <- tolower(names(params2use))
  finalItems <- params2use$item_id
  
  # any polytymous items?-- this will take care of IRTGPC sp = 2, IRTGPC sp = 3, IRTPCL sp = 2, IRTPCL sp = 3
  any.poly <- any(params2use$scorepoints > 1)
  # identifly all other items-- IRTGPC sp = 1, IRT3pln sp = 1, IRTPCL sp = 1, IRT3pl sp = 1
  # in config file, IRTPCL scorepoints = 1 items have b parameter stored in p0, unlike other 3pl items
  # for example, parameters for IRT3pl sp = 1 are stored as p0 = a(1), p1 = b, p2 = c(0)
  pos <- which(params2use$scorepoints == 1 & grepl('3pl', tolower(params2use$model))) # 3pl
  pos1 <- which(params2use$scorepoints == 1 & tolower(params2use$model) %in% 'irtgpc') # gpc
  pos <- c(pos, pos1)
  pos2 <- which(params2use$scorepoints == 1 & tolower(params2use$model) %in% 'irtpcl') # pcl
  
  mcParams <- params2use[pos, paste0('p', 0:2)]
  mcParams$p2 <- with(mcParams, ifelse(is.na(p2), 0, p2)) # make c-parameter if it is missing, this is needed for IRTGPC sp = 1
  
  mcParams2 <- params2use[pos2, paste0('p', 0:2)]
  if(nrow(mcParams2) > 0){
    mcParams2$p1 <- mcParams2$p0 # b is stored in p0 for IRTPCL sp = 1 items
    mcParams2$p0 <- 1 # set a to 1
    mcParams2$p2 <- 0 # set c to 0
  }
  
  mcParams <- rbind(mcParams, mcParams2)
  binary <- c(pos, pos2)
  
  ## inserted 4 lines for modifying
  binaryID <- params2use$item_id[binary]  
  polyID <- params2use$item_id[-binary]
  
  if(any.poly){
    ### This is the code to get the CR items organized for the params list
    numSteps <- length( grep("[p][0-9]", colnames(params2use)) )
    if(pModel == 'PCM') steps <- paste('p', 0:(numSteps-1), sep='')
    if(pModel %in% c('GRM', 'GPCM')) steps <- paste('p', 1:(numSteps-1), sep='')
    if(length(binary) > 0){
      crParams <- params2use[-binary, steps]
    } else {
      crParams <- params2use[, steps]
    } 
    if(pModel == 'PCM') aVar <- rep(1, nrow(crParams))
    if(length(binary) > 0){
      if(pModel %in% c('GRM', 'GPCM')) aVar <- params2use[-binary, 'p0']
    } else{
      if(pModel %in% c('GRM', 'GPCM')) aVar <- params2use[, 'p0']
    }
    crParams$base <- 0
    if(pModel == 'PCM') crParams <- crParams[, c('base', steps)]
    if(pModel == 'GPCM') crParams <- crParams[, c('base', steps)]
    if(pModel == 'GRM') crParams <- crParams[, steps]
    crParams <- as.data.frame(t(crParams))
    crList <- as.list(crParams)
    crList <- lapply(crList, function(x) x[!is.na(x)])
    if(pModel == 'PCM') pModel <- 'GPCM' # so that model can be automatically called by irt.ability
    if(length(binary) > 0){ # test has both CM and polytomous items
      params <- list('3pl' = list(a = mcParams$p0, b = mcParams$p1, c = mcParams$p2), 
                     gpcm = list(a = aVar , d = crList, polyID=polyID), model = pModel, items = finalItems, mcID=binaryID, pyID=polyID, mcPos = binary, 
                     score.pts = params2use$scorepoints)
    } else { # test has only polytomous items
      params <- list('3pl' = NULL, 
                     gpcm = list(a = aVar , d = crList), model = pModel, items = finalItems, mcID=NULL, pyID=polyID, mcPos = NULL, score.pts = params2use$scorepoints)
    }
  } else { # test has only MC items
    if(pModel == 'PCM') pModel <- 'GPCM' # so that model can be automatically called by irt.ability
    params <- list('3pl' = list(a = mcParams$p0, b = mcParams$p1, c = mcParams$p2), 
                   gpcm = NULL, items = finalItems, mcID=binaryID, pyID=NULL, mcPos = binary, score.pts = params2use$scorepoints, model = pModel)
  }
  class(params) <- 'param.list'
  params     
}
