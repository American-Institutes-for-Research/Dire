#' @importFrom tidyr pivot_longer
#' @importFrom stats dnorm
#' @importFrom dplyr all_of
scoreTest <- function(stuDat,
                      dichotParamTab = NULL,
                      polyParamTab = NULL,
                      testScale = NULL,
                      idVar = NULL,
                      Q = 34,
                      minNode = -4,
                      maxNode = 4,
                      polyModel = c('GPCM', 'GRM'),
                      multiCore = FALSE,
                      fast = TRUE,
                      strataVar = NULL,
                      PSUVar = NULL,
                      weightVar = NULL,
                      verbose = 0) {
  call <- match.call()
  if(verbose >= 1) {
    message("Prep data")
  }
  polyModel <- match.arg(polyModel)
  polyModel <- tolower(polyModel)
  if(is.null(idVar)) {
    idVar <- "id"
    stuDat$id <- 1:nrow(stuDat)
  }
  if(missing(weightVar)) {
    paste0("No `weightVar` specified, defaulting to unweighted.")
    weightVar <- "one"
    stuDat$one <- 1
  }
  if(!inherits(stuDat, "data.frame")) {
    stop(paste0("Argument ", dQuote("stuDat"), " must be a data frame."))
  } else {
    if(length(class(stuDat))>1) {
      # not all data.frames are the same, recast
      stuDat <- as.data.frame(stuDat)
    }
  }
  if(!idVar %in% names(stuDat)) {
    stop(paste0(dQuote("idVar"), ", ", dQuote(idVar), ", must be a variable on stuDat."))
  }
  if(length(unique(stuDat[,idVar])) != nrow(stuDat)) {
    stop(paste0(dQuote("idVar"), " must be unique on ", dQuote("stuDat"),"."))
  }
  if(missing(polyParamTab) & missing(dichotParamTab)) {
    stop(paste0("At least one of the arguments ", dQuote("polyParamTab"), " or ", dQuote("dichotParamTab"), " must be defined."))
  }
  if(is.null(polyParamTab) & is.null(dichotParamTab)) {
    stop(paste0("At least one of the arguments ", dQuote("polyParamTab"), " or ", dQuote("dichotParamTab"), " must be defined."))
  }
  dichotParamTab <- cleanDichotParamTab(dichotParamTab) 
  # clean polyParamTab, get the list of d variables on it
  pptd <- cleanPolyParamTab(polyParamTab)
  paramTab <- condenseParamTab(dichotParamTab, pptd$polyParamTab, pptd$dvars)
  items <- paramTab$ItemID
  missingItems <- items[!items %in% colnames(stuDat)]
  items <- items[items %in% colnames(stuDat)]
  stuItems <- stuDat[ , items, drop=FALSE]
  stuItems[ , idVar] <- stuDat[ , idVar, drop=FALSE]
  # wide to long, an item- idVar pair is now a key for a row. Score has the score for that item.
  stuItems <- stuItems %>%
              pivot_longer(
                cols = all_of(items),
                names_to = "key",
                values_to = "score"
              )
  rownames(stuItems) <- NULL
  stuItems <- cleanStuItems(stuItems, stuDat, idVar, paramTab)
  gc()
  # check response ranges
  responseRanges <- generateResponseRanges(paramTab = paramTab, stu = stuItems, polyModel=polyModel, already_warned=FALSE)

  stuDat <- cleanStuDat(stuDat, stuItems, idVar)
  if(!is.null(testScale)) {
    testScale <- checkTestScale(testScale)
  }
  ### These are the exported functions passed to the dopar function
  nodes <- seq(from = minNode, to = maxNode, length.out = Q)

  ### This portion of the code computes all the likelihood evaluations
  ### and does so outside of the function that is maximized
  ### This saves a lot of overhead by using fixed quadrature points
  if(verbose >= 1) {
    message("Calculating likelihood function")
  }
  subtests <- c(unique(paramTab$subtest))
  rr1list <- list()
  for(subtesti in seq_along(subtests)) {
    st <- subtests[subtesti]
    pt <- paramTab[paramTab$subtest %in% st,]
    si <- stuItems[stuItems$key %in% pt$ItemID,]
    # split stuItems into one data frame per person
    rr1 <- buildrr1(pt, si, nodes, idVar, polyModel)
    gc()
    # add names to rr1
    rr1list <- c(rr1list, list(rr1))
    names(rr1list)[length(rr1list)]  <- subtests[subtesti]
  }
  res <- list(stuDat = stuDat,
              idVar = idVar,
              likelihood_list = rr1list,
              paramTab = paramTab,
              testScale = testScale,
              strataVar = strataVar,
              PSUVar = PSUVar,
              nodes = nodes,
              responseRanges = responseRanges,
              weightVar = weightVar,
              polyModel = polyModel,
              fast = fast
              )
  class(res) <- "scoredTest"
  return(res)
}

#'@importFrom data.table rollup data.table
buildrr1 <- function(paramTab, stuItems, nodes, idVar, polyModel) {
  scored_paramTab <- scoreParamTab(paramTab, nodes, polyModel)
  # This will throw an error if the scoring introduced an invalid score
  responseRanges <- generateResponseRanges(paramTab = scored_paramTab, 
                                           stu = stuItems, polyModel=polyModel, already_warned=TRUE)
  if(any(is.na(scored_paramTab$Quadpoint1))) {
    scored_paramTab <- scored_paramTab[is.na(scored_paramTab$Quadpoint1),]
    bad_keys <- unique(scored_paramTab$key)
    stop("There is an issue scoring the following items: ", paste(dQuote(bad_keys), collapse=" "))
  }
  # without all=FALSE unoccupied score points would be included in the final mg
  scoredItems <- stuItems %>%
    left_join(scored_paramTab, by = c("key" = "ItemID", "score" = "score"))
  rm(stuItems)
  # This will throw an error if our merging introduced an invalid score
  responseRanges <- generateResponseRanges(paramTab = paramTab, 
                                           stu = scoredItems, polyModel=polyModel, already_warned=TRUE)
  
  quadCols <- grep("Quadpoint", colnames(scoredItems))
  f2 <- function(scoredItems, idVar, quadCols) {
    ids <- scoredItems[[idVar]]
    cn <- colnames(scoredItems)[quadCols]
    scoreMat <- as.matrix(scoredItems[, cn, drop = FALSE])
    summed <- rowsum(scoreMat, group = ids, na.rm = TRUE)
    pre_lnl <- as.data.frame(summed)
    pre_lnl[[idVar]] <- rownames(summed)
    pre_lnl <- pre_lnl[order(pre_lnl[[idVar]], method="radix"), ]
    
    # Apply exponentiation to the numeric columns and bind ID column at front
    lnl <- cbind(pre_lnl[[idVar]], exp(pre_lnl[, cn, drop = FALSE]))
    colnames(lnl)[1] <- idVar
    lnl
  }
  stuLikelihood <- f2(scoredItems, idVar, quadCols)
  rm(scoredItems)
  if(any(is.na(stuLikelihood$Quadpoint1))) {
    stop("Unknown error aggregating likelihoods.\n")
  }
  # add idVar, and move it to first among all columns
  return(stuLikelihood)
}

# take the paramTab, as it gets built by `Dire:::condenseParamTab`
# at nodes
scoreParamTab <- function(paramTab, nodes, polyModel) {
  res <- paramTab
  # make sure these are in order, d0, d1, ..., d9 (does not select past 9)
  # if they are not in order, this fixes that
  dcols <- sort(colnames(paramTab)[grepl('d[[:digit:]]', colnames(paramTab))])
  ii <- 1
  qp <- paste0("Quadpoint",1:length(nodes))
  res$score <- 0
  for(i in 1:length(nodes)) {
    res[ , paste0("Quadpoint",i)] <- 0
  }
  for(i in 1:nrow(paramTab)) {
    dParams <- paramTab[i,]
    if(dParams$scorePoints == 1) {
      # dichotimous item
      res[ii,] <- c(dParams, 0, mi0 <- multItems(0, dParams$guessing, dParams$D, dParams$slope, nodes, dParams$difficulty))
      ii <- ii + 1
      res[ii,] <- c(dParams, 1, mi1 <- multItems(1, dParams$guessing, dParams$D, dParams$slope, nodes, dParams$difficulty))
      ii <- ii + 1
      if(!is.na(dParams$missingCode)) {
        if(dParams$missingValue == "c") {
          # note that this could be achieved with the following call as well
          #val <- multItems(dParams$guessing, dParams$guessing, dParams$D, dParams$slope, nodes, dParams$difficulty)
          val <- dParams$guessing*mi1 + (1-dParams$guessing)*mi0
        } else {
          mv <- as.numeric(dParams$missingValue) # if some are "c" this will be a character
          # note that this could be achieved with the following call as well
          #val <- multItems(mv, dParams$guessing, dParams$D, dParams$slope, nodes, dParams$difficulty)
          val <- mv*mi1 + (1-mv)*mi0
        }
        res[ii,] <- c(dParams, dParams$missingCode, val)
        ii <- ii + 1
      }
    } else {
      # polytomous item
      # BW: this sets range to e.g. 0,1,2,3,4
      range <- seq(0, dParams$scorePoints)
      if (polyModel == 'gpcm') {
        # BW: what's the reason for this again? range is now 0,2,2,3,4,
        # and this is what will get passed into ansItems
        #range[2] <- range[2] + 1
      } else {
        stop("Cannot score non-GPCM. Ask Dire developers for help.")
      }
      aa <- dParams$slope
      pre_dv <- dParams[dcols]
      dv <- pre_dv[!is.na(pre_dv)]
      D <- dParams$D
      # if fast we use fastApply.cpp 
      Da <- D * aa
      gpcm2 <- function (theta, d, score, Da) {
        if(is.na(score)){
          return(NA)
        } 
        if (score >= length(d)) {
          stop (paste0("Score of ", score," higher than maximum (", length(d),")"))
        }
        tmp <- exp(sum(Da * (theta - d[1:(score+1)]))) / sum(exp(cumsum(Da * (theta - d))))
        return(sum(log(tmp)))
      }
      for(j in seq_along(range)) {
        rg <- range[j]
        Q <- length(nodes)
        v1 <- sapply(1:Q, function(k) gpcm2(d = dv, Da = Da, theta = nodes[k], score = rg))
        res[ii,] <- c(dParams, rg, v1)
        ii <- ii + 1
      }
      if(!is.na(dParams$missingCode)) {
        rg <- range[1]
        Q <- length(nodes)
        v1 <- sapply(1:Q, function(k) gpcm2(d = dv, Da = Da, theta = nodes[k], score = rg))
        res[ii,] <- c(dParams, dParams$missingCode, v1)
        ii <- ii + 1
      }
    }
  }
  return(res)
}

# helper that checks if score points are occupied
generateResponseRanges <- function(paramTab, stu, polyModel, already_warned=FALSE) {
  agg <- data.frame(key = unique(stu$key))
  for(i in 1:nrow(agg)) {
    agg$scorePoints[i] <- subset(paramTab, ItemID == agg$key[i], "scorePoints")[[1]][1]
    vals <- stu$score[stu$key == agg$key[i]]
    # drop missing code
    if(agg$scorePoints[i] == 1) {
      expectedOccupied <- c(0:1)
    } else {
      if(polyModel == "gpcm") {
        expectedOccupied <- seq(0, agg$scorePoints[i], by=1)
      } else {
        expectedOccupied <- seq(1, agg$scorePoints[i]+1, by=1)
      }
    }
    expectedOccupied_ <- expectedOccupied
    if(!is.na(paramTab$missingCode[i])) {
      expectedOccupied <- c(expectedOccupied, paramTab$missingCode[i])
    }
    agg$expectedOccupied[i] <- paste( expectedOccupied, collapse=":")
    agg$occupied[i] <- paste(names(table(vals)), collapse=":")
    # U02713 is a green checkmark, this indicates that every score point is occupied
    # "!" inicates that some are unoccupied
    # NOTE uses expectedOccupied_ (without missing code) so as to apply a green checkmark
    # even if a missing code is not occupied.
    agg$Result[i] <- ifelse( any(!expectedOccupied_ %in% names(table(vals)) ), "!", "\U02713")
    # U274C (a red X) overrules other states and indicates a level is unoccupied
    agg$Result[i] <- ifelse( any(!names(table(vals)) %in% expectedOccupied), "\U274C", agg$Result[i])
  }
  rownames(agg) <- paste0(agg$Result, " " , agg$key)
  agg$score <- NULL
  agg$key <- NULL
  if( any(agg$Result %in% "\U274C") ) {
    agg2 <- subset(agg, Result %in% "\U274C")
    agg2$Result <- NULL
    print(agg2)
    stop("Some items score points inconsistent with expectations.")
  }
  if( any(agg$Result %in% "!") ) {
    agg2 <- subset(agg, Result %in% "!")
    agg2$Result <- NULL
    if(!already_warned){
      warning("Some items score points inconsistent with expectations.")
      print(agg2)
      assign("already_warned", TRUE, envir = parent.frame())
    }
  }
  agg$Result <- NULL
  return(ranges=agg)
}
