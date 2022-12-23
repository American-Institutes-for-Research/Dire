#' @importFrom foreach %dopar% foreach
#' @importFrom iterators iter
calcRR1_dopar <- function(stu, Q, model, paramTab, nodes, fast) {
  itx1 <- iter(stu)
  dcols <- colnames(paramTab)[grepl('d[[:digit:]]', colnames(paramTab))]
  # check if missing present
  doMissing <- all(c("missingCode", "missingValue") %in% colnames(paramTab))
  rr1 <- foreach(dopari = itx1, .export = c('ldbinom3', "gpcm", "grm", "pasteItems")) %dopar% {
    x <- dopari$score
    if(length(x) == 0 || all(is.na(x))) {
      return(rep(1, Q))
    }
    # subset to items this student was shown
    paramTabi <- paramTab[paramTab$ItemID %in% dopari$key,] # rename so object isn't confused wih paramTab argument
    ind.dichot <- which(paramTabi$scorePoints == 1)
    dParams <- paramTabi[ind.dichot, ]
    # get parameter values from items
    slope <- dParams$slope
    difficulty <- dParams$difficulty
    guessing <- dParams$guessing
    D <- dParams$D
    x1 <- x[ind.dichot]
    if(doMissing) {
      missingVal <- dParams$missingValue
      missingVal[missingVal == "c"] <- dParams$guessing
      missingVal <- as.numeric(missingVal)
      missingCode <- dParams$missingCode
      x1 <- ifelse(x1 == missingCode, missingVal, x1)
    }
    # fix multinomial items
    pVec <- !1:length(x) %in% ind.dichot
    if (model == 'gpcm') {
      x2 <- x[pVec] + 1
    } else {
      x2 <- x[pVec]
    }
    # get multiple choice part
    log.mc.part <- if (length(x1) > 0) {
      if(fast){
        multItems(x1, guessing, D, slope, nodes, difficulty)
      } else {
        sapply(1:Q, function(j) sum(ldbinom3(x1, guessing + (1 - guessing) / (1 + exp(-D * slope * (nodes[j] - difficulty)))), na.rm=TRUE))
      }
    } else {
      0
    }
    # get polytimous part
    log.poly.part <- if (length(x2) > 0) {
      pParams <- paramTabi[pVec,]
      aa <- pParams$slope
      d <- list()
      for(i in 1:nrow(pParams)) {
        pre_dv <- pParams[i, dcols]
        dv <- pre_dv[!is.na(pre_dv)]
        d <- c(d, list(dv))
        if(x2[i] > length(dv)) {
          stop(paste0("score of ", dQuote(x2[i] - (model=='gpcm')), " higher than expected maximum of ", dQuote(length(dv) - (model=='gpcm')), ". on item ", dQuote(pParams$ItemID[i]), " person ", dQuote(dopari[1,1])))
        }
      }
      D <- pParams$D
      # if fast we use fastApply.cpp 
      if(fast){
        ansItems(d,aa,nodes,x2, D)
      } else {
        # otherwise R solution 
        sapply(1:Q, function(j) sum(log(mapply(get(tolower(model)), d = d, a = aa, theta = nodes[j], score = x2, D=D)), na.rm=TRUE))
      }
    } else {
      0
    }
    # summping with apply allows for NAs to be removed, + does not
    log.parts <- rbind(log.mc.part, log.poly.part)
    return(exp(apply(log.parts, 2, sum, na.rm=TRUE)))
  }
  return(do.call(cbind, rr1))
}

calcRR1 <- function(stu, Q, model, paramTab, nodes, fast) {
  dcols <- colnames(paramTab)[grepl('d[[:digit:]]', colnames(paramTab))]
  # check if missing present
  doMissing <- all(c("missingCode", "missingValue") %in% colnames(paramTab))
  rr1 <- lapply(stu, function(dopari) {
    x <- dopari$score
    if(length(x) == 0 || all(is.na(x))) {
      return(rep(1, Q))
    }
    # subset to items this student was shown
    if(any(duplicated(dopari$key))) {
      stop(paste0("person ", dQuote(dopari[1,1]), " has duplicated items: ", paste(dopari$key[duplicated(dopari$key)], collapse=", ")))
    }
    paramTab <- paramTab[paramTab$ItemID %in% dopari$key,]
    ind.dichot <- which(paramTab$scorePoints == 1)
    dParams <- paramTab[ind.dichot,]
    # get parameter values from items
    slope <- dParams$slope
    difficulty <- dParams$difficulty
    guessing <- dParams$guessing
    D <- dParams$D
    x1 <- x[ind.dichot]
    if(doMissing) {
      missingVal <- dParams$missingValue
      missingVal[missingVal == "c"] <- dParams$guessing
      missingVal <- as.numeric(missingVal)
      missingCode <- dParams$missingCode
      x1 <- ifelse(x1 == missingCode, missingVal, x1)
    }
    # fix multinomial items
    pVec <- !1:length(x) %in% ind.dichot
    if (model == 'gpcm') {
      x2 <- x[pVec] + 1
    } else {
      x2 <- x[pVec]
    }
    # get multiple choice part
    log.mc.part <- if (length(x1) > 0) {
      if(fast){
        multItems(x1, guessing, D, slope, nodes, difficulty)
      } else {
        sapply(1:Q, function(j) sum(ldbinom3(x1, guessing + (1 - guessing) / (1 + exp(-D * slope * (nodes[j] - difficulty)))), na.rm=TRUE))
      }
    } else {
      0
    }
    # get polytimous part
    log.poly.part <- if (length(x2) > 0) {
      pParams <- paramTab[pVec,]
      aa <- pParams$slope
      d <- list()
      for(i in 1:nrow(pParams)) {
        pre_dv <- pParams[i, dcols]
        dv <- pre_dv[!is.na(pre_dv)]
        d <- c(d, list(dv))
        if(x2[i] > length(dv)) {
          stop(paste0("score of ", dQuote(x2[i] - (model=='gpcm')), " higher than expected maximum of ", dQuote(length(dv) - (model=='gpcm')), ". on item ", dQuote(pParams$ItemID[i]), " person ", dQuote(dopari[1,1])))
        }
      }
      D <- pParams$D
      # if fast we use fastApply.cpp 
      if(fast){
        ansItems(d,aa,nodes,x2, D)
      } else {
        # otherwise R solution 
        sapply(1:Q, function(j) sum(log(mapply(get(tolower(model)), d = d, a = aa, theta = nodes[j], score = x2, D=D)), na.rm=TRUE))
      }
    } else {
      0
    }
    # summping with apply allows for NAs to be removed, + does not
    log.parts <- rbind(log.mc.part, log.poly.part)
    return(exp(apply(log.parts, 2, sum, na.rm=TRUE)))
  })
  return(unname(do.call(cbind, rr1)))
}

