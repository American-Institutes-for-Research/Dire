# stuItems is a data set with student responses on it.
scr <- function(stuItems, Q, paramTab, polyModel, idVar) {
  ids <- unique(stuItems[[idVar]])
  N <- length(ids)

  if(any(!c(idVar, "key", "score") %in% colnames(stuItems))) {
    stop(paste0("Argument ", dQuote("stuItems"), " must have variables ", pasteItems(c(idVar, "key", "score")), "."))
  }
  
  stuItems[[idVar]] <- as.character(stuItems[[idVar]])
  stuItems$key <- as.character(stuItems$key)
  stuItems <- stuItems[order(stuItems[[idVar]]),]

  stuItems <- stuItems[,c(idVar, "key", "score")]
  stuItems <- stuItems[!is.na(stuItems$score),]
  stuItems <- stuItems[stuItems$key %in% paramTab$ItemID,]
  stu <- stuItems[order(stuItems$key),]
  stu <- split(stu, stu[[idVar]])
  paramTab <- paramTab[order(paramTab$ItemID),]

  
  stu <- lapply(stuItems, function(x) {
    # subset to items on paramTab after subsetting to items on this test or subtest
    res <- as.data.frame(x)[c(idVar, "key", "score")]
    res <- res[!is.na(res$score),]
    res <- res[res$key %in% paramTab$ItemID,]
    # sort by key/ItemID
    res <- res[order(as.character(res$key)),]
    return(res)
  })


    rr1 <- calcRR1_dopar(stu, Q, polyModel, paramTab, nodes)
  calcRR1(stu, Q, polyModel, paramTab, nodes)

}
