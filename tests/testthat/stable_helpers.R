library(lsasim)
library(TAM)
library(mice)
library(fastDummies)
library(varhandle)
library(LaplacesDemon)

mice_pv <- function(amp_obj, pv_obj){
  covars <- paste0("q",1:5)
  temp_df <- data.frame(amp_obj,pv_obj$pv[,"PV1.Dim1"]) #this is our temp data with only one PV value, it will help us to determine where to impute
  names(temp_df) <- c(covars, "pv") #we want names to be the same during the entire imputation process
  starter <- mice(temp_df, m=1, maxit = 0, method = "pmm") #this finds us where to impute
  where_to_impute <- starter$where #save where to impute, without this we can't iteratively do imputation
  
  imputed_df_list <- list()
  for(i in 1:5){
    pvi <- c(paste0("PV",i,".Dim1"))
    if(i!=1){ #for the pv1 don't change temp_df, because it already has pv1 and coplete function won't work yet, bec imp didn't start yet 
      temp_df <- cbind(df_imputed, pv_obj$pv[, pvi])
      names(temp_df) <- c(covars, "pv")
    }
    temp_df_imputed <- mice(temp_df, m=1, maxit = 10, method = "pmm", where = where_to_impute) #do the imp with the corresponding pv
    #save the imputed data frame
    df_imputed <- complete(temp_df_imputed,1)[,covars]
    imputed_df_list <- c(imputed_df_list, list(df_imputed))
    names(imputed_df_list)[length(imputed_df_list)] <- paste0("imputed_", pvi)
  }
  return(imputed_df_list)
}



nested_pv_gen <- function(imp_obj){
  nested_pv <- list()
  EAP_res <- matrix(NA, N, 5)
  EAP_sd <-  matrix(NA, N, 5)
  EAP_rel <- vector()
  for(i in 1:5){
    Y_imp <- model.matrix(~ q1+q2+q3+q4+q5, imp_obj[[i]])
    model_imp <- tam.mml.3pl(resp_matrix[,item_names], guess=c,  Y = Y_imp,
                             xsi.fixed = xsi,
                             gammaslope.fixed = gammaslope,
                             ndim = 1, pid = resp_matrix$subject,
                             control=list(snodes=1000, maxiter=5, 
                                          Msteps=50, fac.oldxsi=0.5, 
                                          nodes=seq(-5,5,len=34), 
                                          increment.factor=1.01))
    pv_imp <- tam.pv.mcmc(model_imp, nplausible = 5)
    nested_pv[[i]] <- pv_imp$pv[,2:6]
    EAP_res[,i] <- model_imp$person$EAP
    EAP_sd[,i] <- model_imp$person$SD.EAP
    EAP_rel[i] <-  model_imp$EAP.rel
  }
  return(list(pvs=nested_pv, eaps=EAP_res, sd=EAP_sd, eap_rel = EAP_rel))
}



pooled_mean_se <- function(ym, m=5){ #ym is the mean of each PV
  y.hat <- mean(ym) #eq 1.1
  sigma2.y_imp <- (m+1)/(m*(m-1))*sum((ym-y.hat)^2) #eq 1.2
  return(c(y.hat, sigma2.y_imp^.5))
}



pooled_mean_se_nested <- function(tam.data, k = 5, m = 5) {
  #eq 2.1 
  ym <- colMeans(tam.data)
  
  #eq. 2.2
  y.hat <- mean(ym) 
  
  #eq. 2.4 btw
  a <- rep(NA, m)
  for(i in 1:m){
    a[i] <- (ym[i]-y.hat)^2
  }
  sigma2yb <- (k/(m-1))*sum(a)
  
  #eq. 2.5 within
  b <- rep(NA, m)
  for(i in 1:m){
    b[i] <- sum((tam.data[,i] - ym[i])^2)
  }
  sigma2yw <- sum(b)/(m*(k-1))
  
  #eq. 2.3 
  sigma2yimp <- ((1+1/m)*sigma2yb)/k + (1-1/k)*sigma2yw
  return(c(y.hat,sqrt(sigma2yimp)))
}

rmspe <- function(y, y.pred){
  sqrt(mean((y-y.pred)/y)^2)*100
}

mse <- function(y, y.pred){
  mean((y-y.pred)^2)
}



response_gen2 <- function (subject, item, theta, a_par = NULL, b_par, c_par = NULL, 
                           d_par = NULL, item_no = NULL, ogive = "Logistic") 
{
  if (length(subject) != length(item)) 
    stop("subject and item vectors must be equal length.", 
         call. = FALSE)
  if (is.null(a_par)) 
    a_par <- rep(1, length(unique(item)))
  if (is.null(c_par)) 
    c_par <- rep(0, length(unique(item)))
  if (ogive == "Logistic") 
    DD <- 1
  if (ogive == "Normal") 
    DD <- 1.7
  if (is.null(item_no)) 
    item_no <- seq(length(unique(item)))
  if (is.null(d_par)) 
    b_pars <- split(b_par, seq(length(b_par)))
  if (!is.null(d_par)) {
    d_mat <- do.call("cbind", d_par)
    b_pars <- list()
    for (i in 1:length(b_par)) {
      if (sum(abs(d_mat[i, ])) != 0) {
        b_list <- list()
        for (j in 1:length(d_mat[i, ])) b_list[[j]] <- b_par[i] + 
            d_mat[i, j]
        b_pars[[i]] <- unlist(b_list)
      }
      if (sum(abs(d_mat[i, ])) == 0) 
        b_pars[[i]] <- b_par[i]
    }
  }
  names(b_pars) <- item_no
  
  y <- numeric(length(subject))
  for (n in 1:length(subject)) {
    y[n] <- irt_gen2(theta = theta[subject[n]], a_par = a_par[which(item_no == 
                                                                      item[n])], b_par = b_pars[[which(item_no == item[n])]], 
                     c_par = c_par[which(item_no == item[n])], D = DD)
  }
  df_l <- data.frame(item = item, subject = subject, response = y)
  df_w <- reshape(df_l, timevar = "item", idvar = "subject", 
                  direction = "wide")
  df_item_old <- colnames(df_w)[2:length(df_w)]
  df_item_num <- gsub("[^[:digit:]]", "", df_item_old)
  df_item_new <- ifelse(nchar(df_item_num) == 1, paste0("i00", 
                                                        df_item_num), ifelse(nchar(df_item_num) == 2, paste0("i0", 
                                                                                                             df_item_num), paste0("i", df_item_num)))
  colnames(df_w)[2:length(df_w)] <- df_item_new
  df_w <- df_w[, order(names(df_w))]
  rownames(df_w) <- NULL
  return(y = df_w)
}



irt_gen2 <- function (theta, a_par = 1, b_par, c_par = 0, D = 1)
{
  response_pr <- c_par + (1 - c_par) * (exp(a_par*(theta - b_par)) / (1 + exp(a_par*(theta - b_par))))
  y <- rbinom(1, size = 1, prob = response_pr)
  return(y)
}


