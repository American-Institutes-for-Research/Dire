# takes data and a formula and returns data with new columns on it and
# a new formula that nearly spands the same space as the original vectors
# this renders the coefficients very difficult to interpret but may ease fitting
# by (often drastically) reducing the number of columns used to fit the data
# while also using orthagonal columsn which may make the fitting easier for
# the optimizer
#
# parameter data      a data frame used to form the design matrix 
# parameter formula   the formula for the regression
# retainedInformation a number larger than 0 and no more than 1 that
#                     is the proportion of the information to maintain
# returns a list with elements:
# data: the updated data with columsn svd1, svd2, ..., svdi
# i: the value of the maximum SVD variable
# V: the first i columns of the V matrix from the SVD
# diagS: the first i diagonal elements of the Sigma matrix from the SVD
# formula: the updated formula
reduceInformation <- function(data, formula, retainedInformation=0.9, verbose=TRUE) { 
  if(retainedInformation > 1) {
    stop("retainedInformation must not be more than 1")
  }
  if(retainedInformation <= 0) {
    stop("retainedInformation must be larger than 0")
  }
  trms0 <- terms(formula, data=data)
  trms <- delete.response(trms0)
  X <- model.matrix(trms, data=data)
  if(!attributes(trms)$intercept %in% 1) {
    skip_columns <- c()
  } else {
    skip_columns <- 1
  }
  # drop the intercept
  X <- X[ , -1, drop=FALSE]
  svdX <- svd(X)
  # clean apparenlty zero values to allow the sum correctly
  dvec <- svdX$d
  dvec[dvec <= dvec[1] * .Machine$double.eps] <- 0
  i <- 1
  data[ , paste0("PC", i)] <- svdX$u[ , i, drop=FALSE]
  # reverse to get the most stable sum
  while( sum(rev(dvec[1:i]))/sum(rev(dvec)) < retainedInformation) {
    i <- i + 1
    data[ , paste0("PC", i)] <- svdX$u[ , i, drop=FALSE]
  }
  outcome <- ""
  if(attributes(trms0)$response != 0) {
    vars <- attributes(trms0)$variables
    mode(vars) <- "list"
    outcome <- vars[[attributes(trms0)$response + 1]]
  }
  new_formula <- formula(paste0(outcome, " ~ 1 + ", paste(paste0("PC", 1:i), collapse=" + ")))
  V <- svdX$v[ , 1:i, drop=FALSE]
  diagS <- svdX$d[1:i]
  if(verbose) {
    message(paste0("maintaining ", format(100*sum(rev(dvec[1:i]))/sum(rev(dvec)), digits=4), "% of the variance"))
  }
  # to prooject this to a new subset: Uprime <- newX %*% svdx$v %*% diag(1/svdx$d)
  return(list(data=data, i=i, V=V, diagS=diagS, formula=new_formula, skip_columns = skip_columns))
}

