#' Marginal Maximum Likelihood Estimation of Linear Models
#' @description
#' Implements a survey-weighted marginal maximum estimation, a type of
#' regression where the outcome is a latent trait (such as student ability).
#' Instead of using an estimate, the likelihood function marginalizes student
#' ability. Includes a variety of variance estimation strategies.
#' 
#' @param formula  a \ifelse{latex}{\code{formula}}{\code{\link[stats]{formula}}}
#'                 object in the style of \ifelse{latex}{\code{lm}}{\code{\link[stats]{lm}}}
#' @param stuItems deprecated. Simply put items in \code{stuDat} now.
#' @param stuDat   a \code{data.frame} with a single row per student. Predictors in
#'                 the \code{formula} must be in \code{stuDat} as well as items.
#' @param dichotParamTab a \code{data.frame} of dichotomous item information, see Details
#' @param polyParamTab a \code{data.frame} of polytomous item information, see Details
#' @param testScale a \code{data.frame} of scaling information, see Details
#' @param Q        an integer; the number of integration points
#' @param minNode a numeric; the smallest integration point for the latent variable
#' @param maxNode a numeric; the largest integration point for the latent variable
#' @param polyModel polytomous response model;
#'                  one of \code{GPCM} for the Graded Partial Credit Model
#'                  or \code{GRM} for the Graded Response Model
#' @param weightVar a variable name on \code{stuDat} that is the full sample weight
#' @param idVar a variable name on \code{stuDat} that is the identifier. Every 
#'              ID from \code{stuDat} must appear on \code{stuItems} and vice versa.
#' @param multiCore a logical indicating whether to use parallel processing
#' @param bobyqaControl deprecated. A list that gets passed to the \code{bobyqa} optimizer in \code{minqa}
#' @param composite a logical indicating if an overall test should be treated as
#'                  a composite score; a composite is a weighted average of the
#'                  subscales in it.
#' @param strataVar character naming a variable on \code{stuDat}, the variable indicating the
#'                  stratum for each row. Used in post-hoc robust variance estimation.
#' @param PSUVar character naming a variable on \code{stuDat}; the primary sampling unit
#'               (PSU) variable. Used in post-hoc robust variance estimation. The values
#'               do not need to be unique across strata.
#' @param fast deprecated. Always TRUE now.
#' @param calcCor deprecated. Always TRUE now.
#' @param verbose integer, negative or zero for no details, increasingly verbose messages at one and two
#' @param retainedInformation set to a value of 1 to fit the model as is typically fit. If the value is less than one, a principal component analysis is performed and columns associated with principal components totaling at least \code{retainedInformation} are retained while all other data is discarded. The estimated latent regression coefficients will accordingly be named as PC1, ..., PCN.
#' @param optimizer character naming the optimization method to use;
#'                  one of \code{EM} for Expectation Maximization
#'                  or \code{QN} for Quasi-Newton
#' 
#' @details
#' 
#' The \code{mml} function models a latent outcome conditioning on
#' item response data, covariate data, item parameter information,
#' and scaling information.
#' These four parts are broken up into at least one argument each.
#' Student item response data go into \code{stuItems}; whereas student
#' covariates, weights, and sampling information go into \code{stuDat}.
#' The \code{dichotParamTab} and \code{polyParamTab}
#' contain item parameter information for dichotomous and polytomous items,
#' respectively---the item parameter data is the result of an existing
#' item parameter scaling. In the case of 
#' the National Assessment of Educational Progress (NAEP),
#' they can be found online, for example, at
#' \href{https://nces.ed.gov/nationsreportcard/tdw/analysis/scaling_irt.aspx}{NAEP technical documentation}.
#' Finally, information about scaling and subscale weights for composites are put in \code{testScale}.
#' 
#' The model for dichotomous responses data is, by default, three Parameter Logit
#' (3PL), unless the item parameter information provided by users suggests
#' otherwise. For example, if the scaling used a two Parameter Logit (2PL) model,
#' then the guessing parameter can simply be set to zero. For polytomous
#' responses data, the model is dictated by the \code{polyModel} argument.
#'
#' The \code{dichotParamTab} argument is a \code{data.frame} with a column named
#' \code{ItemID} that identifies the items and agrees with
#' the \code{key} column in the \code{stuItems} argument,
#' and, for  a 3PL item, columns \code{slope},
#' \code{difficulty}, and \code{guessing} for the \dQuote{a}, \dQuote{d}, and
#' \dQuote{g} parameters, respectively; see the vignette for details of
#' the 3PL model. Users can also use the column names directly from the
#' vignette notation (\dQuote{a}, \dQuote{d}, and \dQuote{g}) if they prefer.
#' Items that are missing (\code{NA}) are not used in the likelihood function. 
#' Users wishing to apply a special behavior for a subset of items can use
#' set those items to an invalid score and put that in the \code{dichotParamTab}
#' column \code{missingCode}. They are then scored as if they are \code{missingValue}
#' proportion correct. To use the guessing parameter for the proportion correct
#' set \code{missingValue} to \dQuote{c}.
#' 
#' The \code{polyParamTab} has columns \code{ItemID} that must match with the
#' \code{key} from \code{stuItems}, as well as \code{slope}
#' (which can also be called \code{a}) that corresponds to the \code{a}
#' parameter in the vignette.
#' Users must also specify the location of the cut points (\eqn{d_{cj}} in the vignette)
#' which are named \code{d1}, \code{d2}, ..., up to \code{dn} where \code{n} is
#' one less than the number of score points. Some people prefer to also apply a 
#' shift to all of these and this shift is applied when there is a column named
#' \code{itemLocation} by simply adding that to every \code{d*} column. Items
#' are not included in the likelihood for an individual when their value on \code{stuItems}
#' is \code{NA}, but no provision is made for guessing, nor special provision for 
#' missing codes in polytomous items.
#' 
#' For both \code{dichotParamTab} and \code{polyParamTab} users wishing
#' to use a \code{D} paramter of 1.7 (or any other value) may specify that, per item,
#' in a column named \code{D}. 
#'
#' When there are multiple constructs, subscales, or the user wants a composite
#' score, additional, optional, columns \code{test} and \code{subtest} can be used. 
#' these columns can be numeric or text, they must agree with the same 
#' columns in \code{testScale} to scale the results. 
#' 
#' Student data are broken up into two parts. The item response data goes
#' into \code{stuItems}, and the student covariates for the formula go into
#' \code{stuDat}. Information about items, such as item difficulties, is in 
#' \code{paramTab}. All dichotomous items are assumed to be 
#' 3PL, though by setting the guessing parameter to zero, the user
#' can use a 2PL or the one Parameter Logit (1PL) or Rasch models.
#' The model for polytomous responses data is dictated by the \code{polyModel}
#' argument.
#' 
#' The marginal maximum likelihood then integrates the product of the student
#' ability from the assessment data, and the estimate from the linear model
#' estimates each student's ability based on the \code{formula} provided
#' and a residual standard error term. This integration happens from the
#' minimum node to the maximum node in the \code{control} argument (described
#' later in this section) with \code{Q} quadrature points. 
#' 
#' The \code{stuItems} argument has the scored student data. It is a list where
#' each element is named with student ID and contains
#' a \code{data.frame} with at least two columns.
#' The first required column is named
#' \code{key} and shows the item name as it appears in \code{paramTab};
#' the second column in named
#' \code{score} and shows the score for that item. For dichotomous
#' items, the \code{score} is 0 or 1. For \code{GPCM} items, the scores
#' start at zero as well. For \code{GRM}, the scores start at 1.
#' 
#' For a \code{GPCM} model, \code{P0} is the \dQuote{a} parameter, and the other 
#' columns are the \dQuote{d} parameters; see the vignette for details
#' of the GPCM model.
#' 
#' The quadrature points then are a range from \code{minNode} to \code{maxNode}
#' with a total of \code{Q} nodes. 
#' 
#' @return When called for a single subscale or overall score, returns object of class \code{mmlMeans}. 
#' This is a list with elements: 
#' \itemize{
#' \item{\code{call} the call used to generate this \code{mml.means} object}
#' \item{\code{coefficients} the unscaled marginal maximum likelihood regression coefficients}
#' \item{\code{LogLik} the log-likelihood of the fit model}
#' \item{\code{getX} a function that, when called on \code{stuDat}, returns the design matrix of the marginal maximum likelihood regression}
#' \item{\code{Convergence} a convergence note from the optimizer}
#' \item{\code{location} used for scaling the estimates}
#' \item{\code{scale} used for scaling the estimates}
#' \item{\code{lnlf} the log-likelihood function of the unscaled parameters}
#' \item{\code{funcs} a list of functions that may be evaluated on a vector of parameters. The first returns the log-likelihood of the marginal maximum likelihood regression evaluated at the given parameters. The second and third correspond to the first and second derivatives, respectively, of the MML regression evaluated at the parameters.}  
#' \item{\code{rr1} the density function of each individual, conditional only on item responses in \code{stuItems}}
#' \item{\code{stuDat} the \code{stuDat} argument}
#' \item{\code{likelihood_stus} \code{rr1}, but pivoted long so that each row corresponds to a student.}
#' \item{\code{weightVar} the name of the weight variable on \code{stuDat}}
#' \item{\code{nodes} the nodes the likelihood was evaluated on}
#' \item{\code{iterations} the number of iterations required to reach convergence}
#' \item{\code{obs} the number of observations used}
#' \item{\code{testScale} the \code{testScale} used to scale the data}
#' \item{\code{weightedObs} the weighted N for the observations}
#' \item{\code{strataVar} the column name of the stratum variable on stuDat; potentially used for variance estimation}
#' \item{\code{PSUVar} the column name of the PSU variable on stuDat; potentially used for variance estimation}
#' \item{\code{itemScorePoints} a data frame that shows item IDs, the number of score points, expected scores (both from the paramTab arguments), as well as the occupied score points}
#' \item{\code{formula} the formula passed to \code{mml}}
#' \item{\code{contrasts} the contrasts used in forming the design matrix}
#' \item{\code{xlevels} the levels of the covariates used in forming the design matrix}
#' \item{\code{polyModel} the value of the argument of the same name passed to \code{mml}}
#' \item{\code{paramTab} a data frame that condenses \code{dichotParamTab} and \code{polyParamTab}}
#' \item{\code{fast} the value of the argument of the same name passed to \code{mml}}
#' \item{\code{idVar} the value of the argument of the same name passed to \code{mml}}
#' \item{\code{posteriorEsts} the posterior estimates for the people in \code{stuDat} included in the model}
#' \item{\code{pred} the predicted outcome based on the estimated coefficients}
#' \item{\code{values} a list of values used for posterior estimation}
#' \item{\code{V} a diagonal matrix of the number of columns in the design matrix}
#' }
#' 
#' When a composite score is computed there are several subscales run and the return is a \code{mmlCompositeMeans}. Many elements are themselves list with one element per construct.
#' this is a list with elements:
#' \itemize{
#' \item{\code{call} the call used to generate this \code{mml.means} object}
#' \item{\code{coefficients} matrix of the unscaled marginal maximum likelihood regression coefficients, each row represents a subscale, each column represents a coefficient}
#' \item{\code{X} the design matrix of the marginal maximum likelihood regression}
#' \item{\code{rr1} a list of elements, each the rr1 object for a subscale (see \code{mmlMeans} output)}
#' \item{\code{ids} The ID variable used for each row of \code{stuDat}}
#' \item{\code{Convergence} a vector of convergence notes from the optimizer}
#' \item{\code{lnlfl} a list of log-likelihood functions of the unscaled parameters, by construct}
#' \item{\code{stuDat} a list of \code{stuDat} data frames, as used when fitting each construct, filtered to just relevant student records}
#' \item{\code{weightVar} the name of the weight variable on \code{stuDat}}
#' \item{\code{nodes} the nodes the likelihood was evaluated on}
#' \item{\code{iterations} a vector of the number of iterations required to reach convergence on each construct}
#' \item{\code{obs} a vector of the the number of observations used on each construct}
#' \item{\code{testScale} the \code{testScale} used to scale the data}
#' \item{\code{weightedObs} a vector of the weighted N for the observations}
#' \item{\code{SubscaleVC} the covariance matrix of subscales. The residuals are assumed to be multivariate normal with this covariance matrix}
#' \item{\code{idVar} the name of the identifier used on \code{stuDat} and \code{stuItems} data}
#' \item{\code{resl} list of mmlMeans objects, one per construct}
#' \item{\code{strataVar} the column name of the stratum variable on \code{stuDat}; potentially used for variance estimation}
#' \item{\code{PSUVar} the column name of the PSU variable on \code{stuDat}; potentially used for variance estimation}
#' \item{\code{stuItems} the data frame passed to \code{mml} reformatted for use in mml}
#' \item{\code{formula} the formula passed to \code{mml}}
#' \item{\code{contrasts} the contrasts used in forming the design matrix}
#' \item{\code{xlevels} the levels of the covariates used in forming the design matrix}
#' \item{\code{polyModel} the value of the argument of the same name passed to \code{mml}}
#' \item{\code{posteriorEsts} the list of posterior estimates for the people in \code{stuDat} included in the model}
#' \item{\code{SubscaleVC} the matrix of latent correlations across constructs}
#' \item{\code{modelFrameFull} the full model frame}
#' \item{\code{stuDatSubsets} a logical matrix with a column per subscale and a row for each student; used for subsetting students within subscales}
#' \item{\code{values} a list of values used for posterior estimation}
#' }
#'
#' \code{LogLik} is not returned because there is no likelihood for a composite model.
#' 
#' @example man/examples/de.R
#' @author Harold Doran, Paul Bailey, Claire Kelley, Sun-joo Lee, and Eric Buehler 
#' @export
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr left_join %>%
#' @importFrom methods as 
#' @importFrom stats as.formula model.matrix coef dbinom sd terms reformulate complete.cases delete.response .getXlevels na.omit
#' @importFrom utils head
#' @import Rcpp
#' @import haven 
#' @useDynLib Dire, .registration = TRUE
mml <- function(formula,
                stuItems = NULL,
                stuDat,
                idVar,
                dichotParamTab = NULL,
                polyParamTab = NULL,
                testScale = NULL,
                Q = 66,
                minNode = -5,
                maxNode = 5,
                polyModel = c('GPCM', 'GRM'),
                weightVar = NULL,
                multiCore = FALSE,
                bobyqaControl = NULL,
                composite = TRUE,
                strataVar = NULL,
                PSUVar = NULL,
                fast = NULL,
                calcCor = NULL,
                verbose=0,
                retainedInformation=1,
                optimizer = c("EM", "QN")) {
  call <- match.call()
  optimizer <- match.arg(optimizer)
  if(maxNode < minNode) {
    stop("the maxNode must be larger than the minNode.")
  }
  if(!missing(bobyqaControl)) {
    message("bobyqaControl is deprecated and will be ignored. It was replace with Newton's method.")
  }
  if(!missing(fast)) {
    message("fast is deprecated and will be ignored. The fast methods are always used now.")
  }
  if(!missing(calcCor)) {
    message("calcCor is deprecated and will be ignored. The correlations are always calculated now.")
  }
  
  polyModel <- match.arg(polyModel)
  if(polyModel != 'GPCM') {
    stop(paste0(dQuote(polyModel), " has not been implemented."))
  }
  polyModel <- tolower(polyModel)
  if(!idVar %in% names(stuDat)) {
    stop(paste0(dQuote("idVar"), ", ", dQuote(idVar), ", must be a variable on stuDat."))
  }
  if(!inherits(stuDat, "data.frame")) {
    stop(paste0("Argument ", dQuote("stuDat"), " must be a data frame."))
  } else {
    if(length(class(stuDat))>1) {
      # not all data.frames are the same, recast
      stuDat <- as.data.frame(stuDat)
    }
  }
  if(missing(polyParamTab) & missing(dichotParamTab)) {
    stop(paste0("At least one of the arguments ", dQuote("polyParamTab"), " or ", dQuote("dichotParamTab"), " must be defined."))
  }

  if(missing(weightVar)) {
    weightVar <- "one"
    stuDat$one <- 1
  }
  if(!is.null(stuItems)) {
    items <- unique(stuItems$key)
    cols_to_drop <- colnames(stuDat)[colnames(stuDat) %in% items]
    if(length(cols_to_drop) > 0) {
      warning("'stuItems' is deprecated. Please include the item responses in `stuDat`; dropping colums ", paste(cols_to_drop, collapse=", "), " from stuDat and using data from stuItems instead.")
      stuDat <- stuDat[, !names(stuDat) %in% cols_to_drop]
    } else {
      warning("'stuItems' is deprecated. Please include the item responses in `stuDat`; adding columns ", paste(items, collapse = ", "), " to stuDat.")
    }
    stuItems_wide <- stuItems %>%
      pivot_wider(names_from = key, values_from = score)
    stuDat <- stuDat %>%
      left_join(stuItems_wide, by = idVar)
  }
  # formula may be missing, then define a null model
  if(missing(formula)) {
    formula <- ~ 1
  }
  outcome <- all.vars(formula)[1]
  if(outcome %in% testScale$test) {
    testScale <- testScale[testScale$test == outcome,]
  }
  else if(outcome %in% testScale$subtest) {
    testScale <- testScale[!is.na(testScale$subtest) & testScale$subtest == outcome,]
  }
  if(!composite & nrow(testScale) > 1) {
    testScale <- testScale[is.na(testScale$subtest), ]
  }
  if(nrow(testScale) > 1) {
    composite <- TRUE
    if("subtest" %in% colnames(testScale)) {
      testScale <- testScale[!is.na(testScale$subtest),]
    }
  } else {
    composite <- FALSE
  }
  if(verbose > 0) {
    message("Scoring item responses.")
  }
  scoredTest <- scoreTest(stuDat=stuDat, dichotParamTab=dichotParamTab,
                          polyParamTab=polyParamTab, testScale=testScale,
                          weightVar = weightVar, idVar = idVar,
                          Q=Q, minNode = minNode, maxNode = maxNode,
                          strataVar = strataVar, PSUVar = PSUVar)
  
  if(retainedInformation <= 0 | retainedInformation > 1) {
    stop("argument ", dQuote("Retained information"), " must be larger than 0 and not larger than 1.")
  }
  information_reduced <- FALSE
  if(retainedInformation < 1) {
    information_reduced <- TRUE
    svd_inf <- reduceInformation(data = stuDat, formula = formula, retainedInformation=retainedInformation, verbose=verbose>0)
    formula <- svd_inf$formula
    scoredTest$stuDat <- svd_inf$data
  }
  # overall: boolean for if this is a test (multiple subtest)
  # assume not, check if TRUE later
  overall <- FALSE 

  # break down to subtests and optimize those
  if(composite) {
    trms <- delete.response(terms(formula))
    mFull <- model.frame(trms, scoredTest$stuDat, na.action=na.omit)
    subtests <- sort(unique(scoredTest$paramTab$subtest))
    termLabels <- attr(terms(formula), "term.labels")
    resL <- list()
    for(sti in 1:length(subtests)) {
      if(length(termLabels) == 0) {
        termLabels <- "1"
      }
      formulai <- reformulate(termLabels, response = subtests[sti])
      if(verbose > 0) {
        message("Estimating subscale:", subtests[sti])
      }
      mmli <- mmlu(formula = formulai,
                   scored_data = scoredTest,
                   optimizer = optimizer, 
                   verbose = verbose)
      gc()
      resL <- c(resL, list(mmli))
    }
    if(verbose > 0) {
      message("Estimating posterior correlations.")
    }
    mmlc <- mmlcomp(resL, optimizer="MLE")
    # nodes are all the same
    nodes <- resL[[1]]$nodes
    k <- length(resL[[1]]$coefficients)
    for(sti in 1:length(resL)) {
      # remove, condensation complete
      resL[[sti]]$posteriorEsts <- NULL
      co <- resL[[sti]]$coefficients
      resL[[sti]]$X <- resL[[sti]]$getX(resL[[sti]]$stuDat)
      # remove standard deviation
      Xbi <- as.vector(resL[[sti]]$X %*% co[-length(co)])
      Xbdf <- data.frame(id=rownames(resL[[sti]]$X), xb=Xbi, stringsAsFactors=FALSE)
      colnames(Xbdf)[2] <- paste0("Xb", sti)
      if(sti > 1) {
        Xb <- merge(Xb, Xbdf, by="id", all=TRUE)
      } else {
        Xb <- Xbdf
      }
      # when made into a list, make it Xl
      resL[[sti]]$Xl <- resL[[sti]]$X
      resL[[sti]]$X <- NULL
      resL[[sti]]$rr1l <- resL[[sti]]$rr1
      resL[[sti]]$rr1 <- NULL
      resL[[sti]]$lnlfl <- resL[[sti]]$lnlf
      resL[[sti]]$lnlf <- NULL
      resL[[sti]]$contrastsl <- get("contrasts", envir=environment(resL[[sti]]$getX))
      resL[[sti]]$xlevelsl <- get("xlevels", envir=environment(resL[[sti]]$getX))
      resL[[sti]]$stdev <- resL[[sti]]$coefficients[length(resL[[sti]]$coefficients)]
      resL[[sti]]$stuDat <- NULL
      resL[[sti]]$iter <- resL[[sti]]$iterations
      resL[[sti]]$iterations <- NULL
      resL[[sti]]$coef <- resL[[sti]]$coefficients
      resL[[sti]]$coefficients <- NULL
      resL[[sti]]$wobs <- resL[[sti]]$weightedObs
      resL[[sti]]$weightedObs <- NULL
      resL[[sti]]$call <- NULL
      resL[[sti]]$strataVar <- NULL
      resL[[sti]]$PSUVar <- NULL
    }
    # store full res list
    resl <- resL
    for(sti in 1:length(resL)) {
      # drop nodes
      resL[[sti]]$nodes <- NULL
    }
    # exclude is for variables we don't want to set globally,
    # but need resL to retain for post hoc

    # assign variables from the many calls back onto the current frame (the env for this call to mml)
    setRes(resL, uncollapsed=c("Xl", "rr1l"),
           exclude=c("idVar","formula","polyModel","paramTab","fast","scale",
                     "location", "weightVar","testScale"))

    # add names to everything
    #names(stuDatl) <- 
    names(lnlfl) <- names(rr1l) <- names(Xl) <- names(resl) <- subtests
    names(iter) <- names(obs) <- names(wobs) <- subtests
    res <- structure(list(call = call,
                          coefficients = mmlc$coefficients,
                          X = Xl,
                          rr1 = rr1l,
                          ids = mmlc$stuDat[[idVar]],
                          Convergence = Convergence,
                          lnlfl = lnlfl,
                          stuDat = mmlc$stuDat,
                          weightVar = weightVar,
                          nodes = nodes,
                          iterations = iter,
                          obs = obs,
                          testScale = testScale,
                          weightedObs = wobs,
                          idVar = idVar,
                          resl = resl,
                          strataVar = strataVar,
                          PSUVar = PSUVar,
                          modelFrameFull = mFull,
                          posteriorEsts = mmlc$posteriorEsts,
                          formula = formula,
                          polyModel = polyModel,
                          stuDatSubsets = mmlc$stuDatSubsets,
                          values = mmlc$values,
                          SubscaleVC = mmlc$SubscaleVC
                          ),
                     class = "mmlCompositeMeans")
    if (exists("contrastsl", inherits = TRUE)) {
      res$contrasts <- get("contrastsl", inherits = TRUE)
    }
    if (exists("xlevelsl", inherits = TRUE)) {
      res$xlevels <- get("xlevelsl", inherits = TRUE)
    }
    return(res)
  } else{ # end of composite
    mmli <- mmlu(formula = formula,
                 scored_data = scoredTest,
                 optimizer = optimizer,
                 verbose = verbose)
    mmli$call <- call
    return(mmli)
  }# end if(composite) {
  
  return(res)
}


#' @importFrom dplyr inner_join
cleanStuItems <- function(si, stuDat, idVar, paramTab) {
  # make sure every stuItems student has a stuDat student, and the other way around
  if(!inherits(si, "data.frame")) {
    si <- do.call(rbind, si)
  }
  if(length(class(si))>1) {
    # not all data.frames are the same, recast
    si <- as.data.frame(si)
  }
  expectedSIVars <- c(idVar, "key", "score")
  missingVarsSI <- expectedSIVars[!expectedSIVars %in% colnames(si)]
  if(length(missingVarsSI) > 0) {
    stop("stuItems missing column(s): ", paste(missingVarsSI, collapse=", "))
  }
  # these are the only columns we need, in this order
  si <- si[ , c(idVar, "key", "score")]
  if(any(!si[[idVar]] %in% stuDat[[idVar]])) {
    missing <- (si[[idVar]])[!(si[[idVar]]) %in% stuDat[[idVar]]]
    stop(paste0("The ", dQuote("stuItems"), " argument must be a list with names that correspond to every " , dQuote("idVar"), " in ", dQuote("stuDat"), ". some missing IDs ", pasteItems(dQuote(head(missing,5))), "."))
  }
  # make sure si and stuDat are in the same order
  si[[idVar]] <- as.character(si[[idVar]])
  si$key <- as.character(si$key)
  si <- si[order(si[[idVar]], method="radix"), ]
  si <- si %>%
    inner_join(paramTab, by = c("key" = "ItemID"))
  si <- si[ , c(idVar, "key", "score")]
  # drop rows with no test data on them, they get dropped from the likelihood function anyways
  #si <- si[!is.na(si$score),]
  si <- si[order(si$key, method="radix"), ]
  return(si)
}

cleanDichotParamTab <- function(dpt) {
  if(is.null(dpt)) {
    return(NULL)
  }
  if(!inherits(dpt, "data.frame")) {
    stop(paste0("Argument ", dQuote("dichotParamTab"), " must be a data frame."))
  }
  if(length(class(dpt))>1) {
    # not all data.frames are the same, recast
    dpt <- as.data.frame(dpt)
  }
  if(!"ItemID" %in% colnames(dpt)) {
    stop(paste0("Argument ", dQuote("dichotParamTab"), " must have column ", dQuote("ItemID"), "."))
  }
  if(length(unique(dpt$ItemID)) < nrow(dpt)) {
    stop(paste0("The ", dQuote("ItemID"), " column on the ", dQuote("dichotParamTab"), " argument must be unique."))
  }
  colnames(dpt)[colnames(dpt)=="a"] <- "slope"
  colnames(dpt)[colnames(dpt) %in% c("asymptote", "g") ] <- "guessing"
  colnames(dpt)[colnames(dpt)=="d"] <- "difficulty"
  reqVar <- c("ItemID", "test", "slope", "difficulty", "guessing")
  if(any(!reqVar %in% colnames(dpt))) {
    stop(paste0("dichotParamTab must have all of the columns ", pasteItems(reqVar), "."))
  }
  if("subtest" %in% colnames(dpt)) {
    dpt$subtest <- as.character(dpt$subtest)
  }
  if("test" %in% colnames(dpt)) {
    dpt$test <- as.character(dpt$test)
  }
  return(dpt)
}

cleanPolyParamTab <- function(ppt) {
  dvars <- c()
  if(!is.null(ppt)) {
    if(!inherits(ppt, "data.frame")) {
      stop(paste0("Argument ", dQuote("polyParamTab"), " must be a data frame."))
    } else {
      if(length(class(ppt))>1) {
        # not all data.frames are the same, recast
        ppt <- as.data.frame(ppt)
      }
    }
    if(length(unique(ppt$ItemID)) < nrow(ppt)) {
      stop(paste0("The ", dQuote("ItemID"), " column on the ", dQuote("polyParamTab"), " argument must be unique."))
    }
    colnames(ppt)[colnames(ppt)=="a"] <- "slope"
    colnames(ppt)[colnames(ppt)=="d"] <- "itemLocation"
    if(!"itemLocation" %in% colnames(ppt)) {
      ppt$itemLocation <- NA
    }
    i <- 1
    while(paste0("d",i) %in% colnames(ppt)) {
      var <- paste0("d",i)
      ppt[,var] <- ifelse(is.na(ppt$itemLocation), ppt[,var], ppt$itemLocation - ppt[,var])
      dvars <- c(dvars, var)
      i <- i + 1
    }
    # prevent this from being applied twice for composite
    ppt$itemLocation <- NULL 
    if("d0" %in% colnames(ppt)) {
      if(any(ppt$d0 != 0)) {
        stop("per the vignette, in polyParamTab argument, the lowest level is named d0 and must always be zero")
      }
    } else {
      warning("polyParamTab missing d0 parameter, setting it to zero.")
    }
    ppt$d0 <- 0
    # make a column that is always NA for subsequnt program to check
    i <- i - 1
    if(any(!is.na(ppt[,paste0("d",i)]))) {
      ppt[,paste0("d",i+1)] <- NA
    }
    reqVar <- c("ItemID", "test", "slope", "scorePoints")
    if(any(!reqVar %in% colnames(ppt))) {
      stop(paste0("polyParamTab must have all of the columns ", pasteItems(reqVar), "."))
    }
  }
  # sort the d variables in hte ouput
  cn <- colnames(ppt)
  dvars <- unique(c("d0",dvars))
  cn <- cn[!cn %in% dvars]
  ppt <- ppt[,c(cn,dvars)]
  for(i in seq_along(ppt$scorePoints)) {
    for(j in 1:(ppt$scorePoints[i] + 1)) { # so that varj covers d0 through d(scorePoints)
      varj <- paste0("d",j-1)
      if(!varj %in% colnames(ppt)) {
        print(ppt[i,])
        stop(paste0("in polyParamTab argument, item ", ppt$ItemID, " is an item with ", ppt$scorePoints[i], " score points and so polyParamTab must have (at least) columns d0 to d",ppt$scorePoints[i]))
      }
      if(is.na(ppt[i,varj])) {
        print(ppt[i,])
        stop(paste0("in polyParamTab argument, item ", ppt$ItemID, " is an item with ", ppt$scorePoints[i], " score points and so polyParamTab must have non-missing columns d0 to d",ppt$scorePoints[i], " in the associated row."))
      }
    }
  }
  return(list(polyParamTab=ppt, dvars=dvars))
}

cleanStuDat <- function(stuDat, stuItems, idVar) {
  if(any(!stuDat[[idVar]] %in% stuItems[[idVar]])) {
    stuDat <- stuDat[stuDat[[idVar]] %in% stuItems[[idVar]], ]
  }
  stuDat[[idVar]] <- as.character(stuDat[[idVar]])
  stuDat <- stuDat[order(stuDat[[idVar]]), ]
  return(stuDat)
}

condenseParamTab <- function(dpt, ppt, dvars) {
  # build paramTab with info from both
  if(!is.null(dpt)) {
    paramTab <- dpt
    paramTab$scorePoints <- 1L
    if(!is.null(ppt)) {
      for(dvi in c("d0",dvars)) {
        paramTab[,dvi] <- NA
      }
      # make sure col names agree
      for(ci in colnames(paramTab)) {
        if(!ci %in% colnames(ppt)) {
          ppt[,ci] <- NA
        }
      }
      if("D" %in% colnames(ppt)) {
        if(!"D" %in% colnames(dpt)) {
          if(sd(ppt$D) < sqrt(.Machine$double.eps)) {
            paramTab$D <- ppt$D
          } else {
            paramTab$D <- 1
          }
        }
      }
      for(ci in colnames(ppt)) {
        if(!ci %in% colnames(paramTab)) {
          paramTab[,ci] <- NA
        }
      }
      paramTab <- rbind(paramTab, ppt)
    }
  } else {
    paramTab <- ppt
    paramTab$guessing <- NA
    paramTab$difficulty <- NA
  }
  paramTab$ItemID <- as.character(paramTab$ItemID)
  if(!"subtest" %in% colnames(paramTab)) {
    paramTab$subtest <- "unified"
  }
  if(!"missingCode" %in% colnames(paramTab)) {
    paramTab$missingCode <- NA
  }
  paramTab <- paramTab[order(paramTab$ItemID), ]
  return(paramTab)
}

# optimize fn.regression function fn
# verbose makes the function say more when verbose gets larger.
# needed because we use several methods in sucession
#' @importFrom stats optim
robustOptim <- function(fn, X, par0=NULL, verbose=0, fixedSD=NULL) {
  fnDerivs <- getDerivs(fn, fixedSD)
  if(verbose >= 1) {
    message("Initial optimization with optim using the L-BFGS-B method in optim.")
  }
  if(!is.null(fixedSD)) {
    fn_final_factory <- function(fixedSD, fn) {
      function(x) {
        return(fn(x,fixedSD))
      }
    }
    fn_final <- fn_final_factory(fixedSD, fn)
    par0 
  } else {
    fn_final <- fn
  }
  opt <- optim(par0, fn=fn_final, gr=fnDerivs$grad, method="L-BFGS-B", control=list(maxit=1e5, factr = 1e-10, trace = verbose >= 2, lmm=10, parscale=1/c(pmax(1,apply(X,2,sd)),1)))
  if(is.null(fixedSD)) {
    opt$par[length(opt$par)] <- max(log(1e-6)+2, opt$par[length(opt$par)])
  }
  fng <- fnDerivs$grad
  gr <- fng(opt$par)
  # this is the minimum value. Rest if it steps that far out
  if(sqrt(sum(gr^2)) / max(1, sqrt(sum(opt$par^2))) > 1e-5) {
    if(verbose >= 1) {
      message("Second optimization with lbfgs::lbfgs.")
    }
    opt <- lbfgs::lbfgs(call_eval=fn_final,
                        call_grad=fnDerivs$grad,
                        vars=opt$par,
                        invisible=ifelse(verbose > 2, 0, 1), # a backwards verbose
                        epsilon= 1e-5, # gradient convergence criterion
                        max_iterations=1e5,
                        m=max(6, sqrt(length(par0))))
    # these seem to indicate convergence,
    if(opt$convergence %in% c(-1001, -998)) {
      opt$convergence <- 0
    }
  }
  if(is.null(fixedSD)) {
    opt$par[length(opt$par)] <- max(log(1e-6)+2, opt$par[length(opt$par)])
    opt$par[length(opt$par)] <- sqrt(ifelse(opt$par[length(opt$par)] < 1, exp(opt$par[length(opt$par)] - 1), opt$par[length(opt$par)]^2))
  }
  if(opt$convergence %in% c(0)) {
    convergence <- "converged"
  } else {
    if(opt$convergence %in% 1) {
      convergence <- "Iteration limit reached"
    } else {
      convergence <- "Did not converge"
    }
  }
  opt$convergence <- convergence
  if(!"iterations" %in% names(opt)) {
    opt$iterations <- -1
  }
  if(!is.null(fixedSD)) {
    opt <- c(opt, fixedSD)
  }
  return(opt)
}


getDerivs <- function(fn, fixedSD=NULL) {
  fnGrad <- function(par) {
    if(is.null(fixedSD)) {
      return(fn(par, gr=TRUE))
    } else {
      return(fn(c(par, fixedSD), gr=TRUE))
    }
  }
  fnHess <- function(par) {
    if(is.null(fixedSD)) {
      return(fn(par, hess=TRUE))
    } else {
      return(fn(c(par, fixedSD), hess=TRUE))
    }
  }
  return(list(grad=fnGrad, hess=fnHess))
}

# assign variables from the many calls of a composite back onto the current frame
# resList: a list, each element is an mmlMeans
# pos: where to assign the variables
# pos=1 is the current parent frame (Calling environment)
# uncollapsed: these variables placed in a list that is not collapsed; other variables are collapsed with rbind
# exclude: these vriables are not updated on the frame at 'pos'
setRes <- function(resList, pos=1, uncollapsed=NULL, exclude=NULL) {
  res0 <- resList[[1]]
  rn <- names(res0)
  for(i in 1:length(rn)) {
    resi <- list()
    for(j in 1:length(resList)) {
      resi <- c(resi, list(resList[[j]][[rn[i]]]))
    }
    binder <- "c"
    if(inherits(res0[[rn[i]]], "list")) {
      binder <- "gidentity"
    }
    if(inherits(res0[[rn[i]]], "data.frame")) {
      binder <- "rbind"
    }
    if(rn[i] %in% uncollapsed) {
      binder <- "gidentity"
    }
    resi <- do.call(binder, resi)
    if(!(rn[i] %in% exclude)){
      assign(rn[i], resi, envir=parent.frame())
    }
  }

}

gidentity <- function(...) {
  inputs <- list(...)
  return(inputs)
}

checkTestScale <- function(ts) {
  if(!all(c("test", "location", "scale") %in% colnames(ts))) {
    stop("the argument ", dQuote("testScale"), " must have columns: test, location, and scale")
  }
  if("subtest" %in% colnames(ts)) {
    if(!"subtestWeight" %in% colnames(ts)) {
      stop(paste0("when you include the column ", dQuote("subscales"), " in the ", dQuote("testScale"), " argument, you must also include the column ", dQuote("subtestWeight")))
    }
  }
  ts
}


