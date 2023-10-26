#' Marginal Maximum Likelihood Estimation of Linear Models
#' @description
#' Implements a survey-weighted marginal maximum estimation, a type of
#' regression where the outcome is a latent trait (such as student ability).
#' Instead of using an estimate, the likelihood function marginalizes student
#' ability. Includes a variety of variance estimation strategies.
#' 
#' @param formula  a \ifelse{latex}{\code{formula}}{\code{\link[stats]{formula}}}
#'                 object in the style of \ifelse{latex}{\code{lm}}{\code{\link[stats]{lm}}}
#' @param stuItems a \code{data.frame} where each row represents a single student's response to one item.
#'                 The columns must include the \code{idVar} column, a \code{key} column, and a
#'                 \code{score} column. Values in the \code{score} column are checked against expectations
#'                 (based on \code{dichotParamTab} and \code{polyParamTab}) and when
#'                 \code{verbose} is >= 1 a table of expected and actual levels is printed.
#' @param stuDat   a \code{data.frame} with a single row per student. Predictors in
#'                 the \code{formula} must be in \code{stuDat}.
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
#' @param multiCore allows the \code{foreach} package to be used. You should
#'                  have already setup the 
#' \ifelse{latex}{the \code{registerDoParallel} function in the \code{doParallel} package}{\code{\link[doParallel]{registerDoParallel}}}.
#' @param bobyqaControl deprecated. A list that gets passed to the \code{bobyqa} optimizer in \code{minqa}
#' @param composite a logical indicating if an overall test should be treated as
#'                  a composite score; a composite is a weighted average of the
#'                  subscales in it.
#' @param strataVar character naming a variable on \code{stuDat}, the variable indicating the
#'                  stratum for each row. Used in post-hoc robust variance estimation.
#' @param PSUVar character naming a variable on \code{stuDat}; the primary sampling unit
#'               (PSU) variable. Used in post-hoc robust variance estimation. The values
#'               do not need to be unique across strata.
#' @param fast a logical indicating if cpp code should be used in \code{mml} processes. This should 
#'             yield speed-ups to runs. 
#' @param calcCor set to \code{TRUE} to calculate covariances. Needed to estimate variances and form plausible values
#' @param verbose integer, negative or zero for no details, increasingly verbose messages at one and two
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
#' \item{\code{X} the design matrix of the marginal maximum likelihood regression}
#' \item{\code{Convergence} a convergence note from the optimizer}
#' \item{\code{location} used for scaling the estimates}
#' \item{\code{scale} used for scaling the estimates}
#' \item{\code{lnlf} the log-likelihood function of the unscaled parameters} 
#' \item{\code{rr1} the density function of each individual, conditional only on item responses in \code{stuItems}}
#' \item{\code{stuDat} the \code{stuDat} argument}
#' \item{\code{weightVar} the name of the weight variable on \code{stuDat}}
#' \item{\code{nodes} the nodes the likelihood was evaluated on}
#' \item{\code{iterations} the number of iterations required to reach convergence}
#' \item{\code{obs} the number of observations used}
#' \item{\code{weightedObs} the weighted N for the observations}
#' \item{\code{strataVar} the column name of the stratum variable on stuDat; potentially used for variance estimation}
#' \item{\code{PSUVar} the column name of the PSU variable on stuDat; potentially used for variance estimation}
#' \item{\code{itemScorePoints} a data frame that shows item IDs, the number of score points, expected scores (both from the paramTab arguments), as well as the occupied score points}
#' \item{\code{stuItems} the data frame passed to \code{mml} reformatted for use in mml}
#' \item{\code{formula} the formula passed to \code{mml}}
#' \item{\code{contrasts} the contrasts used in forming the design matrix}
#' \item{\code{xlevels} the levels of the covariates used in forming the design matrix}
#' \item{\code{polyModel} the value of the argument of the same name passed to \code{mml}}
#' \item{\code{paramTab} a data frame that condenses \code{dichotParamTab} and \code{polyParamTab}}
#' \item{\code{fast} the value of the argument of the same name passed to \code{mml}}
#' \item{\code{idVar} the value of the argument of the same name passed to \code{mml}}
#' \item{\code{posteriorEsts} the posterior estimates for the people in \code{stuDat} included in the model}
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
#' \item{\code{SubscaleVC} the covariance matrix of subscales. The residuals are assumed to be multivariate normal with this covairiance matrix}
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
#' }
#'
#' \code{LogLik} is not returned because there is no likelihood for a composite model.
#' 
#' @example man/examples/de.R
#' @author Harold Doran, Paul Bailey, Claire Kelley, Sun-joo Lee, and Eric Buehler 
#' @export
#' @importFrom methods as 
#' @importFrom stats as.formula model.matrix coef dbinom sd terms reformulate complete.cases delete.response .getXlevels na.omit
#' @importFrom utils head
#' @import Rcpp
#' @import haven 
#' @useDynLib Dire, .registration = TRUE
mml <- function(formula,
                stuItems,
                stuDat,
                idVar,
                dichotParamTab = NULL,
                polyParamTab = NULL,
                testScale = NULL,
                Q = 30,
                minNode = -4,
                maxNode = 4,
                polyModel = c('GPCM', 'GRM'),
                weightVar = NULL,
                multiCore = FALSE,
                bobyqaControl = NULL,
                composite = TRUE,
                strataVar = NULL,
                PSUVar = NULL,
                fast = TRUE,
                calcCor = TRUE,
                verbose=0) {
  if(!missing(bobyqaControl)) {
    message("bobyqaControl is deprecated and will be ignored. It was replace with Newton's method.")
  }
  # check for parrallel if multiCore True 
  if(multiCore == TRUE){
    # check for doParallel
    if(!requireNamespace("doParallel")) {
      message("Unable to load package doParallel, setting multiCore to FALSE. Install doParallel to use multiCore option.")
      multiCore <- FALSE
    } else {
      # always use C++ version when multicore
      fast <- TRUE
    }
  }
    
  call <- match.call()
  polyModel <- match.arg(polyModel)
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
  dichotParamTab <- cleanDichotParamTab(dichotParamTab) 
  # clean polyParamTab, get the list of d variables on it
  pptd <- cleanPolyParamTab(polyParamTab)
  paramTab <- condenseParamTab(dichotParamTab, pptd$polyParamTab, pptd$dvars)
  
  # save a copy of the raw stuItems for potential composite call
  stuItems0 <- stuItems
  stuItems <- cleanStuItems(stuItems, stuDat, idVar)
  stuDat <- cleanStuDat(stuDat, stuItems, idVar)

  # overall: boolean for if this is a test (multiple subtest)
  # assume not, check if TRUE later
  overall <- FALSE 
  # formula may be missing, then define a null model
  if(missing(formula)) {
    formula <- ~ 1
  } # end if(missing(formula))
  testScale2 <- testScale # potentially update testScale
  if(length(formula) > 2) {
    outcome <- as(formula[[2]], "character")
    if(outcome %in% c(paramTab$test)) {
      testScale2 <- testScale[testScale$test %in% outcome, ]
      if("subtest" %in% colnames(testScale)) {
        testScale2 <- testScale2[is.na(testScale2$subtest), ]
      }
      paramTab <- paramTab[paramTab$test %in% outcome, ]
      overall <- TRUE
    } else {
      if( "subtest" %in% colnames(paramTab) & outcome %in% paramTab$subtest) {
        paramTab <- paramTab[paramTab$subtest %in% outcome,]
        testScale2 <- testScale[testScale$subtest %in% outcome,]
      } else {
        stop(paste0("Outcome variable ", dQuote(outcome)," not found in any of  ", pasteItems(dQuote(c("dichotParamTab$test", "dichotParamTab$subtest", "polyParamTab$test", "polyParamTab$subtest")), "nor"), "."))
      }
    }
  } else { # else for     if(length(formula) > 2) 
    # the user did not specify a test or subtest
    if(length(unique(paramTab$test)) > 1) {
      stop("when there are multiple tests defined the test must be specified as the outcome in the formula.")
    }
    if(!"subtest" %in% colnames(paramTab)) {
      paramTab$subtest <- "default subtest"
      composite <- FALSE
      overall <- FALSE
    }
    if(length(unique(paramTab$subtest)) > 1) {
      overall <- TRUE
    }
  } # end else for if(length(formula) > 2)
  # in this situation, it is not composite/overall
  if(length(unique(paramTab$subtest)) == 1) {
    overall <- FALSE
    composite <- FALSE
  }
  # subset testScale if not composite
  if(!composite) {
    testScale <- testScale2
  }
  # break down to subtests and optimize those
  if(composite) {
    # now form X
    trms <- delete.response(terms(formula))
    mFull <- model.frame(trms, stuDat, na.action=na.omit)
    subtests <- sort(unique(paramTab$subtest))
    testScale <- testScale[order(testScale$subtest), ]
    testScale <- testScale[testScale$subtest %in% subtests, ]
    coef <- c()
    lnlfl <- list()
    rr1l <- list()
    stuDatl <- list()
    Xl <- list()
    contrastsl <- list()
    xlevelsl <- list()
    iter <- c()
    obs <- c()
    wobs <- c()
    Convergence <- c()
    subt <- c()
    resl <- list()
    stdev <- c()
    posteriorEsts<- data.frame()
    calls <- list()
    for(sti in 1:length(subtests)) {
      subt <- c(subt, subtests[sti])
      termLabels <- attr(terms(formula), "term.labels")
      if(length(termLabels) == 0) {
        termLabels <- "1"
      }
      formulai <- reformulate(termLabels, response = subtests[sti])
      calli <- call
      calli$formula <- formulai
      # this resolves any lazy unresolved calls in variables
      calli$stuItems <- stuItems
      calli$stuDat <- stuDat
      calli$dichotParamTab <- dichotParamTab
      calli$polyParamTab <- polyParamTab
      calli$testScale <- testScale
      calli$strataVar <- strataVar
      calli$PSUVar <- PSUVar
      # this subtest will not itself be a composite
      calli$composite <- FALSE
      # we multicore on the outer loop, not inner
      calli$multiCore <- FALSE
      calls <- c(calls, list(calli))
    }
    
    if(multiCore) {
      if(verbose >= 1) {
        cat("\n")
        message(paste0("Estimating constructs in parallel."))
      }
      itc <- iter(calls)
      resL <- foreach(dopari=itc) %dopar% {
        res <- eval(dopari)
        calliName <- call
        calliName$formula <- formulai
        res$call <- calliName
        return(res)
      }
    } else {
      resL <- list()
      for(i in 1:length(calls)) {
        if(verbose >= 1) {
          cat("\n")
          message(paste0("Estimating construct ", dQuote(subtests[i])))
        }
        calli <- calls[[i]]
        res <- eval(calli)
        calliName <- call
        calliName$composite <- FALSE
        calliName$formula <- formulai
        res$call <- calliName
        resL <- c(resL, list(eval(calli)))
      }
    }
    names(resL) <- subtests
    # nodes are all the same
    nodes <- resL[[1]]$nodes
    k <- length(resL[[1]]$coefficients)
    for(sti in 1:length(resL)) {
      # fix posteriors
      #pei <- resL[[sti]]$posteriorEsts
      nidv <- colnames(resL[[sti]]$posteriorEsts) != "id" # non id variable names
      colnames(resL[[sti]]$posteriorEsts)[nidv] <- make.names(paste0(colnames(resL[[sti]]$posteriorEsts)[nidv], "_", subtests[sti]))
      if(sti == 1) {
        posteriorEsts <- resL[[sti]]$posteriorEsts
      } else {
        posteriorEsts <- merge(posteriorEsts, resL[[sti]]$posteriorEsts, by="id", all=TRUE)
      }
      # remove, condensation complete
      resL[[sti]]$posteriorEsts <- NULL
      # get Xb
      co <- resL[[sti]]$coefficients
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
      resL[[sti]]$contrastsl <- resL[[sti]]$contrasts
      resL[[sti]]$contrasts <- NULL
      resL[[sti]]$xlevelsl <- resL[[sti]]$xlevels
      resL[[sti]]$xlevels <- NULL
      resL[[sti]]$stdev <- resL[[sti]]$coefficients[length(resL[[sti]]$coefficients)]
      resL[[sti]]$stuDatl <- resL[[sti]]$stuDat
      resL[[sti]]$stuDat <- NULL
      resL[[sti]]$iter <- resL[[sti]]$iterations
      resL[[sti]]$iterations <- NULL
      resL[[sti]]$coef <- resL[[sti]]$coefficients
      resL[[sti]]$coefficients <- NULL
      resL[[sti]]$call <- NULL
      resL[[sti]]$strataVar <- NULL
      resL[[sti]]$PSUVar <- NULL
      resL[[sti]]$itemScorePoints <- NULL
      resL[[sti]]$wobs <- resL[[sti]]$weightedObs
      resL[[sti]]$weightedObs <- NULL
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
    setRes(resL, uncollapsed=c("Xl", "rr1l", "stuDatl"),
           exclude=c("idVar","formula","polyModel","paramTab","fast","scale",
                     "location", "weightVar"))
    co <- coefficients
    # merge together posterior estimates
    # nodes are the same for every run

    if(is.null(weightVar)) {
      stuDat$one <- 1
      stuDatl <- lapply(stuDatl, function(x) {
        x$one <- 1
        return(x)
      })
      weightVar <- "one"
    }
    vc <- matrix(0, nrow=length(subtests), ncol=length(subtests))

    diag(vc) <- stdev^2
    wgt <- stuDat[ , c(idVar, weightVar)]
    colnames(wgt)[2] <- "w"
    # vcfmat is the correlation functions, just the portion above the diagonal.
    # use a list because we're storing functions
    if(multiCore) {
      corm <- data.frame(i= rep(1:length(subtests), each=length(subtests)),
                         j= rep(1:length(subtests), length(subtests)))
      corm <- corm[corm$i > corm$j, ]
      corms <- split(corm, 1:nrow(corm))
      itc <- iter(corms)
      if(verbose >= 1) {
        message("Estimating correlations in parallel.")
      }
      cori <- foreach(dopari = itc, .packages="Dire", .export = c("optimize", "calcRrij", "mmlCor", "fnCor")) %dopar% {
        i <- dopari$i
        j <- dopari$j
        xbi <- Xb[ , paste0("Xb",i)]
        xbj <- Xb[ , paste0("Xb",j)]
        # just cases in both
        ijss <- !is.na(xbi) & !is.na(xbj)
        xbi <- xbi[ijss]
        xbj <- xbj[ijss]
        # subset rr1 
        rr1i <- rr1l[[i]]
        rr1j <- rr1l[[j]]
        rr1i <- rr1i[ , colnames(rr1i) %in% Xb$id[ijss]]
        rr1j <- rr1j[ , colnames(rr1j) %in% Xb$id[ijss]]
        if( all.equal(colnames(rr1i), Xb$id[ijss])[1] != TRUE) {
          stop("rr1i names do not agree with Xb names.")
        }
        if( all.equal(colnames(rr1j), Xb$id[ijss])[1] != TRUE) {
          stop("rr1j names do not agree with Xb names.")
        }
        w <- wgt[wgt[,idVar] %in% Xb$id[ijss],]
        if( all.equal(w[,idVar], Xb$id[ijss])[1] != TRUE) {
          stop("weigth names do not agree with Xb names.")
        }
        vcij <- mmlCor(Xb1=xbi,
                       Xb2=xbj,
                       s1=stdev[i],
                       s2=stdev[j],
                       rr1=rr1i,
                       rr2=rr1j,
                       weights=w$w,
                       nodes=nodes,
                       fast=fast)
        return(list(i=i,j=j,vc=vcij))
      }
      # make a long data frame with columns i, j, and cor
      vclong <- as.data.frame(do.call(rbind, lapply(cori, function(x) { unlist(c(x$i, x$j, x$vc$rho)) } )))
      colnames(vclong) <- c("i", "j", "cor")
      for(i in 1:length(subtests)) {
        for(j in 1:length(subtests)) {
          if(i > j) {
            vc[i,j] <- vclong[vclong$i == i & vclong$j == j, "cor"] 
            vc[j,i] <- vc[i,j]
            foundSub <- FALSE
            indi <- 1
            while(!foundSub) {
              corii <- cori[[indi]]
              if(corii$i == i & corii$j == j) {
                foundSub <- TRUE
              } 
              indi <- indi + 1
            }
          }
        }
      }
    } else { # end: if(multiCore)
      # single cor correlation calculation
      for(i in 1:length(subtests)) {
        for(j in 1:length(subtests)) {
          if(i > j) {
            # subset Xb
            xbi <- Xb[,paste0("Xb",i)]
            xbj <- Xb[,paste0("Xb",j)]
            # just cases in both
            ijss <- !is.na(xbi) & !is.na(xbj)
            xbi <- xbi[ijss]
            xbj <- xbj[ijss]
            # subset rr1 
            rr1i <- rr1l[[i]]
            rr1j <- rr1l[[j]]
            rr1i <- rr1i[ , colnames(rr1i) %in% Xb$id[ijss]]
            rr1j <- rr1j[ , colnames(rr1j) %in% Xb$id[ijss]]
            if( all.equal(colnames(rr1i), Xb$id[ijss])[1] != TRUE) {
              stop("rr1i names do not agree with Xb names.")
            }
            if( all.equal(colnames(rr1j), Xb$id[ijss])[1] != TRUE) {
              stop("rr1j names do not agree with Xb names.")
            }
            w <- wgt[wgt[ , idVar] %in% Xb$id[ijss],]
            if( all.equal(w[ , idVar], Xb$id[ijss])[1] != TRUE) {
              stop("weigth names do not agree with Xb names.")
            }
            if(verbose >= 1) {
              message(paste0("Estimating construct correlations between ", subtests[i], " and ", subtests[j]))
            }
            vcf <- mmlCor(Xb1=xbi,
                          Xb2=xbj,
                          s1=stdev[i],
                          s2=stdev[j],
                          rr1=rr1i,
                          rr2=rr1j,
                          weights=w$w,
                          nodes=nodes,
                          fast=fast)
            vc[j, i] <- vc[i, j] <- vcf$rho
          }
        }
      }
    }
    colnames(vc) <- rownames(vc) <- subtests
    # add names to everything
    names(stuDatl) <- names(lnlfl) <- names(rr1l) <- names(Xl) <- names(resl) <- subtests
    names(iter) <- names(obs) <- names(wobs) <- subtests
    coefM <- matrix(coef, nrow=length(subtests), byrow=TRUE)
    rownames(coefM) <- subtests
    colnames(coefM) <- names(coef)[1:ncol(coefM)]
    res <- structure(list(call = call,
                          coefficients = coefM,
                          X = Xl,
                          contrasts = contrastsl,
                          xlevels = xlevelsl,
                          rr1 = rr1l,
                          ids = stuDat[[idVar]],
                          Convergence = Convergence,
                          lnlfl = lnlfl,
                          stuDat = stuDatl,
                          stuItems = stuItems,
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
                          posteriorEsts = posteriorEsts,
                          formula = formula),
                     class = "mmlCompositeMeans")
    if(calcCor) {
      res$SubscaleVC <- vc
    }
    return(res)
  } # end if(composite) {

  # done with the left hand side of the formula, if it exists, drop it
  if(length(formula) > 2) {
    formula[[2]] <- NULL
  }
  # subset stuDat to valid data
  cc <- complete.cases(stuDat[,c(all.vars(formula), weightVar, strataVar, PSUVar)])
  stuDat <- stuDat[cc, ]
  stuItems <- stuItems[stuItems[[idVar]] %in% stuDat[[idVar]],]
  if(nrow(stuDat) == 0) {
    stop(paste0("no complete cases in ", dQuote("stuDat"), "."))
  }
  if(nrow(stuItems) == 0) {
    stop(paste0("no complete cases in ", dQuote("stuItems"), "."))
  }

  # only keep items in the paramTab
  stuItems <- stuItems[stuItems$key %in% paramTab$ItemID, ]
  if(nrow(stuItems) == 0) {
    stop(paste0("no student scored found for any of the relevant items. Check that dichotParamTab and polyParamTab ", dQuote("ItemID"), " columns agree with stuITems ", dQuote("key"), " column values."))
  }
  paramTab <- paramTab[order(paramTab$ItemID), ]
  stu <- stuItems[order(stuItems$key), ]
  # check response ranges
  agg <- data.frame(key = unique(stu$key))
  for(i in 1:nrow(agg)) {
    agg$scorePoints[i] <- subset(paramTab, ItemID == agg$key[i], "scorePoints")[[1]]
    vals <- stu$score[stu$key == agg$key[i]]
    # drop missing code
    vals <- vals[!vals %in% paramTab$missingCode[paramTab$ItemID == agg$key[i]]]
    if(agg$scorePoints[i] == 1) {
      expectedOccupied <- 0:1
    } else {
      if(polyModel == "gpcm") {
        expectedOccupied <- seq(0, agg$scorePoints[i], by=1)
      } else {
        expectedOccupied <- seq(1, agg$scorePoints[i]+1, by=1)
      }
    }
    agg$expectedOccupied[i] <- paste( expectedOccupied, collapse=":")
    agg$occupied[i] <- paste(names(table(vals)), collapse=":")
    agg$Result[i] <- ifelse( any(!expectedOccupied %in% names(table(vals)) ), "!", "\U02713")
    # red check overrules other states
    agg$Result[i] <- ifelse( any(!names(table(vals)) %in% expectedOccupied), "\U274C", agg$Result[i])
  }
  rownames(agg) <- paste0(agg$Result, " " , agg$key)
  agg$score <- NULL
  agg$key <- NULL
  if( any(agg$Result %in% "\U274C") ) {
    agg$Result <- NULL
    if(verbose >= 1) {
      print(agg)
      stop("Some items score points inconsistent with expectations.")
    }
    stop("Some items score points inconsistent with expectations; increase verbose level to at least 1 see table.")
  }
  agg$Result <- NULL
  stu <- split(stu, stu[[idVar]])

  # subset to students who have at least one valid score
  stuDat <- stuDat[stuDat[[idVar]] %in% names(stu), ]

  # now form X
  trms <- delete.response(terms(formula))
  m <- model.frame(trms, data=stuDat, drop.unused.levels = TRUE)
  X <- model.matrix(formula, m)
  eig <- eigen(crossprod(X))
  ev <- eig$values/eig$values[1]
  if(any(ev <= 100 * .Machine$double.eps)) {
    cs <- colSums(X)
    cs[cs <= which.min(cs)]
    message("Thin levels:")
    print(cs)
    stop("Nearly singular design matrix. Consider adjusting the model to improve the design matrix.")
  }

  #for prediction
  contrasts <- attributes(X)$contrasts
  xlevels <- .getXlevels(trms, m)
  rownames(X) <- stuDat[[idVar]]
  K <- ncol(X) # number of fixed parameters to estimate   
  nms <- c(colnames(X), 'Population SD') 
  startVal <- c(rep(0, K), 1)
  if(nrow(X) != nrow(stuDat)) {
    stop("Missing values not allowed in independent variables.")
  }
  ### These are the exported functions passed to the dopar function
  nodes <- seq(from = minNode, to = maxNode, length.out = Q)

  ### This portion of the code computes all the likelihood evaluations
  ### and does so outside of the function that is maximized
  ### This saves a lot of overhead by using fixed quadrature points
  if(verbose >= 1) {
    message("Calculating likelihood function.")
  }
  if(multiCore) {
    # use multiple cores to calculate rr1
    rr1 <- calcRR1_dopar(stu, Q, polyModel, paramTab, nodes, fast)
  } else {
    # calculate rr1 on a single core
    rr1 <- calcRR1(stu, Q, polyModel, paramTab, nodes, fast)
  }

  # add names to rr1
  colnames(rr1) <- names(stu)
  if(!all.equal(colnames(rr1), stuDat[[idVar]])) {
    stop("Sorting error in mml.")
  }
  eDecomp <- eigen(crossprod(X))
  if(verbose >= 1) {
    message(paste0("design matrix condition number = ", round(eDecomp$values[1] / eDecomp$values[length(eDecomp$values)],4) ))
  }
  if( (min(eDecomp$values) < 0) || eDecomp$values[1] / eDecomp$values[length(eDecomp$values)] > 1/.Machine$double.eps) {
    stop("Design matrix exactly singular. Adjust covariates to avoid perfect multicolinearity.")
  }
  
  fn2 <- fn.regression(X_=X, i=NULL, wv=weightVar, rr1=rr1, nodes=nodes, stuDat=stuDat)
  opt <- robustOptim(fn2, startVal, verbose=verbose, X=X)
  posteriorEsts <- fn2(opt$par, returnPosterior=TRUE)
  names(opt$par) <- c(colnames(X), "s")
  # default location and scale
  location <- NA
  scale <- NA
  if(!is.null(testScale) & nrow(testScale) == 1) {
    location <- testScale$location[1]
    scale <- testScale$scale[1]
  }
  if(is.null(weightVar)) {
    obs <- sum(apply(rr1,2,sum)>0)
    weightedObs <- obs
  } else {
    obs <- sum(stuDat[ , weightVar] > 0 & apply(rr1, 2, sum) > 0) # number with positive weight, at least one response
    weightedObs <- sum(stuDat[apply(rr1, 2, sum) > 0, weightVar]) # sum of weights for those with positive weight, at least one response
  }
  coefficients <- opt$par
  names(coefficients) <- nms
  # report on theta scale if no scale found
  if(is.na(scale) & is.na(location)) {
    scale <- 1
    location <- 0
  }
  if(is.na(scale) | is.na(location)) {
    if(is.na(scale)) {
      warning(paste0("Could not find a valid scale. Resetting to theta scale. Check ", dQuote("testDat"), " argument."))
    } else {
      warning(paste0("Could not find a valid location. Resetting to theta scale. Check ", dQuote("testDat"), " argument."))
    }
    scale <- 1
    location <- 0
  }
  assign("insd", FALSE, envir=environment(fun=fn2))
  res <- structure(list(call = call,
                        coefficients = coefficients,
                        LogLik = -1/2*opt$value,
                        X = X,
                        Convergence = opt$convergence,
                        location = location,
                        scale = scale,
                        lnlf= fn2,
                        rr1= rr1,
                        stuDat = stuDat,
                        stuItems = stuItems,
                        weightVar = weightVar,
                        nodes = nodes,
                        iterations = opt$iter,
                        obs = obs,
                        weightedObs = weightedObs,
                        strataVar = strataVar,
                        PSUVar = PSUVar,
                        formula = formula,
                        contrasts = contrasts,
                        xlevels = xlevels,
                        polyModel = polyModel,
                        paramTab = paramTab,
                        fast = fast,
                        idVar = idVar,
                        posteriorEsts = posteriorEsts,
                        itemScorePoints = agg),
                   class = "mmlMeans")
  return(res)
}

#' @importFrom stats optimize optim reshape rnorm
mmlCor <- function(Xb1,
                   Xb2,
                   s1,
                   s2,
                   rr1,
                   rr2,
                   nodes,
                   weights=NULL,
                   fast=TRUE,
                   verbose=0) { 
  # fn2 uses the Fisher-Z transformation
  # fine rebins rr1 for large correlations where numerical integration becomes difficult
  fn2f <- fnCor(Xb1=Xb1, Xb2=Xb2, s1=s1, s2=s2, w=weights, rr1=rr1, rr2=rr2, nodes=nodes, fine=TRUE, fast = fast)
  ## optimize in Fisher-Z space
  # first ball park
  opt <- optimize(fn2f, c(1, 3), tol=.Machine$double.eps^0.25)$minimum
  # transform back to original space and return result
  return(list(rho=tanh(opt) * s1 * s2, corLnl=fn2f))
}


cleanStuItems <- function(si, stuDat, idVar) {
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
  si <- si[,c(idVar, "key", "score")]
  if(any(!si[[idVar]] %in% stuDat[[idVar]])) {
    missing <- (si[[idVar]])[!(si[[idVar]]) %in% stuDat[[idVar]]]
    stop(paste0("The ", dQuote("stuItems"), " argument must be a list with names that correspond to every " , dQuote("idVar"), " in ", dQuote("stuDat"), ". some missing IDs ", pasteItems(dQuote(head(missing,5))), "."))
  }
  # make sure si and stuDat are in the same order
  si[[idVar]] <- as.character(si[[idVar]])
  si$key <- as.character(si$key)
  si <- si[order(si[[idVar]]), ]
  # drop rows with no test data on them, they get dropped from the likelihood function anyways
  si <- si[!is.na(si$score),]
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
    ppt[,"d0"] <- 0
    i <- i - 1
    if(any(!is.na(ppt[,paste0("d",i)]))) {
      ppt[,paste0("d",i+1)] <- NA
    }
    reqVar <- c("ItemID", "test", "slope", "scorePoints")
    if(any(!reqVar %in% colnames(ppt))) {
      stop(paste0("polyParamTab must have all of the columns ", pasteItems(reqVar), "."))
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
  return(paramTab)
}

# optimize fn.regression function fn
# verbose makes the function say more when verbose gets larger.
# needed because we use several methods in sucession
robustOptim <- function(fn, X, par0=NULL, verbose=0) {
  fnDerivs <- getDerivs(fn)
  if(verbose >= 1) {
    message("Initial optimization with optim using the L-BFGS-B method in optim.")
  }
  opt <- optim(par0, fn=fn, gr=fnDerivs$grad, method="L-BFGS-B", control=list(maxit=1e5, factr = 1e-10, trace = verbose >= 2, lmm=10, parscale=1/c(pmax(1,apply(X,2,sd)),1)))
  if(verbose >= 1) {
    message("Second optimization with lbfgs::lbfgs.")
  }
  opt$par[length(opt$par)] <- max(log(1e-6)+2, opt$par[length(opt$par)])
  fng <- fnDerivs$grad
  gr <- fng(opt$par)
  # this is the minimum value. Rest if it steps that far out
  if(sqrt(sum(gr^2)) / max(1, sqrt(sum(opt$par^2))) > 1e-5) {
    opt <- lbfgs::lbfgs(call_eval=fn,
                        call_grad=fnDerivs$grad,
                        vars=opt$par,
                        invisible=ifelse(verbose > 2, 0, 1), # a backwards verbose
                        epsilon= 1e-5, # gradient convergence criterion
                        max_iterations=1e5,
                        m=max(6, sqrt(length(par0))))
  }
  opt$par[length(opt$par)] <- max(log(1e-6)+2, opt$par[length(opt$par)])
  # map the SD term back
  opt$par[length(opt$par)] <- sqrt(ifelse(opt$par[length(opt$par)] < 1, exp(opt$par[length(opt$par)] - 1), opt$par[length(opt$par)]^2))
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
  return(opt)
}

grad_descent <- function(x0, fn, grf, max_it = 10 * length(x0), verbose=0, c1=0.001, c2=0.1) {
  x <- x0
  fx <- fn(x)
  gx <- grf(x)
  alpha <- 1
  itter <- 0
  while(itter < max_it && sqrt(sum(gx^2))/max(1,sqrt(sum(x^2))) > 1e-5) {
    #it's possible to get stuck in a oscelator in a nice quadratic peak. This stops that.
    alpha <- alpha / 2
    itter <- itter + 1
    wolf_conds_met <- FALSE
    mina <- 0
    maxa <- Inf
    cat("  alpha =",alpha,"\n")
    while(!wolf_conds_met) {
      w1met <- checkW1(fn, grf, x, gx, fx, alpha, c1)
      w2met <- checkW2(fn, grf, x, gx, fx, alpha, c2)
      wolf_conds_met <- w1met && w2met
      if(!wolf_conds_met) {
        if(!w1met) {
          maxa <- min(maxa, alpha)
        }
        if(!w2met) {
          mina <- max(mina, alpha)
        }
        if(maxa < Inf) {
          alpha <- (mina + maxa)/ 2
        } else {
          alpha <- alpha * 2
        }
      }
      
    }
    x <- x + alpha * gx
    fx <- fn(x)
    gx <- grf(x)
    if(verbose > 0) {
      cat("value = ", format(fx, digits=16), " gradient condition =", format(sum(abs(gx))/(length(x0) * 1e-5), digits=16), " > 1 \n")
    }
  }
  return(x)
}

getDerivs <- function(fn) {
  fnGrad <- function(par) {
    fn(par, gr=TRUE)
  }
  fnHess <- function(par) {
    fn(par, hess=TRUE)
  }
  return(list(grad=fnGrad, hess=fnHess))
}

Newton <- function(par0, iter0, verbose, fn, fng, fnh) {
  # push to actual convergence
  opt <- list(par=par0, iter=iter0, convergence=0)
  # get an actual gr, maybe we don't need more refinement
  gr <- fng(opt$par)
  # this is the minimum value. Rest if it steps that far out
  opt$par[length(opt$par)] <- max(log(1e-6)+2, opt$par[length(opt$par)])
  # first try BFGS
  gr_prev <- 2*max(abs(gr)) # allows the first condition alone to decide if it is a good idea to try BFGS
  while(sqrt(sum(gr)) / max(1, sqrt(sum(opt$par^2))) < 1e-5 && opt$iter < 50 + iter0) {
    if(iter0 == opt$iter && verbose >= 1) {
      message("Further refining optimization with Newton's method.")
    }
    gr <- fng(opt$par)
    H <- fnh(opt$par)
    if(any(H %in% c(NA, NaN, Inf, -Inf))) {
      stop("Hessian is not defined.")
    }
    # sometimes H is just not PD, adjust it slightly to be PD in that case
    H <- nearPD2(H)
    update <- qr.solve(qr(H), -1*gr)
    opt$par <- opt$par + update
    opt$iter <- opt$iter + 1
    if(verbose >= 2) {
      message(paste0("  step =", opt$iter - iter0))
      message(paste0("    lnl =", fn(opt$par)))
      message(paste0("    ||gradient||2=", sqrt(sum(gr^2))))
      message(paste0("    max gradient=", max(abs(gr))))
      message(paste0("    ||step||2=", sqrt(sum(update^2))))
      message(paste0("    max step =", max(abs(update))))
    }
    gr <- gr / pmax(1, opt$par)
    # this is the minimum value. Reset if it steps that far out
    opt$par[length(opt$par)] <- max( (log(1e-6)+1), opt$par[length(opt$par)])
  }
  if(opt$iter == 50 + iter0) {
    opt$convergence <- "max iterations reached"
    warning("Convergence not reached. Try scaling parameters.")
  }
  opt$par <- as.vector(opt$par)
  opt$value <- fn(opt$par)
  if(verbose >= 2) {
    message(paste0("final lnl=", opt$value))
  }
  return(opt)
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

checkW1 <- function(fn, grf, x, gx, fx, alpha, c1) {
  fn(x + alpha * gx) <= fx + c1 * alpha * sum(gx^2)
}

checkW2 <- function(fn, grf, x, gx, fx, alpha, c2) {
  -gx * grf(x + alpha * gx) <= -c2 * sum(gx^2)
}
