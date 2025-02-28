\name{NEWS}
\title{NEWS for Dire}

\section{Changes in version 3.0.0}{
  \subsection{METHODS UPDATE}{
    \itemize{
      \item The \code{Dire mml} fitting procedures were sped up, and the memory footprint was reduced to allow for data to be fit simultaneously with more rows and more columns. A significant improvement is moving from the poorly scaling Newton’s method to the expectation maximization (EM) method that scales better for large numbers of covariates.
      \item Additional speed ups for larger datasets were afforded by moving data preparation from base R to \code{data.table}, \code{tidyr}, and \code{dplyr}.
      \item The Quasi-Newton (QN) methods now use the \code{lbfgs} package which offers superior convergence checks and removes the need for a few slow Newton steps to assure convergence.
      \item An additional new feature is built-in information reduction, which provides further speed ups. This allows users to fit a model more similar to large-scale data (LSA) such as NAEP, PISA, and TIMSS, by fitting only enough principal components of the covariate matrix to maintain a user-specified proportion of the variance in the design matrix (usually called the X matrix in a regression specification). Users can use the \code{retainedInformation} argument in \code{mml} to take advantage of this feature.
      \item The code now displays more status outputs using the \code{cli} package. This helps users monitor progress of the fit.
    }
  }
  \subsection{API CHANGES}{
    \itemize{
      \item Previously, \code{mml} required users to generate a wide file with covariates and a long file with one row per student and item. The new API preserves this option but also allows the user to simply pass one wide file that has items on it. However, the \code{stuItems} argument is now deprecated and may be removed at a future date because it is easier to simply pass a wide file.
      \item There is an \code{optimizer} argument that allows users to select between EM and QN methods. In most cases, EM is faster, but users can experiment to see if they find QN faster for their case.
      \item The \code{calcCor} argument used to allow users to calculate a composite without covariance between the elements. This is not used in any LSA and so is deprecated. When a composite is estimated, the correlations between the subscales are always calculated automatedly under the hood.
      \item Previously, the \code{fast} argument allowed the user to use faster C++ code optionally. The faster C++ code is now well tested and always used, so the fast argument is deprecated.
      \item Because the principal component analysis makes the formation of the covariate matrix more complicated, returning the X (design) matrix was substituted with a function \code{getX}. Function \code{getX}  takes student data (stuDat) and returns the design matrix or principal components for that data, mirroring if the principal component analysis was used in the original call. This also reduces the memory footprint for the result. The \code{getX} function also manages relevels of covariates.
      \item All summaries now use the gradient approximation to the Hessian. This is because of substantial performance improvements made to the gradient calculation and the intractability of calculating a full Hessian for system with a large number of covariates.
      \item To save memory, a composite fit no longer returns a list of the datasets used to generate each subscale filtered to only students who saw items in that subscale. It now returns the full dataset only once.
      \item For generalized partial credit model (GPCM) items, the \code{polyParamTab} argument must now have a d0 column that is always zero. This was previously done silently. This is intended to make the connection between the vignette description of the GPCM, the inputs, and the outputs more clearly linked to each other.
      \item Also, for GPCM items, the number of score points is consistent with general psychometric understanding that an item with three possible scores (incorrect, partially correct, and correct) has two score points.
      \item drawPVs now always uses the more defensible \code{stochasticBeta=TRUE} based on simulation results that we will share at NCME 2025.
    }
  }
}

\section{Changes in version 2.2.0}{
  \subsection{Methods update}{
    \itemize{
      \item Functions are a bit clearer about some errors when data sets do not agree. In particular, when a PSU/stratum variable is missing and \code{Taylor} variance is selected it gives a plain language error.
      \item turning on \code{multiCore} now fits the latent regressions with multiple cores too. Previously it would only fit the covariance matrix with multiple cores.
      \item added a nearly singular model check
      \item allow \code{stuDat} to have students that are not on the item data without throwing an error.
      \item optimization tries to avoid Newton's method by using the \code{lbfgs} package which allows for a condition on the gradient to be set. Newton's method can be very slow for large datasets.
      \item when Newton's method is used, the output is more verbose.
      \item the C++ implementation of the Hessian has been sped up a bit.
    }
  }
}

\section{Changes in version 2.1.1}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed bug in degrees of freedom replication in composite. This causes summary to fail in many cases.
      \item Fixed version number error in this file. 2.1.0 changes had previously been named 2.0.0.
    }
  }
}

\section{Changes in version 2.1.0}{
  \subsection{New features}{
    \itemize{
      \item Added degrees of freedom and p-values to \code{mml} results
      \item \code{mml} should be faster now
    }
  }
}

\section{Changes in version 2.0.0}{
  \subsection{New features}{
    \itemize{
      \item Added \code{drawPV} functions that draw plausible values from a normal approximation to the posterior distribution. See the \code{drawPVs} function help for details.
      \item the object returned by \code{mml} now includes an object \code{itemScorePoints} that shows, for each item, the expected and actually occupied score points.
      \item If items have invalid score points an error now shows the \code{itemScorePoints} table.
    }
  }
  \subsection{Methods update}{
    \itemize{
      \item The \code{mml} function used to use the bobyqa optimizer and now uses a combination of the \code{optim} function and then a Newton's method optimizer
      \item The \code{mml} function Taylor series covariance calculation for composite results has been updated so that the correlation is calculated for all subscales simultaneously. This results in covariance matrix that is always positive definite. The old method can be used by requesting "Partial Taylor".
    }
  }
}

