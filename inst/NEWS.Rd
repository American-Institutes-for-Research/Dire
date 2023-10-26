\name{NEWS}
\title{NEWS for Dire}

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

