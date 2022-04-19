\name{NEWS}
\title{NEWS for Dire}

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

