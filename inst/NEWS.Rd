\name{NEWS}
\title{NEWS for Dire}

\section{Changes in version 2.0.0}{
  \subsection{New features}{
    \itemize{
      \item Added \code{drawPV} functions that draw plasuible values from a normal approximation to the posterior distribution. See the \code{drawPVs} function help for details.
      \item If items have invalid score points an error now shows the table of items, expected score points (based on the paramTab) and observed score points. This is returned as \code{occupied} when there is no error and shows when there are unoccupied score points.
    }
  }
  \subsection{Methods update}{
    \itemize{
      \item changed from bobyqa to first using optim and then Newton
      \item update covariance calculation for composite results so that the correlation is constant subscales. This results in covariance matrix that is always positive definite.
    }
  }
}