\name{relimplmboot-class}
\docType{class}
\alias{relimplmboot-class}
\alias{show,relimplmboot-method}
%\alias{print,relimplmboot-method}

\title{Class relimplmboot }
\description{  This is the class for output objects from function \code{boot.relimp}. It is needed as input to function \code{booteval.relimp}. }
\section{Objects from the Class}{
Objects should only be created by calls to the function \code{boot.relimp}. 
}
\section{Slots}{
  \describe{
    \item{\code{boot}:}{Object of class \code{"boot"}, output from bootstrapping }
    \item{\code{type}:}{Object of class \code{"character"}, metrics that have been bootstrapped }
    \item{\code{nboot}:}{Object of class \code{"numeric"}, number of bootstrap runs }
    \item{\code{rank}:}{Object of class \code{"logical"}, have ranks been bootstrapped ? }
    \item{\code{diff}:}{Object of class \code{"logical"}, have differences been bootstrapped ? }
    \item{\code{rela}:}{Object of class \code{"logical"}, have the metrics been forced to add up to 100pct (TRUE) or not? }
    \item{\code{fixed}:}{Object of class \code{"logical"}, has bootstrapping been done for fixed (TRUE) or random (FALSE) regressors? }
    \item{\code{namen}:}{Object of class \code{"character"}: variable names of all variables (y and the regressors) }
    \item{\code{nobs}:}{Object of class \code{"numeric"}: number of valid observations used in calculations }
    \item{\code{always}:}{Object of class \code{"numintnull"}: variables always in the model }
    \item{\code{alwaysnam}:}{Object of class \code{"charnull"}: variable names of variables always in model }
  }
}
\section{Methods}{
  This documentation section documents S4 methods only. 
  There is an analogous S3-method for printing.
  \describe{
%    \item{print}{S3 method: meaningful error message, since objects of this class should not be printed }
    \item{show}{S4 method: identical to \code{print}}
  }
}
\author{ Ulrike Groemping, BHT Berlin }
\seealso{
  \code{\link{boot.relimp}}, \code{\link{booteval.relimp}},
  \code{\link{relimplmbooteval-class}},\code{\link{relimplmbootMI-class}}
}
\keyword{classes}
