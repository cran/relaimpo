\name{relimplmboot-class}
\docType{class}
\alias{relimplmboot-class}

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
  }
}
\section{Methods}{
No methods defined with class "relimplmboot" in the signature.
}
\author{ Ulrike Groemping, TFH Berlin }
\seealso{
  See also \code{\link{boot.relimp}}, \code{\link{booteval.relimp}},\code{\link{relimplmbooteval-class}}
}
\keyword{classes}
