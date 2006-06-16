\name{relaimpo-package}
\alias{relaimpo-package}
\alias{relaimpo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Package to calculate relative importance metrics for linear models  }
\docType{package}
\description{
  \pkg{relaimpo} calculates several relative importance metrics for the linear model. 
The recommended metrics are \code{lmg} (\eqn{R^2} partitioned by averaging over orders, like in Lindemann, Merenda and Gold (1980, p.119ff)) 
and \code{pmvd} (a newly proposed metric by Feldman (2005), non-US version only). 
For completeness, several other metrics are also on offer. Other packages with related topics: \code{hier.part}, \code{relimp}.
}
\details{
  \pkg{relaimpo} calculates the metrics and also offers the possibility of bootstrapping them and of displaying results in print and graphically. 

  It is possible to designate a subset of variables as adjustment variables that always stay in the model so that relative importance is only 
  assessed among the remaining variables.
}
\references{ 
Chevan, A. and Sutherland, M. (1991) Hierarchical Partitioning. \emph{The American Statistician} \bold{45}, 90--96.

Darlington, R.B. (1968) Multiple regression in psychological research and practice. \emph{Psychological Bulletin}  \bold{69}, 161--182.

Feldman, B. (2005) Relative Importance and Value. Manuscript (Version 1.1, March 19 2005), downloadable at \url{http://www.prismanalytics.com/docs/RelativeImportance050319.pdf}

Lindeman, R.H., Merenda, P.F. and Gold, R.Z. (1980) \emph{Introduction to Bivariate and Multivariate Analysis}, Glenview IL: Scott, Foresman.

Go to \url{http://www.tfh-berlin.de/~groemp} for further information and references.}
\section{Acknowlegment }{ This package uses as an internal function the function \code{nchoosek} from \pkg{vsn}, authored by Wolfgang Huber, available under LGPL. }
\author{ Ulrike Groemping, TFH Berlin }
\note{ There are two versions of this package. The version on CRAN is globally licensed under GPL version 2 (or later). 
There is an extended version with the interesting additional metric \bold{pmvd} that is licensed according to GPL version 2
under the geographical restriction "outside of the US" because of potential issues with US patent 6,640,204. This version can be obtained 
from Ulrike Groempings website (cf. references section). Whenever you load the package, a display tells you, which version you are loading. }

\section{Warning }{\code{lmg} and \code{pmvd} are computer-intensive. Although they are calculated based on the 
covariance matrix, which saves substantial computing time in comparison to carrying out actual regressions, 
these methods still take quite long for problems with many regressors. Obviously, this is particularly relevant in combination with bootstrapping.
} 

\seealso{ See Also \code{\link{calc.relimp}}, \code{\link{booteval.relimp}}, \code{\link{classesmethods.relaimpo}}, package \code{hier.part} }
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ htest }% __ONLY ONE__ keyword per line

