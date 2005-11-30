\name{calc.relimp}
\alias{calc.relimp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to calculate relative importance metrics for linear models }
\description{
  calc.relimp calculates several relative importance metrics for the linear model. 
The recommended metrics are \code{lmg} (\eqn{R^2} partitioned by averaging over orders, like in Lindemann, Merenda and Gold (1980, p.119ff)) 
and \code{pmvd} (a newly proposed metric by Feldman (2005) that is provided in the non-US version of the package only). 
For completeness and comparison purposes, several other metrics are also on offer (cf. e.g. Darlington (1968)).
}
\usage{
calc.relimp(covg, type = "lmg", diff = FALSE, rank = TRUE, rela = TRUE)
}
\arguments{
  \item{covg   }{ is the covariance matrix of a response y and regressors x, 
	e.g. obtained by cov(cbind(y,x)), if y is a column vector of response values 
	and x a corresponding matrix of regressors }
  \item{type   }{ is the list of selected metrics. On offer are \code{lmg}, \code{pmvd} (non-US version only), \code{last}, \code{first}, \code{betasq} and \code{pratt}. For brief sketches of their meaning cf. Details section.}
  \item{diff   }{ logical; if TRUE, pairwise differences between the relative contributions are calculated; default FALSE }
  \item{rank   }{ logical; if TRUE, ranks of regressors in terms of relative contributions are calculated; default TRUE }
  \item{rela   }{ logical; if TRUE, all metrics are forced to sum to 100pct; if FALSE, details depend on specific method; defaul TRUE }
}
\details{
\describe{
  \item{lmg}{ is the \eqn{R^2} contribution averaged over orders among regressors, cf. e.g. Lindeman, Merenda and Gold 1980, p.119ff or Chevan and Sutherland (1991).} 
  \item{pmvd}{ is the proportional marginal variance decomposition as proposed by Feldman (2005) (non-US version only).} 
  \item{last}{ is each variables contribution when included last, also sometimes called usefulness.} 
  \item{first}{ is each variables contribution when included first, which is just the squared covariance between y and the variable.}
  \item{betasq}{ is the squared standardized coefficient.} 
  \item{pratt}{ is the product of the standardized coefficient and the correlation.}
}
Each metric is calculated using the internal function \dQuote{metric}\code{calc}, e.g. \code{lmgcalc}.

\code{lmg}, \code{pmvd} and \code{pratt} sum to \eqn{R^2}, if \code{rela = FALSE} and to 100pct if \code{rela = TRUE}.

The other metrics are given relative to var(y) but do not sum to \eqn{R^2} if \code{rela = FALSE}. 

If \code{rela = TRUE}, they are artificially forced to sum to 100pct.
}
\value{
  \item{R2 }{the coefficient of determination, \eqn{R^2}}
  \item{lmg }{vector of relative contributions obtained from the \code{lmg} method, if \code{lmg} has been requested in \code{type}}
  \item{lmg.diff }{vector of pairwise differences between relative contributions obtained from the \code{lmg} method, 
if \code{lmg} has been requested in \code{type} and \code{diff=TRUE}}
  \item{lmg.rank }{rank of the regressors relative contributions obtained from the \code{lmg} method, 
if \code{lmg} has been requested in \code{type} and \code{rank=TRUE}}
  \item{metric, metric.diff, metric.rank }{analogous to \code{lmg} for other metrics}
}
\references{ 
Chevan, A. and Sutherland, M. (1991) Hierarchical Partitioning. \emph{The American Statistician} \bold{45}, 90--96.

Darlington, R.B. (1968) Multiple regression in psychological research and practice. \emph{Psychological Bulletin}  \bold{69}, 161--182.

Feldman, B. (2005) Relative Importance and Value. Manuscript (Version 1, March 8 2005), downloadable at \url{http://www.qwafafew.org/?q=filestore/download/268}

Lindeman, R.H., Merenda, P.F. and Gold, R.Z. (1980) \emph{Introduction to Bivariate and Multivariate Analysis}, Glenview IL: Scott, Foresman.

Go to \url{http://www.tfh-berlin.de/~groemp} for further information and references.
}
\author{ Ulrike Groemping, TFH Berlin }
\note{ There are two versions of this package. The version on CRAN is globally licensed under GPL version 2 (or later). 
There is an extended version with the interesting additional metric \bold{pmvd} that is licensed according to GPL version 2
under the geographical restriction "outside of the US" because of potential issues with US patent 6,640,204. This version can be obtained 
from Ulrike Groempings website (cf. references section). Whenever you load the package, a display tells you, which version you are loading. }

\section{Warning }{\code{lmg} and \code{pmvd} are computer-intensive. Although they are calculated based on the 
covariance matrix, which saves substantial computing time in comparison to carrying out actual regressions, 
these methods still take quite long for problems with many regressors.
} 

\seealso{ See Also \code{\link{booteval.relimp}}, \code{\link{classesmethods.relaimpo}} }
\examples{
#####################################################################
### Example: relative importance of various socioeconomic indicators 
###          for Fertility in Switzerland
### Fertility is first column of data set swiss
#####################################################################
data(swiss)
    calc.relimp(cov(swiss), 
       type = c("lmg", "last", "first", "betasq", "pratt"), rela = TRUE )
    # calculation of all available relative importance metrics 
        # non-US version offers the additional metric "pmvd", 
        # i.e. call would be 
        # calc.relimp(cov(swiss), 
        # type = c("lmg", "pmvd", "last", "first", "betasq, "pratt"), 
        # rela = TRUE )
    plot(calc.relimp(cov(swiss), 
        type = c("lmg", "last", "first", "betasq", "pratt"), rela = TRUE ))
    # bar plot of the relative importance metrics

   #of statistical interest in this context: correlation matrix
       cor(swiss)
  }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line