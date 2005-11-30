\name{booteval.relimp}
\alias{boot.relimp}
\alias{booteval.relimp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Functions to Bootstrap Relative Importance Metrics }
\description{
  These functions provide bootstrap confidence intervals for relative importances. 
	\code{boot.relimp} uses the R package boot to do the actual bootstrapping of requested metrics 
      (which may take quite a while), 
	while \code{booteval.relimp} evaluates the results and provides confidence intervals.
      Output from \code{booteval.relimp} is printed with a tailored print method, and a plot method
      produces bar plots with confidence indication of the relative importance metrics. 
}
\usage{
boot.relimp(y, x, b = 1500, type = "lmg", 
    rank = TRUE, diff = TRUE, rela = TRUE)

booteval.relimp(bootrun, bty = "bca", level = 0.95, 
    sort = FALSE, norank = FALSE, nodiff = FALSE, 
    typesel = c("lmg", "pmvd", "last", "first", "betasq", "pratt"))
}
%- end of 'usage'
\arguments{
  \item{y}{ \code{y} is a vector of response data }
  \item{x}{ \code{x} is a matrix or data frame of regressor variables }
  \item{b}{ \code{b} is the number of bootstrap runs requested on boot.relimp (default: b=1500). 
       Make sure to set this to a lower number, if you are simply testing code. }
  \item{type}{ \code{type} is a character vector requesting metrics that are to be calculated. 
       Available metrics: \code{lmg}, \code{pmvd} (non-US version only), \code{last}, \code{first}, \code{betasq}, \code{pratt},
       described in \code{\link{calc.relimp}} }
  \item{rank}{ \code{rank} is a logical requesting bootstrapping of ranks (\code{rank=TRUE}, default) for each metric from type }
  \item{diff}{ \code{diff} is a logical requesting bootstrapping of pairwise differences in relative importance (\code{diff=TRUE}, default) for each metric in type }
  \item{rela}{ \code{rela} is a logical requesting relative importances summing to 100% (\code{rela=TRUE}, default). }
If rela is FALSE, some of the metrics sum to \eqn{R^2} (\code{lmg}, \code{pmvd}, \code{pratt}), others do not have a meaningful sum (\code{last}, \code{first}, \code{betasq}). More detail is given in \code{\link{calc.relimp}}.}
  \item{bootrun}{ \code{bootrun} is an object of class \code{relimplmboot} created by function \code{boot.relimp}. It hands over all relevant information on the bootstrap runs to function \code{booteval.relimp}. }
  \item{bty}{ \code{bty} is the type of bootstrap interval requested, as handed over to the function \code{boot.ci} from package 
       \code{boot}. Possible choices are \code{bca}, \code{perc}, \code{basic} and \code{norm}. \code{student} is not supported.  }
  \item{level}{ \code{level} is a single confidence level or a numeric vector of confidence levels. }
  \item{sort}{ \code{sort} is a logical requesting output sorted by size of relative contribution (\code{sort=TRUE}) or by variable position in list (\code{sort=FALSE}, default). }
  \item{norank}{ \code{norank} is a logical that suppresses of rank letters (\code{norank=TRUE}) even if ranks have been bootstrapped. }
  \item{nodiff}{ \code{nodiff} is a logical that suppresses output of confidence intervals for differences (\code{nodiff=TRUE}), even if differences have been bootstrapped. }
  \item{typesel}{ \code{typesel} provides the metrics that are to be reported. Default: all available ones (intersection of those available in object \code{bootrun} and those requested in \code{typesel})}
}
%- end of arguments

\details{
Calculations of metrics are based on the function \code{calc.relimp}. 
Bootstrapping is done with the R package \code{boot}, 
resampling the full observation vectors (combinations of response and regressors, cf. Fox (2002). 

The output provides results for all selected relative importance metrics. 
The output object can be printed and plotted (description of syntax: \code{\link{classesmethods.relaimpo}}).   

Printed output: In addition to the standard output of \code{calc.relimp} 
(one row for each regressor, one column for each bootstrapped metric), 
there is a table of confidence intervals for each selected metric 
(one row per combination of regressor and metric). 
This table is enhanced by information on rank confidence intervals, 
if ranks have been bootstrapped (\code{rank=TRUE}) and \code{norank=FALSE}. 
In addition, if differences have been bootstrapped (\code{diff=TRUE}) and \code{nodiff=FALSE}, 
there is a table of estimated pairwise differences with confidence intervals.

Graphical output: Application of the plot method to the object created by \code{booteval.relimp} 
yields barplot representations for all bootstrapped metrics (all in one graphics window). Confidence level (\code{lev=}) 
and number of characters in variable names to be used (\code{names.abbrev=}) can be modified. 
Confidence bounds are indicated on the graphs by added vertical lines. 
\code{par()} options can be used for modifying output (exception: mfrow is overridden by the plot method).
}
%- end of details
\value{
  The value of \code{boot.relimp} is of class \code{relimplmboot}. It is designed to be useful as input for \code{booteval.relimp} and is not further described here. 
  \code{booteval.relimp} returns an object of class \code{relimplmbooteval}, the items of which can be accessed by 
  the \code{$} or the \code{@} extractors. 

	In addition to the items described for function \code{calc.relimp}, which are also available here, 
the following items may be of interest for further calculations:
  \item{metric.lower }{matrix of lower confidence bounds for \dQuote{metric}: one row for each confidence level, 
        one column for each element of \dQuote{metric}. \dQuote{metric} can be any of \code{lmg}, \code{lmg.rank}, \code{lmg.diff}, ... 
       (replace \code{lmg} with other available relative importance metrics, cf. \code{\link{calc.relimp}})}
  \item{metric.upper }{matrix of upper confidence bounds for \dQuote{metric}: one row for each confidence level, 
        one column for each element of \dQuote{metric}}
  \item{metric.boot }{matrix of bootstrap results for \dQuote{metric}: one row for each bootstrap run, 
       one column for each element of \dQuote{metric}. 
       Here, \dQuote{metric} can be chosen as any of the above-mentioned and also as \eqn{R^2}}
  \item{nboot }{number of bootstrap runs underlying the evaluations}
  \item{level }{confidence levels}
}
%- end of value
\references{
Chevan, A. and Sutherland, M. (1991) Hierarchical Partitioning. \emph{The American Statistician} \bold{45}, 90--96.

Feldman, B. (2005) Relative Importance and Value. Manuscript (Version 1, March 8 2005), downloadable at \url{http://www.qwafafew.org/?q=filestore/download/268}

Fox, J. (2002) Bootstrapping regression models. \emph{An R and S-PLUS Companion to Applied Regression: A web appendix to the book}.  	
\url{http://cran.r-project.org/doc/contrib/Fox-Companion/appendix-bootstrapping.pdf}.

Lindeman, R.H., Merenda, P.F. and Gold, R.Z. (1980) \emph{Introduction to Bivariate and Multivariate Analysis}, Glenview IL: Scott, Foresman.

Go to \url{http://www.tfh-berlin.de/~groemp} for further information and references.
}
\author{ Ulrike Groemping, TFH Berlin }
\note{ There are two versions of this package. The version on CRAN is globally licensed under GPL version 2 (or later). 
There is an extended version with the interesting additional metric \bold{pmvd} that is licensed according to GPL version 2
under the geographical restriction "outside of the US" because of potential issues with US patent 6,640,204. This version can be obtained 
from Ulrike Groempings website (cf. references section). Whenever you load the package, a display tells you, which version you are loading. }

\section{Warning }{The bootstrap confidence intervals should be used for exploratory purposes only. 
They are somewhat liberal: Simulations have shown that non-coverage probabilities can be twice the 
nominal probabilities. More investigations are needed.

Be aware that the method itself needs some computing time in case of many regressors. 
Hence, bootstrapping should be used with awareness of computing time issues.} 

\seealso{ See also \code{\link{calc.relimp}}, \code{\link{classesmethods.relaimpo}}  }
\examples{
#####################################################################
### Example: relative importance of various socioeconomic indicators 
###          for Fertility in Switzerland
### Fertility is first column of data set swiss
#####################################################################
data(swiss)
   # bootstrapping
       bootswiss <- boot.relimp(swiss[,1], swiss[,2:6], b = 100,  
                type = c("lmg", "last", "first", "pratt"),
                rank = TRUE, diff = TRUE, rela = TRUE)
       # for demonstration purposes only 100 bootstrap replications

   #default output
    booteval.relimp(bootswiss)
         #because of only 100 bootstrap replications, 
         #default bca intervals produce warnings
    plot(booteval.relimp(bootswiss))

    #sorted printout, chosen confidence levels, chosen interval method
    #store as object
        result <- booteval.relimp(bootswiss, bty="perc", 
              sort = TRUE, level=c(0.8,0.9))
    #output driven by print method
        result
    #result plotting with default settings 
    #(largest confidence level, names abbreviated to length 4)
        plot(result)
    #result plotting with modified settings (chosen confidence level, 
    #names abbreviated to chosen length)
        plot(result, level=0.8,names.abbrev=5)
    #result plotting with longer names shown vertically
        par(las=2)
        plot(result, level=0.9,names.abbrev=6)
    #plot does react to options set with par()
    #exception: mfrow-option is set within plot, 
    #depending on the number of metrics to be shown
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ htest }% __ONLY ONE__ keyword per line
