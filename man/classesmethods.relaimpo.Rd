\name{classesmethods.relaimpo}
\alias{classesmethods.relaimpo}
\alias{relimplm}
\alias{relimplmtest}
\alias{relimplmboot}
\alias{relimplmbooteval}
\alias{relimplmbootMI}
\alias{plot.relimplm}
\alias{print.relimplm}
\alias{$.relimplm}
\alias{plot.relimplmbooteval}
\alias{print.relimplmbooteval}
\alias{$.relimplmbooteval}
\alias{plot.relimplmbootMI}
\alias{print.relimplmbootMI}
\alias{$.relimplmbootMI}
\alias{summary.relimplmbootMI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Classes and Methods in Package relaimpo }
\description{
  Output objects from package \pkg{relaimpo} have classes \code{relimplm} (output from calc.relimp), \code{relimplmboot} 
  (output from boot.relimp), \code{relimplmbooteval} (output from booteval.relimp) or \code{relimplmbootMI}. 
  For classes \code{relimplm}, \code{relimplmbooteval} and \code{relimplmbootMI}, there are methods for plotting and printing, 
  usage of which is described below. For class \code{relimplmbootMI}, there is in addition a summary-method,
  which produces a less detailed output than the show / print - method.
  For classes \code{relimplm}, \code{relimplmbooteval} and \code{relimplmbootMI}, 
  there is in addition a method for extracting slots of the class with \code{$}. 

}
\usage{
\method{print}{relimplm}(x,...,show.coeffs = ifelse(any(c("lmg", "pmvd") \%in\% x@type) 
    & is.null(x@always), TRUE, FALSE))

\method{plot}{relimplm}(x,...,names.abbrev=4, ylim=NULL, main=NULL, cex.title=1.5)

\method{print}{relimplmbooteval}(x,...)                                    

\method{print}{relimplmbootMI}(x,...)                                    

\method{plot}{relimplmbooteval}(x, ..., lev=max(x@level), names.abbrev=4, ylim=NULL, 
    main=NULL, cex.title=1.5)   

\method{plot}{relimplmbootMI}(x, ..., lev=max(x@level), names.abbrev=4, ylim=NULL, 
    main=NULL, cex.title=1.5)   

\method{summary}{relimplmbootMI}(object, ..., lev = max(object@level))                                    

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ \code{x} is an output object from package \pkg{relaimpo} of the required class }
  \item{show.coeffs}{ \code{show.coeffs}, if set to \code{FALSE}, suppresses printing of averaged coefficients, 
           which are otherwise printed, if \code{lmg} and/or \code{pmvd} are among the metrics in \code{x} and 
           \code{always} is \code{NULL}. 
           See details for a discussion of the averaged coefficients, if \code{always} is not \code{NULL}.}
  \item{object}{ \code{object} is an object of class \code{relimplmbootMI} 
             (output from function \code{mianalyze.relimp}) }
  \item{...}{ further arguments to functions }
  \item{names.abbrev}{ \code{names.abbrev} is an integer that provides the number of characters to which the bar labels 
             are shortened (default: 4). }
  \item{ylim}{The plot routines try to use appropriate scaling. If adjustments are needed, \code{ylim} can be used 
             like usually on plot. }
  \item{main}{The plot routine uses a default title based on the reponses name. 
             If adjustments are desired, \code{main} can be used for specifying a different title.
             Note that only the first title is affected 
             (in case of plot.relimpbooteval, there is also a sub title that cannot be changed). }
  \item{cex.title}{ \code{cex.title} specifies the text size for the overall title. 
             Thus, the \code{par} option \code{cex.main} can be used for specifying the size of individual plot titles.}
  \item{lev}{ \code{lev} is a numeric that provides the confidence level to be plotted 
             or displayed respectivly (default: maximum available confidence level; 0.5<lev<1  ). }
}
\details{
This documentation part describes S3 methods. In addition there are S4 methods for \code{show} which
coincide with the S3 methods for \code{print} and an S4 method for coercing objects of \code{relimplm} to lists
(of their numeric elements).

Print (and show) methods produce annotated output for \code{calc.relimp}, \code{booteval.relimp}, and 
\code{mianalyze.relimp} (or the objects produced by these functions). Since version 2.1, \code{calc.relimp} 
provides averaged coefficients for different sub model sizes (slot \code{ave.coeffs} of class \code{relimplm}), 
if metrics based on averaging over orderings (\code{lmg} and/or \code{pmvd}) are calculated. 
These are per default printed if the slot \code{always} of \code{x} is \code{NULL}. 
If some variables were forced into all models (non-NULL \code{always}), the averaged coefficients 
refer to the adjusted model after taking residuls from regressions on the \code{always}-columns 
of the X-matrix for both response and the other columns of the X-matrix. The reason is that these could be 
easily and cheaply implemented into the existing code and do correspond to sub models relevant 
for \code{lmg} and \code{pmvd}. Users who are interested in these coefficients, can set option \code{show.coeffs=TRUE} 
in spite of non-NULL \code{always}. 

The plot methods produce barplots of relative contributions, 
either of the metrics alone for output objects of class \code{relimplm} from function \code{calc.relimp}, 
or of the metrics with lines indicating confidence intervals for output objects of class \code{relimplmbooteval} 
from function \code{booteval.relimp} or \code{relimplmbootMI} from function \code{booteval.relimp}.   

Most \code{par()} options can be set and should work on plot. 
Exceptions: mfrow, oma and mar are set by the plot function, 
depending on the number of metrics to plot and the amount of annotating text required. 

The summary-method for class \code{relimplmbootMI} allows to quickly display brief output and to change 
the confidence level versus the level used in the original run (with interval bounds stored in the 
\dQuote{metric}.lower and \dQuote{metric}.upper slots and displayed by print and show methods).

Because of a defined S3-extraction method, 
slots of classes \code{relimplm}, \code{relimplmbooteval} and \code{relimplmbootMI} can be extracted not only 
with the \code{@} extractor but also with \code{$}. 
Hence, output elements from functions \code{calc.relimp}, \code{booteval.relimp}, and \code{mianalyze.relimp}
can be extracted as though the output objects were lists.

Note that there also is an an internally-used class relimplmtest that permits the internal function 
calc.relimp_default.intern to output further detail needed for usage from within other funtions. 

}

\references{ Go to \url{https://prof.bht-berlin.de/groemping/} for further information and references.}
\author{ Ulrike Groemping, BHT Berlin }

\seealso{ \pkg{\link{relaimpo}}, \code{\link{calc.relimp}}, \code{\link{booteval.relimp}}, 
       \code{\link{mianalyze.relimp}}, 
       \code{\link{relimplm-class}}, \code{\link{relimplmboot-class}}, 
       \code{\link{relimplmbooteval-class}}, \code{\link{relimplmbootMI-class}}    }
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ htest }% __ONLY ONE__ keyword per line
