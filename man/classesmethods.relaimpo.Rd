\name{classesmethods.relaimpo}
\alias{classesmethods.relaimpo}
\alias{relimplm}
\alias{relimplmboot}
\alias{relimplmbooteval}
\alias{plot.relimplm}
\alias{print.relimplm}
\alias{$.relimplm}
\alias{$.relimplmbooteval}
\alias{plot.relimplmbooteval}
\alias{print.relimplmbooteval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Classes and Methods in Package relaimpo }
\description{
  Output objects from package \pkg{relaimpo} have classes \code{relimplm} (output from calc.relimp), \code{relimplmboot} 
  (output from boot.relimp) and \code{relimplmbooteval} (output from booteval.relimp). 
  For classes \code{relimplm} and \code{relimplmbooteval}, there are methods for plotting and printing, 
  usage of which is described below. 
  For class \code{relimplm}, there is in addition a method for extracting slots of the class with \$. 
  Since class \code{relimplmbooteval} inherits from \code{relimplm}, this method can also be applied there. 
}
\usage{
\method{print}{relimplm}(x,...)

\method{plot}{relimplm}(x,...,names.abbrev=4)

\method{print}{relimplmbooteval}(x,...)                                    

\method{plot}{relimplmbooteval}(x, ..., lev=max(x@level), names.abbrev=4)   
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ \code{x} is an output object from package \pkg{relaimpo} of the required class }
  \item{...}{ \code{...} further arguments to function \code{plot} }
  \item{names.abbrev}{ \code{names.abbrev} is an integer that provides the number of characters to which the bar labels are shortened (default: 4). }
  \item{lev}{ \code{lev} is a numeric that provides the confidence level to be plotted (default: maximum available confidence level). }
}
\details{
Print methods produce annotated output for \code{calc.relimp} and \code{booteval.relimp}. 
The plot methods produce barplots of relative contributions, 
either of the metrics alone for output objects of class \code{relimplm} from function \code{calc.relimp}, 
or of the metrics with lines indicating confidence intervals for output objects of class \code{relimplmbooteval} from function \code{booteval.relimp}.   

All \code{par()} options can be set and should work on plot. Exception: mfrow is set by the plot function, 
depending on the number of metrics to plot. 

In addition, there is a method that allows using \$ for extracting slots of classes \code{relimplm} 
and \code{relimplmbooteval}, i.e. elements from results of functions \code{calc.relimp} and \code{booteval.relimp} 
can be extracted as though the output objects were lists.
}
\references{ Go to \url{http://www.tfh-berlin.de/~groemp} for further information and references.}
\author{ Ulrike Groemping, TFH Berlin }

\seealso{ See also \code{\link{calc.relimp}}, \code{\link{booteval.relimp}}, \pkg{\link{relaimpo}}    }
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ htest }% __ONLY ONE__ keyword per line
