\name{calc.relimp}
\alias{calc.relimp}
\alias{calc.relimp.default}
\alias{calc.relimp.formula}
\alias{calc.relimp.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to calculate relative importance metrics for linear models }
\description{
  calc.relimp calculates several relative importance metrics for the linear model. 
The recommended metrics are \code{lmg} (\eqn{R^2} partitioned by averaging over orders, like in Lindemann, Merenda and Gold (1980, p.119ff)) 
and \code{pmvd} (a newly proposed metric by Feldman (2005) that is provided in the non-US version of the package only). 
For completeness and comparison purposes, several other metrics are also on offer (cf. e.g. Darlington (1968)).
}
\usage{

## generic function
calc.relimp(object, ...)

## default S3 method, should be called without suffix ".default"
calc.relimp.default(object, x = NULL, ..., 
       type = "lmg", diff = FALSE, rank = TRUE, rela = FALSE, always = NULL)

## S3 method for formula object, should be called without suffix ".formula"
calc.relimp.formula(formula, data, na.action, ..., subset=NULL)

## S3 method for objects of class lm
calc.relimp.lm(object, ...)

}
\arguments{
  \item{object   }{ The class of this object determines which of the methods is used:
        There are special methods for output objects from function \code{\link{lm}} and for formula objects.
        For all other types of object, the default method is used.

        Thus, object can be

        a formula (e.g. y~x1+x2+x3) without interaction terms and without factors

        OR

        the output of a linear model call (class \code{lm}, but not \code{glm} or \code{mlm}); 
        output objects from \code{lm} or \code{aov} (without factors in x variables) work;
        there may be further functions that output objects inheriting from \code{lm} 
        which may or may not work reasonably with \code{calc.relimp}; 
        for \code{calc.relimp} to be appropriate, the underlying model must at least be linear!

        OR 

        the covariance matrix of a response y and regressors x, 
	(e.g. obtained by cov(cbind(y,x)), if y is a column vector of response values 
	and x a corresponding matrix of regressors) 

        OR 
        
        a (raw) data matrix or data frame with the response variable in the first column 
        (numeric variables only, no factors)


        OR

        a response vector or one-column matrix, 
        if \code{x} contains the corresponding matrix or data frame of regressors.}

  \item{formula   }{ The first object, if a formula is to be given; one response, no factors, and no interaction terms }
  \item{x   }{ a (raw) data matrix or data frame containing the regressors (no factors), 
       if \code{object} is a response vector or one-column matrix

        OR 
        
        NULL, if \code{object} is anything else }
  \item{type   }{ can be a character string, character vector or list of character strings.
       It is the collection of metrics that are to be calculated. 
       Available metrics: \code{lmg}, \code{pmvd} (non-US version only), \code{last}, \code{first}, 
       \code{betasq}, \code{pratt}. For brief sketches of their meaning cf. details section. }
  \item{diff   }{ logical; if TRUE, pairwise differences between the relative contributions are calculated; default FALSE }
  \item{rank   }{ logical; if TRUE, ranks of regressors in terms of relative contributions are calculated; default TRUE }
  \item{rela   }{ logical; if TRUE, all metrics are forced to sum to 100pct; if FALSE, details depend on specific method; default FALSE }
  \item{always }{ is a vector of column numbers or names of variables to be always in the 
       model (adjusted for). Valid numbers are 2 to (number of regressors + 1) (1 is reserved for the response), 
       valid character strings are all column names of \code{y} or \code{x} respectively that refer to regressor variables.

       Relative importance is only assessed for the variables not selected in \code{always}. }
  \item{data}{ if first object is of class formula: 
       an optional matrix or data frame that the variables in formula and subset come from; 
       if it is omitted, all names must be meaningful in the environment from which calc.relimp is called}
  \item{subset}{ if first object is of class formula: 
        an optional expression indicating the subset of the rows of \code{data} that should be used in the fit. 
          This can be a logical vector, or a numeric vector indicating which observation numbers are to be included, 
          or a  character  vector of the row names to be included.  All (non-missing) observations are included by
          default. }
  \item{na.action}{ if first object is of class formula: 
        an optional function that indicates what should happen when the data contain 'NA's. 
        The default is first, any na.action attribute of data, second the setting given in the call to calc.relimp,
        third the na.action setting of options. Possible choices are "na.fail",  
        (print an error message and terminate if there are any incomplete observations), 
        "na.omit" or "na.exclude" (equivalent for package \code{relaimpo}, 
        both analyse complete cases only and print a warning, this is also what is done the default method ). }
  \item{...}{ usable for further arguments, particularly all arguments of default method can be given to all other methods }
}
\details{
\describe{
  \item{lmg}{ is the \eqn{R^2} contribution averaged over orderings among regressors, cf. e.g. Lindeman, Merenda and Gold 1980, 
         p.119ff or Chevan and Sutherland (1991).} 
  \item{pmvd}{ is the proportional marginal variance decomposition as proposed by Feldman (2005) (non-US version only).
         It can be interpreted as a weighted average over orderings among regressors, with data-dependent weights. } 
  \item{last}{ is each variables contribution when included last, also sometimes called usefulness.} 
  \item{first}{ is each variables contribution when included first, which is just the squared covariance between y and the variable.}
  \item{betasq}{ is the squared standardized coefficient.} 
  \item{pratt}{ is the product of the standardized coefficient and the correlation.}
}
Each metric is calculated using the internal function \dQuote{metric}\code{calc}, e.g. \code{lmgcalc}.

Three of the metrics in \code{calc.relimp} (\code{lmg}, \code{pmvd} and \code{pratt}) decompose the model \eqn{R^2}. 
If \code{always} requests some variables to be always in the model, these variables are included in the model first. 
Only the remaining \eqn{R^2} that is not explained by these variables is decomposed among the other regressors.
\code{lmg}, \code{pmvd} and \code{pratt} sum to the \eqn{R^2} that is to be decomposed, 
if \code{rela = FALSE} and to 100pct if \code{rela = TRUE}.

The other metrics also (artificially) sum to 100pct if \code{rela = TRUE}. 
If \code{rela = FALSE}, they are given relative to var(y) (or the conditional variance of y after adjusting out the variables 
requested in \code{always}) but do not sum to \eqn{R^2}.

}
\value{
  \item{var.y }{the variance of the response}
  \item{R2 }{the coefficient of determination, \eqn{R^2}}
  \item{R2.decomp }{the part of the coefficient of determination that is decomposed among the 
                variables under investigation }
  \item{lmg }{vector of relative contributions obtained from the \code{lmg} method, if \code{lmg} has been requested in \code{type}}
  \item{lmg.diff }{vector of pairwise differences between relative contributions obtained from the \code{lmg} method, 
if \code{lmg} has been requested in \code{type} and \code{diff=TRUE}}
  \item{lmg.rank }{rank of the regressors relative contributions obtained from the \code{lmg} method, 
if \code{lmg} has been requested in \code{type} and \code{rank=TRUE}}
  \item{metric, metric.diff, metric.rank }{analogous to \code{lmg} for other metrics}
  \item{namen}{names of variables, starting with response}
  \item{type}{character vector of metrics available}
  \item{rela}{Have metrics been normalized to sum 100\% ?}
  \item{always}{column numbers of variables always in the model}
  \item{alwaysnam}{names of variables always in the model}
}
\references{ 
Chevan, A. and Sutherland, M. (1991) Hierarchical Partitioning. \emph{The American Statistician} \bold{45}, 90--96.

Darlington, R.B. (1968) Multiple regression in psychological research and practice. \emph{Psychological Bulletin}  \bold{69}, 161--182.

Feldman, B. (2005) Relative Importance and Value. Manuscript (Version 1.1, March 19 2005), downloadable at \url{http://www.prismanalytics.com/docs/RelativeImportance050319.pdf}

Lindeman, R.H., Merenda, P.F. and Gold, R.Z. (1980) \emph{Introduction to Bivariate and Multivariate Analysis}, Glenview IL: Scott, Foresman.

Go to \url{http://www.tfh-berlin.de/~groemp} for further information and references.
}
\author{ Ulrike Groemping, TFH Berlin }
\note{ There are two versions of this package. The version on CRAN is globally licensed under GPL version 2 (or later). 
There is an extended version with the interesting additional metric \bold{pmvd} that is licensed according to GPL version 2
under the geographical restriction "outside of the US" because of potential issues with US patent 6,640,204. 
This version can be obtained from Ulrike Groempings website (cf. references section). 
Whenever you load the package, a display tells you, which version you are loading. }

\section{Warning }{\code{lmg} and \code{pmvd} are computer-intensive. Although they are calculated based on the 
covariance matrix, which saves substantial computing time in comparison to carrying out actual regressions, 
these methods still take quite long for problems with many regressors.

\code{relaimpo} is a package for univariate linear models. 
Using \code{relaimpo} on objects that inherit from class \code{lm} but are not univariate linear model objects 
may produce nonsensical results without warning. Objects of classes \code{mlm} and \code{glm} lead to an error message. 
} 

\seealso{ See Also \pkg{\link{relaimpo}}, \code{\link{booteval.relimp}}, \code{\link{classesmethods.relaimpo}} }
\examples{
#####################################################################
### Example: relative importance of various socioeconomic indicators 
###          for Fertility in Switzerland
### Fertility is first column of data set swiss
#####################################################################
data(swiss)
    calc.relimp(swiss, 
       type = c("lmg", "last", "first", "betasq", "pratt") )
    # calculation of all available relative importance metrics 
        # non-US version offers the additional metric "pmvd", 
        # i.e. call would be 
        # calc.relimp(cov(swiss), 
        # type = c("lmg", "pmvd", "last", "first", "betasq, "pratt"), 
        # rela = TRUE )
    ## same analysis with formula or lm method and a few modified options
    crf <- calc.relimp(Fertility~Agriculture+Examination+Education+Catholic+Infant.Mortality,swiss, 
        subset = Catholic>40,
        type = c("lmg", "last", "first", "betasq", "pratt"), rela = TRUE )
    crf
    linmod <- lm(Fertility~Agriculture+Examination+Education+Catholic+Infant.Mortality,swiss)
    crlm <- calc.relimp(linmod, 
        type = c("lmg", "last", "first", "betasq", "pratt"), rela = TRUE )
    plot(crlm)
    # bar plot of the relative importance metrics

    #of statistical interest in this context: correlation matrix
       cor(swiss)

    #demonstration of conditioning on one regressor using always
    calc.relimp(swiss, 
       type = c("lmg", "last", "first", "betasq", "pratt"), rela = FALSE,
       always = "Education" )
  }
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
