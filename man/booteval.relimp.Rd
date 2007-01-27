\name{booteval.relimp}
\alias{boot.relimp}
\alias{boot.relimp.default}
\alias{boot.relimp.formula}
\alias{boot.relimp.lm}
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
## generic function
boot.relimp(object, ...)

## default S3 method, should be called without suffix ".default"
boot.relimp.default(object, x = NULL, ..., b = 1000, type = "lmg", 
    rank = TRUE, diff = TRUE, rela = FALSE, always = NULL, 
    groups = NULL, groupnames = NULL, fixed=FALSE)

## S3 method for formula objects, should be called without suffix ".formula"
boot.relimp.formula(formula, data, na.action, ..., subset = NULL)

## S3 method for objects of class lm, should be called without suffix ".lm"
boot.relimp.lm(object, ...)

## function for evaluating bootstrap results
booteval.relimp(bootrun, bty = "bca", level = 0.95, 
    sort = FALSE, norank = FALSE, nodiff = FALSE, 
    typesel = c("lmg", "pmvd", "last", "first", "betasq", "pratt"))
}
%- end of 'usage'
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
  \item{b}{ is the number of bootstrap runs requested on boot.relimp (default: \code{b=1000}). 
       Make sure to set this to a lower number, if you are simply testing code. }
  \item{type}{ can be a character string, character vector or list of character strings.
       It is the collection of metrics that are to be calculated. 
       Available metrics: \code{lmg}, \code{pmvd} (non-US version only), \code{last}, \code{first}, 
       \code{betasq}, \code{pratt}, as
       described in \code{\link{calc.relimp}}.}
  \item{rank}{ is a logical requesting bootstrapping of ranks (\code{rank=TRUE}, default) for each metric from type }
  \item{diff}{ is a logical requesting bootstrapping of pairwise differences in relative importance (\code{diff=TRUE}, 
         default) for each metric in type }
  \item{rela}{ is a logical requesting relative importances summing to 100% (\code{rela=TRUE}). 
         If rela is FALSE (default), some of the metrics sum to \eqn{R^2} (\code{lmg}, \code{pmvd}, \code{pratt}), 
         others do not have a meaningful sum (\code{last}, \code{first}, \code{betasq}). 
         More detail is given in \code{\link{calc.relimp}}.}
  \item{always }{ is a vector of column numbers or names of variables to be always in the 
       model (adjusted for). Valid numbers are 2 to \dQuote{number of regressors + 1} (1 is reserved for the response), 
       valid character strings are all column names of \code{y} or \code{x} respectively that refer to regressor variables.

       Relative importance is only assessed for the variables not selected in \code{always}. }
  \item{groups }{ is a list of vectors of column numbers or names of variables to be combined into groups. 
       If only one group is needed, a vector can be given. The numbers and character strings needed are of the same form 
       as for \code{always}. Numbers and names cannot be mixed. 

       Relative importance is only allocated between groups of regressors, no subdivision within groups is calculated. 
       Regressors that do not occur in any group are included as singletons.
       A regressor must not occur in \code{always} and in \code{groups}.}
  \item{groupnames }{ is a vector of names for the variable groups to be used for annotation of output.}
  \item{fixed}{ is a logical requesting bootstrapping for a fixed design matrix (if TRUE). 
       The default is bootstrapping for randomly drawn samples (fixed = FALSE). }
  \item{data}{ if first object is of class formula: 
       an optional matrix or data frame that the variables in formula and subset come from; 
       if it is omitted, all names must be meaningful in the environment from which boot.relimp is called}
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
  \item{bootrun}{ is an object of class \code{relimplmboot} created by function \code{boot.relimp}. 
         It hands over all relevant information on the bootstrap runs to function \code{booteval.relimp}. }
  \item{bty}{ is the type of bootstrap interval requested (a character string), 
        as handed over to the function \code{boot.ci} from package \code{boot}. 
        Possible choices are \code{bca}, \code{perc}, \code{basic} and \code{norm}. 
        \code{student} is not supported.  }
  \item{level}{ is a single confidence level or a numeric vector of confidence levels. }
  \item{sort}{ is a logical requesting output sorted by size of relative contribution (\code{sort=TRUE}) or by variable position in list (\code{sort=FALSE}, default). }
  \item{norank}{ is a logical that suppresses of rank letters (\code{norank=TRUE}) even if ranks have been bootstrapped. }
  \item{nodiff}{ is a logical that suppresses output of confidence intervals for differences (\code{nodiff=TRUE}), even if differences have been bootstrapped. }
  \item{typesel}{ provides the metrics that are to be reported. Default: all available ones 
       (intersection of those available in object \code{bootrun} and those requested in \code{typesel}). 
       \code{typesel} accepts the same values as \code{type}.}
}
%- end of arguments

\details{
Calculations of metrics are based on the function \code{calc.relimp}. 
Bootstrapping is done with the R package \code{boot}, 
resampling the full observation vectors (combinations of response and regressors, cf. Fox (2002)). 

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
\code{par()} options can be used for modifying output (exceptions: mfrow, oma and mar are overridden by the plot method).
}
%- end of details
\value{
  The value of \code{boot.relimp} is of class \code{relimplmboot}. It is designed to be useful as input for \code{booteval.relimp} and is not further described here. 
  \code{booteval.relimp} returns an object of class \code{relimplmbooteval}, the items of which can be accessed by 
  the \code{$} or the \code{@} extractors. 

	In addition to the items described for function \code{\link{calc.relimp}}, which are also available here, 
the following items may be of interest for further calculations:
  \item{metric.lower }{matrix of lower confidence bounds for \dQuote{metric}: one row for each confidence level, 
        one column for each element of \dQuote{metric}. \dQuote{metric} can be any of \code{lmg}, \code{lmg.rank}, 
        \code{lmg.diff}, ... 
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

Darlington, R.B. (1968) Multiple regression in psychological research and practice. \emph{Psychological Bulletin}  \bold{69}, 161--182.

Feldman, B. (2005) Relative Importance and Value. Manuscript (Version 1.1, March 19 2005), downloadable at \url{http://www.prismanalytics.com/docs/RelativeImportance050319.pdf}

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
They can be somewhat liberal: Limited simulations for percentile intervals have shown that non-coverage probabilities 
can be up to twice the nominal probabilities. More investigations are needed.

Be aware that the method itself needs some computing time in case of many regressors. 
Hence, bootstrapping should be used with awareness of computing time issues.

\code{relaimpo} is a package for univariate linear models. 
Using \code{relaimpo} on objects that inherit from class \code{lm} but are not univariate linear model objects 
may produce nonsensical results without warning. Objects of classes \code{mlm} and \code{glm} lead to an error message. 
} 

\seealso{ See also \pkg{\link{relaimpo}}, \code{\link{calc.relimp}}, \code{\link{classesmethods.relaimpo}}  }
\examples{
#####################################################################
### Example: relative importance of various socioeconomic indicators 
###          for Fertility in Switzerland
### Fertility is first column of data set swiss
#####################################################################
data(swiss)
   # bootstrapping
       bootswiss <- boot.relimp(swiss, b = 100,  
                type = c("lmg", "last", "first", "pratt"),
                rank = TRUE, diff = TRUE, rela = TRUE)
       # for demonstration purposes only 100 bootstrap replications

       #alternatively, use formula interface
       bootsub <- boot.relimp(Fertility~Education+Catholic+Infant.Mortality, swiss, 
              subset=Catholic>40, b = 100, type = c("lmg", "last", "first", "pratt"),
              rank = TRUE, diff = TRUE)
       # for demonstration purposes only 100 bootstrap replications

   #default output
    booteval.relimp(bootswiss)
         #because of only 100 bootstrap replications, 
         #default bca intervals produce warnings
    plot(booteval.relimp(bootswiss))

    #sorted printout, chosen confidence levels, chosen interval method
    #store as object
        result <- booteval.relimp(bootsub, bty="perc", 
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
    #exceptions: mfrow, mar and oma are set within the plot routine itself
}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ htest }% __ONLY ONE__ keyword per line
