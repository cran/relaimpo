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
    groups = NULL, groupnames = NULL, fixed=FALSE, 
    weights = NULL, design = NULL)

## S3 method for formula objects, should be called without suffix ".formula"
boot.relimp.formula(formula, data, weights, na.action, ..., subset = NULL)

## S3 method for objects of class lm, should be called without suffix ".lm"
boot.relimp.lm(object, type = "lmg", groups = NULL, groupnames=NULL, always = NULL, ..., b=1000)

## function for evaluating bootstrap results
booteval.relimp(bootrun, bty = "perc", level = 0.95, 
    sort = FALSE, norank = FALSE, nodiff = FALSE, 
    typesel = c("lmg", "pmvd", "last", "first", "betasq", "pratt"))
}
%- end of 'usage'
\arguments{
  \item{object   }{ cf. \code{\link{calc.relimp}}.}

  \item{formula   }{ cf. \code{\link{calc.relimp}}.
        But note the additional restriction that - in connection with the \code{design=}-option - 
        it is currently not possible to use factors, interactions or calculated quantities in a formula.
        }
  \item{x   }{ cf. \code{\link{calc.relimp}}.}
  \item{b}{ is the number of bootstrap runs requested on boot.relimp (default: \code{b=1000}). 
       Make sure to set this to a lower number, if you are simply testing code. }
  \item{type}{ cf. \code{\link{calc.relimp}}.}
  \item{rank}{ is a logical requesting bootstrapping of ranks (\code{rank=TRUE}, default) for each metric from type }
  \item{diff}{ is a logical requesting bootstrapping of pairwise differences in relative importance (\code{diff=TRUE}, 
         default) for each metric in type }
  \item{rela}{ cf. \code{\link{calc.relimp}}.}
  \item{always }{ cf. \code{\link{calc.relimp}}.}
  \item{groups }{ cf. \code{\link{calc.relimp}}.}
  \item{groupnames }{ cf. \code{\link{calc.relimp}}.}
  \item{weights }{ cf. \code{\link{calc.relimp}} for specification of weights. 
       See also the Details section of this help page for usage of different types of weights. }
  \item{design }{ cf. \code{\link{calc.relimp}}. But note that there are 
       currently some restrictions regarding usability of other possibilities, when using a \code{design} 
       in \code{boot.relimp}: formulae can only be simpler than usual, and factors, interactions or 
       calculated variables in formulae are not permitted.
       
       For a description of the bootstrap method's treatment of designs, see the details section.
       In the current version, using a design in bootstrapping must be considered EXPERIMENTAL.}
  \item{fixed}{ is a logical requesting bootstrapping for a fixed design matrix (if TRUE). 
       The default is bootstrapping for randomly drawn samples (fixed = FALSE). }
  \item{data}{ cf. \code{\link{calc.relimp}}.}
  \item{subset}{ cf. \code{\link{calc.relimp}}.}
  \item{na.action}{ cf. \code{\link{calc.relimp}}.}
  \item{...}{ usable for further arguments, particularly all arguments of the default method can be given to all other methods }
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
resampling the full observation vectors by default (combinations of response, weights and regressors, cf. Fox (2002)). 
If \code{fixed=TRUE} is specified, bootstrapping is based on residuals rather than full observations, 
keeping the X-variables fixed.

If the \code{weights} option is used, weights are resampled together with the full observations, 
and weighted contributions are calculated for each resample (no re-normalization is done within the resamples.)

Please note that usage of weights in linear models can have very different reasons. 
One motivation is a different variability of different observations, where weights are the inverse variances. 
This is the way \code{weights} are treated in function \code{lm} and also in \code{calc.relimp} and in \code{boot.relimp}, 
if a vector of weights is given with the \code{weights} option. Specifically, weights do not affect the 
resampling probability in bootstrapping, i.e. each observation has the same probability to be included in resamples.

If the weights in a data frame represent the multiplicity of each observation (i.e. there are several units with 
identical combination of values in the data, and the weights represent the number of units with exactly this value 
pattern for each row of the data frame), they can also be directly used as weights in \code{calc.relimp} for 
calculating the metrics. However, such frequency weights cannot be appropriately accomodated in \code{boot.relimp}; 
instead, the data frame with frequencies has to be expanded to include one row for each unit before using the 
resampling routine (e.g. using function \code{untable} from package \pkg{reshape} or function \code{expand.table} 
from package \pkg{epitools}.

In survey situations, weights are used to generate a more representative picture of the population: 
an observation's weight is typically the number of units of the population that this single observed unit represents. 
In this situation, there is no reason to consider observations with higher weights as less variable than observations 
with lower weights; thus, while estimates can again be obtained treating the weights in the same way as mentioned before,
their usage in estimation of standard errors and in bootstrapping is different. 
In order to appropriately accomodate survey weights for these purposes, it is not sufficient to only provide the 
weights vector; instead, it is necessary to provide a design generated with package \pkg{survey} or an object of 
class \code{svyglm} (produced by function \code{svyglm}) that includes the appropriate design information.

Clusters are a way to take care of dependency structures like in longitudinal data. 
Thus, while relaimpo does not (currently) cover linear mixed models (e.g. produced by function \code{lme}), it is possible to 
accomodate clustered data by applying function \code{svyglm} with linear link function and gaussian distribution 
to a design that contains clusters. The bootstrapping approach subsequently takes care of the dependence by considering 
clusters as sampling units. Users who want to use this approach can mimic the second example below.

If the \code{design} option is used (experimental), resampling is done within package \pkg{survey}, 
and the resampled contributions are also calculated within package survey. The results from these calculations 
are plugged into an object from package \pkg{boot}, and confidence interval calculation is subsequently handled in 
\pkg{boot}. The approach is considered experimental: so far no simulation studies have been conducted for 
complex survey designs, and because of limited experience (in spite of thorough testing) it is not unlikely that 
bugs will be found by users who are routinely using complex survey designs.

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
\code{par()} options can be used for modifying output (exceptions: \code{mfrow}, \code{oma} and 
\code{mar} are overridden by the plot method).
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
There is an extended version with the interesting additional metric \code{pmvd} that is licensed according to GPL version 2
under the geographical restriction "outside of the US" because of potential issues with US patent 6,640,204. This version can be obtained 
from Ulrike Groempings website (cf. references section). Whenever you load the package, a display tells you, which version you are loading. }

\section{Warning }{The bootstrap confidence intervals should be used for exploratory purposes only. 
They can be somewhat liberal: Limited simulations for percentile intervals have shown that non-coverage probabilities 
can be up to twice the nominal probabilities. More investigations are needed.

Be aware that the method itself needs some computing time in case of many regressors. 
Hence, bootstrapping should be used with awareness of computing time issues.

\code{relaimpo} is a package for univariate linear models. 
Using \code{relaimpo} on objects that inherit from class \code{lm} but are not univariate linear model objects 
may produce nonsensical results without warning. Objects of class \code{mlm} or \code{glm} with link functions other than identity 
or family other than gaussian lead to an error message. 
} 

\seealso{ \pkg{\link{relaimpo}}, \code{\link{calc.relimp}}, \code{\link{mianalyze.relimp}}, 
     \code{\link{classesmethods.relaimpo}}  }
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

   #default output (percentily intervals, as of Version 2 of the package)
    booteval.relimp(bootswiss)
    plot(booteval.relimp(bootswiss))

    #sorted printout, chosen confidence levels, chosen interval method
    #store as object
        result <- booteval.relimp(bootsub, bty="bca", 
              sort = TRUE, level=c(0.8,0.9))
         #because of only 100 bootstrap replications, 
         #default bca intervals produce warnings
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

#####################################################################
### Example: bootstrapping clustered data                            
###          data taken from example.lmm of package lmm
### y is change in pulse (heart beats per minute) 
###    15 min (occ 1 to 3) and 90 min (occ 4 to 6) after 
###    treatment with Placebo (occ 1 or 4) low (occ 2 or 5) 
###    or high (occ 3 or 6) dose of marihuana
### each of 9 subjects is observed under most or all 
### of the 6 possible conditions
#####################################################################
## create example data from package lmm
y <- c(16,20,16,20,-6,-4,
    12,24,12,-6,4,-8,
    8,8,26,-4,4,8,
    20,8,20,-4,
    8,4,-8,22,-8,
    10,20,28,-20,-4,-4,
    4,28,24,12,8,18,
    -8,20,24,-3,8,-24,
    20,24,8,12)
occ <- c(1,2,3,4,5,6,
      1,2,3,4,5,6,
      1,2,3,4,5,6,
      1,2,5,6,
      1,2,3,5,6,
      1,2,3,4,5,6,
      1,2,3,4,5,6,
      1,2,3,4,5,6,
      2,3,4,5)
subj <- c(1,1,1,1,1,1,
       2,2,2,2,2,2,
       3,3,3,3,3,3,
       4,4,4,4,
       5,5,5,5,5,
       6,6,6,6,6,6,
       7,7,7,7,7,7,
       8,8,8,8,8,8,
       9,9,9,9)
### manual creation of dummies
### reference category placebo after 90min (occ=4)
dumm1 <- as.numeric(occ==1)
dumm2 <- as.numeric(occ==2)
dumm3 <- as.numeric(occ==3)
dumm5 <- as.numeric(occ==5)
dumm6 <- as.numeric(occ==6)

## create data frame
dat <- data.frame(y,dumm1,dumm2,dumm3,dumm5,dumm6,subj)

### create design with clusters
des <- svydesign(id=~subj,data=dat)

### apply bootstrapping
### using the design with subjects as clusters implies that 
###     clusters are generally kept in or excluded as a whole
### of course, b=100 is too small, only chosen for speed of package creation 
bt <- boot.relimp(y~dumm1+dumm2+dumm3+dumm5+dumm6,data=dat,
   design=des,b=100,type=c("lmg","first","last","betasq","pratt"))

### calculate and display results
booteval.relimp(bt,lev=0.9,nodiff=TRUE)

}
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ htest }% __ONLY ONE__ keyword per line
