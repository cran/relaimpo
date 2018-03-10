\name{mianalyze.relimp}
\alias{mianalyze.relimp}

\title{Function to do relative importance calculations based on multiply imputed datasets}
\description{
The function mianalyze.relimp takes a list of imputed data frames (or matrices), 
calculates relative importance metrics for each of these and combines these metrics into 
overall estimates with estimated variances according to the method by Rubin (1987). 
The output object can be summarized, printed and plotted.
}
\usage{
mianalyze.relimp(implist, level = 0.95, sort = FALSE, ..., b = 50, type = "lmg", 
    diff = TRUE, no.CI = FALSE, rela = FALSE, always = NULL, groups = NULL, 
    groupnames = NULL, deslist = NULL, bootlist.out = FALSE, formula = NULL, 
    weights = NULL, strata=NULL, ids=NULL)
}
\arguments{
  \item{implist}{list of data frames or matrices containing multiply-imputed datasets, 
       or object of class \code{imputationList}
       
       If no formula is given, the first column of each data frame/matrix is assumed to 
       be the response variable, the other columns are regressors.
       
       If a list of designs is also given, the \code{variables} component of each design must 
       consist of the necessary columns from the respective entry in implist; 
       if no formula is given, the \code{variables} component of each design must coincide 
       (except for the order of columns) with the respective entry in implist.}
  \item{level}{is a single confidence level (between 0.5 and 1)}
  \item{sort}{ is a logical requesting output sorted by size of relative contribution 
       (\code{sort=TRUE}) or by variable position in list (\code{sort=FALSE}, default). }
  \item{\dots}{Further arguments, currently none available}
  \item{b}{ is the number of bootstrap runs requested on boot.relimp (default: \code{b=50}). 
       Make sure to set this to a higher number, if you want to subsequently use the \code{bootlist} slot
       for calculating further confidence intervals with function \code{booteval.relimp}. }
  \item{type}{ cf. \code{\link{calc.relimp}}.}
  \item{diff}{ is a logical requesting bootstrapping of pairwise differences in relative importance (\code{diff=TRUE}, 
         default) for each metric in type }
  \item{no.CI}{ if set to TRUE, suppresses calculation of confidence intervals and only averages estimated metrics
       from all imputed data sets in implist. Currently, \code{no.CI = TRUE} is the only setting for which 
       \code{mianalyze.relimp} works when using models with factors, groups or interactions.}
  \item{rela}{ cf. \code{\link{calc.relimp}}.}
  \item{always}{ cf. \code{\link{calc.relimp}}.}
  \item{groups}{ cf. \code{\link{calc.relimp}}.}
  \item{groupnames}{ cf. \code{\link{calc.relimp}}.}
  \item{deslist}{ is a list of design object of class \code{survey.design} (cf. package \code{survey}).
       You can EITHER specify a \code{deslist} OR \code{weights} and/or \code{strata} and/or \code{ids}. 
       Note that the design list must contain the same data objects (in the \dQuote{variables} element) 
       that are listed in \code{implist}, so that a lot of storage space is needed in case of large datasets.
       
       If deslist is not given, the function creates a list of designs using \code{weights}, 
       \code{strata}, and \code{ids} information. Whenever the required designs are simple enough to be covered by
       specifying \code{weights}, \code{strata}, and \code{ids}, this is by far preferrable in terms of 
       storage.  }
  \item{bootlist.out}{ If TRUE, the individual bootstrap results for each multiply imputed data set 
       are stored in the bootlist slot of the output object (may be storage-intensive).
       }
  \item{formula}{ cf. \code{\link{boot.relimp}}; 
       NOTE: If no.CI = FALSE, i.e. confidence intervals are not suppressed, 
       \code{formula} has to follow the same restrictions as mentioned under item \code{design} for 
       \code{\link{boot.relimp}} (no calculated variables, no interaction terms, no factors),
       since confidence interval calculations in \code{mianalyze.relimp} are design-based,
       even if no \code{deslist-}option is given.}
  \item{weights}{ is a vector of case weights for the observations in the data frame (or matrix).
       You can EITHER specify \code{weights} OR a \code{deslist}. 
       If \code{weights} is NULL, equal weights are assumed, unless otherwise specified in \code{deslist}.
       For the different types of weights and their appropriate treatment for obtaining confidence intervals,
       cf. the \dQuote{Details} section of \code{\link{boot.relimp}}.}
  \item{strata}{ is a strata request that will be handed to function \code{svydesign} for defining the 
       strata in a survey design (to be given to mianalyze without the \dQuote{\~}).
       You can EITHER specify \code{strata} OR a \code{deslist}. 
       If \code{strata} is NULL, one stratum is assumed, unless otherwise specified in \code{deslist}.}
  \item{ids}{ is an id-request that will be handed to function \code{svydesign} for defining the clusters 
       in a survey design (to be given to mianalyze without the \dQuote{\~}).
       You can EITHER specify \code{ids} OR a \code{deslist}. 
       If \code{ids} is NULL, it is assumed that there are no clusters, 
       unless otherwise specified in \code{deslist}.}
}

\details{
  Multiple imputation is a contemporary method for handling data with a substantial missing value problem. 
  It produces a number of completed data sets (e.g. 10) the inference from which is subsequently combined. 
  The most frequently used way of combination is the one by Rubin: estimates from the 
  different completed data sets are averaged, and the variance is estimated by combining the average over the 
  estimated variances (within imputation variance) with an appropriately-scaled variance between estimates,
  and confidence intervals are obtained by using a t-distribution with appropriately 
  chosen degrees of freedom.
  
  The variance-covariance matrix of the vector of estimates for each individual completed data set is obtained 
  from function withReplicates in package \pkg{survey} based on \pkg{survey}'s bootstrap replication weights. 
  On request (\code{bootlist.out=TRUE}), the underlying bootstrap resamples are also stored in the 
  \code{bootlist}-slot of the output object. In this case, list elements of the \code{bootlist}-slot 
  are objects of class \code{relimplmboot} and can be processed by function \code{booteval.relimp}. 
  This can help in getting an impression whether the overall aggregated confidence intervals are heavily distorted 
  towards symmetry. If such sanity-checking is intended, the default value for \code{b} should be substantially 
  increased.
  
  Function \code{mianalyze.relimp} needs a list of multiply-imputed data sets or an object of class 
  \code{imputationList} for input. 
  Multiply imputed data sets can - within R - be obtained from various packages. Hints for creating lists 
  of the form needed for \code{mianalyze.relimp} are given below for users of functions \code{aregImpute}, 
  \code{mice}, and \code{amelia}. Users of packages \pkg{norm}, \pkg{cat}, \pkg{mix}, or \pkg{pan} 
  (who have managed to operate these extremely uncomfortable packages) can of course also produce lists of 
  imputed data sets (only less comfortably). 
  
  For an object \code{imp} of class \code{mids} obtained from function \code{mice} in package \pkg{mice}, 
  the code 
  
  \code{lapply(as.list(1:imp$m),function(obj) complete(imp,action=obj))} 
  
  produces a list of multiply-imputed data sets as needed for function \code{mianalyze.relimp}.
  For an object \code{f} of class \code{aregImpute} produced by function \code{aregImpute} in package \pkg{Hmisc},
  
  \code{lapply(as.list(1:f$m),function(obj) impute(imp,imputation=obj))} 
  
  produces the required list of multiply-imputed data sets.
  For an object \code{output} produced by function \code{amelia} in package \pkg{Amelia}, the code 
  
  \code{output[1:output$amelia.args$m]} 
  
  produces the list of multiply-imputed data sets as needed for function \code{mianalyze.relimp}.
  
  For multiple imputation, practice is in many cases ahead of theory; this is no different with function 
  \code{mianalyze.relimp}. Users should note that the validity of confidence intervals has only been proven for 
  likelihood-based analyses. Since the metrics calculated in \pkg{relaimpo} are not strictly likelihood-based, 
  the confidence intervals from function \code{mianalyze.relimp} must be considered approximate and experimental.
}
\value{
  The value returned by function \code{mianalyze.relimp} is an object of class \code{relimplmbootMI} 
  (if \code{no.CI = FALSE}, default) or an object of class \code{relimplm} (if no.CI=TRUE). 
  It can be printed, plotted and summarized using special methods. 
  For extracting its items, the \code{@} or \code{$} extractors can be used. 

In addition to the items described for function \code{\link{calc.relimp}}, which are also available here, 
the following items from class \code{relimplmbootMI} may be of interest for further calculations:
  \item{metric.lower}{matrix of lower confidence bounds for \dQuote{metric}: one row for each confidence level, 
        one column for each element of \dQuote{metric}. \dQuote{metric} can be any of \code{lmg}, \code{lmg.rank}, 
        \code{lmg.diff}, ... 
       (replace \code{lmg} with other available relative importance metrics, cf. \code{\link{calc.relimp}})}
  \item{metric.upper}{matrix of upper confidence bounds for \dQuote{metric}: one row for each confidence level, 
        one column for each element of \dQuote{metric}}
  \item{nboot}{number of bootstrap runs underlying the evaluations}
  \item{level}{confidence level}
  \item{MIresult}{object of class \code{MIresult} that can be processed by the function summary.MIresult 
       from package survey}
  \item{bootlist}{only available if bootlist.out=TRUE has been chosen; 
       list of objects of class boot.relimp; each list element can be input to function 
       \code{\link{booteval.relimp}} }
}

\references{
Chevan, A. and Sutherland, M. (1991) Hierarchical Partitioning. \emph{The American Statistician} \bold{45}, 90--96.

Darlington, R.B. (1968) Multiple regression in psychological research and practice. \emph{Psychological Bulletin}  \bold{69}, 161--182.

Feldman, B. (2005) Relative Importance and Value. Manuscript (Version 1.1, March 19 2005), downloadable at \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2255827}

Genizi, A. (1993) Decomposition of R2 in multiple regression with correlated regressors. \emph{Statistica Sinica} \bold{3}, 407--420. 
Downloadable at \url{http://www3.stat.sinica.edu.tw/statistica/password.asp?vol=3&num=2&art=10}

Groemping, U. (2006) Relative Importance for Linear Regression in R: The Package relaimpo 
         \emph{Journal of Statistical Software} \bold{17}, Issue 1. 
         Downloadable at \url{http://www.jstatsoft.org/v17/i01}

Lindeman, R.H., Merenda, P.F. and Gold, R.Z. (1980) \emph{Introduction to Bivariate and Multivariate Analysis}, Glenview IL: Scott, Foresman.

Little, R.J.A. and Rubin, D.B. (2002) \emph{Statistical Analysis with Missing Data}, Wiley, New York.

Zuber, V. and Strimmer, K. (2010) \emph{Variable importance and model selection by decorrelation}. Preprint, downloadable at \url{http://arxiv.org/abs/1007.5516}

Go to \url{http://prof.beuth-hochschule.de/groemping/} for further information and references.
}
\author{ Ulrike Groemping, BHT Berlin }

\note{ There are two versions of this package. The version on CRAN is globally licensed under GPL version 2 (or later). 
There is an extended version with the interesting additional metric \code{pmvd} that is licensed according to GPL version 2
under the geographical restriction "outside of the US" because of potential issues with US patent 6,640,204. This version can be obtained 
from Ulrike Groempings website (cf. references section). Whenever you load the package, a display tells you, which version you are loading. }

\section{Warning }{The confidence intervals produced here should be used for exploratory purposes only. 
They can be somewhat liberal and are likely to be too symmetric particularly for small data sets.
The confidence intervals produced by function \code{mianalyze.relimp} need further research into their behaviour 
and are currently considered experimental.

Be aware that the methods themselves (\code{lmg} and even more \code{pmvd}) need some computing time in case of 
many regressors. Hence, bootstrapping of multiple data sets should be used with awareness of computing time issues.

} 

\seealso{ \pkg{\link{relaimpo}}, \code{\link{calc.relimp}}, \code{\link{booteval.relimp}}, 
  \code{\link{classesmethods.relaimpo}}}
\examples{
  ## smi contains a list of 5 imputed datasets (class imputationList) from package mitools
  ## (first element of smi is list of data frames)
  ## it is not a well-suited example for relative importance but easily available for demonstrating 
  ##         multiple imputation-related functionality
  
  data(smi)
  
  ## obtain averaged estimates only, without confidence intervals
  ## works with factors and interactions
  mianalyze.relimp(smi[[1]], formula = cistot ~ drkfre+sex+wave, no.CI = TRUE)
  ## for obtaining all individual estimates, use lapply:
  smi.cr.list <- lapply(smi[[1]], function(obj) calc.relimp(cistot ~ drkfre+sex+wave, data=obj))
  ## display result for first individual imputed data set
  smi.cr.list[[1]]
  
  ## obtain confidence intervals, 
  ## currently only usable for models without calculated variables, factors, groups, interactions
  
  ## call without using weights, strata, clusters or a design list
  mianalyze.relimp(smi[[1]], formula = cistot ~ mdrkfre+sex+wave)  
  ## call using the id column (identical in all smi data sets) for cluster structure
  ident <- smi[[1]][[1]]$id
  mitest <- mianalyze.relimp(smi[[1]], formula = cistot ~ mdrkfre+sex+wave, ids=ident)  
  mitest
      ## postprocess: look at intervals with different confidence level
      summary(mitest,lev=0.8)
  ## call with design list
  deslist = lapply(smi[[1]], function(obj) svydesign(~id,strata=~sex,weights=~cistot,data=obj))
  mitest <- mianalyze.relimp(smi[[1]], formula = cistot ~ mdrkfre+sex+wave, deslist=deslist, 
       level=c(0.8))
  mitest  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
\keyword{ htest }% __ONLY ONE__ keyword per line
