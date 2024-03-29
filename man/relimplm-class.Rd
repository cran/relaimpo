\name{relimplm-class}
\docType{class}
\alias{relimplm-class}
\alias{coerce,relimplm,list-method}
\alias{show,relimplm-method}
%document S4 methods only, do not make index too long
%\alias{print,relimplm-method}
%\alias{plot,relimplm-method}
%\alias{$,relimplm-method}

\title{Class relimplm }
\description{  This is the class of output objects from the function \code{calc.relimp} 
in package \pkg{relaimpo}. Its elements are described in the documentation of \code{\link{calc.relimp}}.
}
\section{Objects from the Class}{
Objects should only be created by calls to the function \code{calc.relimp}.  
}
\section{Slots}{
  \describe{
    \item{\code{var.y}:}{Object of class \code{"numeric"} }
    \item{\code{R2}:}{Object of class \code{"numeric"} }
    \item{\code{R2.decomp}:}{Object of class \code{"numeric"} }
    \item{\code{lmg}:}{Object of class \code{"numeric"}  }
    \item{\code{pmvd}:}{Object of class \code{"numeric"}  }
    \item{\code{first}:}{Object of class \code{"numeric"}  }
    \item{\code{last}:}{Object of class \code{"numeric"}  }
    \item{\code{betasq}:}{Object of class \code{"numeric"}  }
    \item{\code{pratt}:}{Object of class \code{"numeric"}  }
    \item{\code{genizi}:}{Object of class \code{"numeric"}  }
    \item{\code{car}:}{Object of class \code{"numeric"}  }
    \item{\code{lmg.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{pmvd.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{first.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{last.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{betasq.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{pratt.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{genizi.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{car.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{lmg.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{pmvd.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{first.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{last.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{betasq.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{pratt.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{genizi.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{car.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{namen}:}{Object of class \code{"character"}: variable names, starting with name for response }
    \item{\code{nobs}:}{Object of class \code{"numeric"}: number of valid observations used in analysis 
                                            (available only if calculations are based on raw data) }
    \item{\code{ave.coeffs}:}{Object of class \code{"matrix"}: average coefficients for models of different sizes }
    \item{\code{type}:}{Object of class \code{"character"}: metrics that have been calculated }
    \item{\code{rela}:}{Object of class \code{"logical"}: TRUE means that metrics have been normalized to sum 100pct }
    \item{\code{always}:}{Object of class \code{"numintnull"}: variables always in the model (adjusted for),
       in terms of columns in the model matrix (i.e. in case of factors in the model, there are more entries 
       than in the call to function \code{calc.relimp}; first possible column is 2 (1 reserved for response)) }
    \item{\code{alwaysnam}:}{Object of class \code{"charnull"}: variable names of variables always in model }
    \item{\code{groupdocu}:}{Object of class \code{"list"}: information on group definitions }
    \item{\code{call}:}{Object of class \code{"langnull"}: the call that created the object }
  }
}
\section{Methods}{
  This documentation section documents S4 methods only. 
  There are also S3-methods for printing, plotting and extracting slots with \code{$}. 
  For their usage, see \code{\link{classesmethods.relaimpo}}.
  \describe{
    \item{coerce}{S4-method: If an object of this class is coerced to list by \code{as(object, "list")}, 
    the slots are output to the list in the following order: \code{var.y}, \code{R2}, \code{lmg}, \code{lmg.rank}, \code{lmg.diff}, 
next three in non-US version only: \code{pmvd}, \code{pmvd.rank}, \code{pmvd.diff}, \code{last},  
\code{last.rank}, \code{last.diff}, ..., \code{pratt}, \code{pratt.rank}, \code{pratt.diff}, \code{namen}, \code{type}, 
\code{rela}, \code{always}, \code{alwaysnam}
  }
 %   \item{print}{S3-method: meaningful printout of object, see \code{\link{print.relimplm}} }
    \item{show}{S4-method: identical to \code{print}}
 %   \item{plot}{S3-method: barplots of percentages, see \code{\link{plot.relimplm}}}
 %   \item{\code{$}}{S3 method: extracts slots with \code{$} extractor }
}
}
\author{ Ulrike Groemping, BHT Berlin }


\seealso{
  \code{\link{calc.relimp}}, \code{\link{print.relimplm}}, \code{\link{plot.relimplm}}, \pkg{\link{relaimpo}}
}
\keyword{classes}
\keyword{methods}
