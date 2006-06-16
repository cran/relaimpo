\name{relimplm-class}
\docType{class}
\alias{relimplm-class}
\alias{coerce,relimplm,list-method}

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
    \item{\code{lmg}:}{Object of class \code{"numeric"} }
    \item{\code{pmvd}:}{Object of class \code{"numeric"} }
    \item{\code{last}:}{Object of class \code{"numeric"} }
    \item{\code{first}:}{Object of class \code{"numeric"} }
    \item{\code{betasq}:}{Object of class \code{"numeric"} }
    \item{\code{pratt}:}{Object of class \code{"numeric"}  }
    \item{\code{lmg.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{pmvd.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{last.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{first.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{betasq.rank}:}{Object of class \code{"numeric"}  }
    \item{\code{pratt.rank}:}{Object of class \code{"numeric"} }
    \item{\code{lmg.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{pmvd.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{last.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{first.diff}:}{Object of class \code{"numeric"} }
    \item{\code{betasq.diff}:}{Object of class \code{"numeric"}  }
    \item{\code{pratt.diff}:}{Object of class \code{"numeric"} }
    \item{\code{namen}:}{Object of class \code{"character"}: variable names, starting with name for response }
    \item{\code{nobs}:}{Object of class \code{"numeric"}: number of valid observations used in analysis 
                                            (available only if calculations are based on raw data) }
    \item{\code{type}:}{Object of class \code{"character"}: metrics that have been calculated }
    \item{\code{rela}:}{Object of class \code{"logical"}: TRUE means that metrics have been normalized to sum 100pct }
    \item{\code{always}:}{Object of class \code{"numintnull"}: variables always in the model (adjusted for) }
    \item{\code{alwaysnam}:}{Object of class \code{"charnull"}: variable names of variables always in model }
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{If an object of this class is coerced to list by \code{as(object, "list")}, 
    the slots are output to the list in the following order: \code{var.y}, \code{R2}, \code{lmg}, \code{lmg.rank}, \code{lmg.diff}, 
next three in non-US version only: \code{pmvd}, \code{pmvd.rank}, \code{pmvd.diff}, \code{last},  
\code{last.rank}, \code{last.diff}, ..., \code{pratt}, \code{pratt.rank}, \code{pratt.diff}, \code{namen}, \code{type}, 
\code{rela}, \code{always}, \code{alwaysnam}
}
}
}
\author{ Ulrike Groemping, TFH Berlin }


\seealso{
  See also \code{\link{calc.relimp}}, \pkg{\link{relaimpo}}
}
\keyword{classes}
\keyword{methods}
