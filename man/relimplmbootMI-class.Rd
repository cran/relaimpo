\name{relimplmbootMI-class}
\docType{class}
\alias{relimplmbootMI-class}
\alias{show,relimplmbootMI-method}
%document S4 methods only, do not make index too long
%\alias{print,relimplmbootMI-method}
%\alias{plot,relimplmbootMI-method}
%\alias{summary,relimplmbootMI-method}
%\alias{$,relimplmbootMI-method}

\title{Class relimplmbootMI }
\description{  Output object from function \code{\link{mianalyze.relimp}}, described there. }
\section{Objects from the Class}{
Objects should only be created by calls to function \code{mianalyze.relimp}.
}
\section{Slots}{
  \describe{
    \item{\code{lmg.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{lmg.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{lmg.rank.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{lmg.rank.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{lmg.diff.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{lmg.diff.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{pmvd.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{pmvd.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{pmvd.rank.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{pmvd.rank.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{pmvd.diff.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{pmvd.diff.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{last.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{last.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{last.rank.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{last.rank.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{last.diff.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{last.diff.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{first.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{first.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{first.rank.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{first.rank.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{first.diff.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{first.diff.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{betasq.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{betasq.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{betasq.rank.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{betasq.rank.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{betasq.diff.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{betasq.diff.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{pratt.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{pratt.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{pratt.rank.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{pratt.rank.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{pratt.diff.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{pratt.diff.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{genizi.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{genizi.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{genizi.rank.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{genizi.rank.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{genizi.diff.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{genizi.diff.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{car.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{car.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{car.rank.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{car.rank.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{car.diff.lower}:}{Object of class \code{"matrix"}  }
    \item{\code{car.diff.upper}:}{Object of class \code{"matrix"}  }
    \item{\code{MIresult}:}{Object of class \code{"MIresult"}  }
    \item{\code{est}:}{Object of class \code{"numintnull"}  }
    \item{\code{vcov}:}{Object of class \code{"nummatnull"}  }
    \item{\code{bootlist}:}{Object of class \code{"listnull"}, 
          list of objects of class relimplmboot or null }
    \item{\code{level}:}{Object of class \code{"numeric"}  }
    \item{\code{nboot}:}{Object of class \code{"numeric"}  }
    \item{\code{diffnam}:}{Object of class \code{"character"}  }
    \item{\code{rank}:}{Object of class \code{"logical"}  }
    \item{\code{diff}:}{Object of class \code{"logical"}  }
    \item{\code{rela}:}{Object of class \code{"logical"}  }
    \item{\code{fixed}:}{Object of class \code{"logical"}  }
    \item{\code{type}:}{Object of class \code{"character"}  }
    \item{\code{sort}:}{Object of class \code{"logical"}  }
    \item{\code{bty}:}{Object of class \code{"character"}  }
    \item{\code{mark}:}{Object of class \code{"matrix"}  }
    \item{\code{markdiff}:}{Object of class \code{"matrix"}  }
    \item{\code{var.y}:}{Object of class \code{"numeric"}  }
    \item{\code{R2}:}{Object of class \code{"numeric"}  }
    \item{\code{R2.decomp}:}{Object of class \code{"numeric"}  }
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
    \item{\code{namen}:}{Object of class \code{"character"}  }
    \item{\code{nobs}:}{Object of class \code{"numeric"}: number of valid observations used in analysis  }
    \item{\code{always}:}{Object of class \code{"numintnull"}: variables always in the model }
    \item{\code{alwaysnam}:}{Object of class \code{"charnull"}: variable names of variables always in model }
  }
}
\section{Extends}{
  Class \code{"relimplm"}, directly.
    }
\section{Methods}{
  This documentation section documents S4 methods only. 
  There are also S3-methods for printing, plotting, extracting slots with \$ and producing 
  a brief summary. For their usage, see \code{\link{classesmethods.relaimpo}}.
    \describe{
%    \item{print}{ S3-method: meaningful printout of object }
    \item{show}{ S4-method: identical to \code{print}}
%    \item{summary}{ S3-method: applies the function \code{summary} from package \pkg{mitools} 
         to the slot \code{MIresult} for a briefer output than the standard printout }
%    \item{plot}{ S3-method: plots estimated proportions with indication of confidence intervals }
%    \item{\$}{S3 method: extracts slots with \$ extractor }
    }
\author{ Ulrike Groemping, BHT Berlin}
\seealso{
    \code{\link{relimplm-class}}, \code{\link{relimplmbooteval-class}}, 
    \code{\link{booteval.relimp}}, \code{\link{print.relimplmbootMI}}, 
    \code{\link{plot.relimplmbootMI}}, \code{\link{summary.relimplmbootMI}}
}
\keyword{classes}
