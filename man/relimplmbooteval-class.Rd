\name{relimplmbooteval-class}
\docType{class}
\alias{relimplmbooteval-class}
\alias{show,relimplmbooteval-method}
%\alias{print,relimplmbooteval-method}
%\alias{plot,relimplmbooteval-method}
%\alias{$,relimplmbooteval-method}

\title{Class relimplmbooteval }
\description{  Output object from function \code{\link{booteval.relimp}}, described there. }
\section{Objects from the Class}{
Objects should only be created by calls to function \code{booteval.relimp}.
}
\section{Slots}{
  \describe{
    \item{\code{lmg.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{lmg.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{lmg.rank.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{lmg.rank.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{lmg.diff.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{lmg.diff.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pmvd.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pmvd.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pmvd.rank.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pmvd.rank.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pmvd.diff.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pmvd.diff.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{last.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{last.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{last.rank.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{last.rank.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{last.diff.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{last.diff.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{first.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{first.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{first.rank.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{first.rank.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{first.diff.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{first.diff.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{betasq.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{betasq.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{betasq.rank.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{betasq.rank.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{betasq.diff.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{betasq.diff.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pratt.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pratt.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pratt.rank.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pratt.rank.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pratt.diff.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pratt.diff.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genizi.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genizi.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genizi.rank.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genizi.rank.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genizi.diff.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genizi.diff.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{car.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{car.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{car.rank.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{car.rank.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{car.diff.lower}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{car.diff.upper}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{var.y.boot}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{R2.boot}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{R2.decomp.boot}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lmg.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pmvd.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{last.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{first.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{betasq.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pratt.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genizi.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{car.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{lmg.rank.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pmvd.rank.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{last.rank.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{first.rank.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{betasq.rank.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pratt.rank.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genizi.rank.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{car.rank.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{lmg.diff.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pmvd.diff.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{last.diff.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{first.diff.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{betasq.diff.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{pratt.diff.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{genizi.diff.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{car.diff.boot}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{est}:}{Object of class \code{"numintnull"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"numintmatnull"} ~~ }
    \item{\code{level}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nboot}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{diffnam}:}{Object of class \code{"character"} ~~ }
    \item{\code{rank}:}{Object of class \code{"logical"} ~~ }
    \item{\code{diff}:}{Object of class \code{"logical"} ~~ }
    \item{\code{rela}:}{Object of class \code{"logical"} ~~ }
    \item{\code{fixed}:}{Object of class \code{"logical"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{sort}:}{Object of class \code{"logical"} ~~ }
    \item{\code{bty}:}{Object of class \code{"character"} ~~ }
    \item{\code{mark}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{markdiff}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{var.y}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{R2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{R2.decomp}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lmg}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pmvd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{first}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{last}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{betasq}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pratt}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{genizi}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{car}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lmg.rank}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pmvd.rank}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{first.rank}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{last.rank}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{betasq.rank}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pratt.rank}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{genizi.rank}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{car.rank}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{lmg.diff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pmvd.diff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{first.diff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{last.diff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{betasq.diff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pratt.diff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{genizi.diff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{car.diff}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{namen}:}{Object of class \code{"character"} ~~ }
    \item{\code{nobs}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ave.coeffs}:}{Object of class \code{"numintmatnull"} ~~ }
    \item{\code{always}:}{Object of class \code{"numintnull"} ~~ }
    \item{\code{alwaysnam}:}{Object of class \code{"charnull"} ~~ }
    \item{\code{groupdocu}:}{Object of class \code{"list"} ~~ }
    \item{\code{call}:}{Object of class \code{"langnull"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{relimplm}"}, directly.
}
\section{Methods}{
  This documentation section documents S4 methods only. 
  There are also S3-methods for printing, plotting and extracting slots with \code{$}. 
  For their usage, see \code{\link{classesmethods.relaimpo}}.
    \describe{
%    \item{print}{S3 method: meaningful printout of object }
    \item{show}{S4 method: identical to \code{print}}
%    \item{plot}{S3 method: plots estimated proportions with indication of confidence intervals }
%    \item{\code{$}}{S3 method: extracts slots with \code{$} extractor }
    }
}
\author{ Ulrike Groemping, BHT Berlin}
\seealso{
  \code{\link{relimplm-class}}, \code{\link{relimplmboot-class}}, \code{\link{relimplmbootMI-class}}, 
    \code{\link{booteval.relimp}}, \code{\link{print.relimplmbooteval}}, 
    \code{\link{plot.relimplmbooteval}}
}
\keyword{classes}
