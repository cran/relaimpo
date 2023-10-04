calcrelimp.forsurvey <- function(wt,dat,...){
  ergeb <- list2vec(as(calc.relimp_default.intern(cov.wt(dat, wt=wt)$cov, ...),"list"))
  }
