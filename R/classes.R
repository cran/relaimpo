
#classes and methods in the following

setClass("relimplm",representation=representation(
    var.y="numeric",R2="numeric",
    lmg="numeric",pmvd="numeric",first="numeric",last="numeric",
    betasq="numeric",pratt="numeric",
    lmg.rank="numeric",pmvd.rank="numeric",
    first.rank="numeric",last.rank="numeric",
    betasq.rank="numeric",pratt.rank="numeric",
    lmg.diff="numeric",pmvd.diff="numeric",first.diff="numeric",
    last.diff="numeric",
    betasq.diff="numeric",pratt.diff="numeric",
    namen="character",type="character"))
setValidity("relimplm",function(object){
    p<-length(slot(object,"namen"))-1
    var.y <- is(slot(object,"var.y"),"numeric") && slot(object,"var.y")>0
    R2<-is(object@R2,"numeric") && object@R2>=0 && object@R2<=1
    lmg<-is(object@lmg,"numeric") && (length(object@lmg) %in% c(0,p))
    lmg.rank<-is(object@lmg.rank,"numeric") && (length(object@lmg.rank) %in% c(0,p))
    lmg.diff<-is(object@lmg.diff,"numeric") && (length(object@lmg.diff) %in% c(0,p*(p-1)/2))
    pmvd<-is(object@pmvd,"numeric") && (length(object@pmvd) %in% c(0,p))
    pmvd.rank<-is(object@pmvd.rank,"numeric") && (length(object@pmvd.rank) %in% c(0,p))
    pmvd.diff<-is(object@pmvd.diff,"numeric") && (length(object@pmvd.diff) %in% c(0,p*(p-1)/2))
    last<-is(object@last,"numeric") && (length(object@last) %in% c(0,p))
    last.rank<-is(object@last.rank,"numeric") && (length(object@last.rank) %in% c(0,p))
    last.diff<-is(object@last.diff,"numeric") && (length(object@last.diff) %in% c(0,p*(p-1)/2))
    first<-is(object@first,"numeric") && (length(object@first) %in% c(0,p))
    first.rank<-is(object@first.rank,"numeric") && (length(object@first.rank) %in% c(0,p))
    first.diff<-is(object@first.diff,"numeric") && (length(object@first.diff) %in% c(0,p*(p-1)/2))
    betasq<-is(object@betasq,"numeric") && (length(object@betasq) %in% c(0,p))
    betasq.rank<-is(object@betasq.rank,"numeric") && (length(object@betasq.rank) %in% c(0,p))
    betasq.diff<-is(object@betasq.diff,"numeric") && (length(object@betasq.diff) %in% c(0,p*(p-1)/2))
    pratt<-is(object@pratt,"numeric") && (length(object@pratt) %in% c(0,p))
    pratt.rank<-is(object@pratt.rank,"numeric") && (length(object@pratt.rank) %in% c(0,p))
    pratt.diff<-is(object@pratt.diff,"numeric") && (length(object@pratt.diff) %in% c(0,p*(p-1)/2))
    namen <- is(slot(object,"namen"),"character")
    return(var.y && R2 && lmg && lmg.rank && lmg.diff 
        && pmvd && pmvd.rank && pmvd.diff 
        && last && last.rank && last.diff && first && first.rank && first.diff
        && betasq && betasq.rank && betasq.diff 
        && pratt && pratt.rank && pratt.diff && namen)
    })
setAs("relimplm","list",function(from,to){
    to<-slot(from,"var.y")
    to<-append(to,list(R2=slot(from,"R2")))
    if (length(from@lmg)>0) to<-append(to,list(lmg=as.vector(from@lmg)))
    if (length(from@lmg.rank)>0) to<-append(to,list(lmg.rank=as.vector(from@lmg.rank)))
    if (length(from@lmg.diff)>0) to<-append(to,list(lmg.diff=as.vector(from@lmg.diff)))
    if (length(from@pmvd)>0) to<-append(to,list(pmvd=as.vector(from@pmvd)))
    if (length(from@pmvd.rank)>0) to<-append(to,list(pmvd.rank=as.vector(from@pmvd.rank)))
    if (length(from@pmvd.diff)>0) to<-append(to,list(pmvd.diff=as.vector(from@pmvd.diff)))
    if (length(from@last)>0) to<-append(to,list(last=as.vector(from@last)))
    if (length(from@last.rank)>0) to<-append(to,list(last.rank=as.vector(from@last.rank)))
    if (length(from@last.diff)>0) to<-append(to,list(last.diff=as.vector(from@last.diff)))
    if (length(from@first)>0) to<-append(to,list(first=as.vector(from@first)))
    if (length(from@first.rank)>0) to<-append(to,list(first.rank=as.vector(from@first.rank)))
    if (length(from@first.diff)>0) to<-append(to,list(first.diff=as.vector(from@first.diff)))
    if (length(from@betasq)>0) to<-append(to,list(betasq=as.vector(from@betasq)))
    if (length(from@betasq.rank)>0) to<-append(to,list(betasq.rank=as.vector(from@betasq.rank)))
    if (length(from@betasq.diff)>0) to<-append(to,list(betasq.diff=as.vector(from@betasq.diff)))
    if (length(from@pratt)>0) to<-append(to,list(pratt=as.vector(from@pratt)))
    if (length(from@pratt.rank)>0) to<-append(to,list(pratt.rank=as.vector(from@pratt.rank)))
    if (length(from@pratt.diff)>0) to<-append(to,list(pratt.diff=as.vector(from@pratt.diff)))
    if (length(from@namen)>0) to<-append(to,list(namen=as.vector(from@namen)))
    return(to)})


setOldClass("boot") 

setClass("relimplmboot",representation=representation(
    boot="boot",type="character",nboot="numeric",
    rank="logical",diff="logical",rela="logical"))

setClass("relimplmbooteval",representation=representation(
    lmg.lower="matrix",lmg.upper="matrix",
    lmg.rank.lower="matrix",lmg.rank.upper="matrix",
    lmg.diff.lower="matrix",lmg.diff.upper="matrix",
    pmvd.lower="matrix",pmvd.upper="matrix",
    pmvd.rank.lower="matrix",pmvd.rank.upper="matrix",
    pmvd.diff.lower="matrix",pmvd.diff.upper="matrix",
    last.lower="matrix",last.upper="matrix",
    last.rank.lower="matrix",last.rank.upper="matrix",
    last.diff.lower="matrix",last.diff.upper="matrix",
    first.lower="matrix",first.upper="matrix",
    first.rank.lower="matrix",first.rank.upper="matrix",
    first.diff.lower="matrix",first.diff.upper="matrix",
    betasq.lower="matrix",betasq.upper="matrix",
    betasq.rank.lower="matrix",betasq.rank.upper="matrix",
    betasq.diff.lower="matrix",betasq.diff.upper="matrix",
    pratt.lower="matrix",pratt.upper="matrix",
    pratt.rank.lower="matrix",pratt.rank.upper="matrix",
    pratt.diff.lower="matrix",pratt.diff.upper="matrix",
    var.y.boot="numeric",R2.boot="numeric",
    lmg.boot="matrix",pmvd.boot="matrix",last.boot="matrix",
    first.boot="matrix",betasq.boot="matrix",pratt.boot="matrix",
    lmg.rank.boot="matrix",pmvd.rank.boot="matrix",last.rank.boot="matrix",
    first.rank.boot="matrix",betasq.rank.boot="matrix",pratt.rank.boot="matrix",
    lmg.diff.boot="matrix",pmvd.diff.boot="matrix",last.diff.boot="matrix",
    first.diff.boot="matrix",betasq.diff.boot="matrix",pratt.diff.boot="matrix",
    level="numeric",nboot="numeric",diffnam="character",rank="logical",diff="logical",
    rela="logical",type="character",sort="logical",bty="character",mark="matrix",
    markdiff="matrix"),contains="relimplm")


