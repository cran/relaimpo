
#classes and methods in the following
setClassUnion("numintnull", c("numeric", "integer", "NULL"))
setClassUnion("charnull", c("character", "NULL"))

setClass("relimplm",representation=representation(
    var.y="numeric",R2="numeric",R2.decomp="numeric",
    lmg="numeric",pmvd="numeric",first="numeric",last="numeric",
    betasq="numeric",pratt="numeric",
    lmg.rank="numeric",pmvd.rank="numeric",
    first.rank="numeric",last.rank="numeric",
    betasq.rank="numeric",pratt.rank="numeric",
    lmg.diff="numeric",pmvd.diff="numeric",first.diff="numeric",
    last.diff="numeric",
    betasq.diff="numeric",pratt.diff="numeric",
    namen="character",nobs="numeric",type="character",rela="logical",
    always="numintnull",alwaysnam="charnull", groupdocu="list"))
setValidity("relimplm",function(object){
    p<-length(slot(object,"namen"))-1
    var.y <- is(slot(object,"var.y"),"numeric") && slot(object,"var.y")>0
    R2<-is(object@R2,"numeric") && object@R2>=0 && object@R2<=1
    R2.decomp<-is(object@R2.decomp,"numeric") && object@R2.decomp>0 && 
             object@R2.decomp<=object@R2
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

setOldClass("boot") 

setClass("relimplmboot",representation=representation(
    boot="boot",type="character",nboot="numeric",
    rank="logical",diff="logical",rela="logical",fixed="logical", 
    namen="character", nobs="numeric", always="numintnull",alwaysnam="charnull",
    groupdocu="list"))

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
    var.y.boot="numeric",R2.boot="numeric",R2.decomp.boot="numeric",
    lmg.boot="matrix",pmvd.boot="matrix",last.boot="matrix",
    first.boot="matrix",betasq.boot="matrix",pratt.boot="matrix",
    lmg.rank.boot="matrix",pmvd.rank.boot="matrix",last.rank.boot="matrix",
    first.rank.boot="matrix",betasq.rank.boot="matrix",pratt.rank.boot="matrix",
    lmg.diff.boot="matrix",pmvd.diff.boot="matrix",last.diff.boot="matrix",
    first.diff.boot="matrix",betasq.diff.boot="matrix",pratt.diff.boot="matrix",
    level="numeric",nboot="numeric",diffnam="character",rank="logical",diff="logical",
    rela="logical",fixed="logical",type="character",sort="logical",bty="character",mark="matrix",
    markdiff="matrix"),contains="relimplm")


