"plot.relimplm" <-
function(x,...,names.abbrev=4){

if(!(is(x,"relimplm"))) stop("x must be the output from function calc.relimp")
if(!(is.numeric(names.abbrev))) stop("names.abbrev must be a number")

    p<-length(x@namen)-1
    yname<-x@namen[1]
    xnames<-substr(x@namen[2:(p+1)],1,names.abbrev)
    ylab<-paste("Relative Importances for ", yname,sep="")
    max<-max(x@lmg,x@pmvd,x@last,x@first,x@betasq,x@pratt)
    min<-min(0,x@lmg,x@pmvd,x@last,x@first,x@betasq,x@pratt)
    axmax<-ceiling(10*max)/10
    axmin<-floor(10*min)/10
    type<-x@type
    if (length(type)==0) 
        print("Nothing to plot")
    else 
    {
    ntype<-length(type)
    op <- par(no.readonly = TRUE)
    if (ntype==2) par(mfrow=c(1,2))
    if (ntype>2 && ntype<=4) par(mfrow=c(2,2))
    if (ntype>2 && ntype<=4) par(mfrow=c(2,2))
    if (ntype>4) par(mfrow=c(2,3))
    if ("lmg" %in% x@type) barplot(x@lmg,ylab=ylab,main="Method LMG",names.arg=xnames,ylim=c(axmin,axmax))
    if ("pmvd" %in% x@type) barplot(x@pmvd,ylab=ylab,main="Method PMVD",names.arg=xnames,ylim=c(axmin,axmax))
    if ("last" %in% x@type) barplot(x@last,ylab=ylab,main="Method Last",names.arg=xnames,ylim=c(axmin,axmax))
    if ("first" %in% x@type) barplot(x@first,ylab=ylab,main="Method First",names.arg=xnames,ylim=c(axmin,axmax))
    if ("betasq" %in% x@type) barplot(x@betasq,ylab=ylab,main="Method Betasq",names.arg=xnames,ylim=c(axmin,axmax))
    if ("pratt" %in% x@type) barplot(x@pratt,ylab=ylab,main="Method Pratt",names.arg=xnames,ylim=c(axmin,axmax))
    par(op)
    }
}

