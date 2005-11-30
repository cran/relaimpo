"print.relimplm" <-
function(x,...){
if(!(is(x,"relimplm"))) stop("x must be the output from function calc.relimp")

    p<-length(slot(x,"namen"))-1
    cat("Response variable:",slot(x,"namen")[1],"\n")
    cat("Total response variance:",x@var.y,"\n")
    cat(p,"Regressors:",paste(slot(x,"namen")[2:(p+1)],collapse=" "),"\n")
    cat("Proportion of variance explained by model: ",round(100*x@R2,2),"%","\n",sep="")
    cat("\n")
    cat("Relative importance metrics:","\n")
    type<-c("lmg","pmvd","last","first","betasq","pratt")
    if (length(x@lmg)==0) type<-setdiff(type,"lmg")
    if (length(x@pmvd)==0) type<-setdiff(type,"pmvd")
    if (length(x@last)==0) type<-setdiff(type,"last")
    if (length(x@first)==0) type<-setdiff(type,"first")
    if (length(x@betasq)==0) type<-setdiff(type,"betasq")
    if (length(x@pratt)==0) type<-setdiff(type,"pratt")
    if (length(type)==0) 
        cat("No relative importance metrics available. Spelling error in function call?")
    else 
        print(matrix(cbind(x@lmg,x@pmvd,x@last,x@first,x@betasq,x@pratt),p,length(type),
            dimnames=list(x@namen[2:(p+1)],type)))
}

