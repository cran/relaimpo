calc.relimp.lm <- function(object, ...){
    lm<-object
    if ( missing(lm) ) stop("object missing or incorrect")
    if (is.null(lm$terms)) stop ("object does not contain a terms component")
    if (!is(lm,"lm")) stop("object is not of class lm")
    if (length(lm$xlevels)>0) stop("The model must not contain any factors!")
    if (is(lm,"mlm")) stop("relaimpo does not work on multiresponse models")
    if (is(lm,"glm")) stop("relaimpo works on linear models only (not glm, but lm)")
    if (!is.null(lm$weights)) stop("relaimpo currently does not work on weighted models")
    terms <- lm$terms
    resp <- attr(terms,"response")
    if (max(attr(terms,"order")) != 1) stop ("model contains higher order terms")
    if (attr(terms,"intercept") != 1) stop ("model must contain intercept")

    ## selection of columns from model needed because of e.g. lm(y~x1+x2+x3-x2)
    ## selection of columns from model based on formula below 
    ##               works even in case of multi-column terms such as poly(x2,3)

    ## as.matrix makes sure that things work correctly even for multi-column effects like polynomials
 
    DATA <- as.matrix(lm$model[,c(resp,which(rowSums(attr(terms,"factors"))>0))])

    y <- do.call("calc.relimp", list(DATA, ...))
    y
}
