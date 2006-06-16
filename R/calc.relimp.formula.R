calc.relimp.formula <- function(formula, data, na.action, ..., subset=NULL){
    if ( missing(formula) || (max(attr(terms(formula),"order")) != 1) || 
                 (attr(terms(formula),"response") != 1) ) 
        stop("formula missing or incorrect")
    if (missing(na.action)) 
        na.action <- getOption("na.action")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())

    terms <- attr(mf,"terms")
    resp <- attr(terms,"response")
    if (max(attr(terms,"order")) != 1) stop ("formula contains higher order terms")
    if (attr(terms,"intercept") != 1) stop ("model must contain intercept")

    if (any(as.logical(lapply(mf,is.factor)))) stop("The model must not contain any factors!")
    if (!is.null(attr(mf,"na.action"))) warning(naprint(attr(mf,"na.action")))
    if (!is.null(dim(model.response(mf)))){
          if (ncol(model.response(mf))>1) stop("too many response variables")
          }
    ## selection of columns from model needed because of e.g. lm(y~x1+x2+x3-x2)
    ## selection of columns from model based on formula below 
    ##               works even in case of multi-column terms such as poly(x2,3)

    COVA <- cov(data.frame(mf[,c(resp,which(rowSums(attr(terms,"factors"))>0))]))

    y <- do.call("calc.relimp", list(COVA, ...))
    y
}

