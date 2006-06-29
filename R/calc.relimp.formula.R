calc.relimp.formula <- function (formula, data, na.action, ..., subset = NULL) 
{
    if (missing(formula)) stop("formula missing")
    if (missing(na.action)) 
        na.action <- getOption("na.action")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    terms <- attr(mf, "terms")
    resp <- attr(terms, "response")

    if (resp != 1 ) stop("incorrect formula") 
    if (max(attr(terms, "order")) != 1) 
        stop("formula contains higher order terms")
    if (attr(terms, "intercept") != 1) 
        stop("model must contain intercept")
    if (any(as.logical(lapply(mf, is.factor)))) 
        stop("The model must not contain any factors!")
    if (!is.null(attr(mf, "na.action"))) 
        warning(naprint(attr(mf, "na.action")))
    if (!is.null(dim(model.response(mf)))) {
        if (ncol(model.response(mf)) > 1) 
            stop("too many response variables")
    }
    y <- do.call("calc.relimp", list(data.frame(mf[, c(resp, which(rowSums(attr(terms, 
        "factors")) > 0))]), ...))
    ## calling with data matrix better than with covariance matrix, because number of observations is displayed
    y
}