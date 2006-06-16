"boot.relimp.default" <-
function (object, x = NULL, ..., b = 1000, type = "lmg", rank = TRUE, diff = TRUE, 
    rela = FALSE, always = NULL, fixed = FALSE) 
{
   y <- object
    # Author and copyright holder: Ulrike Groemping

    #This routine is distributed under GPL version 2 or newer.
    #The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    # function for simulating percentage contribution and rank with CIs
    # result (ergebnis) contains estimated percentages, estimated ranks, CIs for percentages, CIs for ranks

    # object can be one of the following:
    #  - the variance-covariance matrix of y with all regressor variables (y in first position),
    #    which needs to be square, at least 3x3 (= 2 regressors), and positive definite
    #   - a response variable, 
    #   - a data frame with the response in the first column (no factors allowed)
    #   - a matrix with the response in the first column 
    # x is
    #   - a matrix or data frame of x-variables if y is a response variable 
    #   - NULL otherwise
    # b is number of bootstrap replicates
    # rank (if TRUE) requests bootstrapping ranks
    # diff (if TRUE) requests bootstrapping differences
    # rela (if TRUE) requests forcing to sum 100%
    # type is a character vector or list (or single character string) of relative importance 
    #    types for which bootstrapping is requested
    # always gives regressors to adjust for
    # fixed (if TRUE) requests fixed design matrix bootstrapping

    xcall <- x    

    #error control and initialization
    if (is.null(xcall)) {
       if (!is.data.frame(y) && !is.matrix(y)) 
           stop("If x is NULL, then object must be a data frame or a matrix.")
       if (is.data.frame(y) && any(as.logical(lapply(y,is.factor))))
           stop("Data frame must not contain factors.")
       if (is.matrix(y) && !is.numeric(y))
           stop("Matrix object must be numeric.")
        names <- colnames(y)
        n <- nrow(y)     
       y <- y[complete.cases(y),] 
       nobs <- nrow(y) 
             if (!nobs==n) 
                 warning(paste((n-nobs), "observations deleted due to missing"))
             if (!(nobs > ncol(y)+2)) 
                stop("Too few complete observations for estimating this model")
       x <- y[,2:ncol(y)]
       y <- matrix(y[,1])
       if (is.null(names)) names <- c("y", paste("X",1:(ncol(object)-1),sep=""))
       colnames(x) <- names[2:length(names)]
       colnames(y) <- names[1]
       ### !length(names)==ncol(object) prevents errors for cases with e.g. polynomial regression,
       ### where the colnames vector is shorter than the number of columns (need to find out how to treat this)
    }

    if (!is.null(xcall)) {
        if (!(length(y) == nrow(x))) 
            stop("number of rows in object and x MUST be identical")
       if (is.data.frame(x) && any(as.logical(lapply(x,is.factor))))
           stop("x must not contain factors.")
       if (is.matrix(x) && !is.numeric(x))
           stop("Matrix x must be numeric.")
       if ( (is.matrix(y) && ncol(y)>1) || !is.numeric(y) )
           stop(paste("object must be a numeric vector or one-column matrix,",
              "\n", "if x is also given.",sep=""))
        n <- length(y)
        nomiss <- complete.cases(y,x)  
        nobs <- sum(nomiss)
        if (!nobs==n) 
              warning(paste((n-nobs), "observations deleted due to missing"))
        if (!(nobs > ncol(x)+3)) 
              stop("Too few complete observations for estimating this model")
         ## name must be caught here, because it is lost after manipulating y
         ynam <- deparse(substitute(object))
         if (is.matrix(y)  & !is.null(colnames(y))) ynam <- colnames(y)
         y <- y[nomiss]
         x <- x[nomiss,]
        if (is.null(colnames(x))) 
            colnames(x) <- paste("X", 1:ncol(x), sep = "")
        #provide names for identifying statistics and columns
        names <- c(ynam, colnames(x))
        }

    if (!is.logical(rank)) 
        stop("rank must be a logical")
    if (!is.logical(diff)) 
        stop("diff must be a logical")
    if (!is.logical(rela)) 
        stop("rela must be a logical")

    alltype <- alltype()
    if (!all(type %in% alltype) && (length(alltype) == 6 || !("pmvd" %in% 
        type))) 
        stop("invalid type requested")
    if (!all(type %in% alltype) && length(alltype) == 5 && "pmvd" %in% 
        type) 
        stop("pmvd is not a valid type in this version of relaimpo, obtain the non-US version, if you want to use pmvd")

       #combine y and x and calulate covariance
    daten <- as.matrix(cbind(y, x))
    colnames(daten) <- names
    empcov <- cov(daten)

    p <- ncol(x)

    if (!is.null(always) && !length(always) == length(unique(always)))
        stop("duplicate elements in always")
    if (!is.null(always) && length(always) > p - 2)
        stop("always has too many elements, less than two regressors left")
    if (is.numeric(always) && !always == as.integer(always))
        stop("always must contain integer numbers (numeric or integer) or variable names")
    if ((is.numeric(always) || is.integer(always)) && (min(always)<2 || max(always)>p+1))
        stop(paste("Numbers in always must be between 2 and ", p+1, 
           " (1=response variable, 2=first regressor, ...).", sep=""))
    if (is.character(always) && !all(always %in% names[2:(p+1)]))
        stop("Names in always must be names of the regressor variables.")

    if (!is.null(always) && is.character(always)) always <- which(names %in% always)
      ## now always is numeric

    # prepare output object
    ausgabe <- new("relimplmboot")
    ausgabe@type <- alltype[which(alltype %in% type)]
    ausgabe@nobs <- nobs
    ausgabe@nboot <- b
    ausgabe@rank <- rank
    ausgabe@diff <- diff
    ausgabe@rela <- rela
    ausgabe@always <- always
    ausgabe@fixed <- fixed
    ausgabe@namen <- colnames(daten)

    #provide names for identifying statistics
    names <- colnames(x)[setdiff(1:p,always)]
    preduced <- p - length(always)
    diffnam <- paste(names[nchoosek(preduced, 2)[1, ]], names[nchoosek(preduced, 
        2)[2, ]], sep = "-")

    #run bootstrap
    #options for calc.relimp handed over after b
    if (!fixed) {
    booterg <- boot(daten, calcrelimp.forboot, b, type = type, 
        diff = diff, rank = rank, rela = rela, always = always)
    ausgabe@boot <- booterg
    }
    else {
       linmod <- lm(data.frame(daten),qr=FALSE,model=FALSE)
       e <- linmod$residuals
       fit <- linmod$fitted.values
       booterg <- boot(data.frame(x,fit=fit,e=e), calcrelimp.forboot.fixed, b, type = type, 
             diff = diff, rank = rank, rela = rela, always = always)
       booterg$data <- daten
       slot(ausgabe,"boot") <- booterg
    }
    return(ausgabe)
}

