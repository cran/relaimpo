"calc.relimp.default" <-
function (object, x = NULL, ..., type = "lmg", diff = FALSE, rank = TRUE, rela = FALSE, always = NULL) 
{
    y<-object
    # Author and copyright holder: Ulrike Groemping

    #This routine is distributed under GPL version 2 or newer.
    #The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    # object can be one of the following:
    #  - the variance-covariance matrix of y with all regressor variables (y in first position),
    #    which needs to be square, at least 3x3 (= 2 regressors), and positive definite
    #   - a response variable, 
    #   - a data frame with the response in the first column (no factors allowed)
    #   - a matrix with the response in the first column 
    # x is
    #   - a matrix or data frame of x-variables if y is a response variable 
    #   - NULL otherwise

    #error control and initialization

    nobs <- NULL

    if (is.null(x)) {
       if (!is.data.frame(y) && !is.matrix(y)) 
           stop("If x is NULL, then object must be a data frame or a matrix.")
       if (is.data.frame(y) && any(as.logical(lapply(y,is.factor))))
           stop("Data frame object must not contain factors.")
       if (is.matrix(y) && !is.numeric(y))
           stop("Matrix object must be numeric.")
       # not a covariance matrix
       if (!(ncol(y) == nrow(y))) {
             names <- colnames(y)
             n <- nrow(y)     
             y <- y[complete.cases(y),] 
             nobs <- nrow(y)
             if (!nobs==n) 
                 warning(paste((n-nobs), "observations deleted due to missing"))
             if (!(nobs > ncol(y)+2)) 
                stop("Too few complete observations for estimating this model")
             covg <- cov(y)
       }
       else covg <- y
    }

    if (!is.null(x)) {
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
        covg <- cov(cbind(y,x)[nomiss,])
        if (is.null(colnames(covg))) colnames(covg) <- c("y", paste("X", 1:ncol(x), sep = ""))
        if (colnames(covg)[1]=="y") colnames(covg)[1] <- deparse(substitute(object))       
        if (is.null(colnames(x))) colnames(covg)[2:ncol(covg)] <- paste("X", 1:ncol(x), sep = "")
        }


    # type is a character listing (vector incl. single character string, list, matrix) 
    #    that chooses one or more of the available relative importance metrics (in any order)
    #    The selected metrics will always appear in the result in the order given below:
    #    "lmg","pmvd","last","first","betasq","pratt" (without "pmvd" for the global version)


    # diff is a logical requesting differences between metrics if TRUE
    # rank is a logical requesting ranks of metrics if TRUE (largest=1, smallest=p)
    # rela is a logical requesting normalization to a percentage scale (sum 100%) if TRUE
        # if FALSE, all metrics are relative to var(Y), lmg, pmvd and pratt sum to R^2, 
        # first, last and betasq do not sum to R^2 but neither to 100%
    # always is a vector of integer numbers or a character vector given column positions 
    #    of regressors that are always kept in the model, i.e. are adjusted for in all calculations
    #    numbers refer to the position of the regressors in the variable list (1 refers to response),
    #    character strings to variable names

    # now, after treatment of the data frame, all covg need to be square matrices
    if (!is.matrix(covg)) stop(paste("If covg is square, ", "\n", 
          "it must be a covariance matrix (of type matrix)."))

    # check covariance matrix properties
    hilf <- eigen(covg, only.values = T)$values
    if (is.complex(hilf) || min(hilf) <= 0)
        stop(paste("covg must be", "\n", 
          "a positive definite covariance matrix", "\n",
          "or a data matrix / data frame with linearly independent columns."))

    if (!is.logical(rank)) 
        stop("rank must be a logical")
    if (!is.logical(diff)) 
        stop("diff must be a logical")
    if (!is.logical(rela)) 
        stop("rela must be a logical")

    # alltype is set in zzz.R
    alltype <- alltype()
    if (!all(type %in% alltype) && (length(alltype) == 6 || !("pmvd" %in% 
        type))) 
        stop("invalid type requested")
    if (!all(type %in% alltype) && length(alltype) == 5 && "pmvd" %in% 
        type) 
        stop("pmvd is not a valid type in this version of relaimpo, obtain the non-US version, if you want to use pmvd")

    # no of regressors
    p <- ncol(covg) - 1
    names <- colnames(covg)
    if (is.null(names)) 
        names <- c("y", paste("X", 1:p, sep = ""))

    if (!is.null(always) && !length(always) == length(unique(always)))
        stop("Duplicate elements in always")
    if (!is.null(always) && length(always) > p - 2)
        stop("always has too many elements, less than two regressors left")
    if (is.numeric(always) && !always == as.integer(always))
        stop("always must contain integer numbers (numeric or integer) or variable names.")
    if ((is.numeric(always) || is.integer(always)) && (min(always)<2 || max(always)>p+1))
        stop("Numbers in always must be between 2 and p+1, corresponding to the position of regressors in covg.")
    if (is.character(always) && !all(always %in% names[2:(p+1)]))
        stop("Names in always must come from the names of columns 2 to p+1 of covg.")
 
    #initialise output matrices
    wahr <- matrix(0, 1, p)
    #vector for setdiff
    alle <- 1:(p + 1)
    var.y <- covg[1,1]
    covall <- covg
    andere <- alle
    if ((is.numeric(always) || is.integer(always)) && !all(always %in% 2:(p+1)))
        stop(paste("Numbers in always must refer to columns 2 to ", p+1, " in covg", sep=""))


    # adjust out variables that are supposed to always stay in the model
    if (!is.null(always))
         {if (is.character(always)) always <- which(names %in% always)
             covall <- covg
             andere <- setdiff(alle, always)
             if (length(always)==1)
             covg <- covg[andere, andere] - covg[andere,always]%o%covg[always,andere]/covg[always,always]
             if (length(always)>1)
             covg <- covg[andere, andere] - covg[andere,always]%*%solve(covg[always,always],covg[always,andere])
                 var.y.distrib <- covg[1,1]  
             alwaysnam <- names[always]
             names <- names[andere]
             R2.always <- (var.y - covg[1,1])/var.y
             p <- p - length(always)
             alle <- 1:(p+1)
         }

    # start by calculating all conditional variances of y
    # construct all subsets with function nchoosek (included in relaimpo, taken from package vsn),
    # calculate conditional variances as top left corner of appropriate conditional covariance matrix
    # and write them to vectors in same position as with index matrices

    #initialise lists of index matrices and variance rows

    #first element of list unconditional, list initialized in full length
    indices <- rep(list(0), p + 1)
    variances <- rep(list(covg[1, 1]), p + 1)
    #conditioning on all variables, i.e. var=s^2
    indices[[p + 1]] <- matrix(1:p, p, 1)
    variances[[p + 1]] <- covg[1:1] - covg[1, 2:(p + 1)] %*% 
        solve(covg[2:(p + 1), 2:(p + 1)], covg[2:(p + 1), 1])
    hilf <- varicalc(type, alle, covg, p, indices, variances)

    indices <- hilf$indices
    variances <- hilf$variances

    #output R-squared in order to show the total that is subdivided
    if (!is.null(always)) ausgabe <- new("relimplm", var.y = var.y, 
        R2 = as.numeric(1 - variances[[p + 1]]/var.y), 
        R2.decomp = as.numeric(variances[[1]] - variances[[p + 1]])/var.y)
    if (is.null(always)) ausgabe <- new("relimplm", var.y = as.numeric(variances[[1]]), 
        R2 = as.numeric(1 - variances[[p + 1]]/variances[[1]]), 
        R2.decomp = as.numeric(1 - variances[[p + 1]]/variances[[1]]))

    if ("lmg" %in% type) 
        {
        ausgabe <- lmgcalc(ausgabe, p, indices, variances, rank, 
            diff, rela, var.y)
        names(ausgabe@lmg)<-names[2:(p+1)]
        if (rank) names(ausgabe@lmg.rank)<-names[2:(p+1)]
        }
    if ("pmvd" %in% type) 
        {
        ausgabe <- pmvdcalc(ausgabe, p, indices, variances, rank, 
            diff, rela)
        names(ausgabe@pmvd)<-names[2:(p+1)]
        if (rank) names(ausgabe@pmvd.rank)<-names[2:(p+1)]
        }
   if ("last" %in% type) 
        {
        ausgabe <- lastcalc(ausgabe, p, variances, rank, diff, 
            rela, var.y)
        names(ausgabe@last)<-names[2:(p+1)] 
        if (rank) names(ausgabe@last.rank)<-names[2:(p+1)]
        }
   if ("first" %in% type) 
        {
        ausgabe <- firstcalc(ausgabe, p, variances, rank, diff, 
            rela, var.y)
        names(ausgabe@first)<-names[2:(p+1)]
        if (rank) names(ausgabe@first.rank)<-names[2:(p+1)]
        }
    if ("betasq" %in% type) 
        {
        ausgabe <- betasqcalc(ausgabe, covg, p, variances, rank, 
            diff, rela, var.y)
        names(ausgabe@betasq)<-names[2:(p+1)]
        if (rank) names(ausgabe@betasq.rank)<-names[2:(p+1)]
        }
    if ("pratt" %in% type) 
        {
        ausgabe <- prattcalc(ausgabe, covg, p, rank, diff, rela, var.y)
        names(ausgabe@pratt)<-names[2:(p+1)]
        if (rank) names(ausgabe@pratt.rank)<-names[2:(p+1)]
        }

    #ausgabe contains (in this order) var.y, R2, lmg, rank.lmg, diff.lmg, 
    #                        pmvd, rank.pmvd, diff.pmvd,  (non-US version only)
    #                        last, rank.last, diff.last, first, rank.first, diff.first,
    #                                 betasq, rank.betasq, diff.betasq, pratt, rank.pratt, diff.pratt
    # as far as requested by the call
    # default: R2, lmg, rank.lmg
    # in addition, some logicals and names are included

    slot(ausgabe, "rela") <- rela
    slot(ausgabe, "namen") <- names
    if (!is.null(nobs)) slot(ausgabe, "nobs") <- nobs
     slot(ausgabe, "always") <- always
    if (!is.null(always)) slot(ausgabe, "alwaysnam") <- alwaysnam
    slot(ausgabe, "type") <- alltype[which(alltype %in% type)]
           # this cryptic approach ensures the correct order of types
           # and makes it possible for type to be a list without generating an error
    return(ausgabe)
}

