"boot.relimp.default" <-
function (object, x = NULL, ..., b = 1000, type = "lmg", rank = TRUE, diff = TRUE, 
    rela = FALSE, always = NULL, groups = NULL, groupnames = NULL, fixed = FALSE) 
{
   y <- object
   ogroupnames <- groupnames
   ogroups <- groups
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

    alle <- 1:(p+1)
    alwaysnam <- NULL
    if (!is.null(always)) {
        if (is.character(always)) always <- which(names %in% always) 
        alwaysnam <- names[always]
        ## now always is numeric
        andere <- setdiff(alle, always)
        }
    g <- p - length(always)
    if (!is.null(groups)) 
    {
    if (any(c("betasq","pratt") %in% type)) stop("Metrics betasq and pratt do not work with groups.") 
    if (!is.list(groups)) {
        if (is.character(groups)) 
            groups <- which(names %in% groups)
        if ((is.numeric(groups) || is.integer(groups)) && !all(groups %in% 2:(p+1)))
            stop(paste("Numbers in groups must refer to columns 2 to ", p+1, " in cov(Y,X1,...,Xp)", sep=""))
       if (length(groups) <= 1) 
            stop("groups must list groups of more than one regressor.")
       groups=list(groups)
        }
    else {
       groups <- lapply(groups, function(obj){        
        if (is.character(obj)) 
            obj <- which(names %in% obj)
        if ((is.numeric(obj) || is.integer(obj)) && !all(obj %in% 2:(p+1)))
            stop(paste("Numbers in elements of groups must refer to columns 2 to ", p+1, " in cov(Y,X1,...,Xp).", sep=""))
        if (length(obj) <= 1) 
            stop("Each element of groups must contain more than one regressor.")  
        obj } )

       if (!length(list2vec(groups))==length(unique(list2vec(groups))))
            stop("Overlapping groups are not permitted!")
       }

    if (any(always %in% list2vec(groups))) 
       stop("groups must not refer to regressors that also occur in always.")

    if (!is.null(groupnames)) { 
       if (!length(groups)==length(groupnames)) 
       stop(paste("groupnames must have one entry for each group.", "\n", 
           "There are", length(groups), "groups and", length(groupnames), "group names.")) }
    else groupnames <- paste("G", 1:length(groups), sep="") 
       ## groupdocu will support printing the meaning of groups
       groupdocu <- list(groupnames, lapply(groups, function(obj){names[obj]}))
                 ### use correct columns of x-matrix
       groupdocu[[1]] <- append(groupdocu[[1]],  as.list(names[setdiff(2:(p+1), c(list2vec(groups),always))]))
       groupdocu[[2]] <- append(groupdocu[[2]],  as.list(names[setdiff(2:(p+1), c(list2vec(groups),always))]))

       ## groups will be used for picking appropriate elements from covariance matrix
       ## g ist number of groups
       groupnames <- c(groupnames, names[setdiff(alle, c(1,list2vec(groups),always))]) 
       groups <- append(groups, as.list(setdiff(alle, c(1,list2vec(groups),always))))

## oder
##       groupdocu <- list(as.list(groupnames), lapply(groups, function(obj){names[obj]}))
##       names <- cbind(groupnames, names[setdiff(alle, list2vec(groups))]) 
##       groups <- append(groups, as.list(setdiff(2:(p+1), list2vec(groups))))

       g <- length(groups)
       }    ## end if !is.null(groups)
##    ## always-Manipulationen, die nach groups-Abschnitt erfolgen sollen
    if (!is.null(always)) {
       names <- names[andere]
       p <- p - length(always)
       alle <- 1:(p+1)
##    ### ??? !!! hier groups bekommt zu viel abgezogen, wenn mehr als ein Element!!!
    if (!is.null(groups)) {
         groups <- lapply(groups, function(obj){
             obj - rowSums(matrix(obj,length(obj),length(always),byrow=F)>matrix(always,length(obj),
             length(always),byrow=T))
         } 
         )
         names <- c(names[1],as.character(groupnames))
         }
         }
##    ## always==obj impossible because of error checking

    
    
    # prepare output object
    ausgabe <- new("relimplmboot")
    ausgabe@type <- alltype[which(alltype %in% type)]
    ausgabe@nobs <- nobs
    ausgabe@nboot <- b
    ausgabe@rank <- rank
    ausgabe@diff <- diff
    ausgabe@rela <- rela
    ausgabe@always <- always
    ausgabe@alwaysnam <- alwaysnam
    ausgabe@fixed <- fixed
    ausgabe@namen <- names    # variable names of y and decomposition variables
    if (!is.null(groups)) ausgabe@groupdocu <- groupdocu

    #provide names for identifying statistics
##    names <- colnames(x)[setdiff(1:p,union(always, groups))]
##    if (!is.null(groups)) names=c(groupnames,names)  # groupnames needs to have g elements

##    preduced <- p - length(always)      # no. of variables used in calculations
##    g <- length(names)                  # no. of groups, potentially lower than preduced (not higher) 

##    diffnam <- paste(names[2:(g+1)][nchoosek(g, 2)[1, ]], names[2:(g+1)][nchoosek(g, 
##        2)[2, ]], sep = "-")

    #run bootstrap
    #options for calc.relimp handed over after b
    if (!fixed) {
    booterg <- boot(daten, calcrelimp.forboot, b, type = type, 
        diff = diff, rank = rank, rela = rela, always = always, groups=ogroups, groupnames=ogroupnames)
    ausgabe@boot <- booterg
    }
    else {
       linmod <- lm(data.frame(daten),qr=FALSE,model=FALSE)
       e <- linmod$residuals
       fit <- linmod$fitted.values
       booterg <- boot(data.frame(x,fit=fit,e=e), calcrelimp.forboot.fixed, b, type = type, 
             diff = diff, rank = rank, rela = rela, always = always, groups=ogroups, groupnames=ogroupnames)
       booterg$data <- daten
       slot(ausgabe,"boot") <- booterg
    }
##       ## groupdocu will support printing the meaning of groups
##       if (!is.null(groups)) {
##       groupdocu <- list(groupnames, lapply(groups, function(obj){ausgabe@namen[obj]}))
##       groupdocu[[1]] <- append(groupdocu[[1]],  as.list(ausgabe@namen[setdiff(2:(p+1), 
##                 union(list2vec(groups),always))]))
##       groupdocu[[2]] <- append(groupdocu[[2]],  as.list(ausgabe@namen[setdiff(2:(p+1), 
##                 union(list2vec(groups),always))]))
##       }

##       names <- c(ausgabe@namen[1],groupnames, ausgabe@namen[setdiff(2:(p+1), union(list2vec(groups,always))])
##       groups <- append(groups, as.list(setdiff(2:(p+1), union(list2vec(groups),always) ) ) )
##       g <- length(groupdocu[[2]])


##    ausgabe@namen <- names
##    ausgabe@covg <- cov(daten)

    return(ausgabe)
}

