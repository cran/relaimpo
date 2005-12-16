"varicalc" <- 
function (type, alle, covg, p, indices, variances) 
{
    # Author and copyright holder: Ulrike Groemping

    #This routine is distributed under GPL version 2 or newer.
    #The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    #routine for calculating all needed residual variances
    #together with indices denoting the variables included in the model
    liste <- c(1, p - 1)
    if (any(c("lmg", "pmvd") %in% type)) 
        liste <- 1:(p - 1)
    for (k in liste) {
        jetzt <- nchoosek(p, k)
        indices[[k + 1]] <- jetzt
        varjetzt <- matrix(0, 1, choose(p, k))
        for (j in 1:(choose(p, k))) {
            diese <- jetzt[, j] + 1
            andere <- setdiff(alle, diese)
            varjetzt[j] <- (covg[andere, andere] - covg[andere, 
                diese] %*% solve(covg[diese, diese], matrix(covg[diese, 
                andere], k, p + 1 - k)))[1, 1]
        }
        variances[[k + 1]] <- varjetzt
    }
    return(list(indices = indices, variances = variances))
}

