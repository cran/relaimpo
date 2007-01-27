"varicalc" <- 
function (type, alle, covg, p, indices, variances, g, groups) 
{
    # Author and copyright holder: Ulrike Groemping

    #This routine is distributed under GPL version 2 or newer.
    #The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    #routine for calculating all needed residual variances
    #together with indices denoting the variables included in the model
    liste <- c(1, g - 1)
    if (any(c("lmg", "pmvd") %in% type)) 
        liste <- 1:(g - 1)
    for (k in liste) {
        jetzt <- nchoosek(g, k)
        indices[[k + 1]] <- jetzt
###bevor man mit den Gruppen arbeiten kann, muss man zunächst eine neue Reihenfolge festlegen
###der Einfachheit halber die Gruppen zuerst, dann die ungruppierten
###groups kann nach Erstellung der Dokumentationsmatrix durch die Einzelspalten ergänzt werden
        varjetzt <- matrix(0, 1, choose(g, k))
        for (j in 1:(choose(g, k))) {
            if (is.null(groups)) diese <- jetzt[,j] + 1
            else diese <- list2vec(groups[jetzt[,j]])
               
#### ist das folgende neuer oder älter ???
##            diese <- jetzt[,j] + 1
##            cat(diese, "\n")
##            cat(groups[[1]],"\n")
##            if (!is.null(groups)) diese <- list2vec(groups[diese - 1])

            andere <- setdiff(alle, diese)
            varjetzt[j] <- (covg[andere, andere] - covg[andere, 
                diese] %*% solve(covg[diese, diese], matrix(covg[diese, 
                andere], length(diese), p + 1 - length(diese))))[1, 1]
        }
        variances[[k + 1]] <- varjetzt
    }
    return(list(indices = indices, variances = variances))
}

