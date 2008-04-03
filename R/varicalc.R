"varicalc" <- 
function (type, alle, covg, p, indices, variances, g, groups, ngroups=NULL, WW=NULL) 
{
if (!(is.null(ngroups) || length(groups)==length(ngroups))) stop ("unexpected error occurred in varicalc")
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
###bevor man mit den Gruppen arbeiten kann, muss man zunaechst eine neue Reihenfolge festlegen
###der Einfachheit halber die Gruppen zuerst, dann die ungruppierten
###groups kann nach Erstellung der Dokumentationsmatrix durch die Einzelspalten ergaenzt werden
        varjetzt <- matrix(0, 1, choose(g, k))
        for (j in 1:(choose(g, k))) {
            if (is.null(groups)) diese <- jetzt[,j] + 1
            else diese <- list2vec(groups[jetzt[,j]])

            andere <- setdiff(alle, diese)
            varjetzt[j] <- (covg[andere, andere] - covg[andere, 
                diese] %*% solve(covg[diese, diese], matrix(covg[diese, 
                andere], length(diese), p + 1 - length(diese))))[1, 1]
        }
        variances[[k + 1]] <- varjetzt
    }
    ## change by ML
    if (length(WW[[1]]) > 0) {
      for (j in 2:(length(indices)-1)) {
        if (!is.null(ngroups)) WWc <- apply(indices[[j]],2,checkWW,WW[[1]],ngroups)
        else WWc <- apply(indices[[j]],2,checkWW,WW[[1]])
        variances[[j]] <- variances[[j]][,WWc]
        if (is.vector(indices[[j]][,WWc])) {
           indices[[j]] <- matrix(indices[[j]][,WWc], ncol = 1)
        } else {
          indices[[j]] <- indices[[j]][,WWc]
      }
      }
          indices[[2]] <- matrix(indices[[2]], nrow = 1)
    }
    ## change by ML end
    return(list(indices = indices, variances = variances))
}

