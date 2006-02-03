"pmvdcalc" <-
function (ausgabe, p, indices, variances, rank, diff, rela) 
{
    # Author and copyright holder: Ulrike Groemping

    # This routine is distributed under GPL version 2, with the following geographical restriction:
    #   Distribution  o u t s i d e   t h e   U S   o n l y   !
    # The text of GPL an be found at http://www.gnu.org/copyleft/gpl.html.
    # The geographical restriction is declared in accordance with paragraph 8 of GPL version 2, 
    # which is cited here for convenience:

    #  8. If the distribution and/or use of the Program is restricted in certain countries 
    #  either by patents or by copyrighted interfaces, the original copyright holder who places 
    #  the Program under this License may add an explicit geographical distribution limitation 
    #  excluding those countries, so that distribution is permitted only in or among countries 
    #  not thus excluded. In such case, this License incorporates the limitation as if written 
    #  in the body of this License. 


    # program that calculates pmvd
    # based on the recursive potential definition of the proportional value
    # worth of coalition defined as variance reduction when coalition withdraws from full model

    # structured similarly to lmg
    prec <- rep(list(0), p)
    anz = 0
    indicesneu <- indices
    variancesneu <- variances
    # contributions of each variable in last position
    # w({i})
    prec[[1]] <- rev(variances[[p]]) - as.numeric(variances[[p + 
        1]])
    if (any(prec[[1]] == 0)) {
        nullen <- which(prec[[1]] == 0)
        anz = length(nullen)
        for (i in 2:(p - anz + 1)) {
            checkmat <- matrix(is.element(indicesneu[[i]], nullen), 
                i - 1, NCOL(indicesneu[[i]]))
            spalten <- which(colSums(checkmat) == 0)
            indicesneu[[i]] <- matrix(indicesneu[[i]][, spalten], 
                i - 1, length(spalten))
            variancesneu[[i]] <- variancesneu[[i]][, spalten]
        }
        # shorten prec[[1]] to nonnull contributors
        prec[[1]] <- prec[[1]][setdiff(1:p, nullen)]
    }
    for (i in 2:(p - anz)) {
        # indices[[2]] contains all individual variables, corresponds to prec[[1]]
        nset <- NCOL(indicesneu[[i + 1]])
        if (i == p - anz) 
            indicesneu[[i + 1]] <- matrix(indicesneu[[i + 1]], 
                p - anz, 1)
        # first indices for k together with last indices for p-k elements make up all p elements
        # determine w(S)
        prec[[i]] <- rev(variancesneu[[p - anz - i + 1]]) - as.numeric(variancesneu[[p + 
            1 - anz]])
        for (j in 1:nset) {
            # determine sum of relevant prec entries of preceding row
            spalte <- which(colSums(matrix(is.element(indicesneu[[i]], 
                indicesneu[[i + 1]][, j]), i - 1, NCOL(indicesneu[[i]]))) == 
                i - 1)
            # previous models leaving out one of the current variables
            prec[[i]][j] <- prec[[i]][j]/sum(1/(prec[[i - 1]][spalte]))
        }
    }
    # matrix P is finished
    pmvd <- rev(as.numeric(prec[[p - anz]])/prec[[p - anz - 1]])
    if (anz > 0) {
        index <- 1
        pmvdnonnull <- pmvd
        pmvd <- rep(0, p)
        for (i in 1:p) {
            if (!is.element(i, nullen)) {
                pmvd[i] <- pmvdnonnull[index]
                index <- index + 1
            }
        }
    }

    # normalize
    if (rela) 
        pmvd <- pmvd/sum(pmvd)
    else pmvd <- pmvd/variances[[1]]

    # ranking
    raengepmvd <- p + 1 - rank(pmvd)

    # pairwise differences
    if (diff & p > 2) 
        diffpmvd <- pmvd[nchoosek(p, 2)[1, ]] - pmvd[nchoosek(p, 
            2)[2, ]]
    if (diff & p == 2) 
        diffpmvd <- pmvd[1] - pmvd[2]

    # output results
    slot(ausgabe, "pmvd") <- pmvd
    if (rank) 
        slot(ausgabe, "pmvd.rank") <- raengepmvd
    if (diff) 
        slot(ausgabe, "pmvd.diff") <- diffpmvd
    return(ausgabe)
}

