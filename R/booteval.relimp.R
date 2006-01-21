"booteval.relimp" <-
function (bootrun, bty = "bca", level = 0.95, sort = FALSE, norank = FALSE, 
    nodiff = FALSE, typesel = c("lmg", "pmvd", "last", "first", 
        "betasq", "pratt")) 
{
    # Author and copyright holder: Ulrike Groemping

    # This routine is distributed under GPL version 2 or newer.
    # The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    # function for simulating percentage contribution and rank with CIs
    # result (ergebnis) contains estimated percentages, estimated ranks, CIs for percentages, CIs for ranks

    # bootrun is the result of bootstrap runs from boot.relimp
    # bty is bty of bootstrap intervals, default BCa, 
    #       eventually intended that list can be given like in package boot
    #       currently bty only takes one single bty
    #      (rank intervals always with bty percentile, 
    #       since BCa does not work properly with ranks and normal not reasonable)
    # level is the confidence level, list can be given like in package boot
    # sort (if TRUE) requests that output is sorted by size of relative importances
    # norank (if TRUE) requests suppression of results for ranks (although ranks have been bootstrapped)
    # nodiff (if TRUE) requests suppression of results for differences (although differences have been bootstrapped)
    # typesel is the list of metric types for which evaluation is requested

    # error control
    if (!(is(bootrun, "relimplmboot"))) 
        stop("bootrun must be output from routine boot.relimp")
    if (!bty %in% c("perc", "bca", "norm", "basic")) 
        stop("bty type MUST be one of ", "perc ", "bca ", "norm ", 
            "basic")
    if (!(min(level) >= 0.5 && max(level) < 1)) 
        stop("invalid confidence levels requested: ", paste(level, 
            collapse = " "))
    type <- bootrun@type
    rank <- bootrun@rank
    diff <- bootrun@diff
    rela <- bootrun@rela

    #combine y and x and calulate covariance
    empcov <- cov(bootrun@boot$data)

    p <- ncol(empcov) - 1
    nlev <- length(level)

    # prepare output object by first providing estimates themselves
    ausgabe <- calc.relimp(empcov, type = type, diff = diff, 
        rank = rank, rela = rela)
    # extend output object
    class(ausgabe) <- "relimplmbooteval"
    ausgabe@level <- level
    ausgabe@nboot <- bootrun@nboot
    ausgabe@bty <- bty
    ausgabe@rank <- rank
    ausgabe@diff <- diff
    #provide names for identifying statistics
    names <- colnames(bootrun@boot$data)[2:(p + 1)]
    diffnam <- paste(names[nchoosek(p, 2)[1, ]], names[nchoosek(p, 
        2)[2, ]], sep = "-")
    ausgabe@var.y.boot <- bootrun@boot$t[, 1]
    ausgabe@R2.boot <- bootrun@boot$t[, 2]


    #assign names to elements from bootrun in order to be able to refer to them later
    #columns of matrices can be referred to by their colnames (in quotes in square brackets instead of index)
    #elements of vectors analogously
    zaehl <- 3
    bootnames <- c("var.y", "R2")
    typname <- ""
    for (a in c("lmg", "pmvd", "last", "first", "betasq", "pratt")) {
        if (a %in% type) {
            bootnames <- c(bootnames, paste(names, ".", a, sep = ""))
            if (a %in% typesel) 
                slot(ausgabe, paste(a, "boot", sep = ".")) <- bootrun@boot$t[, 
                  zaehl:(zaehl + p - 1)]
            zaehl <- zaehl + p
            if (rank) {
                bootnames <- c(bootnames, paste(names, ".", a, 
                  "rank", sep = ""))
                if (a %in% typesel) 
                  slot(ausgabe, paste(a, "rank", "boot", sep = ".")) <- bootrun@boot$t[, 
                    zaehl:(zaehl + p - 1)]
                zaehl <- zaehl + p
            }
            if (diff) {
                bootnames <- c(bootnames, paste(diffnam, ".", 
                  a, "diff", sep = ""))
                if (a %in% typesel) 
                  slot(ausgabe, paste(a, "diff", "boot", sep = ".")) <- matrix(bootrun@boot$t[, 
                    zaehl:(zaehl + p * (p - 1)/2 - 1)],1,p * (p - 1)/2)
                zaehl <- zaehl + p * (p - 1)/2
            }
            if (a %in% typesel && !typname[1] == "") 
                typname <- c(typname, a)
            if (a %in% typesel && typname[1] == "") 
                typname <- a
        }
    }
    colnames(bootrun@boot$t) <- bootnames
    names(bootrun@boot$t0) <- bootnames
    ntype <- length(typname)
    percentages <- bootrun@boot$t0[paste(names, ".", matrix(typname, 
        p, ntype, byrow = T), sep = "")]
    if (rank) 
        ranks <- bootrun@boot$t0[paste(names, ".", matrix(typname, 
            p, ntype, byrow = T), "rank", sep = "")]
    if (diff) 
        diffs <- bootrun@boot$t0[paste(diffnam, ".", matrix(typname, 
            p * (p - 1)/2, ntype, byrow = T), "diff", sep = "")]
    perclower <- matrix(0, nlev, ntype * p, dimnames = list(level, 
        names(percentages)))
    if (rank) 
        ranklower <- matrix(0, nlev, ntype * p, dimnames = list(level, 
            names(percentages)))
    if (diff) 
        difflower <- matrix(0, nlev, ntype * p * (p - 1)/2, dimnames = list(level, 
            names(diffs)))
    percupper <- matrix(0, nlev, ntype * p, dimnames = list(level, 
        names(percentages)))
    if (rank) 
        rankupper <- matrix(0, nlev, ntype * p, dimnames = list(level, 
            names(percentages)))
    if (diff) 
        diffupper <- matrix(0, nlev, ntype * p * (p - 1)/2, dimnames = list(level, 
            names(diffs)))

    #strategy: if all bootstrap samples have same value-> lower=upper=that value
    #otherwise: boot.ci
        #determine confidence intervals
        #percentages
     for (j in 1:length(percentages)) {
        var <- var(bootrun@boot$t[, names(percentages)[j]])
        perclower[, j] <- percentages[j]
        percupper[, j] <- percentages[j]
        if (var > 0) {
            temp <- boot.ci(bootrun@boot, index = names(percentages)[j], 
                type = bty, conf = level)
            if (bty %in% c("perc", "bca", "basic")) {
                eval(parse(text = paste("perclower[,j]<-temp$", 
                  bty, "[,4]", sep = ""), n = 1))
                eval(parse(text = paste("percupper[,j]<-temp$", 
                  bty, "[,5]", sep = ""), n = 1))
            }
            if (bty == "norm") {
                eval(parse(text = paste("perclower[,j]<-temp$", 
                  "normal", "[,2]", sep = ""), n = 1))
                eval(parse(text = paste("percupper[,j]<-temp$", 
                  "normal", "[,3]", sep = ""), n = 1))
            }
        }
    }
        #determine confidence intervals
        #ranks
     if (rank && !norank) {
        for (j in 1:length(ranks)) {
            var <- var(bootrun@boot$t[, names(ranks)[j]])
            ranklower[, j] <- ranks[j]
            rankupper[, j] <- ranks[j]
            if (var > 0) {
                temp <- boot.ci(bootrun@boot, index = names(ranks)[j], 
                  type = "perc", conf = level)
                ranklower[, j] <- temp$percent[, 4]
                rankupper[, j] <- temp$percent[, 5]
            }
        }
    }
        #determine confidence intervals
        #diffs
     if (diff) {
        for (j in 1:length(diffs)) {
            var <- var(bootrun@boot$t[, names(diffs)[j]])
            difflower[, j] <- diffs[j]
            diffupper[, j] <- diffs[j]
            if (var > 0) {
                temp <- boot.ci(bootrun@boot, index = names(diffs)[j], 
                  type = bty, conf = level)
                if (bty %in% c("perc", "bca", "student", "basic")) {
                  eval(parse(text = paste("difflower[,j]<-temp$", 
                    bty, "[,4]", sep = ""), n = 1))
                  eval(parse(text = paste("diffupper[,j]<-temp$", 
                    bty, "[,5]", sep = ""), n = 1))
                }
                if (bty == "norm") {
                  eval(parse(text = paste("difflower[,j]<-temp$", 
                    "normal", "[,2]", sep = ""), n = 1))
                  eval(parse(text = paste("diffupper[,j]<-temp$", 
                    "normal", "[,3]", sep = ""), n = 1))
                }
            }
        }
    }
    ausgabe@type <- typname

    #write confidence bounds to output object
    for (a in typname) {
        slot(ausgabe, paste(a, "lower", sep = ".")) <- matrix(perclower[, 
            paste(names, a, sep = ".")], nlev, p)
        if (rank && !norank) 
            slot(ausgabe, paste(a, "rank", "lower", sep = ".")) <- matrix(ranklower[, 
                paste(names, a, sep = ".")], nlev, p)
        if (diff && !nodiff) 
            slot(ausgabe, paste(a, "diff", "lower", sep = ".")) <- matrix(difflower[, 
                paste(diffnam, paste(a, "diff", sep = ""), sep = ".")], 
                nlev, p * (p - 1)/2)
        slot(ausgabe, paste(a, "upper", sep = ".")) <- matrix(percupper[, 
            paste(names, a, sep = ".")], nlev, p)
        if (rank && !norank) 
            slot(ausgabe, paste(a, "rank", "upper", sep = ".")) <- matrix(rankupper[, 
                paste(names, a, sep = ".")], nlev, p)
        if (diff && !nodiff) 
            slot(ausgabe, paste(a, "diff", "upper", sep = ".")) <- matrix(diffupper[, 
                paste(diffnam, paste(a, "diff", sep = ""), sep = ".")], 
                nlev, p * (p - 1)/2)
    }
    #show confidence intervals with rank marks
    #only possible, if ranks are bootstrapped
       #initialize character matrix for showing sorted results with confidence info
    if (rank && !norank) 
        mark <- matrix(rep("", (p * ntype + ntype - 1) * (3 * 
            nlev + 1)), p * ntype + ntype - 1, 3 * nlev + 1, 
            dimnames = list(rep("", p * ntype + ntype - 1), c("percentage", 
                rep(level, 3))))
    else mark <- matrix(rep(0, (p * ntype + ntype - 1) * (2 * 
        nlev + 1)), p * ntype + ntype - 1, 2 * nlev + 1, dimnames = list(rep("", 
        p * ntype + ntype - 1), c("percentage", rep(level, 2))))
    if (sort) 
        marksort <- mark
    for (aa in 1:ntype) {
        a <- typname[aa]
        percent <- slot(ausgabe, a)
        names(percent) <- names
        sortiert <- sort(percent, decreasing = T, index = T)
        percsort <- sortiert$x
        names(percsort) <- names[sortiert$ix]
        cilower <- matrix(slot(ausgabe, paste(a, "lower", sep = ".")), 
            nlev, p)
        ciupper <- matrix(slot(ausgabe, paste(a, "upper", sep = ".")), 
            nlev, p)
        if (rank && !norank) {
            lower <- matrix(slot(ausgabe, paste(a, "rank", "lower", 
                sep = ".")), nlev, p)
            upper <- matrix(slot(ausgabe, paste(a, "rank", "upper", 
                sep = ".")), nlev, p)
            for (j in 1:p) {
                # j is the rank that might or might not be in the confidence interval 
                # for the k-th sorted X
             hilf <- matrix(rep("", p * nlev), p, nlev)
                for (k in 1:p) {
                    #k is the k-th X, i.e. the k-th entry in percent for this type
                   for (i in 1:nlev) {
                         #i is the confidence level index
                      if (j < lower[i, k] | j > upper[i, k]) 
                      hilf[k, i] <- "_"
                    else hilf[k, i] <- LETTERS[j]
                  } # loop i
                } # loop k
                # append latest letter to mark
                mark[((aa - 1) * (p + 1) + 1):(aa * (p + 1) - 
                  1), 2:(1 + nlev)] <- matrix(paste(mark[((aa - 
                  1) * (p + 1) + 1):(aa * (p + 1) - 1), 2:(1 + 
                  nlev)], matrix(hilf, p, nlev), sep = ""), p, 
                  nlev)
            } # loop j
            mark[((aa - 1) * (p + 1) + 1):(aa * (p + 1) - 1), 
                ] <- matrix(cbind(percent, mark[((aa - 1) * (p + 
                1) + 1):(aa * (p + 1) - 1), 2:(1 + nlev)], t(cilower), 
                t(ciupper)), p, 3 * nlev + 1)
        } #if rank and !norank

        if (!rank || norank) {
            mark[((aa - 1) * (p + 1) + 1):(aa * (p + 1) - 1), 
                ] <- matrix(cbind(percent, t(cilower), t(ciupper)), 
                p, 2 * nlev + 1)
        }
        rownames(mark)[((aa - 1) * (p + 1) + 1):(aa * (p + 1) - 
            1)] <- paste(names(percent), a, sep = ".")
        if (sort) {
            marksort[((aa - 1) * (p + 1) + 1):(aa * (p + 1) - 
                1), ] <- mark[((aa - 1) * (p + 1) + 1):(aa * 
                (p + 1) - 1), ][sortiert$ix, ]
            rownames(marksort)[((aa - 1) * (p + 1) + 1):(aa * 
                (p + 1) - 1)] <- rownames(mark[((aa - 1) * (p + 
                1) + 1):(aa * (p + 1) - 1), ])[sortiert$ix]
        }
    } # loop aa

    # reduce number of displayed digits in percentages
    if (rank && !norank) 
        mark[, c(1, (2 + nlev):(3 * nlev + 1))] <- substr(mark[, 
            c(1, (2 + nlev):(3 * nlev + 1))], 1, 6)
    else mark <- round(mark, digits = 4)
    if (sort && rank && !norank) 
        marksort[, c(1, (2 + nlev):(3 * nlev + 1))] <- substr(marksort[, 
            c(1, (2 + nlev):(3 * nlev + 1))], 1, 6)
    if (sort && (!rank || norank)) 
        marksort <- round(marksort, digits = 4)
    if (sort) 
        ausgabe@mark <- marksort
    else ausgabe@mark <- mark

    # differences
    if (diff && !nodiff) {
        mark <- matrix(rep("", (p * (p - 1) * ntype/2 + ntype - 
            1) * (3 * nlev + 1)), p * (p - 1) * ntype/2 + ntype - 
            1, 3 * nlev + 1, dimnames = list(rep("", p * (p - 
            1) * ntype/2 + ntype - 1), c("difference", rep(level, 
            3))))
        if (sort) 
            marksort <- matrix(rep("", (p * (p - 1) * ntype/2 + 
                ntype - 1) * (3 * nlev + 1)), p * (p - 1) * ntype/2 + 
                ntype - 1, 3 * nlev + 1, dimnames = list(rep("", 
                p * (p - 1) * ntype/2 + ntype - 1), c("difference", 
                rep(level, 3))))
        for (aa in 1:ntype) {
            a <- typname[aa]
            differ <- slot(ausgabe, paste(a, "diff", sep = "."))
            sortiert <- sort(abs(differ), decreasing = T, index = T)
            names(differ) <- paste(diffnam, a, sep = ".")
            difflower <- matrix(slot(ausgabe, paste(a, "diff", 
                "lower", sep = ".")), nlev, p * (p - 1)/2)
            diffupper <- matrix(slot(ausgabe, paste(a, "diff", 
                "upper", sep = ".")), nlev, p * (p - 1)/2)
            hilf <- matrix(rep("", p * (p - 1) * nlev/2), p * 
                (p - 1)/2, nlev)
            for (k in 1:(p * (p - 1)/2)) {
                   # k refers to the k-th difference
                 for (i in 1:nlev) {
                       #i is the confidence level index
                   if (0 < difflower[i, k] | 0 > diffupper[i, 
                    k]) 
                    hilf[k, i] <- "*"
                  else hilf[k, i] <- " "
                } # loop i
            } # loop k

           #append latest letter to mark 
            mark[((aa - 1) * (p * (p - 1)/2 + 1) + 1):(aa * (p * 
                (p - 1)/2 + 1) - 1), 2:(1 + nlev)] <- matrix(paste(mark[((aa - 
                1) * (p * (p - 1)/2 + 1) + 1):(aa * (p * (p - 
                1)/2 + 1) - 1), 2:(1 + nlev)], matrix(hilf, p * 
                (p - 1)/2, nlev), sep = ""), p * (p - 1)/2, nlev)
            mark[((aa - 1) * (p * (p - 1)/2 + 1) + 1):(aa * (p * 
                (p - 1)/2 + 1) - 1), ] <- matrix(cbind(differ, 
                mark[((aa - 1) * (p * (p - 1)/2 + 1) + 1):(aa * 
                  (p * (p - 1)/2 + 1) - 1), 2:(1 + nlev)], t(difflower), 
                t(diffupper)), p * (p - 1)/2, 3 * nlev + 1)
            rownames(mark)[((aa - 1) * (p * (p - 1)/2 + 1) + 
                1):(aa * (p * (p - 1)/2 + 1) - 1)] <- names(differ)
            if (sort) {
                marksort[((aa - 1) * (p * (p - 1)/2 + 1) + 1):(aa * 
                  (p * (p - 1)/2 + 1) - 1), ] <- mark[((aa - 
                  1) * (p * (p - 1)/2 + 1) + 1):(aa * (p * (p - 
                  1)/2 + 1) - 1), ][sortiert$ix, ]
                rownames(marksort)[((aa - 1) * (p * (p - 1)/2 + 
                  1) + 1):(aa * (p * (p - 1)/2 + 1) - 1)] <- names(differ)[sortiert$ix]
            }
        } # loop aa

        ##reduce number of displayed digits in percentage differences
        mark[, c(1, (2 + nlev):(3 * nlev + 1))] <- substr(mark[, 
            c(1, (2 + nlev):(3 * nlev + 1))], 1, 6)
        if (!sort) 
            ausgabe@markdiff <- mark
        if (sort) {
            marksort[, c(1, (2 + nlev):(3 * nlev + 1))] <- substr(marksort[, 
                c(1, (2 + nlev):(3 * nlev + 1))], 1, 6)
            ausgabe@markdiff <- marksort
        }
    } # if diff and !nodiff

    # set correct options for printing in output object
    ausgabe@diff <- diff && !nodiff
    ausgabe@rank <- rank && !norank
    ausgabe@sort <- sort
    return(ausgabe)
}

