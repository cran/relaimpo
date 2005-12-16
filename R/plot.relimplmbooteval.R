"plot.relimplmbooteval" <- 
function (x, ..., lev = max(x@level), names.abbrev = 4) 
{
    # function shows barplots with error bars indicating confidence interval for chosen level
    # if chosen level not available, confidence interval for largest available level is produced

    #error control
    if (!(is(x, "relimplmbooteval"))) 
        stop("x must be the output from function booteval.relimp")
    if (!(is.numeric(names.abbrev))) 
        stop("names.abbrev must be a number")
    if (!(is.numeric(lev))) 
        stop("lev must be a number")

    #no good way of tilting available
    #vertical labels do not work as desired but overlap with sub text
    #horizontal plotting does not work either
    #current solution: abbreviation of names to at most names.abbrev characters, default 4

    p <- length(x@namen) - 1
    yname <- x@namen[1]
    xnames <- substr(x@namen[2:(p + 1)], 1, names.abbrev)
    level <- x@level
    pick <- which(level == lev)
    if (length(pick) == 0) {
        pick <- which.max(level)
        cat("Chosen confidence level ", 100 * lev, "% not available,", 
            "\n")
        cat("largest available level (=default) plotted instead.", 
            "\n")
        cat("Available levels: ", level, "\n", sep = " ")
    }
    ylab <- paste("Relative Importances for ", yname, sep = "")
    subtext <- paste("Bootstrap ", 100 * level[pick], "% CIs shown by vertical lines", 
        sep = "")
    type <- x@type
    if (length(type) == 0) 
        print("Nothing to plot")
    else {
        maxi <- 0
        mini <- 0
        for (a in type) {
            maxi <- max(maxi, slot(x, paste(a, "upper", sep = "."))[pick, 
                ])
            mini <- min(mini, slot(x, paste(a, "lower", sep = "."))[pick, 
                ])
        }
        axmax <- ceiling(10 * maxi)/10
        axmin <- floor(10 * mini)/10
        ntype <- length(type)
        op <- par(no.readonly = TRUE)
        if (ntype == 2) 
            par(mfrow = c(1, 2))
        if (ntype > 2 && ntype <= 4) 
            par(mfrow = c(2, 2))
        if (ntype > 2 && ntype <= 4) 
            par(mfrow = c(2, 2))
        if (ntype > 4) 
            par(mfrow = c(2, 3))
        if ("lmg" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@lmg, decreasing = T, index = T)$ix
            plt <- barplot(x@lmg[index], sub = subtext, main = "Method LMG", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab)
            segments(plt, x@lmg.lower[pick, ][index], plt, x@lmg.upper[pick, 
                ][index])
        }
        if ("pmvd" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@pmvd, decreasing = T, index = T)$ix
            plt <- barplot(x@pmvd[index], sub = subtext, main = "Method PMVD", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab)
            segments(plt, x@pmvd.lower[pick, ][index], plt, x@pmvd.upper[pick, 
                ][index])
        }
        if ("last" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@last, decreasing = T, index = T)$ix
            plt <- barplot(x@last[index], sub = subtext, main = "Method Last", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab)
            segments(plt, x@last.lower[pick, ][index], plt, x@last.upper[pick, 
                ][index])
        }
        if ("first" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@first, decreasing = T, index = T)$ix
            plt <- barplot(x@first[index], sub = subtext, main = "Method First", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab)
            segments(plt, x@first.lower[pick, ][index], plt, 
                x@first.upper[pick, ][index])
        }
        if ("betasq" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@betasq, decreasing = T, index = T)$ix
            plt <- barplot(x@betasq[index], sub = subtext, main = "Method Betasq", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab)
            segments(plt, x@betasq.lower[pick, ][index], plt, 
                x@betasq.upper[pick, ][index])
        }
        if ("pratt" %in% type) {
            index <- 1:p
            if (x@sort) 
                index <- sort(x@pratt, decreasing = T, index = T)$ix
            plt <- barplot(x@pratt[index], sub = subtext, main = "Method Pratt", 
                names.arg = xnames[index], ylim = c(axmin, axmax), 
                ylab = ylab)
            segments(plt, x@pratt.lower[pick, ][index], plt, 
                x@pratt.upper[pick, ][index])
        }
        par(op)
    }
}

