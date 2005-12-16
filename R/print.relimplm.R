"print.relimplm" <-
function (x, ...) 
{
    if (!(is(x, "relimplm"))) 
        stop("x must be the output from function calc.relimp")
    p <- length(slot(x, "namen")) - 1
    cat("Response variable:", slot(x, "namen")[1], "\n")
    cat("Total response variance:", x@var.y, "\n")
    cat(p, "Regressors:", paste(slot(x, "namen")[2:(p + 1)], 
        collapse = " "), "\n")
    cat("Proportion of variance explained by model: ", round(100 * 
        x@R2, 2), "%", "\n", sep = "")
    cat("\n")
    cat("Relative importance metrics:", "\n")
    type <- slot(x, "type")
    print(matrix(cbind(x@lmg, x@pmvd, x@last, x@first, x@betasq, 
        x@pratt), p, length(type), dimnames = list(x@namen[2:(p + 
        1)], type)))
}

