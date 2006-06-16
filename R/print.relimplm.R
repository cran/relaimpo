"print.relimplm" <-
function (x, ...) 
{
    if (!(is(x, "relimplm"))) 
        stop("x must be the output from function calc.relimp")
    p <- length(slot(x, "namen")) - 1
    cat("Response variable:", slot(x, "namen")[1], "\n")
    cat("Total response variance:", x@var.y, "\n")
    if (length(x@nobs)>0) cat("Analysis based on", x@nobs, "observations", "\n")
    if (!is.null(x@always)) {
      cat("\n")
      cat(p + length(slot(x, "always")), "Regressors:", "\n", sep = " ")
      cat("Proportion of variance explained: ", round(100 * 
          x@R2, 2), "%", "\n", "\n", sep = "")
      if (length(slot(x, "always")) == 1) {
         cat("One Regressor always included in model:", "\n", 
             paste(slot(x, "alwaysnam"), collapse = " "), "\n")
         cat(round(100 * (x@R2 - x@R2.decomp), 2), "%", "of variance explained by this regressor", 
          "\n", "\n", sep = " ") 
      } else {
         cat(length(slot(x, "always")), "Regressors always included in model:", "\n", 
             paste(slot(x, "alwaysnam"), collapse = " "), "\n")     
         cat(round(100 * (x@R2 - x@R2.decomp), 2), "%", "of variance explained by these", 
          "\n", "\n", sep = " ")
      } 
    cat("Relative importance of", p, "regressors assessed:", "\n",
        paste(slot(x, "namen")[2:(p + 1)], collapse = " "), "\n") 
    cat(round(100 * x@R2.decomp, 2), "%", "of variance decomposed among these", "\n",
        sep = " ")
    cat("\n")
    } 
    if (is.null(x@always))
    {
      cat("\n") 
      cat(p, "Regressors:", paste(slot(x, "namen")[2:(p + 1)], 
          collapse = " "), "\n")
      cat("Proportion of variance explained by model: ", round(100 * 
          x@R2, 2), "%", "\n", sep = "")
    }
    if (x@rela) cat("Metrics are normalized to sum to 100% (rela=TRUE).", "\n")
    else cat("Metrics are not normalized (rela=FALSE).", "\n")
    cat("\n")

    cat("Relative importance metrics:", "\n")
    cat("\n")
    type <- slot(x, "type")
    print(matrix(cbind(x@lmg, x@pmvd, x@last, x@first, x@betasq, 
        x@pratt), p, length(type), dimnames = list(x@namen[2:(p + 
        1)], type)))
}

