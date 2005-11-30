"plot" <-
structure(function (x, y, ...) 
standardGeneric("plot"), class = structure("standardGeneric", package = "methods"), generic = structure("plot", package = "graphics"), package = "graphics", group = list(), valueClass = character(0), signature = c("x", 
"y"), default = structure(list(), methods = structure(list(ANY = structure(function (x, 
    y, ...) 
{
    if (is.null(attr(x, "class")) && is.function(x)) {
        nms <- names(list(...))
        if (missing(y)) 
            y <- {
                if (!"from" %in% nms) 
                  0
                else if (!"to" %in% nms) 
                  1
                else if (!"xlim" %in% nms) 
                  NULL
            }
        if ("ylab" %in% nms) 
            plot.function(x, y, ...)
        else plot.function(x, y, ylab = paste(deparse(substitute(x)), 
            "(x)"), ...)
    }
    else UseMethod("plot")
}, class = structure("derivedDefaultMethod", package = "methods"), target = structure(character(0), .Names = character(0), class = structure("signature", package = "methods")), defined = structure(character(0), .Names = character(0), class = structure("signature", package = "methods")))), .Names = "ANY"), argument = quote(x), allMethods = list(), class = structure("MethodsList", package = "methods")), skeleton = quote(function (x, 
    y, ...) 
{
    if (is.null(attr(x, "class")) && is.function(x)) {
        nms <- names(list(...))
        if (missing(y)) 
            y <- {
                if (!"from" %in% nms) 
                  0
                else if (!"to" %in% nms) 
                  1
                else if (!"xlim" %in% nms) 
                  NULL
            }
        if ("ylab" %in% nms) 
            plot.function(x, y, ...)
        else plot.function(x, y, ylab = paste(deparse(substitute(x)), 
            "(x)"), ...)
    }
    else UseMethod("plot")
}(x, y, ...)))
