.onLoad <- function(lib, pkg){
   if(!require(MASS, quietly = TRUE))
        stop("Could not load package MASS")
   if(!require(boot, quietly = TRUE))
        stop("Could not load package boot")
   if(!require(methods, quietly = TRUE))
        stop("Could not load package methods")

#cat("This is the non-US version of package relaimpo including the metric pmvd.","\n")
#cat("Please make sure that you are entitled to using it.","\n")
#cat("If you are a US-user, please use the global version (without pmvd) that is available on CRAN.","\n")


cat("This is the global version of package relaimpo.","\n")
cat( "If you are a non-US user, a version with the interesting additional metric pmvd is available","\n")
cat("from Ulrike Groempings web site at prof.tfh-berlin.de/groemping/.","\n")

}

