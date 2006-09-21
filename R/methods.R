
#classes and methods in the following
#Method for coercing x of class relimplm to list 
#with as(x,"list")
setAs("relimplm","list",as.relimplm)
setMethod("show",signature(object="relimplm"),function(object) print.relimplm(object))
setMethod("show",signature(object="relimplmbooteval"),function(object) print.relimplmbooteval(object))
setMethod("show",signature(object="relimplmboot"),function(object) 
    {
    cat("Objects of class relimplmboot should not be printed.", "\n", 
        "To see their structure, use the function str().", "\n")
    })
setMethod("print",signature(x="relimplm"),function(x) print.relimplm(x))
setMethod("print",signature(x="relimplmbooteval"),function(x) print.relimplmbooteval(x))
setMethod("print",signature(x="relimplmboot"),function(x) 
    {
    cat("Objects of class relimplmboot should not be printed.", "\n", 
        "To see their structure, use the function str().", "\n")
    })

