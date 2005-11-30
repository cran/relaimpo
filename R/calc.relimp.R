calc.relimp <- function(covg, type="lmg", diff=FALSE, rank=TRUE, rela=TRUE)
{
# Author and copyright holder: Ulrike Groemping

#This routine is distributed under GPL version 2 or newer.
#The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

### error control
if(!(ncol(covg)==nrow(covg))) stop("covg must be a SQUARE covariance matrix")
hilf<-eigen(covg,only.values=T)$values
if(is.complex(hilf) || min(hilf)<=0 ) stop("covg must be a positive definite covariance matrix")
if(!is.logical(rank)) stop("rank must be a logical")
if(!is.logical(diff)) stop("diff must be a logical")
if (!is.logical(rela)) stop("rela must be a logical")

#alltype is set in zzz.R
alltype <- alltype()
if(!all(type %in% alltype) && (length(alltype)==6 || !("pmvd" %in% type)) ) stop ("invalid type requested")
if(!all(type %in% alltype) && length(alltype)==5 && "pmvd" %in% type ) stop ("pmvd is not a valid type in this version of relaimpo, obtain the non-US version, if you want to use pmvd")

#covg is the variance-covariance matrix of y with all regressor variables
#it needs to be square and at least 3x3 (= 2 regressors)

#type is a character listing (vector incl. single character string, list, matrix) 
#    that chooses one or more of the available relative importance metrics
#    The metrics can be listed in any order (and even wrong text can be in there),
#    but the selected metrics will always appear in the result in the order given below:
#    "lmg","pmvd","last","first","betasq","pratt"


#diff is a logical requesting differences between metrics if TRUE
#rank is a logical requesting ranks of metrics if TRUE (largest=1, smallest=p)
#rela is a logical requesting normalization to a percentage scale (sum 100%) if TRUE
    # if FALSE, all metrics are relative to var(Y), lmg, wlmg, pmvdalt and pratt sum to R^2, 
    # first and last do not sum to R^2 but neither to 100%, 
    # betasq (the squared standardized correlation coefficients) do not sum to R^2 either

#no of regressors
p <- ncol(covg)-1
#initialise output matrices
wahr<-matrix(0,1,p)
#vector for setdiff
alle<-1:(p+1)
names<-colnames(covg)
if (is.null(names)) names<-c("y",paste("X",1:p,sep=""))

# covg is the covariance matrix of y with all p regressors, dimension (p+1)x(p+1), y corresponds to first position
# start by calculating all conditional variances of y

# construct all subsets with function nchoosek from package vsn,
# calculate conditional variances as top left corner of appropriate conditional covariance matrix
# and write them to vectors in same position as with index matrices

#initialise lists of index matrices and variance rows

#first element of list unconditional, list initialized in full length
indices <-rep(list(0),p+1)
variances <- rep(list(covg[1,1]),p+1)
#conditioning on all variables, i.e. var=s^2
indices[[p+1]]<-matrix(1:p,p,1)
variances[[p+1]]<-covg[1:1]-covg[1,2:(p+1)]%*%solve(covg[2:(p+1),2:(p+1)],covg[2:(p+1),1])

hilf<-varicalc(type,alle,covg,p,indices,variances)
indices<-hilf$indices
variances<-hilf$variances

#output R-squared in order to show the total that is subdivided
ausgabe<-new("relimplm",var.y=as.numeric(variances[[1]]),
    R2=as.numeric(1-variances[[p+1]]/variances[[1]]))

if ("lmg"  %in% type) ausgabe<-lmgcalc(ausgabe,p,indices,variances,rank,diff,rela)

if ("pmvd" %in% type) ausgabe<-pmvdcalc(ausgabe,p,indices,variances,rank,diff,rela)

# variance improvement when entering each regressor last into the model
if ("last" %in% type) ausgabe<-lastcalc(ausgabe,p,variances,rank,diff,rela)

# variance improvement when entering each regressor first into the model
if ("first" %in% type) ausgabe<-firstcalc(ausgabe,p,variances,rank,diff,rela)

if ("betasq"  %in% type) ausgabe<-betasqcalc(ausgabe,covg,p,variances,rank,diff,rela)

if ("pratt"  %in% type) ausgabe<-prattcalc(ausgabe,covg,p,rank,diff,rela)

#ausgabe contains (in this order) var.y, R2, lmg, rank.lmg, diff.lmg, 
#                        pmvd, rank.pmvd, diff.pmvd, 
#                        last, rank.last, diff.last, first, rank.first, diff.first,
#                                 betasq, rank.betasq, diff.betasq, pratt, rank.pratt, diff.pratt
# as far as requested by the call
# default: R2, lmg, rank.lmg
slot(ausgabe,"namen")<-names
slot(ausgabe,"type")<-type
#ausgabe<-as(ausgabe,"relimplm")

return(ausgabe)
}

