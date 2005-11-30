"boot.relimp" <-
function(y,x,b=1500,type="lmg",rank=TRUE,diff=TRUE,rela=TRUE)
{
# Author and copyright holder: Ulrike Groemping

#This routine is distributed under GPL version 2 or newer.
#The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

#--------- start of function

#function for simulating percentage contribution and rank with CIs
#result (ergebnis) contains estimated percentages, estimated ranks, CIs for percentages, CIs for ranks
# y is criterion variable, 
# x is matrix of x-variables, 
# b is number of bootstrap replicates
# bty is bty of bootstrap intervals, default BCa, 
#       eventually intended that list can be given like in package boot
#       currently bty only takes one single bty
#      (rank intervals always with bty percentile, 
#       since BCa does not work properly with ranks and normal not reasonable)
# level is the confidence level, list can be given like in package boot
# diff (if TRUE) requests bootstrap intervals for differences
# rank (if TRUE) requests bootstrap intervals for ranks
# type is a list of relative importance types for which bootstrapping is requested

if(!(length(y)==nrow(x))) stop("number of rows in y and x MUST be identical")
if(!is.logical(rank)) stop("rank must be a logical")
if(!is.logical(diff)) stop("diff must be a logical")
if (!is.logical(rela)) stop("rela must be a logical")
#alltype is set in zzz.R
alltype <- alltype()
if(!all(type %in% alltype) && (length(alltype)==6 || !("pmvd" %in% type)) ) stop ("invalid type requested")
if(!all(type %in% alltype) && length(alltype)==5 && "pmvd" %in% type ) stop ("pmvd is not a valid type in this version of relaimpo, obtain the non-US version, if you want to use pmvd")



#combine y and x and calulate covariance
if (is.null(colnames(x))) colnames(x)<-paste("X",1:ncol(x),sep="")
daten<-as.matrix(cbind(y,x))
empcov <- cov(daten)
p<-ncol(x)

ausgabe<-new("relimplmboot")
ausgabe@type<-type
ausgabe@nboot<-b
ausgabe@rank<-rank
ausgabe@diff<-diff
ausgabe@rela<-rela

#provide names for identifying statistics
names<-colnames(x)
#names<-paste("X",1:ncol(x),sep="")
diffnam<-paste(names[nchoosek(p,2)[1,]],names[nchoosek(p,2)[2,]],sep="-")

   #run bootstrap
   #options for calc.relimp to be handed over after b
   booterg<-boot(daten,calcrelimp.forboot,b,type=type,diff=diff,rank=rank,rela=rela)
    #calcrelimp.forboot returns the results from analysis as a vector with the following elements
    #first var.y and R^2, then in the order lmg lmgrank lmgdiff pmvd pmvdrank pmvddiff 
    #                    last lastrank lastdiff first firstrank firstdiff
    #                    betasq betasqrank betasqdiff pratt prattrank prattdiff
    #                    namen
    #those statistics that have been requested (type relimplm)
ausgabe@boot<-booterg
return(ausgabe)
#----------end of function
}

