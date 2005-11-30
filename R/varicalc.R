"varicalc" <-
function(type,alle,covg,p,indices,variances)
{
# Author and copyright holder: Ulrike Groemping

#This routine is distributed under GPL version 2 or newer.
#The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

#routine for calculating all needed residual variances
#together with indices denoting the variables included in the model

liste<-c(1,p-1)
if (any(c("lmg","pmvd")  %in% type)) liste<-1:(p-1)
    #condition also asks for old/potentially new types, does not do any harm, if not needed

for (k in liste)
{
    # k is no. of regressors that are conditioned upon
    # jetzt becomes matrix with k rows and choose(p,k) columns
    jetzt<-nchoosek(p,k)
    #insert in list of indices
    indices[[k+1]]<-jetzt
    # calculation of conditional variances
    # initialise row vector of variances
    varjetzt<-matrix(0,1,choose(p,k))

for (j in 1:(choose(p,k)))
{
    #column of index-matrix, one added for picking appropriate elements from covg
    diese<-jetzt[,j]+1
    #indices for other elements of covg that are not conditioned upon
    andere<-setdiff(alle,diese)

    varjetzt[j]<-(covg[andere,andere]-covg[andere,diese]%*%solve(covg[diese,diese],matrix(covg[diese,andere],k,p+1-k)))[1,1]
}
    #insert in list of variances
    variances[[k+1]]<-varjetzt
}
return(list(indices=indices,variances=variances))
}

