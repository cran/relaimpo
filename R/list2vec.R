"list2vec" <-
function(liste)
{
#function that transforms numerical list to vector
#taken from bioconductor package vsn (author Wolfgang Huber, license LGPL)
    ausgabe<-liste[[1]]
    if (length(liste)>1)
    {
    for (a in 2:length(liste))
    {
    if (is.numeric(liste[[a]])) ausgabe<-c(ausgabe,liste[[a]])
    }
    }
    return(ausgabe)
}

