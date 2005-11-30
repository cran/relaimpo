"calcrelimp.forboot" <-
function(data,indices,...){
# Author and copyright holder: Ulrike Groemping

#This routine is distributed under GPL version 2 or newer.
#The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    data <- data[indices,]
    cova <- cov(data)
    ausgabe <- list2vec(as(calc.relimp(cova,...),"list"))
    return(ausgabe)
}

