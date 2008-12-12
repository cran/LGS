`haplo2string`<- function(x= NULL){

xstr<- x[,-ncol(x)]
xstr<- apply(xstr, 1, paste, collapse="")
xnclones<- x[,ncol(x)]
    return(data.frame(haplotype=xstr, nclones=xnclones))
}
