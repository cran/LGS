`topClones`<- function(
     x= NULL,  ## A vector of counts of each haplotype (e.g. output of selClones$nhaplos)
     np= c(0.68, 0.95, 0.99)  ## Percentage of population
     ){
xv<- cumsum(sort(x, decreasing=TRUE))  ## Cumulative sum ordered from most represented haplotype
xp<- sum(x)*np   ## Value in cumsum that would correspond to the desired percentage of pop.
nc<- vector(length=0)

for(i in 1:length(np)){
     nc[i]<- which(abs(xv-xp[i]) == min(abs(xv-xp[i])))
    }
cat("", np, "\n", nc, "\n", sep="\t")
return(nc=nc) # Number of clones necessary to make up the proportions in np
}

### Find the number of clones that account for n percent of the whole population ###
### There must be an R function that does this! ###
