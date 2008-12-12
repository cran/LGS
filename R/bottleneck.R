`bottleneck`<- function(population=NULL, neck=NULL, keep.zero.count=FALSE){

if(is.matrix(population)){
    ## If a matrix is supplied take the last column with hap counts
    nhaplos<- population[,ncol(population)] 
    }
if(is.vector(population)){
    ## If a vector is supplied take it as a haplotype count
    nhaplos<- population
    }

sh<- sum(nhaplos)

if(sh<neck){
    stop("Less individuals in the current population than individuals after bottleneck")
    }

cullingF<- neck/sh ## Culling factor
expIndiv<- cullingF * nhaplos
newIndiv<- sapply(expIndiv, function(x) rpois(1,x))
    
if(is.matrix(population)){
    ## If a matrix was supplied, replace the previuos haplotype count
    ## with the one after bottleneck
    botMat<- population
    botMat[,ncol(botMat)]<- newIndiv
    if(keep.zero.count==FALSE){
        botMat<- botMat[which(botMat[,ncol(botMat)]>0),]
        }
   }

if(is.vector(population)){
    botMat<- newIndiv
    if(keep.zero.count==FALSE){
        botMat<- botMat[botMat>0]
        }
    }
return(botMat)
}