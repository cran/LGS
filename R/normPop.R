`normPop`<- function(hapList=NULL, proportions= NULL, popSize=1, onePop= TRUE){

if(!is.list(hapList)){
    stop("hapCount should be a list of vectors or matrices: x= list(...).")
    }
if(is.null(proportions)){
    proportions<- rep(1/length(hapList), length(hapList))   ## If a vector of proportions is not given, default to equal weight each population
    } else {
        proportions<- proportions/sum(proportions)  ## Rescale the given proportions to sum up to 1
        }
if(length(hapList)!=length(proportions)){
    stop("Number of populations (hapList) and vector of proportions are not of the same length.")
    }  
hapCount<- lapply(hapList, 
            function(x) if(is.matrix(x)){   ## if there are matrices in hapList, extract the last column which should contain the counts
                            x[, ncol(x)]
                            }else{
                                x
                                }
            )  
            
sx<- lapply(hapCount, sum) ## sum of individuals in each population

normx<- vector(mode="list", length=length(hapCount)) ## List where normalized counts will be stored

for(i in seq(length(hapCount)) ){
    normx[[i]]<- unlist(hapCount[i])/unlist(sx[i])  * (popSize) * proportions[i]   ## Normalize haplotypes
    } 

normHaps<- hapList

for(i in seq(length(hapList))){
    normHaps[[i]]<- if(is.matrix(hapList[[i]])){   ## if there are matrices in x, extract the last column which should contain the counts
                        hm1<- normHaps[[i]]       ## Temporary matrix
                        hm1[, ncol(hm1)]<- unlist(normx[[i]]) # Replace original count with the normalized one
                        normHaps[[i]]<- hm1       ## Replace original matrix with normalized one
                        }else{
                            unlist(normx[[i]])
                            }
    }  

if(onePop==TRUE){   # If only one population is to be returned...
    if(all(sapply(normHaps, is.vector))){   # If only counts (vectors) are supplied), return a single vetctor
        return(unlist(normHaps))
        }
    if(all(sapply(normHaps, is.matrix))){    # If only populations (matrices of haplotypes and counts) return a single population
        mh1<- matrix(data=NA, nrow=0, ncol=ncol(normHaps[[1]]))  # Empty matrix to store the merged list
        for(i in seq(length(normHaps))){
            mh1<- rbind(mh1, normHaps[[i]])    # Append each matrix in normHap to a single matrix
            colnames(mh1)<- colnames(hapList[[1]]) # Give column names from first matrix
            }
        return(mh1)
        }else{
            return(normHaps)    # If a mixture of counts and pops are given, return a normalized list
            }
    }else{         # If onePop==FALSE (no mixing is required), return a list
        return(normHaps)
        }
}
