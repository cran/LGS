`selClones`<- function(haplos= NULL, nhaplos= haplos[,ncol(haplos)], 
                     selCoeff= rep(0, (ncol(haplos)-1)), nindiv= NULL)
{     

if(nrow(haplos)!=length(nhaplos)){
    stop("Number of rows (haplotypes) in 'haplos' and length of 'nhaplos' (haplotype counts) differ!")
    }
if((ncol(haplos)-1)!=length(selCoeff)){
    stop("Number of loci in 'haplos' (ncol-1) and length of 'selCoeff' differ!")
    }
if(nindiv <= sum(nhaplos)){
    warning("Final number of individuals is less than or equal to the starting one (nindiv <= sum(nhaplos))!")
    }
if(min(nhaplos)<0){
    stop("Negative count found in nhaplos.")
    }
if(max(selCoeff)>1 | min(selCoeff)< -1){
    stop("Found selection coefficient below 0 or above 1.")
    }

SelCoef<- abs(selCoeff[selCoeff!=0])
selLoc<- which(selCoeff!=0) # position of the markers under selection

SA<- na.omit(
	 ifelse(selCoeff>0, 0,
     ifelse(selCoeff<0, 1, NA))
	 )

selMat<- matrix(nrow=nrow(haplos), ncol=length(selLoc))  # matrix to store selection coefficients

for(i in seq(length=length(selLoc))  ){
    selMat[,i]<- ifelse(haplos[,selLoc[i]]==SA[i], 1, 1-SelCoef[i])
    }
        

survProb<- apply(selMat, 1, prod) ## Vector of survival probabilities (dead or alive after selection)

    
while(sum(nhaplos) < nindiv){      ## check if the current number of individual is less than the desired
    expClones<- nhaplos*survProb 
    ## New clones are produced by a rpois with expectation = existing number 
    newclones<- sapply(expClones, function(x) rpois(1,x)) 
    prevgen<- nhaplos              ## Store the previous generation to check if it is closer to nindiv then nhaplos
    nhaplos<- nhaplos+newclones    ## Increase the number of each haplotype each round of divisions (while loop)
    }

## Use the previous generation if near to nindiv
if(abs(sum(prevgen)-nindiv)<abs(sum(nhaplos)-nindiv)){
    nhaplos<-prevgen
    }
 
asexP<- nhaplos*haplos   ### Matrix where the 1 alleles in each haplotype are multiplied by the corresponding number of haplotypes 

LGS<- 1-(apply(asexP[,-ncol(asexP)], 2, sum)/sum(nhaplos)) ## Calculate allele ratios

return(list(LGS= LGS,          ### Vector of allele ratios
            survPr= survProb,  ### Survival probability of each gamete
            nhaplos= nhaplos     ### Number of times each haplotype is represented
			)) 
}
