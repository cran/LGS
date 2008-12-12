`zbins`<- function(allRatio=NULL, bins=NULL){

# Test if allRatio is has class "LGS" (if yes, it should have allRatio as "LGS")
if(class(allRatio)=="LGS"){
    lgsar<- allRatio$ChrMap$LGS
	} else{
    if(is.null(allRatio)){
	    stop("Vector of allele ratios not available.")
        }
	lgsar<- allRatio
    }
		
# Use "Bins" in simlGS list if supplied
	
if(!is.null(bins)){
	bins<- bins
	}
if(is.null(bins) & class(allRatio)=="LGS"){
    bins<- allRatio$ChrMap$Bins
	} 

if(length(bins)!=length(lgsar)){
    stop("Length of 'allRatio' and length of 'bins' differ.")
	}

## define z-score function ##
zscore<- function(x) (x - mean(x) )/sd(x)

#z1<- aggregate(lgsar, list(bins), mean)$x ## Mean allele ratio of the bins 
z1<- aggregate(lgsar, list(Bin= bins), mean)
z1<- z1[match(unique(bins), z1$Bin), ]$x   ## Order dataframe according to original order of bins (i.e. not necessarily in alphabeltical order)

z2<- zscore(z1) ## Zscore of bins
return(z2)   ### Vector of bin z-scores
}