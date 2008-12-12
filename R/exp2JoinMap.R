`exp2JoinMap`<- function(popMatrix= NULL, markerNames= NULL, recombinantOnly= TRUE, outfile= "JMPtmp.loc")
{

if(is.null(markerNames)){
    if("nclones" %in% colnames(popMatrix)){
        markerNames<- paste("m", 1:(ncol(popMatrix)-1), sep="")
	} else {
	    markerNames<- paste("m", 1:ncol(popMatrix), sep="")
		}
	}
	
popMatrix<- as.data.frame(popMatrix)  # Convert input matrix to dataframe

if("nclones" %in% colnames(popMatrix)){ # Test if it has the "nclones" column which represents the numebr of times each haplotype is represented
    popMatrix<- sapply(popMatrix[ ,-ncol(popMatrix)], rep, popMatrix$nclones) # Multiply each haplotype by nclones
    }

if(recombinantOnly==TRUE){
	selProg<- apply(popMatrix, 1, sum) # Select recombinant gametes (no parentals) WARNING: This function discharges gametes that have only 0 or only 1 alleles. If a genome is too small, some truly recombinant gametes will not have crossing overs and will be erroneously rejected
    popMatrix<- t(popMatrix[selProg>0 & selProg<dim(popMatrix)[2], ]) # Subset and transpose progeny
	} else {
	    popMatrix<- t(popMatrix)
        }
popMatrix[popMatrix==0]<- "a"  # convert to accomodate JoinMap allele coding
popMatrix[popMatrix==1]<- "b"  # 

mapPop2<- matrix(nrow=dim(popMatrix)[1], ncol=1)  # Empty matrix to store formatted genotypes

for(i in 1:dim(popMatrix)[1]){                    # Alleles collapsed to a single string
   xv<- paste(popMatrix[i,], collapse="")
   mapPop2[i,]<- xv
   }

mapPop2<- data.frame(markerNames, mapPop2)  # Marker names

write.table(x= c(  
                "; Output of LGS simulation", 
                "",
                "name = LGS",
                "popt = HAP1",
                paste("nloc =", dim(popMatrix)[1]),
                paste("nind =", dim(popMatrix)[2]),
                ""
               ),
           file=outfile, append=FALSE, quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)	

write.table(x= mapPop2, file= outfile, append=TRUE, quote=FALSE, row.names=FALSE, sep="\n", col.names=FALSE)  # write output
cat("Population matrix successfully exported!\n")
}
