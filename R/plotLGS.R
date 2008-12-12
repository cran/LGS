`plotLGS`<- function(
	 simLGS=NULL, # output of simLGS function or list with similar structure
     cumDist= simLGS$ChrMap$CumDist,  # cumDist, chr, allRatio can be specified to plot data that doesn't come from simLGS()
     chr= simLGS$ChrMap$Chr,
	 allRatio= simLGS$ChrMap$LGS,
	 width.cm= 23/2.54,  ## width  and heigth of the graphic window
     height.cm= 11/2.54,
     barp=FALSE, # should barplots of clone frequency be showed?     
	 add= FALSE, # should data be plotted on current graph?
	 show.range= TRUE, # should horizontal lines for mean and sd be showed? (better to set this to FALSE if more than one scan is plotted)
	 new.window= TRUE,
	 ...   # further arguments to be passed to plot or points
)
{
if(is.data.frame(simLGS)){     # if simLGS is a dataframe, use it assuming it has LGS, CumDist, Chr columns
		if(is.null(simLGS$LGS)) {stop("Column 'LGS' not found")}
		if(is.null(simLGS$CumDist)) {stop("Column 'CumDist' not found")}
		if(is.null(simLGS$Chr)) {stop("Column 'Chr' not found")}
		# if(is.null(simLGS$SA) | is.null(simLGS$SelCoef)) {warning("Column SA or SelCoef not found (necessary to highlight selected loci)")}
	simLGS$ChrMap<- simLGS
	}

        
if(is.null(simLGS)){     # if no simLGS list is supplied, use cumDist, allRatio, chr as specified in the input
    cumDist<- as.numeric(as.character(cumDist))   # remove anything that is not numeric in cumulative distances and allele ratios
    allRatio<- as.numeric(as.character(allRatio))
	tt<- as.data.frame(na.omit(cbind(cumDist, chr, allRatio)))   # remove rows with NAs
	simLGS$ChrMap$CumDist<- tt$cumDist
	simLGS$ChrMap$Chr<- tt$chr
	simLGS$ChrMap$LGS<- tt$allRatio
	}
	
### Set graphical window ###
if(add==FALSE & new.window==TRUE){ # set window only if a new plot is desired
	x11(width=width.cm, height=height.cm)
	} 
### Collect populations from LGS object before and after selection ###
if(barp==TRUE & add==FALSE){    ## Draw barplots only if asked an if add=FALSE

allPop<- rbind(simLGS$ParPop, simLGS$recPop, simLGS$selPop) #

index<- c(rep("Par.", dim(simLGS$ParPop)[1]), rep("Prog", dim(simLGS$recPop)[1]), rep("Sel", dim(simLGS$selPop)[1]))

mar1<- apply(allPop,1, function(x) length(x[x==0])/length(x)) #For each row (gamete) in AllPop, find if recombinant or parental
mar2<- sapply(mar1, function(x) ifelse(x==1,"P0",ifelse(x==0,"P1","Rec")))

barData<- cbind(Pop=rep(unique(index),each=3), GamType=rep(c("P0", "P1", "Rec"), 3))
ratio<- aggregate(mar2, list(GamType=mar2, Pop=index), length); colnames(ratio)[3]<- "countGam"
sumIn<- aggregate(index, list(Pop=index), length); colnames(sumIn)[2]<- "countSum"

merged<- merge(ratio, sumIn, by= "Pop")
merged2<- merge(merged, barData, by=c("Pop", "GamType"), all=T); merged2<- cbind(merged2, gamRatio=merged2$countGam/merged2$countSum)
merged2[is.na(merged2)]<- 0

#### Set graphic; layoutlayout.show(3) ####

layout(matrix(c(1,2,3,3,3), nrow=1, ncol=5, byrow=T))
par(xpd=T)

##### Barplot of ratios #####

barplot(matrix(merged2$gamRatio, nrow=3), 
   names.arg=unique(index), main="Gamete proportion",
   xlab="Population", ylim=c(0,1.05),
   col = c("lightblue", "mistyrose", "cyan"))
legend(x=-0.8, bty="n", y=1.1, horiz=T, legend=c("Susc", "Res", "Rec"), fill=c("lightblue", "mistyrose", "cyan"))
par(xpd=F)

##### Barplot of absolute number of individuals #####

barplot(matrix(merged2$countGam, nrow=3), 
   names.arg=unique(index), main="Gamete count",
   xlab="Population",
   col = c("lightblue", "mistyrose", "cyan"))

} ## Close if(barp==TRUE) 

#### LGS profile ####
par(mar=c(5,4,6,2)+0.1)

if(is.null(simLGS$selPop)==FALSE){      ## Number of recombinants after selection
    ncl<- nrow(simLGS$selPop)
    }

if(add==FALSE){   # make a new plot if add=FALSE
    plot(y= simLGS$ChrMap$LGS, x= simLGS$ChrMap$CumDist, ylim=c(0, 1), xlab="Cumulative distance", ylab="Allele ratio", 
	    type="o", cex=0.75, lwd= 1.75, main=paste("LGS profile", ifelse(exists("ncl")==TRUE, paste("\n(", ncl, " clones)", sep=""), ""), sep=""), ...)   ## xlim=c(min(simLGS$ChrMap$CumDist), max(simLGS$ChrMap$CumDist))
    }
if(add==TRUE){    # if add=T use points() instead of plot() above
	points(y= simLGS$ChrMap$LGS, x= simLGS$ChrMap$CumDist, type="o", cex=0.75, lwd= 1.75, ...)
	}

chrpos<- as.vector(tapply(simLGS$ChrMap$CumDist, simLGS$ChrMap$Chr, max)) # End and start of each chromosome
chrpos0<- as.vector(tapply(simLGS$ChrMap$CumDist, simLGS$ChrMap$Chr, min))



chrpos1<- c(0,chrpos) # Ending position of each chromosome (leading 0 for compatibility)
chrsize<- diff(chrpos1) # Size of each chr

if(length(unique(simLGS$ChrMap$Chr))>1){    # Draw lines shading chromosomes only if there is more than 1 chr.

    m<- seq(1, length(simLGS$ChrMap$LGS))  # marker index
    m1<- tapply(m, simLGS$ChrMap$Chr, max) # marker index at the end and start of each chr
    m0<- tapply(m, simLGS$ChrMap$Chr, min)
    ar1<- simLGS$ChrMap$LGS[m1[-length(m1)]]  # allele ratio at the end (excluding last) and start (excluding first) of each chr
    ar0<- simLGS$ChrMap$LGS[m0[-1]]  # 
    
    segments(x0= chrpos[-length(chrpos)], ## end of each chr excluding last
             y0= ar1, #allele ratio in correspondence of the last marker of each chomosome
             x1= chrpos0[-1], ## start of each chr
             y1= ar0, # allele ratio at the first marker of each chr)
             col= "white", lty="solid")
}

if(add==FALSE){ 	  ### Draw chromosome baundaries, annotation etc only if a new plot is called
    if(show.range==TRUE){  # plot lines of mean and sd if required
        abline(h=c(mean(simLGS$ChrMap$LGS),  
          (mean(simLGS$ChrMap$LGS)-2*sd(simLGS$ChrMap$LGS)), (mean(simLGS$ChrMap$LGS)+2*sd(simLGS$ChrMap$LGS)), 
          (mean(simLGS$ChrMap$LGS)-3*sd(simLGS$ChrMap$LGS)), (mean(simLGS$ChrMap$LGS)+3*sd(simLGS$ChrMap$LGS))),
            lty=c("solid", rep("dashed",4)), col=c(rep("grey",3), rep("red",2)))
        }


    abline(v=c(0, chrpos), col="white", lty="solid")
    abline(v=c(0, chrpos), col="grey", lty="dotted")  # Chromosome boundaries

    ## Highlight the selected loci ##
    points(x=simLGS$ChrMap$CumDist,    ## Extract from map dataframe (simLGS$ChrMap) the loci that are under selection (SelCoef>0)
           y=simLGS$ChrMap$LGS, ## Get the allele ratio corresponding to the loci in x
           type="p", pch=19, cex=1, 
           col=ifelse(simLGS$ChrMap$SelCoef<0, "red", 
		       ifelse(simLGS$ChrMap$SelCoef>0, "blue", "transparent")) ## Assign col blue or red depending on direction of selection
		   )   
		   #col=ifelse(simLGS$ChrMap$SA[simLGS$ChrMap$SelCoef>0]==1, "red", "blue"))   ## Assign col blue or red
		   
    mtext(side=3, at=c(min(simLGS$ChrMap$CumDist), chrpos-chrsize/2),   ## chromosome names on top
          text=c("Chr", unique(simLGS$ChrMap$Chr)), line=0.5, cex=0.7)

    ##legend(x=1, y=1, cex=1, bg= "white", legend=c("Allele ratio", "Mean ratio", "Mean ratio +/- SD", "Selected locus"), lty=c("solid", "solid", "dashed", NA), col=c("black","grey","grey","blue"), pch=c(21,NA,NA,19))
    
    par(mar=c(5,4,4,2)+0.1)

    rug(simLGS$ChrMap$CumDist)  # Plot on the x axis where the markers are
}
}
