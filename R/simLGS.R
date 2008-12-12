`simLGS`<- function(offspring= 50, fix.before.sel= FALSE, NoP1= 10, NoP2= 10, map= NULL, 
           scenario= NULL, gen.dist= "k", kb2cM= 15)
{
if(is.null(map)){
    stop("Argument to map not supplied.")
    }

if(!is.null(scenario)){
	nam<- names(map)         ### Assign the column containing the selection coefficients to map$SelCoef

	if(is.character(scenario) & is.na(match(scenario,nam))){   ## Return error if argment passed to scenario is not in the map file
       stop(paste("\"",scenario,"\"", " is not the name of a column in the map file", sep=""))
    }
    if(is.numeric(scenario) & scenario>length(nam)){   ## Return error if argment passed to scenario is not in the map file
       stop(paste("Cannot find column index ", scenario, ": the map file has only ", length(nam), " columns!", sep=""))
    }
        
    sc<- ifelse(is.numeric(scenario), scenario, which(nam==scenario)) ## If scenario is given as column index, 
	map$SelCoef<- map[, sc]                                     ## assign that column to map$SelCoef. If scenario is given as column name, find the column index corresponding to it
	
   }

selLoc<- which(map$SelCoef!=0) # position of the markers under selection
SA<- na.omit(
	 ifelse(map$SelCoef>0, 0,
     ifelse(map$SelCoef<0, 1, NA))
	 )

     
if(max(map$SelCoef)>1 | min(map$SelCoef)< -1){    ## Check if selection coefficients are 0<x<1
       cat("Selection coefficients out of bounds:\n", map$SelCoef[map$SelCoef<0 | map$SelCoef>1], "\n")
       stop("Some selection coefficients are below 0 or above 1")
       }

SelCoef<- abs(map$SelCoef[map$SelCoef!=0])
   

intMarkSpa<- c(diff(map$uCumDist), 0)

haldane2recFreq<- function(x) (0.5*(1 - exp(1)^(-2*  (x/100)  ) ))*100        # function to convert Haldane distances to recombination frequency
kosambi2recFreq<- function(x) ((0.5 * ((exp(1)^(4*(x/100)))-1) )/ ((exp(1)^(4*(x/100))+1)) )* 100# function to convert Kosambi cM to recombination frequency

options(warn=-1) # turn off warning message about NA introduced

intMar<- as.vector(as.numeric(as.character(intMarkSpa)))

options(warn=0) # reset warnings to default

if(gen.dist=="h"){
   intMar<- round(haldane2recFreq(intMar), 2) # convert user's distances (Haldane cM) to recombination frequency
   intMar[aggregate(seq(1,length(intMarkSpa)), list(map$Chr), max)$x]<- 50 # replace the last interval of each chr with 50 rec. freq.
   }
if(gen.dist=="k"){
   intMar<- round(kosambi2recFreq(intMar), 2) # convert user's distances (Kosambi cM) to recombination frequency
   intMar[aggregate(seq(1,length(intMarkSpa)), list(map$Chr), max)$x]<- 50 # replace the last interval of each chr with 50 rec. freq.
   }
if(gen.dist=="phy"){
   intMar<- intMar/kb2cM
   intMar[aggregate(seq(1,length(intMarkSpa)), list(map$Chr), max)$x]<- 50 # replace the last interval of each chr with 50 rec. freq.
   }
if(gen.dist=="rf"){
   intMar[aggregate(seq(1,length(intMarkSpa)), list(map$Chr), max)$x]<- 50 # replace the last interval of each chr with 50 rec. freq.
   }

intMar<- as.numeric(as.character(intMar))

#cat("Preparing map and parents... ")

ChrMap<- as.data.frame(map)     # Chromosomes and markers
nchr<- length(unique(ChrMap$Chr))  # count number of chromosomes from map dataframe
nloc<- length(ChrMap$Locus)        # count number of markers

if(is.matrix(NoP1)){          # Parental populations 
   ParMat1<- NoP1[,-ncol(NoP1)] # Strip column with haplotype count
   nhaplo1<- NoP1[,ncol(NoP1)]  # Store haplotype count
   } else {
      ParMat1<- matrix(ncol=nloc, nrow=NoP1, data=0)  # If parent 1 is given as number of haplotypes, make a matrix 
      nhaplo1<- rep(1, nrow(ParMat1))  # Each haplotype counted only 1
	  }
   
if(is.matrix(NoP2)){
   ParMat2<- NoP2[,-ncol(NoP2)]
   nhaplo2<- NoP2[,ncol(NoP2)]
   } else {
      ParMat2<- matrix(ncol=nloc, nrow=NoP2, data=1)
      nhaplo2<- rep(1, nrow(ParMat2))
	  }
   
#cat("Ok\n")

### GENERATING RECOMBINANT POPULATION ###
#cat("Generating recombinant population... \n")

selPop<- matrix(nrow=0, ncol=nloc) #Initialize selected pop., recombinant pop before selection, and rec. pop. before shuffling chr. 
recPop<- matrix(nrow=0, ncol=nloc)
oriRecPop<- matrix(nrow=0, ncol=nloc)
fixpop<- 0  ## Counter to keep track of the number of clones generated.

ParPop<- rbind(ParMat1,ParMat2) # Parental population
hapwt<- c(nhaplo1, nhaplo2)   ## Weight of each parental haplotype (default to equal frequency if not supplied)
rownames(ParPop)<- c(rep("P1", nrow(ParMat1)), rep("P2", nrow(ParMat2)))

while(fixpop<offspring){ # Start of the while loop to generate progeny. Evaluate if the number of clones required has been reached. If not, keep crossing and selecting 

    mate1<- as.vector(t(ParPop[sample(dim(ParPop)[1], size=offspring, replace=T, prob=hapwt),]))  # "First" parents sampled at random from the parental pop. Transformed to vector for ease
    mate2<- as.vector(t(ParPop[sample(dim(ParPop)[1], size=offspring, replace=T, prob=hapwt),]))  # "Second" 
                                                                                 # Each mating gives 1 haploid individual (This is obviously not true in reality)
                                                                                 
    vecCO<- rep(intMar, offspring)    # vector of recombination frequencies extracted from ChrMap. Repeated offspring times


    vecCO1<- (sample(seq(0:10000), size= length(vecCO), replace=TRUE))/100  # Generate crossing overs.
    vecCO2<- ifelse(vecCO-vecCO1<0, "0", "1")                       # If the element of vecCO1 is below the corresponding genetic distance, a crossing over will occur in that marker interval

    noCO<- sum(vecCO)    # Total number of crossing overs (probably unnecessary with some rehandling of the code below)
    posCO<- unique(c(0, which(vecCO2==1), nloc*offspring))  # Index position of the marker intervals with cross over

    if(noCO==0)               {x0<-0; x1<-0} else
      if(length(posCO)%%2==0) x0<- posCO[seq(from=2, to=(length(posCO)-1), by=2)]+1 else
                              x0<- posCO[seq(from=2, to=length(posCO), by=2)] +1
    
    if(noCO==0)               {x0<-0; x1<-0} else
      if(length(posCO)%%2==0) x1<- posCO[seq(from=3, to=(length(posCO)-1), by=2)] else
                              x1<- posCO[seq(from=3, to=length(posCO), by=2)]
    
    tmpRecP<- mate1 #vector(length=offspring*nloc)
    if(noCO==0) tmpRecP<- mate1 else
      for(i in seq(length.out=length(x0))){
      
          tmpRecP[x0[i]:x1[i]]<- mate2[x0[i]:x1[i]]

         }

    tmpRecP<- matrix(tmpRecP, ncol=nloc, byrow=T)
    tmpRP<- tmpRecP #Recombinant population before shuffling chr (needed for summary stats)

recGam<- ifelse(apply(tmpRecP, 1, mean)==0 | apply(tmpRecP, 1, mean)==1, FALSE, TRUE) ## Vector of recombinant gametes (not sure I need this)      


   ### SELECTION ###

    selMat<- matrix(nrow=dim(tmpRecP)[1], ncol=length(selLoc))  # matrix to host selection coefficients

    for(i in seq(length=length(selLoc))  ){
        selMat[,i]<- ifelse(tmpRecP[,selLoc[i]]==SA[i], 1, 1-SelCoef[i])
        }
        
    survProb<- apply(selMat, 1, prod) ## Vector of surviving gametes (dead or alive after selection)
    fate<- sapply(survProb, function(x) sample(c("Alive", "Dead"), size=1, prob=c(x, 1-x), replace=TRUE)) ## Decide whether gamete i will survive
    tmpSelPop<- tmpRecP[which(fate=="Alive"),] # From the recombinant population keep only "alive" gametes

 selPop<- na.omit(rbind(tmpSelPop, selPop)) # na.omit removes NA from the initialize matrix (not sure it's necessary)
 recPop<- na.omit(rbind(tmpRecP, recPop))
 oriRecPop<- na.omit(rbind(tmpRP, oriRecPop)) 

fixpop<- ifelse(fix.before.sel==TRUE, dim(recPop)[1], dim(selPop)[1])   ## Should the counter increase (and stop the while loop) with the number of clones generated before or after selection?

#cat(paste("Clones:", fixpop, "\n"))
} # close while loop

#cat("Ok\n") 
#cat("Recombinant population done. Calculating LGS statistics... ")

if(fix.before.sel==FALSE){     ## These lines set the number of progeny clones exactly to the desired one
    trimPop<- offspring/nrow(selPop) # %of redundant gametes  
    recPop<- matrix(data= recPop[1:(nrow(recPop)*trimPop),], ncol=nloc) # Trim recPop to account for the eccess of gametes.
    selPop<- matrix(data= selPop[1:offspring,], ncol=nloc)
    }

 
if(nrow(selPop)>1){
    LGS<- apply(selPop, 2, function(x) length(x[x==0]))/dim(selPop)[1]  
    } 

if(nrow(selPop)==1){   ## If only one progeny survives, use its genotype as allele ratio (only 0 or 1)
    LGS<- as.vector(selPop)
    LGS[LGS==0]<- 2
    LGS[LGS==1]<- 0
    LGS[LGS==2]<- 1
    } 

if(nrow(selPop)==0){   ## If no individual survive (selPop has zero rows), set LGS to NA
    LGS<- rep(NA, nloc)    
    }

#cat("\nDone! Outputting LGS object\n")

recPop<- cbind(recPop, nclones=1) # Add a column of 1 to the populations

# Add a column of 1 to selPop... If selPop has rows!
if(nrow(selPop)>0){
   selPop<- cbind(selPop, nclones= 1)  # Each offspring is unique even if they might have identical genotype
   } else {
   selPop<- matrix(nrow=0, ncol=(nloc+1))  ## If selPop has no rows, add 1 emtpy column for compatibility with other matrices
     } 
rownames(selPop)<- NULL

ParPop<- cbind(ParPop, nclones=hapwt)


### Cumulative marker distances ###
cumDist<- as.vector(intMarkSpa)
cumDist[aggregate(seq(1,length(intMarkSpa)), list(map$Chr), max)$x]<- 0 # replace the last interval of each chr with 0cM to calculate cumulative distance
cumDist<- c(0, cumDist[1:(length(cumDist)-1)]) # Set first marker at 0cM and remove last element of the vector (i.e. the last end-of-chr)
cumDist<- cumsum(cumDist) # Calculate cumdist

ChrMap$Chr<- as.factor(ChrMap$Chr)    ### Convert Chr and Locus to factor (better for linear models)
ChrMap$Locus<- as.factor(ChrMap$Locus)
ChrMap<- cbind(ChrMap, CumDist= cumDist, LGS= LGS)

if(min(diff(map$uCumDist))<0){
   ng<- which(intMarkSpa<0)+1
   cat("\n")
   warning(paste("Negative intermarker distance(s) found at line(s)", ng, "of the map file", "\n", "Is the map file correct and sorted by the cumulative distance?"))
   }

lgsObj<- list(recPop= recPop, ParPop= ParPop, selPop= selPop, ChrMap= ChrMap)
class(lgsObj)<- "LGS"
return (lgsObj)
}
