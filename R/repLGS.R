`repLGS`<- function(offspring= NULL, scenario= NULL, off.by.scen= TRUE, nsim= 1,
map=NULL, fix.before.sel=FALSE, ...){
	
if(off.by.scen==TRUE){  ## Generate all the combinations of offspring by scenario if required
    noff<- rep(rep(offspring, length(scenario)), each=nsim)
	nscen<- rep(rep(scenario, each=length(offspring)), each=nsim)
	}
if(off.by.scen==FALSE){
    noff<- offspring
	nscen<- scenario
	if(length(noff)!=length(scenario)) {
	    stop("Vector in argument offspring and scenario are not of the same length while off.by.scen==FALSE")
		}
	}
mlgs<- matrix(data=NA, nrow=length(noff), ncol=nrow(map))  # matrix to store allele ratios (each row is a simulation)
npop<- vector(length=length(noff))

## Start simulations ##
for(i in 1:length(noff)){
    
	tmp<- simLGS(offspring=noff[i], scenario=nscen[[i]], map=map, fix.before.sel=fix.before.sel, ...)
	mlgs[i,]<- tmp$ChrMap$LGS
	if(fix.before.sel==TRUE)   ## If the number of progeny is held fixed, store the number of individuals surviving
	   npop[i]<- nrow(tmp$selPop)
	if(fix.before.sel==FALSE)  ## If the number of individuals surviving is held fixed, store the total number of progeny
	   npop[i]<- nrow(tmp$recPop)
}

# Prepare output #
if(fix.before.sel==TRUE){
     nSel<- npop
	 nRec<- noff
	 }
if(fix.before.sel==FALSE){
     nSel<- noff
	 nRec<- npop
	 }
if(is.list(nscen)){      ## If scenario is list (e.g. mixture of column names and column indexes)
     nscen<- unlist(nscen)
	 }


xv<- aggregate(tmp$ChrMap$Chr, list(tmp$ChrMap$Chr), length)   # Number of markers on each chr
m1<- unlist(sapply(xv$x, function(x) seq(1,x)))   # Marker position on chr
	  
colnames(mlgs)<- paste("c", tmp$ChrMap$Chr, "_", m1, sep="")   # Column names for LGS matrix: c1_1, c1_2... 

return( data.frame(scenario= nscen, nRec, nSel, data.frame(mlgs)))
}
