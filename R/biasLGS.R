`biasLGS` <-
function(simLGS=simLGS,  ## lgs list to bias or vector of allele ratios
                   bias= 0.02,     ## Degree of bias, default 5%
                   meanBias= 0     ## bias towards one parent or the other? Default to 0 (random)
          )
{
if(is.list(simLGS)){
    xvlgs<- simLGS$ChrMap$LGS
    } else {
      xvlgs<- simLGS
      }

## Generate a random bias for each marker
if(length(xvlgs) > 1000000) {     ## This condition prevents my laptop to exceed memory when the input dataset is too big
    biasLGS<- vector(length=0)
	for(i in 1:ceiling(length(xvlgs)/1000000 )){      # Concatenate vectors of random numbers. The numebr of vectors to concatenate is a multiple of the input matrix (transformed into vector)/1000000
	     biasLGS<- c(biasLGS, rnorm(n= length(xvlgs)*(1000000/length(xvlgs)), mean= meanBias, sd= bias) )
         biasLGS<- biasLGS[1:length(xvlgs)]  ### Trim vector to fit input
	     }
} else {
  nr<- length(xvlgs)
  biasLGS<- rnorm(n= nr, mean= meanBias, sd= bias)    
  }

xvlgs<- xvlgs+((1+xvlgs)*biasLGS)                              ## Bias the original LGS
xvlgs[xvlgs<0]<- 0 # Reset markers to 0 or 1 where allele ratio is <0 or >1 
xvlgs[xvlgs>1]<- 1

if(is.list(simLGS)){     ### If an simLGS() list is supplied, make a new one equal to the input but with ChrMap$LGS replaced by the biased LGS. Original LGS is in ChrMap$unbiased
   simLGS$ChrMap$unbiasedLGS<- simLGS$ChrMap$LGS
   simLGS$ChrMap$LGS<- xvlgs
   list(recPop= simLGS$recPop, ParPop= simLGS$ParPop, selPop= simLGS$selPop, ChrMap= simLGS$ChrMap, depth= simLGS$depth)
   } 
   else {                ### If instead input is a vector, return a vector
     as.vector(xvlgs)
    }   
}

