`crossLGS`<- function(
     cross= "self", generations= 0, ...)
{
lgs1<- simLGS(...)   ## First generation
  
if(generations==0){     ## If generations is 0, just return the output  of simLGS()
   return(lgs1)
   }     

lgsScans<- lgs1$ChrMap$LGS ## Here the LGS from each generation is archived

for(i in 1:generations){

if(cross=="self"){      ## Set parental populations for advanced crosses
   NoP1<- lgs1$selPop
   NoP2<- lgs1$selPop
   }
if(cross=="bcR"){
   NoP1<- lgs1$selPop
   NoP2<- nrow(lgs1$selPop)
   }
if(cross=="bcS"){
   NoP1<- nrow(lgs1$selPop)
   NoP2<- lgs1$selPop
   }

   lgs1<- simLGS(NoP1=NoP1, NoP2=NoP2, ...)     ## simLGS using populations above
   
   i<- i+1 ## Increase counter of generations
   lgsScans<- cbind(lgsScans, lgs1$ChrMap$LGS)
}
colnames(lgsScans)<- paste("gen", c(0, seq(1,generations)), sep="")
lgs1$ChrMap$LGS<- NULL
lgs1$ChrMap<- cbind(lgs1$ChrMap, lgsScans)
return(lgs1)
}