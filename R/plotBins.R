`plotBins`<- function(lgs= NULL, bins= NULL, bandw= 30, shading= c("white", "lightyellow2"), 
sigv= c(-1.8, -2.5, 1.8, 2.5), ...)
{
if(class(lgs)=="LGS"){
  LGS<- lgs$ChrMap$LGS
  if("Bins" %in% names(lgs$ChrMap)){
      bin2plot<- lgs$ChrMap$Bins
      }
}
if(!is.null(lgs)){
    LGS<- lgs
    }

if(!is.null(bins)){
    bin2plot<- bins
    }

### Calculate zscores ###

bin1<- zbins(LGS, bin2plot)

plot(bin1, xaxt="n", xlab="Marker bins (chr_bin)", ylab="z-score", ...)

### Bin shadings ###
bsh<- cbind(aggregate(1:length(bin2plot), list(bin2plot), min), coor=seq(1, length(unique(bin2plot))))
with(bsh, 
  abline(v=coor+0.1, lwd=bandw, col=ifelse((1:length(x))%%2!=0, shading[1], shading[2]))
  )
### Re-plot data ###
if(!is.null(sigv)){
    abline(h=sigv, col="darkblue", lty="dotted")
    }
points(bin1, ...)
axis(labels=unique(bin2plot), at=seq(1,length(bin1)), side=1, las=2, cex.axis=0.75)
box()
}
