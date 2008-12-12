`summaryLGS`<- function(lgs= NULL){
if(class(lgs)!="LGS"){
    stop("Argument to 'summaryLGS' must be of class 'LGS'.")
    }

summLGS<- summary(lgs$ChrMap$LGS)
collMatrix<- function(x){  # Function to collapse haplotypes to strings. x is population matrix
    xstr<- x[,-ncol(x)]
    xstr<- apply(xstr, 1, paste, collapse="")
    xnclones<- x[,ncol(x)]
    return(data.frame(haplotype=xstr, nclones=xnclones))
    }
p1p2<- row.names(lgs$ParPop)
sParPop<- collMatrix(lgs$ParPop)

p1<- sParPop[which(p1p2=="P1"),]
p2<- sParPop[which(p1p2=="P2"),]

dfRec<- collMatrix(lgs$recPop)
dfSel<- collMatrix(lgs$selPop)

p1inRecPop<- dfRec[dfRec$haplotype %in% p1$haplotype, ]
p1inSelPop<- dfSel[dfSel$haplotype %in% p1$haplotype, ]

p2inRecPop<- dfRec[dfRec$haplotype %in% p2$haplotype, ]
p2inSelPop<- dfSel[dfSel$haplotype %in% p2$haplotype, ]

nsize<- c(
    sum(p1$nclones),
    sum(p2$nclones),
    sum(dfRec$nclones),   ## Size of progeny, selected, parental pops
    sum(dfSel$nclones)
    )

nhaplotypes<- c(       ## Number of distinct haplotypes 
    length(unique(p1$haplotype)),
    length(unique(p2$haplotype)),
    length(unique(dfRec$haplotype)),
    length(unique(dfSel$haplotype))
    )

nParent1<- c(           ## Number of parental individuals found in the progeny and selected progeny
    NA,
    NA,
    sum(p1inRecPop$nclones),
    sum(p1inSelPop$nclones)
    )
    
nParent2<- c(           ## Number of parental strains found in the progeny and selected progeny
    NA,
    NA,
    sum(p2inRecPop$nclones),
    sum(p2inSelPop$nclones)
    )

nRecombinants<- nsize - (nParent1 + nParent2)  ## Number of recombinants (non parental)
    
popdf<- as.data.frame(
    rbind(nsize, nhaplotypes, nParent1, nParent2, nRecombinants), 
    row.names=c("nsize", "nhaplotypes", "nParent1", "nParent2", "nRecombinants")
    )   

colnames(popdf)<- c("Parent1", "Parent2", "Progeny", "SelProgeny")

popFreq<- as.data.frame(
    as.matrix(popdf[c(-1,-2),c(-1,-2)])/ rep(as.matrix(popdf[1,c(-1,-2)]), each=nrow(popdf)-2)
    )


return(list(
        LGS= summLGS, 
        pops= popdf,
        popFreq= popFreq
        )
       )       
}
