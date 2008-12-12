\name{selClones}
\alias{selClones}
\title{ Reproduce clonally and select a population }
\description{
\code{selClones} clones by binary division each individual haplotype of a population matrix. Cloning will continue until a specified number of total individuals is reached. This function is similar to \code{\link{simLGS}} in that it reproduces and selects a starting population. In contrast to \code{\link{simLGS}}, \code{selClones} does not perform sexual reproduction and therefore no recombination occurs.

Note that \code{selClones} never eliminates individuals even if their genotype is strongly counter-selected. Rather, unfavoured genotypes have little chance to reproduce so that their frequency will dicrease as cloning progresses. (This behaviuor might or might not be biologically realistic).

}
\usage{
selClones(haplos= NULL, nhaplos= haplos[,ncol(haplos)], 
          selCoeff= rep(0, (ncol(haplos)-1)), nindiv= NULL)
}
\arguments{
  \item{haplos}{ A matrix representing a population (rows= haplotypes, columns= loci plus an additional column of haplotype counts). }
  \item{nhaplos}{ An optional vector of haplotype counts. By default, the count given in \code{haplos} (last column). }
  \item{selCoeff}{ A vector of selection coeafficients applied on each marker. Defaults to a vector of 0s (no selection). }
  \item{nindiv}{ The final number of individuals required. }
  }
\details{
  \code{nindiv} can be lower than the number of individuals of the starting population (i.e. nindiv <= sum(nhaplos)). In this case, a warning will be issued since the result might not be very meaningful as haplotypes are never removed from the population.
   
   Leaving \code{selCoeff} as no selection (vector of 0s) will simply cause the starting population to double in number until \code{nindiv} is reached. Individuals in excess will be removed.
   
   The cloning algorithm multiplies the exisiting number of haplotype counts by the selection coefficient established for each haplotype to produce the expected number of progeny individuals from each haplotype (e.g. if selection coefficient of a haplotype is 1, its count is expected to double each generation. If the coefficient is 0.1, is expected to increase by 10% each generation). The realized numebr of progeny clones is given by a random number from a poisson distribution with \code{lambda} equal to the expected clones.
   
   Cloning will stop at the generation that produce the total number of individuals nearest to \code{nindiv}. This should typically give a final number equal to \code{nindiv}+/-10%
  }
\value{
A list with components:
\item{LGS}{ A vector of allele ratios at each marker. }
\item{survPr}{ A vector of survival probabilities for each haplotype in \code{haplos}. }
\item{nhaplos}{ A vector of final counts of each haplotype. }
}

%\references{ ~put references to the literature/web site here ~ }
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
%\seealso{ \code{\link{zbins}} to simply extract z-scores. }
\examples{
data(pcMap)
pcMap$NoSel<- 0 # A scenario of no selection
lgs1<- simLGS(offspring= 2000, scenario="NoSel", map=pcMap, gen.dist="phy")
plotLGS(lgs1) # Just 1000 random gametes, no selection

# Clone and select. One population up to 100000, another up to 10000
sel1<- selClones(lgs1$selPop, nindiv=100000, selCoeff= pcMap$HybridAdvantage80)
sel2<- selClones(lgs1$selPop, nindiv=10000, selCoeff= pcMap$HybridAdvantage80)
# Allele ratios after cloning
plotLGS(allRatio= sel1$LGS, cumDist=pcMap$uCumDist, chr=pcMap$Chr, show.range=FALSE)
plotLGS(allRatio= sel2$LGS, cumDist=pcMap$uCumDist, chr=pcMap$Chr, add=TRUE, col="blue")

# 2000 haplotypes present but few of them account for most of the population!
topClones(sel1$nhaplo)
topClones(sel2$nhaplo)
plot(sort(sel1$nhaplos, decreasing= TRUE), 
     xlab="Haplotype ID", ylab="Haplotype count (log scale)", log="y",
     main="After cloning:\nFew haplotypes make most of the population")

}
\keyword{ datagen }

