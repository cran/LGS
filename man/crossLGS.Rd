\name{crossLGS}
\alias{crossLGS}
\title{ Perform advanced crosses }
\description{
  \code{crossLGS} uses \code{\link{simLGS}} to simulate several generations of backcrossing and selfing. Different mating ststems are available. 
}
\usage{
crossLGS(cross= "self", generations= 0, \dots)
}
\arguments{
  \item{cross}{ The type of experimental cross to simulate (see \code{Details})}
  \item{generations}{ Number of advanced generations to be simulated in addition to the first cross }
  \item{\dots}{ Arguments to be passed to \code{\link{simLGS}}}
}
\details{
\code{crossLGS} collects the selected population of \code{\link{simLGS}} (\code{simLGS$selPop}) and mates it with itself or with a parental pure strain. In each case, mating occurs at random.
  
The type of advanced cross is specified in \code{cross} and should be one of: \code{"self"} to perform random mating of the population; \code{"bcR"} or \code{"bcS"} to backcross the selected population to a population of resistant or susceptible parents respectively.

\code{generations} specifies how many rounds of mating and selection should be performed after the first mating. Note that the default, zero, implies a single event of mating and selection and it is therefore equivalent to simply using \code{\link{simLGS}}.

The arguments to be passed to \code{\link{simLGS}} should include at least the map file (\code{map} argument), the number of offspring desired before or after selection (\code{offspring}), and a scenario of selection (\code{scenario}, a column name or index in the map file).
}
\value{
A list with components:
\item{recPop}{ The recombinant population at the end of the cycles of crosses but before last selection is applied. }
\item{ParPop}{ Parental population used for the last cross. }
\item{selPop}{The \code{recPop} after the last round of selection}
\item{ChrMap}{The original map file with additional columns. Each additional column contains the allele ratios after each cycle of cross and selection. These are named: \code{'gen0'} for the first mating, then \code{'gen1'}, \code{'gen2'} etc. for the advanced crosses}
}
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
\examples{
#Load the example map file
data(pcMap)
#Let's imagine a cross bewteen two strains, one resistant to a drug 
#and the other susceptible. Resistance is confered by a single locus 
#and no other loci are under selection.
#The susceptible allele has 99% chance to be counter-selected.

#Five generations of backcrossing to the susceptible strain are simulated:
exp1<- crossLGS(cross="bcS", generations=5, map=pcMap, 
       offspring= 100, scenario= "EqualFitness99",
       gen.dist="phy")
bc<- exp1$ChrMap # Put map and allele ratios in 'bc' for convenience

#Plot the allele ratios of the first cross (gen0), third backcross (gen3),
#and fifth backcross (gen5)
plotLGS(cumDist=bc$CumDist, chr=bc$Chr, allRatio=bc$gen0, show.range=FALSE)
plotLGS(cumDist=bc$CumDist, chr=bc$Chr, allRatio=bc$gen3, col="red", add=TRUE)
plotLGS(cumDist=bc$CumDist, chr=bc$Chr, allRatio=bc$gen5, col="blue", add=TRUE)

#Calculate the baseline (genome-wide mean of the allele ratios), for each generation,
#i.e. first mating plus five backcrosses (allele ratios are in columns 13 to 18)
(gwar<- apply(bc[, 13:18],2, mean))

#Find the allele ratio at the selected locus and the difference with the baseline
#i.e. the depth of the selection valley
(m46<- bc[bc$Locus=="m46", 13:18])
(depth<- 1-(gwar-m46))

#...and plot it
plot(1:6, depth, xaxt="n", xlab="", ylab="Depth of selection valley", ylim=c(0,1),
     main="Effect of advanced backcrossing")
axis(at= 1:6, labels= c("F 1", paste("bc", 1:5)), side=1)
}
\keyword{ datagen }

