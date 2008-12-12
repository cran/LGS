\name{simLGS}
\alias{simLGS}
\title{ Simulate crosses and apply selection }
\description{
\code{simLGS} generates a recombinant population with given genome and recombination map and applies selection at specified loci (see \code{\link{pcMap}} for a template of a map file). This is the most important function of the LGS package.
}
\usage{
simLGS(offspring = 50, fix.before.sel = FALSE, NoP1 = 10, NoP2 = 10, map = NULL, scenario = NULL, gen.dist = "k", kb2cM = 15)
}
\arguments{
  \item{offspring}{ Number of progeny clones to be generated. }
  \item{fix.before.sel}{ Should the number of clones be reached before or after selection? }
  \item{NoP1}{ If numeric: number of parental parasites. If a matrix, the population to be used as parental population. }
  \item{NoP2}{ Same as \code{NoP1} (i.e. the other parent). }
  \item{map}{ Table with map specification and selection coefficients. }
  \item{scenario}{ A column (name or index) specifying the selection coefficients at each locus. }
  \item{gen.dist}{ Metrics of the user's genetic distances (see details). }
  \item{kb2cM}{ Ratio kilobases/cM to convert physical distance to recombination frequency (if \code{gen.dist="phy"}). }
}
\details{
\code{simLGS} simulates random mating and recombiantion between two parental population that, by default, are polymorphic at all the loci specified in the map file \code{map}. By default, each parental individual in \code{NoP1} and  \code{NoP2} is represented by a vector of 0 and 1 alleles, respectively. 

Note that the absolute number in \code{NoP1} and  \code{NoP2} is irrelevant, what matters is their ratio.

\code{NoP1} and  \code{NoP2} can be matrices (rows= individuals and columns= loci) either used specified or obtained by other functions.

If \code{fix.before.sel} is set to \code{TRUE}, simLGS will generate \code{offspring} and then apply selection. If \code{fix.before.sel}=\code{FALSE}, simLGS will keep generate and select until \code{offspring} number is reached.
   
\code{scenario} must be a column present in \code{map}. See example map file \code{\link{pcMap}}.

\code{gen.dist} specifies which measure of marker distances is being used in \code{map}. It must be one of the following: \code{"k"} for Kosambi centi-morgan (the default), \code{"h"} for Haldane cM, \code{"rf"} for recombination frequency, \code{"phy"} for physical distance. If \code{"phy"} is chosen, \code{kb2cM} must be specified. 

The default of kb2cm is 15kb/cM and corresponds approximately to the ratio found in P. chabaudi.
}
\value{
  A list of class \code{LGS} with components:
  \item{recPop }{A dataframe of recombinant (progeny) population before selection.}
  \item{ParPop }{Dataframe of the parental population.}
  \item{selPop }{Dataframe of progeny population (\code{recPop}) after selection.}
  \item{ChrMap}{Dataframe corresponding to \code{map} plus a column for cumulative distances (\code{CumDist}) and a column of allele ratios (\code{LGS}).}
  Datafrmaes of populations have: rows= individuals (haplotypes), columns= loci plus an additional column \code{nclones} specifiying the number of clones represented by each haplotype.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }

\seealso{ \code{\link{plotLGS}} to plot the LGS profile; \code{\link{biasLGS}}, etc. }
\examples{
# Generate a map dataframe with 3 chromosomes, 5 markers/chromosomes, 10cM marker spacing.
# Third marker on chromosome 2 is under selection
SelCoef<- rep(0,15)
SelCoef[8]<- -0.99
(genMap<- data.frame(Chr= rep(1:3, each=5), Locus= seq(1,15),uCumDist= seq(0, 140, length.out=15), SelCoef= SelCoef))
lgs1<- simLGS(offspring =100, map= genMap)
# Plot LGS profile
plotLGS(lgs1)

# Import P. chabaudi map
data(pcMap)
# One locus under selection, sensitive allele is counterselected in 99% of the individuals
lgs99<- simLGS(map=pcMap, scenario= "EqualFitness99", gen.dist="phy") # Physical distances
plotLGS(lgs99)
# Same as above but selection is 80%
lgs80<- simLGS(map=pcMap, scenario= "EqualFitness80", gen.dist="phy") # Physical distances
plotLGS(lgs80, add=TRUE, col="blue")
}
\keyword{ datagen }
