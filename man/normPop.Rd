\name{normPop}
\alias{normPop}
\title{ Normalize populations }
\description{
  \code{normPop} takes a list of populations or a list of haplotype counts and adjusts the proportion of each population to a specified one. The frequency of each haplotype within populations remains unchanged.
}
\usage{
normPop(hapList= NULL, proportions= NULL, popSize= 1, onePop= TRUE)
}
\arguments{
  \item{hapList}{ A list of population matrices or a list of vectors of haplotype counts. }
  \item{proportions}{ A vector of proportions of the same length as \code{hapList}. }
  \item{popSize}{ Final number of total individuals. }
  \item{onePop}{ Should the populations be merged into a single one? }
}
\details{
    Use \code{\link{list}} to return \code{hapList} as a list of matrices or vectors.
	
	\code{proportions} does not have to sum up to one. The default vector of proportions will make each population equally represented.

    \code{popSize} defaults to 1 although this is not biologically meaningful. However, in most cases what matters is the relative amount of each haplotype rather than its absolute number
}
\value{
  If \code{onePop=FALSE}, \code{normPop} will return a list of matrices and/or vectors if \code{hapList} is a list of matrices and/or vectors.
  
  If \code{onePop=TRUE}, \code{normPop} will return a single matrix (or single vector) with number of rows equal to the number of haplotypes in \code{hapList} and column of counts adjusted according to \code{proportions}.
  
  If \code{hapList} contains both matrices and vectors, \code{onePop} will be ignored and the value returned will be a list.
  
}
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
\examples{
# Two populations wih different number of haplotype counts
(pop1<- cbind(matrix(data=0, nrow=5, ncol=10), nclones=1:5))
(pop2<- cbind(matrix(data=1, nrow=5, ncol=10), nclones=10:14))

# Normalize to have each population equally represented (500 individuals each)
(n1<- normPop(list(pop1, pop2), onePop=FALSE, proportions=c(1,1), popSize=1000))
sum(n1[[1]][,11])
sum(n1[[2]][,11])

data(pcMap)
pop1<- cbind(matrix(data=0, nrow=10, ncol=80), nclones=1)
pop2<- cbind(matrix(data=1, nrow=10, ncol=80), nclones=1)

# Effect of varying proportions on LGS profiles
n11<- normPop(list(pop1, pop2), onePop=FALSE, proportions=c(1,1))
n51<- normPop(list(pop1, pop2), onePop=FALSE, proportions=c(5,1))
n201<- normPop(list(pop1, pop2), onePop=FALSE, proportions=c(20,1))

lgs11<- simLGS(offspring=100, NoP1=n11[[1]], NoP2=n11[[2]], 
       scenario="EqualFitness99", map=pcMap, gen.dist="phy")
lgs51<- simLGS(offspring=100, NoP1=n51[[1]], NoP2=n51[[2]], 
       scenario="EqualFitness99", map=pcMap, gen.dist="phy")
lgs201<- simLGS(offspring=100, NoP1=n201[[1]], NoP2=n201[[2]], 
       scenario="EqualFitness99", map=pcMap, gen.dist="phy")

plotLGS(lgs11, show.range=FALSE)
plotLGS(lgs51, add=TRUE, col="blue")
plotLGS(lgs201, add=TRUE, col="red")
legend("topleft", horiz=TRUE, legend=c("1:1", "5:1", "20:1"), 
        col=c("black", "blue", "red"), lty="solid", bg="white", pch=21, cex=0.85)

# Varying proportions in backcrosses
# First generation
lgsG1<- simLGS(offspring=100, scenario="EqualFitness99", map=pcMap, gen.dist="phy")
plotLGS(lgsG1, show.range=FALSE)

# Prepare backcross
pop1<- cbind(matrix(data=0, nrow=10, ncol=80), nclones=1)
popBC<- normPop(list(lgsG1$selPop, pop1), onePop=FALSE, proportions=c(1,2))
lgsBC1<- simLGS(offspring=100, NoP1=popBC[[1]], , NoP2=popBC[[2]], scenario="EqualFitness99", map=pcMap, gen.dist="phy")
plotLGS(lgsBC1, add=TRUE, col="blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

