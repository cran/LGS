\name{zbins}
\alias{zbins}
\title{ z-score: Assessing selection }
\description{
  A criterion to assess the statistical significance of selection acting on a given marker or group of markers is to calculate its deviation from the mean of all the markers assayed. \code{zbins} groups adjacent markers in bins and calculates for each bin the zscore of the allele ratio relative to the whole set of bins. The bin z-score can therefore be considered a measure of signifincae for a putative selection valley.
}
\usage{
zbins(allRatio=NULL, bins=NULL)
}
\arguments{
  \item{allRatio}{ A vector of allele ratios or an object of class \code{LGS} (i.e. output of \code{\link{simLGS}}) }
  \item{bins}{ A vector of the same length as \code{allRatio} specifying which markers should be grouped together. If the argument passed to \code{allRatio} is a \code{link{simLGS}} list, \code{zbins} will search for the \code{simLGS$ChrMap$Bins} vector. }
}
\details{
  \code{zbins} first calculates the z-score of each marker, then assigns to each bin the mean z-score based on the markers contained in each bin. 
  If a vector of bins is supplied both within a \code{\link{simLGS}} list and thorugh \code{bins} argument, the vector in \code{bins} will be used.
  }
\value{
 A vector of z-scores, one for each different bin in \code{bins}.
}

%\references{ ~put references to the literature/web site here ~ }
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
\seealso{ \code{\link{plotBins}} to plot bin z-score as a genome scan. }
\examples{
data(pcMap)
# An example of how bins are constructed (show chr 14)
pcMap[pcMap$Chr=="14", 1:4]

# A simulation to work with
lgs1<- simLGS(offspring=100, scenario="EqualFitness99", gen.dist="phy", map=pcMap)
plotLGS(lgs1)

# Using bins provided in pcMap
bin1<- zbins(lgs1, pcMap$Bins)
# z-score at the most selected valley
min(bin1) # (Assuming we are interested in valleys, not peaks. 
          # Use max(bin1) or max(abs(bin1)) otherwise)
# Bin containing the lowest score
unique(pcMap$Bins)[which(min(bin1)==bin1)]

# z-score of each chromosome
binChr<- zbins(lgs1, pcMap$Chr)
# z-score of the most selected chromosome
min(binChr)
unique(pcMap$Chr)[which(binChr == min(binChr))]

# z-score of each marker
zMarker<- zbins(lgs1, pcMap$Locus)
# z-score of the most selected marker
min(zMarker)
# Marker at the bottom of the valley
pcMap$Locus[which(zMarker==min(zMarker))]

# In this case, maybe better
plotLGS(allRatio= zMarker, cumDist= lgs1$ChrMap$CumDist, chr= lgs1$ChrMap$Chr)

}
\keyword{ htest }

