\name{plotBins}
\alias{plotBins}
\title{ Plot z-scores of bins of markers }
\description{
\code{plotBins} merges \code{\link{plotLGS}} and \code{\link{zbins}} in that it calculates the z-score of groups of markers (bins) and plots the score to produce an equivalent of a genome scan. Like other plots of genome scans, \code{plotBins} can be useful to visually inspect the outcome of a LGS experiment.
}
\usage{
plotBins(lgs= NULL, bins= NULL, bandw= 30, shading= c("white", "lightyellow2"), sigv= c(-1.8, -2.5, 1.8, 2.5), ...)
}
\arguments{
  \item{lgs}{ A vector of allele ratios or an object of class \code{LGS} (i.e. output of \code{\link{simLGS}}). }
  \item{bins}{ A vector of the same length as \code{allRatio} specifying which markers should be grouped together. If the argument passed to \code{allRatio} is a \code{\link{simLGS}} list, \code{zbins} will search for the \code{simLGS$ChrMap$Bins} vector. }
  \item{bandw}{ Width of the bands produced on the plot to distinguih adjacent bins. }
  \item{shading}{ A pair of colours to use to shade bins. }
  \item{sigv}{ Optional vector of significance thresholds to be plotted as horizontal lines }
  \item{\dots}{ Further arguments to pass to \code{\link{plot}}. }
}
\details{
  \code{plotBins} uses \code{\link{zbins}} to calculate z-scores. As in \code{\link{zbins}}, if a set of bins is provided both with \code{lgs} (as list of class \code{LGS} containing a \code{"Bins"} component) and as a argument to \code{bins}, this latter will have precedence.
  }
\value{
 A plot of z-scores vs bins is produced.
}

%\references{ ~put references to the literature/web site here ~ }
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
\seealso{ \code{\link{zbins}} to simply extract z-scores. }
\examples{
data(pcMap)
lgs1<- simLGS(offspring=100, scenario="GrowthAdvantage99", map=pcMap, gen.dist="phy")
# Plot allele ratios using plotLGS
plotLGS(lgs1)

# Plot z-score of bins in pcMap
# plotBins(lgs1)
# Plot z-scores of chromosomes
# plotBins(lgs1, bins=lgs1$ChrMap$Chr)
}
\keyword{ htest }

