\name{plotLGS}
\alias{plotLGS}
\title{ Plot the output of a LGS experiment as a genome scan }
\description{
  \code{plotLGS} plots marker allele ratios vs genome position to produce an overview of a genome scan. Typically it takes the output from \code{simLGS}. User's specified values to plot real or simulated data can also be supplied. Optional graphs to summarize the progeny population before and after selection are also implemented.
}
\usage{
plotLGS(simLGS = NULL, cumDist = simLGS$ChrMap$CumDist, chr = simLGS$ChrMap$Chr, allRatio = simLGS$ChrMap$LGS, width.cm = 23/2.54, height.cm = 11/2.54, barp = FALSE, add = FALSE, show.range = TRUE, new.window = TRUE, \dots)
}
\arguments{
  \item{simLGS}{ Output list from \code{simLGS}. Leave \code{NULL} if user's values are supplied.}
  \item{cumDist}{ Vector of cumulative marker distances.}
  \item{chr}{ Vector of chromosomes names corresponding to the respective markers.}
  \item{allRatio}{ Vector of allele ratios. }
  \item{width.cm}{ Width of the plotted window.}
  \item{height.cm}{ Height of the plotted window. }
  \item{barp}{ Should a barplot of population composition be shown? }
  \item{add}{ Should the data be overplotted on the current graph? }
  \item{show.range}{ Show mean and sd of the allele ratios. }
  \item{new.window}{ Close the current device before opening the new one. }
  \item{\dots}{ Further arguments to be passed to \code{plot} e.g. \code{col}. }
}
\details{
  If the output of \code{simLGS} is supplied to \code{simLGS}, no other argument needs to be specified. Otherwise, \code{cumDist, chr} and \code{allRatio} have to be given.
  It is useful to set width and height to easily compare different scans using plots of the same size. However, the set size might not correspond exactly with the actual one.
  \code{barp} can be used only if \code{simLGS} is given since it requires the appropriate matrices. It might not be very meaningful if the two parental strains are recombinants.
  \code{show.range} will plot horizontal lines corresponding to: mean of allele ratios (solid, grey), mean +/- 2sd (dashed, grey) and mean +/- 3sd (dashed, red). Not plotted if \code{add = TRUE}.
  Marker positions are shown as a rug. Chromosome boundaries are marked by vertical lines.
  If \code{simLGS} is supplied as input, markers with selection coeffiecients greater than 0 or less than 0 (from \code{simLGS$ChrMap$SelCoef}) will be plotted as blue and red dots, respectively.
}
\value{
A plot is returned on the current device.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Dario Beraldi \email{dario.beraldi@ed.ac.uk} }

\seealso{ \code{\link{simLGS}}}
\examples{
# Simulate a genome with 8 chromosomes, 3 to 10 markers/chromosome
genome<- sort(rnorm(n= 52, mean= 10, sd= 3))
chromos<- rep(1:8, times=3:10)
allr<- biasLGS(rep(0.5, 52)); allr[20]<- 0.01 # Selection on the second marker on chr 2
plotLGS(cumDist=genome, chr= chromos, allRatio= allr)

# Add another profile and change line style
allr2<- biasLGS(rep(0.7, 52)); allr2[30]<- 0.1 # Selection on the second marker on chr 2
plotLGS(cumDist=genome, chr= chromos, allRatio= allr2, add=TRUE, col="blue", lty="dotted")
}
\keyword{ aplot }