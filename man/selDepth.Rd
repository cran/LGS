\name{selDepth}
\alias{selDepth}
\title{ Find depth of a selection valley }
\description{
  \code{selDepth} returns the depth of a selection calculated as difference between baseline (genome-wide mean allele ratio) and allele ratio at the most extreme marker(s). The baseline can be calculated from all the markers or from selected regions.
}
\usage{
selDepth(x= NULL, notAt= NULL, at= NULL, baseline= NULL, out= "df", ties= TRUE)
}
\arguments{
  \item{x}{ A vector of allele frequencies (returns a warning if \code{x} values are not 0<x<1 ). }
  \item{notAt}{ An optional vector of marker positions to exclude from the search of the most divergent value. }
  \item{at}{ An optional vector to subset the whole set of markers. The most extreme markers will be searched in this subset. }
  \item{baseline}{ An optional vector to specify which markers should or should not be included in the calculation of the baseline. }
  \item{out}{ Output desired: "df" a dataframe with \code{maxDepth} and \code{atPos}, "depth" for only the depth, and  "pos" for only the position. }
  \item{ties}{ Should ties be returned? If \code{FALSE}, only the first marker found will be returned. }
}
\details{
   \code{selDepth} is useful to quickly find the deepest selection valley or highest selection peak.
   
   If \code{notAt} is left \code{NULL}, \code{selDepth} will scan all the markers and return the one with largest absolute difference with the baseline. If instead a vector of marker indexes is supplied, \code{selDepth} will return the absolute difference between baseline and the markers not included in \code{notAt}.
   
   \code{baseline} can be used to select regions to be used as baseline, for example chromosomes known to be neutrally selected, or to exclude regions known to be selected, for example chromosomes known to contain selected genes. Use positive or negative indexes to include or exclude regions, respectively.
   
   Option "df" in \code{out} combines both "depth" (maximum depth of selection) and "pos" (index of the marker).
}
\value{
    Options "depth" and "pos" in \code{out} will return a vector of depths (difference between marker and baseline) and positions (vector indexes) respectively. Option "df" combines the two and will return a dataframe with columns:
	    \item{maxDepth}{ Absolute difference between marker and baseline (output of "depth"). }
		\item{atPos}{ Marker index corresponding to where \code{maxDepth} was found (output of "pos"). }
}

%\references{ ~put references to the literature/web site here ~ }
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
\seealso{ \code{\link{which}} and \code{\link{\%in\%}} to include/exclude regions. }
\examples{
data(pcMap)
lgs1<- simLGS(offspring=200, scenario="GrowthAdvantage99", gen.dist="phy", map=pcMap)
plotLGS(lgs1)
scan1<- lgs1$ChrMap$LGS # LGS profile to search

selDepth(scan1) # Strongest selected markers

# Exclude chromosome 11 from the search of selection
d<- selDepth(scan1, notAt=which(pcMap$Chr \%in\% 11))
cbind(pcMap[d$atPos, 1:3], d$maxDepth)

# Don't use chromosomes under selection to compute the baseline
d<- selDepth(scan1, baseline= -which(pcMap$Chr \%in\% c(2,4,7,11,14)))
cbind(pcMap[d$atPos, 1:3], d$maxDepth)

# Investigate the depth of selection in 50 replicates of two scenarios
repli<- repLGS(offspring=c(50), scenario=c("EqualFitness99", "EqualFitness80"), nsim= 50, map=pcMap, gen.dist="phy")
dd<- apply(repli[,-c(1,2,3)], 1, selDepth, ties=FALSE, out="depth", at=46)
boxplot(dd~repli$scenario, ylab="Depth of selection")
}
\keyword{ htest }

