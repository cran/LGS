\name{bottleneck}
\alias{bottleneck}
\title{ Samples a population to simulate a bottleneck }
\description{
Given a population matrix or vector of haplotype counts, \code{bottleneck} reduces each count proportionally. In the context of microorganismsm, this function simulates events of sampling or random mortality.
}
\usage{
bottleneck(population= NULL, neck= NULL, keep.zero.count=FALSE)
}
\arguments{
  \item{population}{ A matrix representing a population (rows= haplotypes, columns= loci plus an additional column of haplotype counts) or a vector of haplotype counts }
  \item{neck}{ Total population size after culling. }
  \item{keep.zero.count}{ Should haplotype with zero count be retained? If \code{FALSE} such haplotypes will be dropped from the returned population }
  }
\details{
\code{bottlenck} calculates the expected number of clones after selection by multiplying the current number by a 'culling factor' code{neck/sum(haplotypeCounts)}. The realized number is given by a random number drawn from a poisoon distribution with expectation (/code{lambda}) given by the expected number of clones. Therefore, the final number \code{neck} is not exactly returned. 
  }
\value{
If \code{population} is a population matrix, the same matrix with haplotype count updated by the bottlenck is returned. If \code{population} is vector of counts, a vector is returned accordingly.
}

%\references{ ~put references to the literature/web site here ~ }
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
%\seealso{ \code{\link{zbins}} to simply extract z-scores. }
\examples{
## An hypothtical vector of haplotype counts
hapCount<- rpois(100, 10000)
bhap<- bottleneck(hapCount, 1000, keep.zero.count=TRUE)
plot(hapCount, bhap)
}
\keyword{ datagen }

