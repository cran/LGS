\name{biasLGS}
\alias{biasLGS}
\title{ Bias a vector of allele ratios (LGS) }
\description{
  \code{biasLGS} reproduces the technical imprecision in the quantification of the allele ratios by adding a random value to each element of a vector of allele ratios.
}
\usage{
biasLGS(simLGS = simLGS, bias = 0.02, meanBias = 0)
}
\arguments{
  \item{simLGS}{ An LGS list output of \code{\link{simLGS}} or a vector of marker allele ratios. }
  \item{bias}{ Average value added or subtracted to the marker allele ratios. }
  \item{meanBias}{ Average bias across all markers. }
}
\details{
  This function uses \code{\link{rnorm}} with \code{mean=meanBias} and \code{sd=bias}. Values above 1 and below 0 are reset to 1 and 0 respectively.
}
\value{
  If \code{simLGS} is an LGS list then \code{biasLGS} will create a new list equal to the input one but with \code{simLGS$ChrMap$LGS} replaced by the biased value and the original \code{simLGS$ChrMap$LGS} renamed to \code{simLGS$ChrMap$unbiasedLGS}.
  If simLGS is a vector, then a biased vector of the same length will be returned.
}
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
\examples{
# Bias a matrix of allele ratios (rows= scans, columns= loci)
(mlgs<- matrix(nrow= 3, ncol= 10, data= runif(n=3*10, min=0, max=1)))
(mlgsB<- t(apply(mlgs, 1, function(x) biasLGS(x, bias=0.02))))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

