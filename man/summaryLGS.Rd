\name{summaryLGS}
\alias{summaryLGS}
\title{ Summary statistics of an LGS object }
\description{
  \code{summaryLGS} provides summary statistics about the allele ratio and populations of an LGS object.
}
\usage{
summaryLGS(lgs= NULL)
}
\arguments{
  \item{lgs}{ An object of class \code{LGS}, typically the output of \code{\link{simLGS}}. }
}
\value{
  A list with components:
  \item{LGS}{ Summary statistics for \code{lgs$ChrMap$LGS} using \code{\link{summary}}. }
  \item{pops}{ A dataframe with four columns and five rows containing summaries about the populations crossed and selected.
  First two columns summaries the two parental strains (\code{Parent1, Parent2}), the second two columns (\code{Progeny, SelProgeny}) relate to the progeny population before and after selection, respectively.
  For each population in the four colummns, rows report: \code{nSize}, the total count of individuals; \code{nhaplotypes}, the count of distinct haplotypes found; \code{nParent1, nParent2, nRecombinants}, count of parental haplotypes and recombinant haplotypes found in the progeny populations (therefore only for \code{Progeny} and \code{SelProgeny}).}
  \item{popFreq}{ A dataframe with two columns and three rows describing the progeny populations in terms of frequency of haplotypes. 
  Columns are \code{Progeny} and \code{SelProgeny} as above. 
  Rows are \code{nParent1, nParent2, nRecombinants} with the same information as above but expressed in terms of frequency.}
}
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
\examples{
data(pcMap)
lgs1<- simLGS(offspring= 100, map=pcMap, scenario= "EqualFitness99", gen.dist="phy")
summaryLGS(lgs1)

lgs2<- simLGS(offspring= 100, map=pcMap, scenario= "HybridAdvantage80", gen.dist="phy")
summaryLGS(lgs2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
