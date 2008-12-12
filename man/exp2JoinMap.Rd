\name{exp2JoinMap}
\alias{exp2JoinMap}
\title{ Make input file for JoinMap }
\description{
JoinMap (Van Oijen et al. ref) is a useful program to estimate genetic distances in populations obtained through line crosses.

\code{exp2JoinMap} takes as input a population matrix obtained, for example, from \code{\link{simLGS}} or \code{\link{crossLGS}} and re-formats it to suite JoinMap 4 input.
}
\usage{
exp2JoinMap(popMatrix= NULL, markerNames= NULL, recombinantOnly= TRUE, outfile= "JMPtmp.loc")
}
\arguments{
  \item{popMatrix}{ Population matrix with rows= individuals and columns= loci. This can be either with or without a column of haplotype counts (\code{nclones}). }
  \item{markerNames}{ Vector of marker names to be exported. It will be created if left \code{NULL}. }
  \item{recombinantOnly}{ Should only the recombianant individuals be exported (see \code{Details})? }
  \item{outfile}{ Name of the destination file which will contain the input for JoinMap. }
}
\details{
JoinMap 4 has been used by the author to check whether the simulations generated recombination distances consistent with those specified in the map file. \code{exp2JoinMap} facilitated the process of re-formatting and exporting the matrices of populations to JoinMap. This function has been included in \code{LGS} even if it might not be required by most of the users.

\code{recombinantOnly= TRUE} is used to submit to JoinMap only those individuals that have originated by two different strains. This is because JoinMap is designed to analyze progenies from line crosses not from random mating. This setting will remove all the clones that present only '0' or '1' alleles. Therefore, individuals that originated from two different strains but had only '1' or '0' alleles would be erroneously removed. However, this is unlikely to happen unless the number of chromosomes and their size is small. (This behaviour could be changed in the future by keeping track of the pedigree that originated each haplotype.)
}
\value{
A file named \code{outfile} is written to the current working directory or to the directory specified in \code{outfile}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }

\examples{
data(pcMap)
map1<- cbind(pcMap, testJMP= rep(0, nrow(pcMap))) # "testJMP": scenario with no selection
lgs1<- simLGS(offspring= 1000, map= map1, scenario= "testJMP", fix.before.sel=TRUE, gen.dist="phy")

# Format and export progeny population, recombinants only 
exp2JoinMap(lgs1$recPop, outfile="testJMP.loc")

#Show first 20 lines of testJMP.loc
readLines("testJMP.loc", n=20)
}
\keyword{ utilities }
