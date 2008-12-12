\name{haplo2string}
\alias{haplo2string}
\title{ Convert haplotypes to strings }
\description{
  Population matrices have one column per locus. \code{haplo2string} converts the columns of loci to a single column where alleles are collapsed to  a string representing haplotypes.
}
\usage{
haplo2string(x= NULL)
}
\arguments{
  \item{x}{ A population matrix with one column per locus and a final column of haplotype counts (\code{nclones}). }
}
\value{
  A data frame with columns:
  \item{haplotype}{ Haplotypes obtained by collapsing loci to string. }
  \item{nclones}{ Haplotype or clone counts as found in the original input matrix. }
}
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
\seealso{ \code{\link{paste}} to collapse vectors to strings. }
\examples{
(m1<- cbind(matrix(data=c(0,1), nrow=10, ncol=10), nclones=1:10))
haplo2string(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
