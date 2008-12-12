\name{topClones}
\alias{topClones}
\title{ Count the most represented clones}
\description{
  \code{topClones} takes a vector of haplotype counts and returns the minimum number of haplotypes (clones) that account for a given percentage of the whole population. It is useful to have an idea of how diverse a population is.
}
\usage{
topClones(x = NULL, np = c(0.68, 0.95, 0.99))
}
\arguments{
  \item{x}{ A vector of haplotype counts. It could be the output of \code{selClones}. }
  \item{np}{ Vector of percentages. \code{topClones} will return the minimal number of clones to reach these proportions (i.e. the most represented haplotypes). }
}
\value{
  A vector of number of haplotypes necessary to make the proportions specified in \code{np} (also printed on standard output).
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Dario Beraldi dario.beraldi@ed.ac.uk}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
# A population where few clones account for a large proportion of total individuals
set.seed(1234)
xv<- round(c(runif(n= 100, min=1, max=10), runif(n= 10, min= 2000, max= 10000))) # Each value is a count of different haplotypes
xy<- topClones(xv) # The six most represented haplotypes make 68% of all the individuals. The 22 top haplotypes make 99%.
}
\keyword{ manip }
