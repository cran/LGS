\name{repLGS}
\alias{repLGS}
\title{ Replicate LGS experiments }
\description{
  \code{repLGS} replicates cross & selection via \code{\link{simLGS}} for a specified number of times. It is therefore useful to gather summary statistics about a given scenario of selection.
}
\usage{
repLGS(offspring= NULL, scenario= NULL, off.by.scen= TRUE, nsim= 1, map= NULL, fix.before.sel= FALSE, \dots)
}
\arguments{
  \item{offspring}{ Vector of number of offspring to be generated. }
  \item{scenario}{ Vector or list of scenarios, i.e. column names or column indexes in the map file that have to be simulated. }
  \item{off.by.scen}{ Should all the combinations of offspring by scenario be simulated? If \code{FALSE}, \code{offspring} and \code{scenario} must be of the same length. }
  \item{nsim}{ Vector of number of simulations for each combination of offspring by scenario. }
  \item{map}{ The map file that will be used by \code{\link{simLGS}}. }
  \item{fix.before.sel}{ Should the number of offspring in \code{offspring} be realized before or after selection? }
  \item{\dots}{ Further arguments to be passed to \code{\link{simLGS}}. }
  }
\details{
  \code{repLGS} can be used to replicate several times different combinations of population parameters and scenarios of selection. It is therefore useful to have an idea of the likely outcome of an LGS experiment or to plan an experiment before hand.
  
  \code{scenario} specifies which scenarios of selection should be used for the simulation. Column names or column indexes are both allowed but if a mixture of the two is used, they should be stored in list by means of \code{\link{list}}.
  
  \code{off.by.scen} determines whether each number offspring in \code{offspring} should be tested for each scenario in \code{scenario}. Therefore, if \code{off.by.scen=TRUE} and, for example, \code{offspring=c(10,50,100)} and \code{scenario=c("scenario1", "scenario2", "scenario3")}, then a total of nine combinations will be simulated and each combination will be replicated \code{nsim} times. If \code{off.by.scen=FALSE}, the  combinations that will be simulated will be simply those specified in the paired elements of \code{offspring} and \code{scenario}. In the example above, \code{offspring= 10} and \code{scenario= "scenario1"}, \code{offspring= 50} and \code{scenario= "scenario2"}, \code{offspring= 100} and \code{scenario= "scenario3"}. Again, each combination replicated \code{nsim} times.  
}
\value{
    A dataframe with a row for each simulation and columns representing:
	    \item{scenario}{ The scenario of selection. }
		\item{nRec}{ Number of offspring generated before selection is applied. This will be equal to the corresponding value in \code{offspring} if \code{sel.before.sel} in \code{\link{simLGS}} is \code{TRUE}}
		\item{nSel}{ Number of offspring after selection. Equal to the value in \code{offspring} if \code{fix.before.sel= FALSE}. }
	The remaining columns are the allele ratios for each marker in the map file. These columns are named:
	    \item{c1_1, c1_2 \\dots cn_n}{ Allele ratio at the first chromosome, first marker (\code{c1_1}, first chromosome, second marker (\code{c1_2} etc.))}
}

%\references{ ~put references to the literature/web site here ~ }
\author{ Dario Beraldi <\email{dario.beraldi@ed.ac.uk}> }
\examples{
data(pcMap)
offv<- c(5, 10, 15) # Investigate the hypothesis that 5, 10, 15 individuals survive selection
scenv<- c("EqualFitness99", "EqualFitness80")

rep1<- repLGS(offspring=offv, scenario=scenv, map=pcMap, gen.dist="phy", nsim= 5)
nrow(rep1) # = 3 offspring * 2 scenarios * 5 replicates = 30 rows
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }

