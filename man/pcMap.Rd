\name{pcMap}
\alias{pcMap}
\docType{data}
\title{ A map file for Plasmodium chabaudi }
\description{
  This is an example of a map file based on the genome of P. chabaudi. This is in part a real data as it is based on real markers. Scenarios of selection are hypothetical but realistic.
}
\usage{data(pcMap)}
\format{
  A dataframe with 80 rows where each row is a marker with characteristics of map position and selection as:
    \enumerate{
        \item{\code{Chr}}{ Chromosome name. }
        \item{\code{Locus}}{ Locus/marker name (a factor). }
        \item{\code{uCumDist}}{ A numeric vector of user's cumulative distances (in Kosambi or Haldane cM or recombination frequency or physical distance). }
        \item{\code{Bins}}{ Bins grouping markers to be tested for selection (a factor). }
        \item{\code{EqualFitness99}}{ This and the remaining columns are numeric vectors represent scenarios of selection (see Details). }
    }
}
\details{
  Users should make sure that their map files contain at least the columns \code{Chr, Locus, uCumDist} (column names should match these names) and one or more scenario column (any name). The \code{Bins} column is instead optional. Additional columns are allowed and will be ignored by \code{LGS} functions. The order of the columns is irrelevant.
  
  The columns of scenarios of selection specify how intensively a locus is counterselected. Each value represents the probability of an allele to be eliminated by selection and it must be between -1 and 1. Therefore, 0 means no selection (neutral marker locus). A negative or positive number is the probability that an individual carrying the 0-allele or 1-allele at that locus will be counterselected, respectively. For example, the scenario \code{EqualFitness99} presents only one non-neutral marker (m46 on chr. 11) and selection on this locus is set to -0.99. This means that individuals with allele '1' at this locus will survive while individuals 
with '0' have 99% chance to be eliminted.
  }
\source{
  Katarzyna Modrzynska
}
\references{
    Experimental evolution, genetic analysis and genome re-sequencing reveal the mutation conferring artemisinin resistance in an isogenic lineage of malaria parasites.
    Hunt P. et al. BMC Genomics 2010.
}
\examples{
data(pcMap)
# Show map for chromosome 11 and the first two scenarios of selection
pcMap[pcMap$Chr=="11", 1:5]
}
\keyword{datasets}
