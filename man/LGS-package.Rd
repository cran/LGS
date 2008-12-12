\name{LGS-package}
\alias{LGS-package}
\alias{LGS}
\docType{package}
\title{
LGS: A package to simulate Linkage Group Selection
}
\description{
Biological background:

LGS is an experimental procedure aimed at mapping loci under selection in organisms where single individuals cannot be genotyped. In a typical LGS experiment, two parental strains with different fitness are crossed and the progeny subjected to selection. By estimating the allele frequencies at molecular markers distributed across the whole genome, it is possible to identify the signatures of selection as skewed allele ratios.

Simulating LGS:

\code{LGS} simulates the biological steps uderlying a linkage group selection experiment. The simulation models are mostly based on the biology of Plasmodium chabaudi. It is assumed that individuals are haploid, mate randomly, and if required, reproduce clonally. Parental strains and progenies are represented by matrices where rows are individuals and columns are marker loci with alleles 0 or 1. Details of recombination frequency, marker distances, parental proportions and cross type can be specified by the user.

\code{LGS} will be useful to biologists and theoreticians to plan and interpret experiments. By comparing the outcome of simulations where all the variables are controlled with real data or real scenarios it will be possible to assess the influence of known variables and to detect unpredected factors.
    }
\details{
\tabular{ll}{
Package: \tab LGS\cr
Type: \tab Package\cr
Version: \tab 0.9\cr
Date: \tab 2008-12-12\cr
License: \tab GPL (version 2 or later)\cr
LazyLoad: \tab yes\cr
}
simLGS is the core function. Crosses strains and applies selection
}
\author{
Dario Beraldi

Maintainer: Dario Beraldi <dario.beraldi@ed.ac.uk>
}
\references{
    Martinelli et al. (Proc Natl Acad Sci U S A, 2005) for the biological and theoretical basis of LGS.

    Hunt et al. (BMC Genomics, 2010) for an application of the simulations implemented here.

    Beraldi et al. (in prep) for a description of the LGS simulations.
}
\keyword{ package }
\examples{
# Example to be added
}
