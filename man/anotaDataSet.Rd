\name{anotaDataSet}
\docType{data}
\alias{anotaDataSet}
\alias{anotaDataP}
\alias{anotaDataT}
\alias{anotaPhenoVec}
\title{Sample data set for anota}
\description{
  6 samples with data from 2 sample categories, both cytosolic
  (anotaDataT) and translational (anotaDataP) together with a sample
  class vector (anotaPhenoVec).
}

\usage{data(anotaDataSet)}


\format{Each data matrix (anotaDataT and anotaDataP) has 1000 rows
  (1000 first identifiers from complete data set) and 6 columns (noAA or
  rich). The anotaPhenoVec vector contains the sample class of each
  sample and anotaDataT, anotaDataP and phenoVec follow the same order.}

\source{Ingolia, NT et al. Genome-wide analysis in vivo of translation
  with nucleotide resolution using ribosome profiling. Science, 2009, 10;324(5924):218-23}

\examples{
##load data set
data(anotaDataSet)
##check dimensions
dim(anotaDataP)
head(anotaDataP)

dim(anotaDataT)
head(anotaDataT)

anotaPhenoVec
}
\keyword{datasets}
