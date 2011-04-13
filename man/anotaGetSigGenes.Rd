\name{anotaGetSigGenes}
\alias{anotaGetSigGenes}

\title{
  Identify genes that are under translational control independent of
  cytosolic mRNA levels
}

\description{
  This function uses analysis of partial variance (APV) to identify genes that are under translational
  regulation independent of cytosolic mRNA levels.
}

\usage{
anotaGetSigGenes(dataT=NULL, dataP=NULL, phenoVec=NULL, anotaQcObj=NULL,
correctionMethod="BH", contrasts=NULL, useRVM=TRUE, useProgBar=TRUE)
}

\arguments{
  \item{dataT}{A matrix with cytosolic mRNA data. Non numerical rownames
  are needed.}
  \item{dataP}{A matrix with translational activity data. Non numerical rownames
  are needed.}
  \item{phenoVec}{A vector describing the sample classes (each class
    should have a unique identifier). Note that dataT, dataP and phenoVec
    have to have the same sample order so that column 1 in dataP is
    the translational data for a sample, column 1 in dataT is the
    cytosolic mRNA data and position 1 in phenoVec describes the sample
    class.}
  \item{anotaQcObj}{The object returned by anotaPerformQc.}
  \item{correctionMethod}{anota corrects p-values for multiple testing using the
    multtest package. Correction method can be "Bonferroni", "Holm",
    "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH" or
    "TSBH" as implemented in the multtest package or "qvalue" as
  implemented in the qvalue package. Default is "BH".}
  \item{contrasts}{When there is more than 2 sample categories it is
    possible to use custom contrasts. The order of the sample
    classes needs to be correct and can be seen in the object generated
    from anotaPerformQc in the phenoClasses slot (see details section).}
  \item{useRVM}{Should the Random Variance Model be applied. Default is
    TRUE.}
  \item{useProgBar}{Should the progress bar be shown. Default is TRUE,
    show progress bar.}
  
}

\details{
  The function performs APV on two or
  more sample categories. When more than two sample classes are
  compared it is possible to set custom contrasts to compare the sample
  classes of interest. Otherwise "treatment" contrasts are used which follow
  the alphabetical order of the sample classes. The order of the sample classes which
  the contrast matrix should follow can be found in the output of the
  anotaPerformQc function in the phenoClasses slot. Contrasts are supplied
  as a matrix where the sample classes are rows (same order as phenoClasses) and
  the columns are the different contrasts used. Contrasts are
  coded by using e.g. -1 for group a, 0 for group b and 1 for
  group c to compare group a and c; -2 for group a, 1 for group b and 1 for
  group c to compare group a to b & c. Each column of the contrast matrix
  should sum to 0 and to analyze orthagonal contrasts the products of
  all pairwise rows should sum to 0. The results will follow the order of the contrasts,
  i.e. the anocovaStats slot in the output object is a list with positions
  1...n where 1 is the first contrast and n is the last.

  A rare error can occur when data within dataT or dataP from any gene
  and any sample class has no variance. This is reported as "ANOVA
  F-TEST on essentially perfect fit...". In this case those genes that
  show no variance for a sample class within either dataT or dataP need
  to be removed before analysis. Trying a different normalization method
  may fix the problem.

}

\value{
  anotaGetSigGenes creates a plot showing the fit of the inverse gamma
  distribution used in RVM (similar output as from anotaPerformQc). anotaGetSigGenes also returns
  a list object with the following slots:
  \item{apvStats}{A list object (each slot named from 1 to the number
    of contrasts) where  each slot contains a matrix with statistics from
    the applied APV for that contrast. Columns are "apvSlope" (the
    common slope used in APV); "apvSlopeP" (if the slope is <0 or >1 a p-value
    for the slope being <0 or >1 is calculated; if the slope is >=0 & <=1
    this value is set to 1); "unadjustedResidError" (the residual error
    before calculating the effective residual error); "apvEff" (the group
    effect); "apvMSerror" (the effective mean square error); "apvF" (the F-value);
    "residDf" (the residual degrees of freedom); "apvP" (the p-value);
    "apvPAdj" (the adjusted p-value).}
  
  \item{apvStatsRvm}{A summary list object (each slot named from 1 to
    the number of contrasts) where each slot contains a matrix with RVM
    statistics from the applied APV. Columns are "apvSlope" (the
    common slope used in APV); "apvSlopeP" (if the slope is <0 or >1 a p-value
    for the slope being <0 or >1 is calculated; if the slope is >=0 & <=1
    this value is set to 1); "apvEff" (the group
    effect); "apvRvmMSerror" (the effective mean square error after RVM); "apvRvmF"
    (the RVM F-value); "residRvmDf" (the residual degrees of freedom
    after RVM); "apvRvmP" (the RVM p-value); "apvRvmPAdj" (the adjusted
    RVM p-value).}
  \item{correctionMethod}{The multiple testing correction  method used to adjust
    the p-values.}
  \item{usedContrasts}{A matrix with the contrasts used. Order is same as in the
    statistical outputs (column wise) so that the first contrast is
    found in the first slot of the apvStats and the apvStatsRvm lists.}
  \item{abList}{A list object containing the a and b parameters from the
    inverse gamma fits. Same order as the contrasts.}
}

\author{Ola Larsson \email{ola.larsson@ki.se}, Nahum Sonenberg
  \email{nahum.sonenberg@mcgill.ca}, Robert Nadon \email{robert.nadon@mcgill.ca}}

\examples{
 ## See example for \code{\link{anotaPlotSigGenes}}
}

\seealso{\code{\link{anotaPerformQc}},
  \code{\link{anotaResidOutlierTest}}, \code{\link{anotaPlotSigGenes}}} 
\keyword{methods}
