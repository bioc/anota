\name{anotaPerformQc}
\alias{anotaPerformQc}

\title{
  Perform quality control to ensure that the supplied data set is
  suitable for Analysis of Partial Variance (APV) within anota.
}

\description{
  Generates a distribution of interaction p-values which are compared to
  the expected NULL distribution. Also assesses the frequency of highly
  influential data points using dfbetas for the regression slope and
  compares the dfbetas to randomly generated simulation data. Calculates
  omnibus class effects.
}

\usage{
anotaPerformQc(dataT=NULL, dataP=NULL, phenoVec=NULL,
generatePlot=FALSE, file="ANOTA_Total_vs_Polysomal_regressions.pdf",
nReg=200, correctionMethod="BH", useDfb=TRUE, useDfbSim=TRUE,
nDfbSimData=2000, useRVM=TRUE, onlyGroup=FALSE, useProgBar=TRUE)
}

\arguments{
  \item{dataT}{A matrix with cytosolic mRNA data. Non numerical rownames
  are needed.}
  \item{dataP}{A matrix with translational activity data. Non numerical rownames
  are needed.}
  \item{phenoVec}{A vector describing the sample classes (each class
    should have a unique identifier). Note that dataT, dataP and phenoVec
    must have the same sample order so that column 1 in dataP is
    the translational activity data for a sample, column 1 in dataT is the
    cytosolic mRNA data and position 1 in phenoVec describes the sample
    class.}
  \item{generatePlot}{anota can plot the regression for each gene. However,
    as there are many genes, this output is normally not
    informative. Default is FALSE, no individual plotting.}
  \item{file}{If generatePlot is set to TRUE use file to set
    desired file name (prints to current directory as a pdf).
    Default is "ANOTA_Total_vs_Polysomal_regressions.pdf"}
  \item{nReg}{If generatePlot is set to TRUE, nReg can be used to
    limit the number of output plots. Default is 200.}
  \item{correctionMethod}{anota adjusts the omnibus interaction and
    sample class p-values for multiple testing. Correction method can be "Bonferroni", "Holm",
    "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH" or
    "TSBH" as implemented in the multtest package or "qvalue" as
    implemented in the qvalue package. Default is "BH".}
  \item{useDfb}{Should anota assess the occurrence of highly influential
    data points (defult is TRUE)?}
  \item{useDfbSim}{The random occurrence of dfbetas can be
    simulated. Default is TRUE. FALSE represses simulation which reduces
    computation time but makes interpretation of the dfbetas difficult.}
  \item{nDfbSimData}{If useDfbSim is TRUE the user can select the
    number of samplings that will be performed per step (10 steps with
    different correlations between the translationally activty and the
    cytosolic mRNA level). Default is 2000.}
  \item{useRVM}{The Random Variance Model (RVM) can be
    used for the omnibus sample class comparison. In this case the effect of RVM
    on the distribution of the interaction significances needs to be
    tested as well. Default (TRUE) leads to calculation of RVM p-values
    for both omnibus interactions and omnibus sample class effects.}
  \item{onlyGroup}{It is possible to suppress the omnibus interaction
    analysis and only perform the omnibus sample class effect analysis. Default
    is FALSE (analyse both interactions and sample class effects.)}
  \item{useProgBar}{Should the progress bar be shown. Default is TRUE,
    show progress bar.}
}

\details{
  The anotaPerformQc performs the basic quality control of the data
  set. Two levels of quality control are assessed, both of which need to show good
  performance for valid application of anota. First, anota assumes that
  there are no interactions (for slopes). The output for this analysis is  both a density plot and a histogram plot of both the raw p-values and
  the p-values adjusted by the selected multiple correction method (if
  RVM was used, the second page shows the same presentation using RMV
  p-values). anota requires a uniform distribution of the raw
  interaction p-values for valid analysis of differential translation. anota also assesses if
  there are more data points with high influence on the regression analyses than
  would be expected by chance. anota identifies influential data points as
  data points that influence the slope of the regression using standardized dfbeta
  (dfbetas). In the literature there are multiple suggestions of what
  should be regarded as an outlier dfbetas (dfbetas>1, dfbetas>2,
  dfbetas>3, dfbetas>(2/sqrt(N)), dfbetas>(3/sqrt(N)),
  dfbetas>(3.5*IQR)). Independent of which threshold is preferred, what is of
  interest is the comparison to the underlying distribution. As this
  distribution is unknown, we simulate random data sets assuming that the
  cytosolic mRNA level and the translationally active mRNA levels are normally
  distributed and that there is a correlation between the cytosolic and the
  translationally active mRNA level. Following such simulation the frequencies of
  outlier dfbetas (using all thresholds) is compared to the frequencies found
  in the simulated data set. The function also performs an omnibus sample class
  effect test if there are more than 2 sample classes. It is possible to use RVM for the omnibus
  sample class statistics. If RVM is used, it
  is necessary to verify that the interaction RVM p-values also follow the
  expected NULL distribution.
  A rare error can occur when data within dataT or dataP from any gene
  and any sample class has no variance. This is reported as "ANOVA
  F-TEST on essentially perfect fit...". In this case those genes that
  show no variance for a sample class within either dataT or dataP need
  to be removed before analysis. Trying a different normalization method
  may fix the problem.
}

\value{
  anotaPerformQc generates several graphical outputs. One output
  ("ANOTA_interaction_p_distribution.pdf") shows the distribution of
  p-values and adjusted p-values for the omnibus interaction (both using
  densities and histograms). The second page of the pdf displays the
  same plots but for the RVM statistics if RVM is used. One output
  ("ANOTA_simulated_vs_obtained_dfbs.pdf") shows  bar graphs of the
  frequencies of outlier dfbetas using different dfbetas thresholds. If
  the simulation was enabled (recommended) these are compared to the
  frequencies from the random data set. One optional graphical output
  shows the gene by gene regressions with the sample classes
  indicated. In the case where RVM is used, a Q-Q plot and a comparison of the
  CDF of the variances to the theoretical CDF of the F-distribution is
  generated (output as "ANOTA_rvm_fit_for_....jpg") for both the omnibus
  sample class and the omnibus interaction test. The function also
  outputs a list object containing the following data:
  
  \item{omniIntStats}{A matrix with a summary of the statistics from the
    omnibus interaction analysis containing the following columns:
    "intMS" (the mean square for the interaction); "intDf" (the degrees
    of freedom for the interaction); "residMS" (the residual error mean
    square); "residDf" (the degrees of freedom for the residual error);
    "residMSRvm" (the mean square for the residual error after applying
    RVM); "residDfRvm"(the degrees of freedom for the residual error after
    applying RVM); "intRvmFval" (the F-value for the RVM statistics);
    "intP" (the p-value for the interaction); "intRvmP" (the p-value for
    the interaction using RVM statistics); "intPAdj" (the adjusted [for
    multiple testing using the selected multiple testing correction method]
    p-value of the interaction); "intRvmPAdj"(the adjusted [for
    multiple testing using the selected multiple testing correction method]
    p-value of the interaction using RVM statistics).}
  
  \item{omniGroupStats}{A matrix with a summary of the statistics from the
    omnibus sample class analysis containing the following
    columns:"groupSlope" (the common slope used in APV); "groupSlopeP" (if the slope is <0 or >1 a p-value
    for the slope being <0 or >1 is calculated; if the slope is >=0 & <=1
    this value is set to 1); "groupMS" (the mean square for sample classes);
    "groupDf" (the degrees  of freedom for the sample classes); "groupResidMS"
    (the residual error mean square); "groupResidDf" (the degrees of freedom
    for the residual error); "residMSRvm" (the mean square for the
    residual error after applying RVM); "groupResidDfRvm"(the degrees of
    freedom for the residual error after applying RVM); "groupRvmFval" (the F-value for the RVM statistics);
    "groupP" (the p-value for the sample class effect); "groupRvmP" (the p-value for
    the sample class effect using RVM statistics); "groupPAdj" (the adjusted [for
    multiple testing using the selected multiple testing correction method]
    p-value of the sample class effect); "groupRvmPAdj"(the adjusted [for
    multiple testing using the selected multiple testing correction method]
    p-value of the sample class effect using RVM statistics).}
  
  
  \item{correctionMethod}{The multiple testing correction method used to adjust the nominal p-values.}

  \item{dsfSummary}{A vector with the obtained frequencies of outlier
    dfbetas without the interaction term in the model.}
  
  \item{dfbetas}{A matrix with the dfbetas from the model without the
  interaction term in the model.}
  
  \item{residuals}{The residuals from the regressions without the
  interaction term in the model.}
  
  \item{fittedValues}{A matrix with the fitted values from the
    regressions without the interaction term in the model.}
  
  \item{phenoClasses}{The sample classes used in the analysis. The
    sample class order
    can be used to create the contrast matrix when identifying
    differential translation using anotaGetSigGenes.}
  
  \item{sampleNames}{A vector with the sample names (taken from the
    translationally active samples).}
  \item{abParametersInt}{The ab parameters for the inverse gamma fit for the interactions within RVM.}
  \item{abParametersGroup}{The ab parameters for the inverse gamma fit
    for sample classes within RVM.}
}

\author{Ola Larsson \email{ola.larsson@ki.se}, Nahum Sonenberg
  \email{nahum.sonenberg@mcgill.ca}, Robert Nadon \email{robert.nadon@mcgill.ca}}

\examples{
 ## See example for \code{\link{anotaPlotSigGenes}}

}

\seealso{\code{\link{anotaResidOutlierTest}}, \code{\link{anotaGetSigGenes}},\code{\link{anotaPlotSigGenes}}} 
\keyword{methods}
