\name{anotaPlotSigGenes}
\alias{anotaPlotSigGenes}

\title{
  Filter and plot genes to visualize/summarize genes that are
  differentially translated.
}

\description{
  This function filters the output from the anotaGetSigGenes function based on
  many user defined thresholds and flags to generate a summary table and
  optional per gene plots.
}

\usage{
anotaPlotSigGenes(anotaSigObj, selIds=NULL, selContr=NULL, minSlope=NULL, maxSlope=NULL, slopeP=NULL, minEff=NULL, maxP=NULL, maxPAdj=NULL, maxRvmP=NULL, maxRvmPAdj=NULL, selDeltaPT=NULL, selDeltaP=NULL, sortBy=NULL, performPlot=TRUE, fileName="ANOTA_selected_significant_genes_plot.pdf", geneNames=NULL)
}

\arguments{
  \item{anotaSigObj}{The output from the anotaGetSigGenes function.}
  \item{selIds}{The function can consider only a subset of the
    identifiers from the input data set (which can be further filtered)
    or used for custom plotting of identifiers of interest (leaving all
    filters as NULL). For custom selection of identifiers, supply a vector of identifiers
    (row names from the original data set) to be
    included. Default is NULL i.e. filtering is performed on all
    identifiers. Minimum length of selIds is currently 2. However, if
    only one identifier is of interest this identifier can be at
    position one and two of the supplied vector which will lead to that
    the data for the identifier of interested will be plotted twice.}
  \item{selContr}{Which contrast should be evaulated during the
    filtering, sorting and plotting? Descriptions of
    the contrasts can be found in the output from the anotaGetSigGenes
    object in the usedContrasts slot. Indicate the contrast by the column
    number.}
  \item{minSlope}{The output can be filtered so that genes whose
    identified slopes are too small can be excluded. Default is NULL
    i.e. no filtering based on lower boundary of the slope. To exclude
    genes with e.g. a slope <(-1) assign -1 to minSlope.}
  \item{maxSlope}{The output can be filtered so that genes whose
    identified slopes are too large can be excluded. Default is NULL
    i.e. no filtering based on upper boundary of the slope. To exclude
    genes with e.g. a slope >2 assign 2 to maxSlope.}
  \item{slopeP}{A p-value for the slope being <0 or >1 is calculated if the estimate
    for the slope is <0 or >1. This p-value can be used to
    filter the output based on unrealistic slopes. When set low fewer
    genes will be disqualified. Default is NULL i.e. no filtering based
    on slope p-value. We recommend setting slopeP between 0.01 and 0.1
    depending on data set characteristics.}
  \item{minEff}{The output can be filtered based on minimum effect for
    inclusion. The value is applied both to negative and positive effects:
    e.g. a value of 1 will evaluate if the effects are >1 OR
    <(-1). Default is NULL i.e. no filtering based on effect.}
  \item{maxP}{The output can be filtered based on raw p-values from the
    anota analysis without RVM  (i.e. smaller compared to assigned
    value). Default is NULL i.e. no filtering.}
  \item{maxPAdj}{The output can be filtered based on adjusted p-values from the
    anota analysis without RVM (i.e. smaller compared to assigned
    value). The adjustment method that was used when running
    anotaGetSigGenes will be evaluated. Default is NULL i.e. no filtering.}
  \item{maxRvmP}{The output can be filtered based on raw p-values from the
    anota analysis with RVM (i.e. smaller compared to assigned
    value). Default is NULL i.e. no filtering.}
  \item{maxRvmPAdj}{The output can be filtered based on adjusted p-values from the
    anota analysis with RVM (i.e. smaller compared to assigned value). The
    adjustment method that was used when running
    anotaGetSigGenes will be evaluated. Default is NULL i.e. no filtering.}
  \item{selDeltaPT}{The output can be filtered based on the mean log2(translational
    activity data / cytosolic mRNA data) between groups difference. The
    groups are defined by the selected contrast. Default is
    NULL i.e. no filtering.}
  \item{selDeltaP}{The output can be filtered based on the translational
    activity data only so that the minimum absolute between groups delta
    translation is used for gene inclusion. The groups are defined by
    the selected contrast. Default is NULL i.e. no filtering.}
  \item{sortBy}{The output can be sorted by effect ("Eff"), raw
    p-value("p") or raw RVM p-value ("apvRvmP"). Default is NULL i.e. no sorting.}
  \item{performPlot}{The function can generate a graphical output per
    gene. Default is TRUE i.e. generate plots.}
  \item{fileName}{The plots are printed to a file whose file name is
    given here. Default is "ANOTA_selected_significant_genes_plot.pdf".}
  \item{geneNames}{When anotaPlotSigGenes plots the individual gene plots
    they will be named by the original row names supplied to the
    anotaGetSigGenes function. geneNames allows the user to add additional
    names when plotting to e.g. include gene symbols. Input is a matrix
    with one column where the original row names match the row names of the input matrix and
    the desired new names are given in column 1.  Default is NULL
    i.e. no additional names.}
}

\details{
  This function allows the user to filter the output generated from the
  anotaGetSigGenes function to derive a reduced selection of genes that
  are considered for further evaluation. This is done by setting one or
  several of the filtering parameters described above. The function also
  generates a graphical output which helps when evaluating a single gene's
  regulation. In the graphical output, the results for each gene is
  displayed on separate rows. The first graph shows all samples and
  per sample class regression lines using the common slope with different colors
  for each sample class. The magnitude of the common slope is indicated. The second graph
  shows key statistics for the gene without the RVM model for all
  contrasts analyzed when running anotaGetSigGenes but any ordering and
  selection of genes is performed on the contrast given by the selContr
  argument as described above. The third graph is similar to the second but
  with RVM statistics instead (if RVM was used in the anotaGetSigGenes analysis).
}

\value{
  anotaPlotSigGenes generates a graphical output as described above and
  a list object containing summary data for those genes that passed the
  selected set of filters. The output list object contains the following
  slots:
  \item{selectedData}{A matrix containing non-RVM data for the filtered
    identifiers. Columns are "apvSlope" (the common slope used in APV);
    "apvSlopeP" (if the slope is <0 or >1 a p-value for the slope being <0
    or >1 is calculated; if the slope is >=0 & <=1 this value is set to
    1); "unadjustedResidError" (the residual error before calculating the
    effective residual error); "apvEff" (the group effect); "apvMSerror"
    (the effective mean square error); "apvF" (the F-value); "residDf"
    (the residual degrees of freedom); "apvP" (the p-value); "apvPAdj"
    (the adjusted p-value).} 
  \item{selectedRvmData}{A matrix containing RVM data for the filtered
    identifiers. Columns are "apvSlope" (the common slope used in APV);
    "apvSlopeP" (if the slope is <0 or >1 a p-value for the slope being <0
    or >1 is calculated; if the slope is >=0 & <=1 this value is set to
    1); "apvEff" (the group effect); "apvRvmMSerror" (the effective mean
    square error after RVM); "apvRvmF" (the RVM F-value); "residRvmDf"
    (the residual degrees of freedom after RVM); "apvRvmP" (the RVM
    p-value); "apvRvmPAdj" (the adjusted RVM p-value).}
  \item{groupIntercepts}{A matrix with the group intercepts, i.e. the
    translational activity for each group independent of cytosolic mRNA level. Can be used for
    e.g. clustering of translational activity. Data for all groups
    defined when using the anotaGetSigGenes function are supplied
    although the filtering is based on the contrast defined under the
    selContr argument.}
  \item{deltaData}{Mean delta translational activity data ("deltaP"), mean
    delta cytosolic mRNA data ("deltaT") and  mean delta log ratio data
  ("deltaPT") comparing the sample classes specified by the selected contrast.}
  \item{usedThresholds}{A list object with the user set values for the filtering.}
}

\author{Ola Larsson \email{ola.larsson@ki.se}, Nahum Sonenberg
  \email{nahum.sonenberg@mcgill.ca}, Robert Nadon \email{robert.nadon@mcgill.ca}}

\examples{
 ## Load the library and dataset (two phenotypes)
  library(anota)
  data(anotaDataSet)
 ## Quality control of the data set.
  anotaQcOut <- anotaPerformQc(dataT= anotaDataT[1:200,], dataP=anotaDataP[1:200,],
phenoVec=anotaPhenoVec, nDfbSimData=500)
 ##Test normality of residuals
  anotaResidOut <- anotaResidOutlierTest(anotaQcObj=anotaQcOut)
 ##Identify differentially translated genes.
  anotaSigGeneOut <- anotaGetSigGenes(dataT= anotaDataT[1:200,],
dataP=anotaDataP[1:200,], phenoVec=anotaPhenoVec, anotaQcObj=anotaQcOut)
 ##Plot some of the differentially expressed mRNAs
  anotSigGeneOutFiltered <-
anotaPlotSigGenes(anotaSigObj=anotaSigGeneOut, selContr=1,
maxP=0.05,slopeP=0.05, maxSlope=1.5, minSlope=(-0.5), selDeltaPT=0.5)
}

\seealso{\code{\link{anotaPerformQc}}, \code{\link{anotaResidOutlierTest}}\code{\link{anotaGetSigGenes}}} 
\keyword{methods}
