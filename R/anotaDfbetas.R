######################################################################
######################################################################
anotaSimDfbs <- function(nData=2000, phenoVec=phenoVec, mode=mode, useProgBar=useProgBar){
    nSamples <- length(phenoVec)
    ##Create matrix for simulated data
    dataPMatrix <- dataTMatrix <- matrix(ncol=nSamples, nrow=nData, data=NA)
    rownames(dataPMatrix) <- rownames(dataTMatrix) <- c(1:nData)
    ##The analysis is performed over a range of realistic correlations which are defined by the sd of the covariate.
    ##We have shown that this setting does not matter for the outcome because we use standardized dfbeta
    sdVec <- c(31:40/10)
    corMat <- matrix(ncol=length(sdVec), nrow=nData)
    ##A few structures to collect outputs
    pDfbCollect <- matrix(ncol=6, nrow=length(sdVec))
    row.names(pDfbCollect) <- sdVec
    corMeanMedCollect <- matrix(ncol=2, nrow=length(sdVec))
    rownames(corMeanMedCollect) <- sdVec
    colnames(corMeanMedCollect) <- c("Mean_correlation", "Median_correlation")
######################################################
    ##Start the analysis
    total <- length(sdVec)
    if(useProgBar==TRUE){
        pb <- txtProgressBar(min=0, max=total, style=3)
    }
    ##Over all the selected sds
    for(k in 1:length(sdVec)){
        if(useProgBar==TRUE){
            setTxtProgressBar(pb, k)
        }
        ##P and T is generated per gene
        for(j in 1:dim(dataPMatrix)[1]){
            ##sd=4 is arbitrary selected but gives a reasonable spread of the data and the setting does not influence the result.
            dataP <- rnorm(n=nSamples, mean=0, sd=4)
            dataT <- c(rep(NA, length(dataP)))
            for(i in 1:length(dataP)){
                dataT[i] <- rnorm(n=1, mean=dataP[i], sd=sdVec[k])
            }
            dataPMatrix[j,] <- dataP
            dataTMatrix[j,] <- dataT
        }
        ##Calculate and store the correlation per gene
        corVec <- c(rep(NA, nData))
        for(l in 1:length(corVec)){
            corVec[l] <- cor(dataPMatrix[l,], dataTMatrix[l,])
        }
        corMat[,k] <- corVec
        ##Store mean and median correlations
        corMeanMedCollect[k,1] <- mean(corVec)
        corMeanMedCollect[k,2] <- median(corVec)
        ##Run the dfb analysis and store the data
        dfbOut <- anotaDfbsSummaryOnly(dataT = dataTMatrix, dataP=dataPMatrix, phenoVec = phenoVec, mode=mode)
        ##Store the obtained dfb summary
        pDfbCollect[k,] <- as.vector(unlist(dfbOut$dfbSummary))
        colnames(pDfbCollect) <- colnames(dfbOut$dfbSummary)
    }
    cat("\n\n")
########################################
    ##Generate output for return
    outputList <- list(
                       "Correlation"=corMeanMedCollect,
                       "SimulatedDfbs" = pDfbCollect)
    return(outputList)
}


#############################################################################
#############################################################################
anotaDfbsSummaryFull <- function(lmDfb, mode, useDfbSim, filename, nDfbSimData, phenoVec, useProgBar=useProgBar){
    dsfSummary <- getSummaryDfb(lmDfb)

    ##Perform dfb simulation to get thresholds
    if(useDfbSim==TRUE){
        cat("\tPerforming dfbetas simulation\n")
        ##perform simulation
        nT <- length(phenoVec)
        ##mode decides how the simulation should be performed
        dfbSimOut <- anotaSimDfbs(nData=nDfbSimData, phenoVec=phenoVec, mode=mode, useProgBar=useProgBar)
        ##Get the output from anotaSimDfbs
        simDfbThreshold <- unlist(dfbSimOut$SimulatedDfbs)
        ##create a plot that compares obtined to simulated
        simDfbThresholdMean <- apply(simDfbThreshold, 2, mean)
        dsfSummary <- as.vector(unlist(dsfSummary))
        names(dsfSummary) <- c("dfb>1", "dfb>2", "dfb>3", "dfb>2/sqrt(N)", "dfb>3/sqrt(N)", "dfb>3.5*iqr")
        dfbDifference <- simDfbThresholdMean-dsfSummary
        dfbDifference <- round(dfbDifference, digits=4)
        tmpDfbMat <- cbind(simDfbThresholdMean, dsfSummary)
        tmpDfbMax <- apply(tmpDfbMat, 1, max)
        tmpVec <- c()
        for(r in 1:6){
            tmpVec <- c(tmpVec, dsfSummary[r], simDfbThresholdMean[r])
        }
        spaceVec <- c(0,0,rep(c(1,0),5))
        colVec <- rep(c(1,2),6)
        dfbsNames = c("dfb>1", "dfb>2", "dfb>3", "dfb>2/sqrt(N)", "dfb>3/sqrt(N)", "dfb>3.5*iqr")
        pdf(file=filename, width=6, height=4, pointsize=1/300)
        plot(x=c(0,17), y=c(0, c(max(tmpDfbMax))+0.01), ylab="Proportion of data points", xlab="Cut off method", main="Proportion of outliers in regression assessment using dfbetas", pch="", xaxt="n")
        barplot(tmpVec, space=spaceVec, col=colVec, add=TRUE)
        legend(x=1, y=max(tmpVec-0.002), legend=c("Obtained", "Simulated"), fill=c(1,2))
        text(y=c(tmpDfbMax+0.002), x=c(1, 4, 7, 10, 13, 16), labels=dfbDifference) 
        dev.off()
    }    
    ##Create a plotting output if no simulation was selected
    if(useDfbSim==FALSE){
        cat("\tNo dfbetas simulation is performed\n")
        dfbsNames <- c("dfb>1", "dfb>2", "dfb>3", "dfb>2/sqrt(N)", "dfb>3/sqrt(N)", "dfb>3.5*iqr")
        pdf(file=filename, width=6, height=4, pointsize=1/300)
        barplot(as.vector(unlist(dsfSummary)), names.arg=dfbsNames, cex.names=1, ylab="Proportion of data points", xlab="Cut off method", main="Proportion of outliers in regression assessment using dfbetas")
        dev.off()
    }
    return(dsfSummary)
}


#################################################################################
#################################################################################
anotaDfbsSummaryOnly <- function(dataT=NULL, dataP=NULL, phenoVec=NULL, mode=mode){
    ##Warnings
    if(is.null(dataT)){
        cat("ERROR: No data for total samples\n")
        stop()
    }
    if(is.null(dataP)){
        cat("ERROR: No data for polysomal samples\n")
        stop()
    }
    if(is.null(phenoVec)){
        cat("ERROR: No phenotypes specified\n")
        stop()
    }
    if(is.null(mode)){
        cat("ERROR: No mode specified\n")
        stop()
    }
    if(identical(rownames(dataT), rownames(dataP))==FALSE){
        cat("ERROR: Polysomal and Total rownames do not follow the same order\n")
        stop()
    }
################################################
    ##calculate n data points and factorise phenoVec
    nData <- dim(dataP)[1]
    phenoVec <- as.factor(phenoVec)
    ##structures to collect outputs
    lmDfb <- matrix(ncol=dim(dataT)[2], nrow=nData)
################################################
    ##Start analysis in a per gene loop
    for(i in 1:nData){
        ##lm model with interaction
        if(mode=="int"){
            tmpLm <- lm(dataP[i,]~dataT[i,]*phenoVec)
        }
        if(mode=="add"){
            tmpLm <- lm(dataP[i,]~dataT[i,]+phenoVec)
        }
        ##get dfbs for the slope i.e. in column 2
        tmpDfb <- dfbetas(tmpLm)
        lmDfb[i,] <- tmpDfb[,2]
    }
#################################################
    ##Evaluate outputs
    dfbSummary <- getSummaryDfb(lmDfb)
    
    ##Create a return object
    dataOut <- list(
                    "dfbSummary"=dfbSummary)
    return(dataOut) 
}

getSummaryDfb <- function(lmDfb){
        ##Evaluate outputs using different thresholds for dfbetas that have been suggested in the litterature.
    dfb1 <- lmDfb>1
    dfb2 <- lmDfb>2
    dfb3 <- lmDfb>3
    dfb2Sqrt <- lmDfb>(2/sqrt(dim(lmDfb)[2]))
    dfb3Sqrt <- lmDfb>(3/sqrt(dim(lmDfb)[2]))
    ##3.5X IQR. IQR is calculated per gene and each gene is tested separately but summarized across all.
    dfbIqr <- apply(lmDfb, 1, IQR)    
    dfbIqrTh <- 3.5*dfbIqr
    dfbIqrLog <- lmDfb>dfbIqrTh
    ##generate a summary output
    ##compare as a proportion of all tests    
    nTests <- ncol(lmDfb) * nrow(lmDfb)
    dfb1P <- sum(dfb1)/nTests
    dfb2P <- sum(dfb2)/nTests
    dfb3P <- sum(dfb3)/nTests
    dfb2SqrtP <- sum(dfb2Sqrt)/nTests
    dfb3SqrtP <- sum(dfb3Sqrt)/nTests
    dfbIqrLogP <- sum(dfbIqrLog)/nTests
    dsfSummary <- list(
                       "Proportion data points with dfb>1"=dfb1P,
                       "Proportion data points with dfb>2"=dfb2P,
                       "Proportion data points with dfb>3"=dfb3P,
                       "Proportion data points with dfb>2/sqrt(N)"=dfb2SqrtP,
                       "Proportion data points with dfb>3/sqrt(N)"=dfb3SqrtP,
                       "Proportion data points with dfb>3.5*iqr"=dfbIqrLogP)
    return(dsfSummary)
}
