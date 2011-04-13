
#####################################################################
#####################################################################
anotaGetSigGenes <- function(dataT=NULL, dataP=NULL, phenoVec=NULL, anotaQcObj=NULL, correctionMethod="BH", contrasts=NULL, useRVM=TRUE, useProgBar=TRUE){
    ##Collect data set info
    phenoVecOrg <- phenoVec
    phenoVec <- as.factor(phenoVec)
    nPheno <- length(levels(phenoVec))
    phenoLev <- levels(phenoVec)
    nData <- dim(dataT)[1]
########################################################
    ##Warnings
    if(identical(rownames(dataT), rownames(dataP))==FALSE){
        cat("ERROR: Polysomal and Total rownames do not follow the same order\n")
        stop()
    }
    if(is.null(dataT)){
        cat ("ERROR: No data for total samples\n")
        stop()
    }
    if(is.null(dataP)){
        cat ("ERROR: No data for polysomal samples\n")
        stop()
    }
    if(is.null(phenoVec)){
        cat("ERROR: No phenotypes specified\n")
        stop()
    }
    ##analysis halts if user gives custom contrasts that are too few or too many.
    if(is.null(contrasts)==FALSE){
        if(dim(contrasts)[2]!=(nPheno-1)){
            if(dim(contrasts)[2]>(nPheno-1)){
                cat("Too many custom contrasts supplied.\nPlease check your contrast matrix\n")
                stop()
            }
            if(dim(contrasts)[2]<(nPheno-1)){
                cat("Too few custom contrasts supplied.\nPlease check your contrast matrix\n")
                stop()
            }
        }
    }
    cat("Running anotaGetSigGenes\n")
    if(is.null(contrasts)){
        cat("\tUsing default \"treatment\" contrasts between (custom contrasts can be set):\n")
        cat("\t", paste(levels(phenoVec)))
        cat("\n")
    }
##########################################################
    ##Create data structures
    residRvmDf <- apvP <- apvRvmP <- apvPAdj <- apvRvmPAdj <- apvEff <- apvMSerror <- apvRvmMSerror<- apvT<- apvF <- apvFP <- apvRvmF <- apvN <- apvRvmDf <- apvDf <- matrix(nrow=nData, ncol=nPheno-1)
    residDf  <- unadjustedResidError <-apvSlope <- apvSlopeP <- c(rep(NA, nData))

    rownames(apvP) <- rownames(apvRvmP) <- rownames(apvRvmDf)<- rownames(apvDf)<- rownames(apvEff) <- rownames(apvRvmF)<- rownames(apvRvmF) <- rownames(apvRvmMSerror) <- names(residDf)<- names(apvSlope) <- names(apvSlopeP) <- rownames(residRvmDf)  <- rownames(dataP)

##########################################################
    ##lm model without interaction
    total <- nData
    for(i in 1:nData){
        ##create a list object with data
        tmpList <- list("PolyRNA"=dataP[i,], "TotalRNA"=dataT[i,], "phenoType"=phenoVec)
        tmpApvSum <- c()
        tmpDim <- c()
        tmpApv <- c()
        ##add custom contrasts to the analysis
        if(is.null(contrasts)==FALSE){
            contrasts(tmpList$phenoType) <- contrasts
        }
        if(is.null(contrasts)==TRUE){
            ##if no specific contrasts are added. Extract the treatment contrasts
            ##change to our coding scheme (-1 for control instead of the standard 0)
            tmpCont <- contrasts(tmpList$phenoType)
            tmpCont[1,] <- (-1)
            contrasts(tmpList$phenoType) <- tmpCont
        }
        if(i==1){
            ##test that contrasts are correct
            tmpOut <- contrasts(tmpList$phenoType)
            tmpColSum <- apply(tmpOut, 2, sum)
            if(sum(tmpColSum==0)!=length(tmpColSum)){
                cat("\nERROR: Incorrect contrasts. Each column of the contrast matrix should sum to 0\n\n")
                stop()
            }
            colnames(tmpOut) <- paste("contrast", c(1:dim(tmpOut)[2]))
            cat("\tThese contrasts will be evaluated:\n\n")
            print(tmpOut)
            cat("\n")
            if(useProgBar==TRUE){
                pb <- txtProgressBar(min=0, max=total, style=3)
            }
        }
        
        attach(tmpList)
        tmpApv <- lm(PolyRNA~TotalRNA+phenoType)
        detach(tmpList)
        tmpApvSum <- summary(tmpApv)

        ##Get the residual error and slope and slope P
        tmpError <- tmpApvSum$sigma
        unadjustedResidError[i] <- tmpError
        tmpSlope <- tmpApvSum$coefficients[2,1]
        apvSlope[i] <- tmpSlope
        
        ##Set slope P to 1. If estimate is <0 or >1 this will be modified
        apvSlopeP[i] <- 1
        if(apvSlope[i]<0 | apvSlope[i]>1){
            apvSlopeP[i] <- anotaSlopeTest(tmpLm=tmpApv, curSlope=apvSlope[i])
        }
        contMat <- contrasts(tmpList$phenoType)
        ##Go through the contrast matrix column by column and perform analysis
        for(j in 1:dim(contMat)[2]){
            ##select all groups that do not have a 0 in the contrast matrix
            ##notice that contMatRed becomes a vector
            contMatRed <-contMat[contMat[,j]!=0,j]
            ##create a summary matrix with the adjusted means, the contrast and the covariate means
            tmpContMat <- matrix(nrow=length(contMatRed), ncol=4)
            rownames(tmpContMat) <- names(contMatRed)
            colnames(tmpContMat) <- c("contrast", "estimate", "groupN", "T_mean")
            for(k in 1:dim(tmpContMat)[1]){
                tmpGroup <- rownames(tmpContMat)[k]
                tmpEst <- mean(dataP[i,phenoVecOrg==tmpGroup]) - (tmpSlope * mean(dataT[i,phenoVecOrg==tmpGroup]))
                tmpCovMean <- mean(dataT[i,phenoVecOrg==tmpGroup])
                tmpContMat[k,"contrast"] <- contMatRed[k]
                tmpContMat[k,"estimate"] <- tmpEst
                tmpContMat[k,"groupN"] <- sum(phenoVecOrg==tmpGroup)
                tmpContMat[k,"T_mean"] <- tmpCovMean
            }
            ##Get the diff effect for the set of contrasts
            ##getpostitive contrsts means
            tmpContMatCopy <- tmpContMat
            ##set contrast sums to 1
            tmpContMatCopy[tmpContMatCopy[,"contrast"]>0,"contrast"] <- tmpContMatCopy[tmpContMatCopy[,"contrast"]>0,"contrast"] /
                sum(abs(tmpContMatCopy[tmpContMatCopy[,"contrast"]>0,"contrast"]))

            tmpContMatCopy[tmpContMatCopy[,"contrast"]<0,"contrast"] <- tmpContMatCopy[tmpContMatCopy[,"contrast"]<0,"contrast"] /
                sum(abs(tmpContMatCopy[tmpContMatCopy[,"contrast"]<0,"contrast"]))
            tmpDiffEff <- sum(tmpContMatCopy[,"contrast"] * tmpContMatCopy[,"estimate"])
            ##Get covariate error sums of squares
            lmT <- lm(dataT[i,]~as.factor(phenoVecOrg))
            lmTAov <- anova(lmT)
            tmpCovSS <- lmTAov[2,2]
            ##Get the covariate mean for all observations
            tmpCovMean <-  mean(dataT[i,])
            ##Get the corrected MSerror and the N term 
            tmpCov <- sum(tmpContMatCopy[,"contrast"]*(tmpContMat[,"T_mean"]-tmpCovMean))
            ##nTerm checked 
            tmpN <- sum(tmpContMatCopy[,"contrast"]*tmpContMatCopy[,"contrast"]/tmpContMat[,"groupN"])
            tmpErrorAdjusted <- tmpError * tmpError * c(tmpN + (tmpCov * tmpCov/tmpCovSS))
            ##Collect data of import
            apvMSerror[i,j] <- tmpErrorAdjusted
            apvEff[i,j] <- tmpDiffEff
            apvDf[i,j] <- 1
        }
        residDf[i] <- tmpApvSum$df[2]
        if(useProgBar==TRUE){
            setTxtProgressBar(pb, i)
        }
    }
    cat("\n\n")
########################################
    ##Statistics
    apvF <- apvEff  * apvEff  / apvMSerror
    ##calculate and correct P
    for(j in 1:(nPheno-1)){
        apvFP[,j] <- 1- pf(apvF[,j], apvDf[,j], residDf)
        if(correctionMethod!="qvalue"){
            apvPAdj[,j] <- anotaAdjustPvals(pVals=apvFP[,j], correctionMethod=correctionMethod)
           }
        if(correctionMethod=="qvalue"){
            apvPAdj[,j] <- anotaAdjustPvalsQ(pVals=apvFP[,j])
        }
    }
    ##Generate an significance LIST object. One table for each contrast
    statsList <- list()
    for(j in 1:dim(contMat)[2]){
        tmpMat <- cbind(apvSlope, apvSlopeP,  unadjustedResidError, apvEff[,j], apvMSerror[,j],  apvF[,j], residDf, apvFP[,j], apvPAdj[,j])
        colnames(tmpMat) <- c("apvSlope", "apvSlopeP",  "unadjustedResidError", "apvEff", "apvMSerror", "apvF", "residDf", "apvP", "apvPAdj")
        statsList[[j]] <- tmpMat
    }
    contrastMat <- contrasts(tmpList$phenoType)
    contrastNames <- colnames(contMat)
###########################################################
    ##RVM to each contrast
    abList <- list()
    statsListRvm <- list()
    if(useRVM==TRUE){
        cat("\tCalculating RVM statistics\n")
        jpeg("ANOTA_rvm_fit_for_all_contrasts_group.jpg", width=800, height=c((nPheno-1)*400), quality=100)
        par(mfrow=c(c(nPheno-1), 2))
        ##Get ab and modify MS and DFs per contrast
        for(j in 1:c(nPheno-1)){
            cat("\tAssessing contrast", j, "\n")
            anotaPlotIGFit(apvMSerror[,j], residDf[1], qqName=paste("Fit for contrast", j))            
            ##adjust the groupResidMSs based on the ab parameters
            tmpRVM <- anotaPerformRVM(MS=apvEff[,j]*apvEff[,j], Df=apvDf[,j], residDf=residDf, residMS=apvMSerror[,j])
            apvRvmMSerror[,j] <- tmpRVM$residMSRvm
            residRvmDf[,j] <-tmpRVM$residDfRvm
            apvRvmF[,j] <- tmpRVM$rvmFval
            apvRvmP[,j] <- tmpRVM$rvmP
            abList[[j]] <- tmpRVM$ab
            ##correct
            if(correctionMethod!="qvalue"){
                apvRvmPAdj[,j] <- anotaAdjustPvals(pVals=apvRvmP[,j], correctionMethod=correctionMethod)
            }
            if(correctionMethod=="qvalue"){
                apvRvmPAdj[,j] <- anotaAdjustPvalsQ(pVals=apvRvmP[,j])
            }
        }
        dev.off()
        
        ##Add new data to summary tables
        for(j in 1:dim(contMat)[2]){
            tmpMat <- cbind(apvSlope, apvSlopeP, apvEff[,j], apvRvmMSerror[,j], apvRvmF[,j], residRvmDf[,j], apvRvmP[,j], apvRvmPAdj[,j])
            tmpColnames <- c("apvSlope", "apvSlopeP", "apvEff", "apvRvmMSerror", "apvRvmF", "residRvmDf", "apvRvmP", "apvRvmPAdj")
            colnames(tmpMat) <- tmpColnames
            statsListRvm[[j]] <- tmpMat
        }
    }
    
########################################
    outputList <- list(
                       "apvStats" = statsList,
                       "apvStatsRvm" = statsListRvm,
                       "correctionMethod"=correctionMethod,
                       "usedContrasts"=contrastMat,
                       "abList"=abList,
                       "groupIntercepts"=anotaQcObj$groupIntercepts,
                       "inputData"=list("dataT"=dataT, "dataP"=dataP, "phenoVec"=phenoVec),
                       "useRVM"=useRVM
                       )
    
    return(outputList)
}

######################################################################
######################################################################
anotaPlotSigGenes <- function(anotaSigObj, selIds=NULL, selContr=NULL, minSlope=NULL, maxSlope=NULL, slopeP=NULL, minEff=NULL, maxP=NULL, maxPAdj=NULL, maxRvmP=NULL, maxRvmPAdj=NULL, selDeltaPT=NULL, selDeltaP=NULL, sortBy=NULL, performPlot=TRUE, fileName="ANOTA_selected_significant_genes_plot.pdf", geneNames=NULL){
    if(is.null(anotaSigObj)==TRUE){
        cat("No significant object")
        stop()
    }

    if(is.null(selContr)==TRUE){
        cat("No contrast selected")
        stop()
    }
    dataT <- anotaSigObj$inputData$dataT
    dataP <- anotaSigObj$inputData$dataP
    phenoVec <- anotaSigObj$inputData$phenoVec

    ##geneNames is a matrix with rownames accoring to the same rownames as in the data matrices and the first column containing an additional identifier which will be used as the main title for the plots.
    useGeneNames <- FALSE
    if(is.null(geneNames)==FALSE){
        useGeneNames <- TRUE
        tmp <- rownames(geneNames)
        geneNames <- as.vector(geneNames)
        names(geneNames) <- tmp
    }
        
    ##stop if no contrast
    if(selContr>dim(anotaSigObj$usedContrasts)[2]){
        cat ("ERROR: Specified contrast does not exist")
        stop()
    }
    tmpData <- anotaSigObj$apvStats[[selContr]]
    ##The analysis may not include RVM. Therefore the RVM needs to be optional
    tmpDataRvm <- NULL
    useRVM <- anotaSigObj$useRVM
    if(useRVM==TRUE){
        tmpDataRvm <- anotaSigObj$apvStatsRvm[[selContr]]
    }

    ##select which Ids that will be used. Can either be custom OR a combination of different thresholds and flags.
    ##user has supplied IDs through the selIds command
    if(is.null(selIds)==FALSE){
        useIds <- selIds
        tmpData <- tmpData[useIds,]
        if(useRVM==TRUE){
            tmpDataRvm <- tmpDataRvm[useIds,]
        }
    }
    if(is.null(selIds)==TRUE){
        ##reduce data to current thresholds
        if(is.null(minSlope)==FALSE){
            tmpData <- tmpData[tmpData[, "apvSlope"]>minSlope,]
            if(useRVM==TRUE){
                tmpDataRvm <- tmpDataRvm[tmpDataRvm[, "apvSlope"]>minSlope,]
            }
        }
        
        if(is.null(maxSlope)==FALSE){
            tmpData <- tmpData[tmpData[, "apvSlope"]<maxSlope,]
            if(useRVM==TRUE){
                tmpDataRvm <- tmpDataRvm[tmpDataRvm[, "apvSlope"]<maxSlope,]
            }
        }
        
        if(is.null(slopeP)==FALSE){
            tmpData <- tmpData[tmpData[, "apvSlopeP"]>slopeP,]
            if(useRVM==TRUE){
                tmpDataRvm <- tmpDataRvm[tmpDataRvm[, "apvSlopeP"]>slopeP,]
            }
        }
        
        if(is.null(minEff)==FALSE){
            tmpData <- tmpData[abs(tmpData[, "apvEff"])>minEff,]
            if(useRVM==TRUE){
                tmpDataRvm <- tmpDataRvm[abs(tmpDataRvm[, "apvEff"])>minEff,]
            }
        }
        
        if(is.null(maxP)==FALSE){
            tmpNames <- rownames(tmpData[tmpData[,"apvP"]<maxP,])
            tmpData <- tmpData[tmpNames,]
            if(useRVM==TRUE){
                tmpDataRvm <- tmpDataRvm[tmpNames,]
            }
        }
        
        if(is.null(maxPAdj)==FALSE){
            tmpNames <- rownames(tmpData[tmpData[,"apvPAdj"]<maxPAdj,])
            tmpData <- tmpData[tmpNames,]
            if(useRVM==TRUE){
                tmpDataRvm <- tmpDataRvm[tmpNames,]
            }
        }
        
        if(is.null(maxRvmP)==FALSE & useRVM==TRUE){
            tmpNames <- rownames(tmpDataRvm[tmpDataRvm[,"apvRvmP"]<maxRvmP,])
            tmpData <- tmpData[tmpNames,]
            tmpDataRvm <- tmpDataRvm[tmpNames,]
        }
        
        if(is.null(maxRvmPAdj)==FALSE & useRVM==TRUE){
            tmpNames <- rownames(tmpDataRvm[tmpDataRvm[,"apvRvmPAdj"]<maxRvmPAdj,])
            tmpData <- tmpData[tmpNames,]
            tmpDataRvm <- tmpDataRvm[tmpNames,]
        }
        
        ##if user selects RVM based filtering but did not run RVM analysis a warning is needed
        if((is.null(maxRvmP)==FALSE | is.null(maxRvmPAdj)==FALSE) & useRVM==FALSE){
            cat("WARNING: You have selected a RVM based filtering but did not perform RVM based analysis using anotaGetSigGenes.\nThis filtering is therefore ignored.\nIf you intend to use RVM based filtering please rerun anotaSigGenes with RVM enabled.\n")
        }
        
        if(dim(tmpData)[1]<1){
            cat("ERROR: No genes pass selected thresholds.\nTry to change the thresholds.\n")
            stop()
        }
        
        useIds <- rownames(tmpData)
    }

    
    ##Get some delta translation and delta transcription values
    phenoLev <- levels(as.factor(phenoVec))
    deltaT <- matrix(ncol=length(phenoLev), nrow=length(useIds))
    deltaP <- matrix(ncol=length(phenoLev), nrow=length(useIds))
    for(i in 1:length(useIds)){
        for(j in 1:length(phenoLev)){
            deltaT[i,j] <- mean(dataT[useIds[i],phenoVec==phenoLev[j]])
            deltaP[i,j] <- mean(dataP[useIds[i],phenoVec==phenoLev[j]])
        }
    }
    deltaPT <- deltaP - deltaT
    rownames(deltaP) <- rownames(deltaT) <- rownames(deltaPT) <- useIds
    colnames(deltaP) <- colnames(deltaT) <- colnames(deltaPT) <- phenoLev
    
    ##allow filtering based on deltaP and deltaPT but also get the delta data and return to user
    ##P
    tmp1 <- t(t(deltaP[useIds,]) * anotaSigObj$usedContrasts[,selContr])
    tmp2 <- apply(tmp1, 1, sum)
    deltaP <- tmp2
    if(is.null(selDeltaP)==FALSE){        
        tmp3 <- (tmp2 > selDeltaP & anotaSigObj$apvStats[[selContr]][useIds,"apvEff"]>0) | (tmp2 < (-selDeltaP) & anotaSigObj$apvStats[[selContr]][useIds,"apvEff"]<0)  
        useIds <- useIds[tmp3==TRUE]
        tmpData <- tmpData[useIds,]
        if(useRVM==TRUE){
            tmpDataRvm <- tmpDataRvm[useIds,]
        }
    }
    ##T
    tmp1 <- t(t(deltaT[useIds,]) * anotaSigObj$usedContrasts[,selContr])
    tmp2 <- apply(tmp1, 1, sum)
    deltaT <- tmp2
    
    ##PT
    tmp1 <- t(t(deltaPT[useIds,]) * anotaSigObj$usedContrasts[,selContr])
    tmp2 <- apply(tmp1, 1, sum)
    deltaPT <- tmp2
    if(is.null(selDeltaPT)==FALSE){
        tmp3 <- (tmp2 > selDeltaPT & anotaSigObj$apvStats[[selContr]][useIds,"apvEff"]>0) | (tmp2 < (-selDeltaPT) & anotaSigObj$apvStats[[selContr]][useIds,"apvEff"]<0)  
        useIds <- useIds[tmp3==TRUE]
        tmpData <- tmpData[useIds,]
        if(useRVM==TRUE){
            tmpDataRvm <- tmpDataRvm[useIds,]
        }
    }
    deltaP <- as.matrix(deltaP[useIds])
    deltaT <- as.matrix(deltaT[useIds])
    deltaPT <- as.matrix(deltaPT[useIds])
    colnames(deltaP) <- "deltaP"
    colnames(deltaT) <- "deltaT"
    colnames(deltaPT) <- "deltaPT"

    ##allow sorting based on user selection
    if(is.null(sortBy)==FALSE){
        if(sortBy=="Eff"){
            tmpSorter <- sort(tmpData[useIds,"apvEff"], decreasing=TRUE)
            useIds <- names(tmpSorter)
        }
        if(sortBy=="p"){
            tmpSorter <- sort(tmpData[useIds,"apvP"], decreasing=FALSE)
            useIds <- names(tmpSorter)
        }
        if(sortBy=="rvmP" & (useRVM==TRUE)){
            tmpSorter <- sort(tmpDataRvm[useIds,"apvRvmP"], decreasing=FALSE)
            useIds <- names(tmpSorter)
            
        }
        if(sortBy=="rvmP" & (useRVM==FALSE)){
            cat("WARNING: You have selected RVM based sorting but did not perform RVM based analysis using anotaGetSigGenes\nThis sorting is therefore ignored.\nIf you intend to use RVM based filtering please rerun anotaSigGenes with RVM enabled.\n")
        }
    }
    
    ##Generate plots for those Ids that are left
    ##plotting is optional
    if(performPlot==TRUE){
        pdf(fileName, width=12, height=12)
        par(mfrow=c(3,3))

        for(i in 1:length(useIds)){
            tmpSlope <- anotaSigObj$apvStats[[1]][useIds[i], "apvSlope"]
            tmpTmin <- min(dataT[useIds,])
            tmpTmax <- max(dataT[useIds,])
            tmpPmin <- min(dataP[useIds,])
            tmpPmax <- max(dataP[useIds,])
            mainTitle <- paste(useIds[i], "Slope:", round(tmpSlope, digits=2))
            if(useGeneNames==TRUE){
                mainTitle <- paste(useIds[i], geneNames[useIds[i]], "Slope:", round(tmpSlope, digits=2))
            }
            plot(y=c(tmpPmin, tmpPmax), x=c(tmpTmin, tmpTmax), pch="", main=mainTitle, xlab="Transcription", ylab="Translation")
            phenoLev <- levels(as.factor(phenoVec))
            for(j in 1:length(phenoLev)){
                text(y=dataP[useIds[i],phenoVec==phenoLev[j]], x=dataT[useIds[i],phenoVec==phenoLev[j]], labels=phenoVec[phenoVec==phenoLev[j]], col=j)   
                tmpT <- mean(dataT[useIds[i],phenoVec==phenoLev[j]])
                tmpP <- mean(dataP[useIds[i],phenoVec==phenoLev[j]])
                tmpInt <- tmpP-(tmpSlope * tmpT)
                lines(x=c(tmpTmin, tmpTmax), y=c((tmpInt + tmpSlope*tmpTmin), (tmpInt + tmpSlope*tmpTmax)),  col=j)
            }
        
            ##add some stats for each contrast
            deltaLine <- 1
            deltaLine2 <- 0.5
            xShift <- 4
            ##no RVM
            plot(y=c(0, 10), x=c(0, 10), main=paste(useIds[i], "APV statistics without RVM"), pch="", xaxt="n", yaxt="n", xlab="", ylab="")
            nCont <- dim(anotaSigObj$usedContrasts)[2]
            lineCount <- 10
            xPos <- 2
            for(j in 1:nCont){
                sampCl1 <- rownames(anotaSigObj$usedContrasts)[anotaSigObj$usedContrasts[,j]<0]
                sampCl2 <- rownames(anotaSigObj$usedContrasts)[anotaSigObj$usedContrasts[,j]>0]
                tmpEff <- anotaSigObj$apvStats[[j]][useIds[i], "apvEff"]
                tmpP <- anotaSigObj$apvStats[[j]][useIds[i], "apvP"]
                tmpPadj <- anotaSigObj$apvStats[[j]][useIds[i], "apvPAdj"]
                text(y=lineCount, x=xPos, labels=paste("Contrast:", j), font=2, cex=1.2)
                lineCount <- lineCount - deltaLine2
                ##disabled because of problems with long contrast names or when several classes are combined to one
                ##text(y=lineCount, x=xPos, labels=paste(sampCl2, "vs.", sampCl1))
                ##lineCount <- lineCount - deltaLine2
                text(y=lineCount, x=xPos, labels=paste("Effect:", round(tmpEff, digits=2)))
                lineCount <- lineCount - deltaLine2
                text(y=lineCount, x=xPos, labels=paste("p-value:", round(tmpP, digits=4)))
                lineCount <- lineCount - deltaLine2
                text(y=lineCount, x=xPos, labels=paste("adjusted p-value:", round(tmpPadj, digits=3)))
                lineCount <- lineCount -deltaLine
                if(lineCount<3){
                    lineCount=10
                    xPos = xPos + xShift
                }
            }
            ##with RVM
            plot(y=c(0, 10), x=c(0, 10), main=paste(useIds[i], "APV statistics with RVM"), pch="", xaxt="n", yaxt="n", xlab="", ylab="")
            if(useRVM==TRUE){
                nCont <- dim(anotaSigObj$usedContrasts)[2]
                lineCount <- 10
                xPos <- 2
                for(j in 1:nCont){
                    sampCl1 <- rownames(anotaSigObj$usedContrasts)[anotaSigObj$usedContrasts[,j]<0]
                    sampCl2 <- rownames(anotaSigObj$usedContrasts)[anotaSigObj$usedContrasts[,j]>0]
                    tmpEff <- anotaSigObj$apvStatsRvm[[j]][useIds[i], "apvEff"]
                    tmpP <- anotaSigObj$apvStatsRvm[[j]][useIds[i], "apvRvmP"]
                    tmpPadj <- anotaSigObj$apvStatsRvm[[j]][useIds[i], "apvRvmPAdj"]
                    text(y=lineCount, x=xPos, labels=paste("Contrast:", j), font=2, cex=1.2)
                    lineCount <- lineCount - deltaLine2
                    ##disabled because of problems with long contrast names or when several classes are combined to one
                    ##text(y=lineCount, x=xPos, labels=paste(sampCl2, "vs.", sampCl1))
                    ##lineCount <- lineCount - deltaLine2
                    text(y=lineCount, x=xPos, labels=paste("Effect:", round(tmpEff, digits=2)))
                    lineCount <- lineCount - deltaLine2
                    text(y=lineCount, x=xPos, labels=paste("p-value:", round(tmpP, digits=4)))
                    lineCount <- lineCount - deltaLine2
                    text(y=lineCount, x=xPos, labels=paste("adjusted p-value:", round(tmpPadj, digits=3)))
                    lineCount <- lineCount -deltaLine
                    if(lineCount<3){
                        lineCount=10
                        xPos = xPos + xShift
                    }
                }
            }
        }
        dev.off()
    }
    tmpDataOut <- tmpData[useIds,]
    tmpDataRvmOut <- NULL
    if(useRVM==TRUE){
        tmpDataRvmOut <- tmpDataRvm[useIds,]
    }
    ##Generate output
    sigGroupIntercepts <- anotaSigObj$groupIntercepts[useIds,]
    tmpOut <- list("selectedData"=tmpDataOut,
                   "selectedRvmData"=tmpDataRvmOut,
                   "groupIntercepts"=sigGroupIntercepts[useIds,],    
                   "deltaData"=cbind("deltaP"=deltaP[useIds,], "deltaT"=deltaT[useIds,], "deltaPT"=deltaPT[useIds,]),                   
                   "usedThresholds" = list(
                   "selContr" = selContr,
                   "minSlope" = minSlope,
                   "maxSlope" = maxSlope,
                   "slopeP" = slopeP,
                   "minEff" = minEff,
                   "maxP" = maxP,
                   "maxPAdj" = maxPAdj,
                   "maxRvmP" = maxRvmP,
                   "maxRvmPAdj" = maxRvmPAdj,
                   "selDeltaPT" = selDeltaPT,
                   "selDeltaP" = selDeltaP)
                   )
    return(tmpOut)
}

#############################################################################
#############################################################################
anotaPlotIGFit <- function(useVar, df, title="Empirical", doQQ=TRUE, qqName=NULL) {
    ##Internal fuctions in anotaPlotIGFit
    ## Generate empirical probability density function for data.
    ## data - vector of data points
    ## q - trimming value.  remove 1-q points as outliers from greater tail.                                    
    my.cum <- function(data, q=.9) {
        len <- getLen(data, quan=q)
        maxi <- sort(data)[len]
        x <- seq(min(data), maxi, length=len)
        p <- rep(0,len-1)
        lenny <- length(data)
        for(i in 1:len)
            p[i] <- (sum(data<x[i]))/lenny
        return(cbind(x, p))
    }
#####################################
    flik <- function(p,y){
        ## log liklihood for a*b*x from an F distribution with m and 2*a degrees of freedom
        ## y is a vector containing data and the m values, p contains a and b.
        x<-y[1:(length(y)/2)]
        m<-y[(length(y)/2+1):length(y)]
        p<-abs(p)
        a<-p[1]
        b<-p[2]
        x<-x*(a*b)
        n<-2*a
        out<-base::log(df(x,m,n))+base::log(a*b)
        sum(-out)
    }
######################################
    ## Get rid of the very large data points so graphs scale better
    ## quan - trimming quantile
    getLen <- function(data, quan=.90) {
        return(trunc(quan*length(data)))
    }    
######################################
    degFreedom1 <- df
    vars <- as.vector(useVar)
    ab <- anotaGetab(vars, rep(degFreedom1, length(vars)))
    cat("\tThe a and b parameters for the inverse gamma distribution are:\n\t", paste(c("a:",ab[1], "b:", ab[2])), "\n")
    adj <- ab[1]*ab[2]
    adjVars <- vars*adj[1]
    scum <- my.cum(adjVars, q=0.9)
    probF <- pf(scum[,1], degFreedom1, 2*ab[1])
    lineWidth <- 2
    if(doQQ) {
        num <- length(adjVars)
        theoF <- rf(num, degFreedom1, df2=2*ab[1])
        qqplot(theoF, adjVars,xlab="Theoretical", ylab="Empirical", main=qqName)
        abline(0,1)
    }
    ##Turn off warnings temporary as there is always an informative warning message
    options(warn=(-1))
    kRes <- ks.test(x=adjVars, y="pf", degFreedom1, 2*ab[1])$p.value
    options(warn=0)
    plot(scum[,1], scum[,2], type="l", lwd=lineWidth, xlab="Var", ylab="cdf",
         main=paste(title, ": KS p-value=", signif(kRes,3), sep=""))
    lines(scum[,1], scum[,2], col=2)
    lines(scum[,1], probF, col=5, lwd=lineWidth)
    legend(x=c(1), y=c(0.3),  legend=c( title, "Theoretical F"), fill=c(2,5))
    options(warn=(1))
}

#############################################################################
#############################################################################
anotaGetab <- function(sig,n){
    flik <- function(p,y){
        ## log liklihood for a*b*x from an F distribution with m and 2*a degrees of freedom
        ## y is a vector containing data and the m values, p contains a and b.
        x<-y[1:(length(y)/2)]
        m<-y[(length(y)/2+1):length(y)]
        p<-abs(p)
        a<-p[1]
        b<-p[2]
        x<-x*(a*b)
        n<-2*a
        out<-base::log(df(x,m,n))+base::log(a*b)
        sum(-out)
    }
    set<-(!is.na(sig)&n>0&sig>0)
    sig<-sig[set]
    n<-n[set]
    set<-n>4
    if (sum(set)>0){
        m1<-(n[set]-2)/((n[set])*sig[set])
        m2<-(n[set]-2)*(n[set]-4)/((n[set])*sig[set])^2
        m1<-mean(m1,na.rm=TRUE)
        m2<-mean(m2,na.rm=TRUE)
        b<-m2/m1-m1
        a<-m1^2/(m2-m1^2)
    }
    else{ a<-b<-1}
    strt<-c(a,b)
###PATCH
    g <- function(p,yunq) flik(p,yunq)
########
    options(warn=(-1))
    a<-nlm(g,strt, yunq=c(sig,n))
    options(warn=0)
    a$estimate<-abs(a$estimate)
}
