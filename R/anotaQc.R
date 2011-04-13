##Perform quality control of the data set
##Perform ominubus analysis
anotaPerformQc<- function(dataT=NULL, dataP=NULL, phenoVec=NULL, generatePlot=FALSE, file="ANOTA_Total_vs_Polysomal_regressions.pdf", nReg=200, correctionMethod="BH", useDfb=TRUE, useDfbSim=TRUE, nDfbSimData=2000, useRVM=TRUE, onlyGroup=FALSE, useProgBar=TRUE){
    ##Get some basic descriptors of the data set
    nData <- dim(dataP)[1]
    phenoVecOrg <- phenoVec
    phenoVec <- as.factor(phenoVec)
    phenoLev <- levels(phenoVec)
    nPheno <- length(phenoLev)
    
    ##Warnings
    ##do the data structures fit with each other?
    anotaPerformQcWarnings1(dataT=dataT, dataP=dataP, phenoVec=phenoVec)
    ##is there sufficient replication?
    anotaPerformQcWarnings2(nPheno=nPheno, phenoLev=phenoLev, phenoVecOrg=phenoVecOrg, onlyGroup=onlyGroup)
##############################
    ##initiation of objects for qc (dfbetas, interactions, slopes) and analysis (omnibus group, intercepts and rvm)
    ##dfbetas structures
    lmFittedValsAdd <- lmResidAdd <- lmDfbAdd <- matrix(ncol=dim(dataT)[2], nrow=nData)
    colnames(lmFittedValsAdd) <- colnames(lmDfbAdd) <- colnames(lmResidAdd) <- colnames(dataP)
    rownames(lmFittedValsAdd) <- rownames(lmDfbAdd) <- rownames(lmResidAdd) <- rownames(dataP)
    
    ##interactions structures
    intP <- intPAdj<- intMS <- intRvmFval <- intDf <- intResidMS <- intResidMSRvm <- intRvmP <- intRvmPAdj <- intResidDf <- intResidDfRvm <- c(rep(NA, nData))
    names(intP) <- names(intPAdj) <- names(intMS) <- names(intRvmFval) <- names(intDf) <- names(intResidMS) <- names(intResidMSRvm) <- names(intRvmP) <- names(intRvmPAdj) <- names(intResidDf) <- names(intResidDfRvm) <- rownames(dataP)

    ##slope structures
    groupSlope <- groupSlopeP <- c(rep(NA, nData))
    names(groupSlope) <- names(groupSlopeP) <- rownames(dataP)
    
    ##group structures
    groupP <- groupPAdj <- groupMS <- groupRvmFval <-  groupDf <-groupResidMS <-  groupResidMSRvm  <- groupRvmP <-groupRvmPAdj<- groupResidDf <- groupResidDfRvm <- c(rep(NA, nData))
    names(groupP) <- names(groupPAdj) <- names(groupMS) <- names(groupRvmFval) <-  names(groupDf) <- names(groupResidMS) <-  names(groupResidMSRvm)  <- names(groupRvmP) <- names(groupRvmPAdj) <- names(groupResidDf) <- names(groupResidDfRvm) <- rownames(dataP)

    ##intercept structure
    groupIntercepts <- matrix(nrow=nData, ncol=nPheno)
    colnames(groupIntercepts) <- phenoLev
    rownames(groupIntercepts) <- rownames(dataP)

    ##rvm structures
    names <- rownames(dataP)
    abInt <- rvmSummary <- rvmSummaryGroup <- abGroup <- dsfSummaryAdd <- NULL
##############################
    ##initiate optional regression plot
    if(generatePlot==1){
        geneNames <- rownames(dataP)
        pdf(file, width=8, height=11, pointsize=1/600)
        par(mfrow=c(4,2))
    }
    cat("Running anotaPerformQc quality control\n")
##############################
    ##Start analysis in a per gene loop
    cat("\tCalculating omnibus interactions & effects and dfbetas\n")
    total <- nData
    if(useProgBar==TRUE){
        pb <- txtProgressBar(min=0, max=total, style=3)
    }
    for(i in 1:nData){
        if(useProgBar==TRUE){
            setTxtProgressBar(pb, i)
        }
        tmpList <- list("PolyRNA"=dataP[i,], "TotalRNA"=dataT[i,], "phenoType"=phenoVec)
        attach(tmpList)
        if(onlyGroup==FALSE){
            ##perform regression with interactions
            tmpLm <- lm(PolyRNA~TotalRNA*phenoType)
            ##Get omnibus interactions and stats
            tmpLmAov <- anova(tmpLm)
            intDf[i] <- tmpLmAov[3,1]
            intMS[i] <- tmpLmAov[3,3]
            intP[i] <- tmpLmAov[3,5]
            intResidDf[i] <- tmpLmAov[4,1]
            intResidMS[i] <- tmpLmAov[4,3]
        }
        ##get omnibus group effects without interaction
        tmpLm <- lm(PolyRNA~TotalRNA+phenoType)
        tmpLmAov <- anova(tmpLm)
        detach(tmpList)
        groupDf[i] <- tmpLmAov[2,1]
        groupMS[i] <- tmpLmAov[2,3]
        groupP[i] <- tmpLmAov[2,5]
        groupResidDf[i] <- tmpLmAov[3,1]
        groupResidMS[i] <- tmpLmAov[3,3]
        
        ##if slope is <0 or >1 test if there is significant difference else set to 1.
        groupSlope[i] <- tmpLm$coefficients[2]
        groupSlopeP[i] <- 1
        if(groupSlope[i]>1 | groupSlope[i]<0){
            groupSlopeP[i] <- anotaSlopeTest(tmpLm=tmpLm, curSlope=groupSlope[i])
        }
        
        ##get dfbetas for the slope i.e. in column 2 from the no interaction model
        if(useDfb==TRUE){
            tmpDfb <- dfbetas(tmpLm)
            lmDfbAdd[i,] <- tmpDfb[,2]
        }
        
        ##Collect residuals and fittedvalues
        lmResidAdd[i,] <- tmpLm$residuals
        lmFittedValsAdd[i,] <- tmpLm$fitted.values
        
        ##get class intercepts
        groupIntercepts[i,] <- anotaGetIntercepts(x=dataT[i,], y=dataP[i,], slope=groupSlope[i], phenoVecOrg=phenoVecOrg, phenoLev=phenoLev)
        
        ##Plot single gene regression
        if(generatePlot==1 & i<=nReg){
            anotaPlotSingleRegression(x=dataT[i,], y=dataP[i,], geneName=geneNames[i], intercepts=groupIntercepts[i,], slope=groupSlope[i], phenoVecOrg=phenoVecOrg, phenoLev=phenoLev)
        }
    }
    cat("\n\n")
    ##End plotting
    if(generatePlot==1){
        dev.off()
    }
    ##done per gene analysis
#############################################
    ##Dfb analysis with or without simulation
    if(useDfb==TRUE){
        cat("\tAssessing dfbetas for model without interaction\n")
        dsfSummaryAdd <- anotaDfbsSummaryFull(lmDfb=lmDfbAdd, mode="add", filename="ANOTA_simulated_vs_obt_dfbetas_without_interaction.pdf", useDfbSim=useDfbSim, nDfbSimData, phenoVec=phenoVecOrg, useProgBar=useProgBar)
    }
    
#############################################
    ##RVM analysis
    ##Interactions
    if(useRVM==TRUE & onlyGroup==FALSE){
        cat("\tUsing RVM for omnibus interaction statistics\n")
        jpeg("ANOTA_rvm_fit_for_interactions.jpg", width=800, height=400, quality=100)
        par(mfrow=c(1,2))
        anotaPlotIGFit(intResidMS, intResidDf[1], qqName="Fit for interactions")
        dev.off()
        ##adjust the intResidMSs based on the ab parameters
        tmpRVM <- anotaPerformRVM(MS=intMS, Df=intDf, residMS=intResidMS, residDf=intResidDf)
        intResidMSRvm <- tmpRVM$residMSRvm
        intResidDfRvm <- tmpRVM$residDfRvm
        intRvmFval <- tmpRVM$rvmFval
        intRvmP <- tmpRVM$rvmP
        abInt <- tmpRVM$ab
    }
    rvmSummary <- cbind(intMS, intDf, intResidMS, intResidDf, intResidMSRvm, intResidDfRvm,  intRvmFval, intP, intRvmP)
    
    ##Group effect
    if(useRVM==TRUE){
        cat("\tUsing RVM for omnibus group statistics\n")
        jpeg("ANOTA_rvm_fit_for_omnibus_group.jpg", width=800, height=400, quality=100)
        par(mfrow=c(1,2))
        anotaPlotIGFit(groupResidMS, groupResidDf[1], qqName="Fit for omnibus group")
        dev.off()
        ##adjust the groupResidMSs based on the ab parameters
        tmpRVM <- anotaPerformRVM(MS=groupMS, Df=groupDf, residMS=groupResidMS, residDf=groupResidDf)
        groupResidMSRvm <- tmpRVM$residMSRvm
        groupResidDfRvm <- tmpRVM$residDfRvm
        groupRvmFval <- tmpRVM$rvmFval
        groupRvmP <- tmpRVM$rvmP
        abGroup <- tmpRVM$ab
    }
    rvmSummaryGroup <- cbind(groupSlope,groupSlopeP, groupMS, groupDf, groupResidMS, groupResidDf, groupResidMSRvm, groupResidDfRvm,  groupRvmFval, groupP, groupRvmP)
    
########################################
    ##Multiple testing adjustments
    cat("\tAdjusting p-values for multiple testing\n\n")    
    ##one set for when multtest is used
    if(correctionMethod!="qvalue"){
        if(onlyGroup==FALSE){
            intPAdj <- anotaAdjustPvals(pVals=intP, correctionMethod=correctionMethod)
        }
        groupPAdj <- anotaAdjustPvals(pVals=groupP, correctionMethod=correctionMethod)
        if(useRVM==TRUE){
            if(onlyGroup==FALSE){
                intRvmPAdj <- anotaAdjustPvals(pVals=intRvmP, correctionMethod=correctionMethod)
            }
            groupRvmPAdj <- anotaAdjustPvals(pVals=groupRvmP, correctionMethod=correctionMethod)
        }
        rvmSummary <- cbind(rvmSummary, intPAdj, intRvmPAdj)
        rvmSummaryGroup <- cbind(rvmSummaryGroup, groupPAdj, groupRvmPAdj) 
    }
    
    ##another set for when it is storey qvalue
    if(correctionMethod=="qvalue"){
        if(onlyGroup==FALSE){
            intPAdj <- anotaAdjustPvalsQ(intP)
        }
        groupPAdj <- anotaAdjustPvalsQ(groupP)
        if(useRVM==TRUE){
            if(onlyGroup==FALSE){
                intRvmPAdj <- anotaAdjustPvalsQ(intRvmP)
            }
            groupRvmPAdj <- anotaAdjustPvalsQ(groupRvmP)
        }
        rvmSummary <- cbind(rvmSummary, intPAdj, intRvmPAdj)
        rvmSummaryGroup <- cbind(rvmSummaryGroup, groupPAdj, groupRvmPAdj) 
    }
    
###############################
    ##Plot for interaction p-values
    if(onlyGroup==FALSE){
        anotaPlotIntPvals(intP=intP, intPAdj=intPAdj, intRvmP=intRvmP, intRvmPAdj=intRvmPAdj, useRVM=useRVM, correctionMethod=correctionMethod)
    }
################################
    ##Create a return object
    dataOut <- list(
                    "omniIntStats"=rvmSummary,
                    "omniGroupStats"=rvmSummaryGroup,
                    "groupIntercepts"=groupIntercepts,
                    "correctionMethod"=correctionMethod,
                    "dsfSummary"=dsfSummaryAdd,
                    "dfbetas"=lmDfbAdd,
                    "residuals"=lmResidAdd,
                    "fittedValues"=lmFittedValsAdd,
                    "phenoClasses"=levels(phenoVec),
                    "sampleNames"=colnames(dataP),
                    "abParametersInt"=abInt,
                    "abParametersGroup"=abGroup)
    
    return(dataOut)
}


anotaPlotIntPvals <- function(intP, intPAdj, intRvmP, intRvmPAdj, useRVM, correctionMethod){
    cAx <- 2
    cLab <- 2
    cMain <- 2
    pdf("ANOTA_interaction_p_distribution.pdf", height=10, width=10, pointsize=1/300)
    par(mar=c(5,6,4,2))
    par(mfrow=c(2,2))
    
    plot(density(intP), main="Omnibus interaction p-values", xlab="p-value", cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
    
    tmpMain <- paste("Omnibus interaction adjusted p-values", "(", correctionMethod, ")")
    plot(density(intPAdj), main=tmpMain, xlab="Adjusted p-value", cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
    hist(intP,  main="Omnibus interaction p-values",breaks=c(0:40)/40, xlab="p-value",
         cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
    
    hist(intPAdj,  main=tmpMain, breaks=c(0:40)/40, xlab="Adjusted p-value",
         cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
    
    if(useRVM==TRUE){
        plot(density(intRvmP), main="Omnibus interaction RVM p-values",
             xlab="RVM p-value", cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
        
        tmpMain <- paste("Omnibus interaction adjusted RVM p-values", "(", correctionMethod, ")")
        plot(density(intRvmPAdj), main=tmpMain, xlab="Adjusted RVM p-value",
             cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
        
        hist(intRvmP,  main="Omnibus interaction RVM p-values",
             breaks=c(0:40)/40, xlab="RVM p-value", cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
        
        hist(intRvmPAdj,  main=tmpMain, breaks=c(0:40)/40, xlab="Adjusted RVM p-value",
             cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
    }
    dev.off()    
}

anotaAdjustPvals <- function(pVals, correctionMethod){
    tmpAdj <- mt.rawp2adjp(as.vector(pVals), proc=correctionMethod)
    pAdj <- tmpAdj$adjp[order(tmpAdj$index),2]
    return(pAdj)
}
anotaAdjustPvalsQ <- function(pVals){
    tmpAdj <- qvalue(as.vector(pVals))
    pAdj <- as.vector(tmpAdj$qvalues)
    return(pAdj)
}
anotaPerformRVM <- function(MS, Df, residDf, residMS){
    ab <- anotaGetab(residMS, residDf)
    names(ab) <- c("a", "b")
    residMSRvm <- ((residDf*residMS)+2/ab[2])/(residDf+(2*ab[1]))
    residDfRvm <- residDf+2*ab[1]
    rvmFval <- MS / residMSRvm
    rvmP <- 1 - pf(rvmFval, Df, residDfRvm)
    return(list(residMSRvm=residMSRvm,
                residDfRvm=residDfRvm,
                rvmFval=rvmFval,
                rvmP=rvmP,
                ab=ab)
           )
}
anotaSlopeTest <- function(tmpLm, curSlope){
    tmpLmSum <- summary(tmpLm)
    if(curSlope<0){
        ##we are doing a one tailed test compared to the 2 tailed test in the output i.e. divide p by 2
        slopeP <- tmpLmSum$coefficients[2,4] / 2
    }
    if(curSlope>1){
        ##compare if slope is sig different to 1
        tmpSlopeEst <- tmpLmSum$coefficients[2,1] - 1
        tmpSlopeErr <- tmpLmSum$coefficients[2,2]
        tmpSlopeT <- tmpSlopeEst/tmpSlopeErr
        ##using resiual dfs to test p-value
        tmpSlopeDf <- tmpLm$df.residual
        slopeP <- 1-pt(tmpSlopeT, tmpSlopeDf)
    }
    return(slopeP)
}

anotaGetIntercepts <- function(x, y, slope, phenoLev, phenoVecOrg){
    tmpInt <- rep(NA, length(phenoLev))
    for(k in 1:length(phenoLev)){
        tmpX <- mean(x[phenoVecOrg==phenoLev[k]])
        tmpY <- mean(y[phenoVecOrg==phenoLev[k]])
        tmpInt[k] <- tmpY-(slope * tmpX)
    }
    return(tmpInt)
}

anotaPlotSingleRegression <- function(x, y, geneName, intercepts, slope, phenoVecOrg, phenoLev){
    plot(x=x, y=y, pch="", main=geneName)
    text(x=x, y=y, labels=phenoVecOrg)
    for(k in 1:length(phenoLev)){
        tmpMinX <- min(x, na.rm=TRUE)
        tmpMaxX <- max(x, na.rm=TRUE)
        lines(x=c(tmpMinX, tmpMaxX), y=c((intercepts[k] + slope*tmpMinX), (intercepts[k] + slope*tmpMaxX)),  lty=k)
    }
}

anotaPerformQcWarnings1 <- function(dataT=dataT, dataP=dataP, phenoVec=phenoVec){
    if(is.null(dataT)){
        cat("ERROR: No total RNA data\n")
        stop()
    }
    if(is.null(dataP)){
        cat("ERROR: No polysomal RNA data\n")
        stop()
    }
    if(is.null(phenoVec)){
        cat("ERROR: No phenotypes specified\n")
        stop()
    }
    if(identical(rownames(dataT), rownames(dataP))==FALSE){
        cat("ERROR: Polysomal and Total rownames do not follow the same order\nMake sure that the rownames are identical\n")
        stop()
    }
}

anotaPerformQcWarnings2 <- function(nPheno, phenoLev, phenoVecOrg, onlyGroup){
    ##require 3 samples per category unless using onlyGroup=TRUE
    for(i in 1:nPheno){
        if(sum(phenoVecOrg==phenoLev[i])<3){
            cat(paste(phenoLev[i], "only has", sum(phenoVecOrg==phenoLev[i]), "sample(s)"), "\n")
            if(onlyGroup==FALSE){
                cat("ANOTA requires 3 samples per sample class (in phenoVec) to run\n")
                cat("\tIf there is at least two samples per group and more than two groups, the onlyGroup mode can be used to assess omnibus group effects\n")
                stop()
            }
            if(onlyGroup==TRUE){
                if(sum(phenoVecOrg==phenoLev[i])==1){
                    cat("\tNo analysis can be run whith only one sample despite using the onlyGroup mode\n")
                    stop()
                }
                if(nPheno<3){
                    cat("\tThree sample classes is required to run the onlyGroup analysis when there is <3 samples per group")
                    stop()
                }
                if(nPheno>2){
                    cat("\tThere is only two samples in the sample class but >2 sample classes so onlyGroup analysis can be performed\n")
                }
            }
        }
    }
}
