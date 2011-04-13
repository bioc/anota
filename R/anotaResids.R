anotaResidOutlierTest <- function(anotaQcObj=NULL, confInt=0.01, iter=5, generateSingleGraph=FALSE, nGraphs=200, generateSummaryGraph=TRUE, residFitPlot=TRUE, useProgBar=TRUE){
    ##Get data
    residualMatrix <- anotaQcObj$residuals
    nData <- dim(residualMatrix)[1]
    geneNames <- rownames(residualMatrix)
######################################################
    ##Warnings
    if(is.null(residualMatrix)){
        cat("ERROR:No residuals as input\n")
        stop()
    }
######################################################
    ##Single gene plotting
    if(generateSingleGraph==TRUE){
        pdf("ANOTA_residual_distributions_single.pdf", width=9, height=12)
        par(mfrow=c(4,3))
    }
######################################################
    ##Data structures
    residualMatrixOutlierSum <- matrix(ncol=iter, nrow=nData, dimnames=list("rownames"=rownames(residualMatrix)))
    residualMatrixOutlier <- matrix(ncol=dim(residualMatrix)[2], nrow=nData, dimnames=list("rownames"=rownames(residualMatrix), "colnames"=colnames(residualMatrix)))
    allSortScaleTrue <- allXs <- matrix(nrow=dim(residualMatrix)[2], ncol=nData)
    residualMatrixOutlierSumP <- rep(NA, iter)
##############################################
    ##Run analysis using the set n iterations
    cat("Running anotaResidOutlierTest\n")
    total <- iter
    if(useProgBar==TRUE){
        pb <- txtProgressBar(min=0, max=total, style=3)
    }
    for(j in 1:iter){
        if(useProgBar==TRUE){
            setTxtProgressBar(pb, (j))
        }
        ##generate the random normally distributed data. 
        rnormMat <- matrix(data=rnorm(dim(residualMatrix)[2]*((1/confInt)-1)),
                           nrow=dim(residualMatrix)[2], ncol=((1/confInt)-1))

        ##Calculate the upper and lower limits of the data and scale
        rnormMat <- apply(scale(rnormMat),2,sort)
        env <- t(apply(rnormMat, 1, range))
        for(i in 1:nData){            
            ##Scale true data per gene sort and cbind to rnorm data set
            trueVec <- sort(scale(residualMatrix[i,]))
            sampMat <- cbind(trueVec, rnormMat)
            ##get real data set qq
            rs <- sampMat[,1]
            xs <- qqnorm(rs, plot=FALSE)$x
            ##get range of the sampled distribution sort position
            ##Calculate if obtained residuals falls outside expected from rnorm
            rsLog <- rs<env[,1] | rs>env[,2]
            residualMatrixOutlierSum[i,j] <- sum(rsLog)
            residualMatrixOutlier[i,] <- rsLog
            ##save true data
            allSortScaleTrue[,i] <- trueVec
            allXs[,i] <- xs
            if(iter==j & generateSingleGraph==TRUE & nGraphs>=i){
                anotaResidOutlierPlot(xs=xs, rs=rs, env=env, geneName=geneNames[i])
            }
            
        }
        ##Collect data for single
        residualMatrixOutlierSumP[j] <- sum(residualMatrixOutlier>0)
    }
    cat("\n\n")
    if(generateSingleGraph==1){
        dev.off()
    }
################################################
    ##calcaulte obtained expected
    ##only create full summary for last iteration
    residualMatrixOutlierLog <- residualMatrixOutlier>0
    residualMatrixOutlierSumP <- residualMatrixOutlierSumP/(nData*dim(residualMatrixOutlier)[2])
    obtVsExpected <- residualMatrixOutlierSumP[j]/confInt
    expected <- nData *dim(residualMatrixOutlier)[2] *confInt
    obtained <- sum(residualMatrixOutlierLog)
################################################
    outputList <- list(
                       "confInt"=confInt,
                       "inputResiduals"=residualMatrix,
                       "rnormIter"=iter,
                       "outlierMatrixLog"=residualMatrixOutlierLog,
                       "meanOutlierPerIteration"=residualMatrixOutlierSumP,
                       "obtainedComparedToExpected"=obtVsExpected,
                       "nExpected"=expected,
                       "nObtained"=obtained)
#################################################
    ##Plotting summary
    if(generateSummaryGraph==TRUE){
        jpeg("ANOTA_residual_distribution_summary.jpeg", width=800, height=800, quality=100)
        anotaResidOutlierPlotAll(all=allSortScaleTrue, xsAll=allXs, env=env, obtained=obtained, expected=expected, obtRelExpected=obtVsExpected, confInt=confInt)
        dev.off()
    }
    ##plot fitted vs residuals
    if(residFitPlot==TRUE){
        jpeg("ANOTA_residual_vs_fitted.jpeg", width=900, height=900, quality=100)
        par(mfrow=c(2,1))
        plot(x=as.vector(anotaQcObj$fittedValues), y=as.vector(anotaQcObj$residuals), ylab="residuals", xlab="Fitted values", main="Residual vs fitted values")
        dev.off()
    }
    return(outputList)   
}


#################################################################################
#################################################################################
anotaResidOutlierPlot <- function(rs=NULL, xs=NULL, env=NULL, geneName=""){  
    matplot(xs, cbind(rs, env), type="pnn", pch=4, mkh=0.06, axes=FALSE, xlab="", ylab="", main=geneName)
    ##gets the limits and calculate lengths for the bars that mark the boundaries of the envelope
    xyul <- par("usr")
    smidge <- min(diff(c(xyul[1], xs, xyul[2])))/2
    segments(xs-smidge, env[,1], xs+smidge, env[,1])
    segments(xs-smidge, env[,2], xs+smidge, env[,2])
    ##get axis with nice ticks
    xul <- trunc(10*xyul[1:2])/10
    axis(side=1, at=seq(xul[1], xul[2], by=0.1), labels=FALSE, tck=0.01)
    xi <- trunc(xyul[1:2])
    axis(side=1, at=seq(xi[1], xi[2], by=0.5), tck=0.02)
    yul <- trunc(5*xyul[3:4])/5
    axis(side=2, at=seq(yul[1], yul[2], by=0.2), labels=FALSE, tck=0.01)
    yi <- trunc(xyul[3:4])
    axis(side=2, at=yi[1]:yi[2], tck=0.02)
    box(bty="l")
    mtext("Quantiles of Standard Normal", side=1, line=2.5, font=3)
    mtext(expression(R[1]), side=2, line=2, at=yul[2])
}

#####################################################################
#####################################################################
anotaResidOutlierPlotAll <- function(all=NULL, xsAll=xsAll,  env=env, obtained, expected, obtRelExpected, confInt){
    allMax=max(all)
    allMin=min(all)
    xsMin=min(xsAll)
    xsMax=max(xsAll)
    ##Get number of outliers per rankposition
    allLog <- all<env[,1] | all>env[,2]
    allLogSum <- apply(allLog, 1, sum)
    allLogSumP <- 100*(allLogSum / dim(allLog)[2])
    allLogSumP <- round(allLogSumP, digits=3)
    ##plot
    plot(x=c(xsMin,xsMax), y=c(allMin-0.2, allMax+0.4), pch="", axes=TRUE, xlab="Quantiles of standard normal", ylab="R", main="Summary of all residuals")
    for(i in 1:dim(all)[2]){
        points(x=xsAll[,i], y=all[,i], pch=16, cex=0.2)
    }
    if(is.null(obtained)==FALSE){
        text(x=xsMin+0.4, y=allMax-0.4-(0.05*allMax), labels=paste("Expected outliers: ", confInt*100, "%", sep=""))
        text(x=xsMin+0.4, y=allMax-0.4-(0.1*allMax), labels=paste("Obtained outliers: ", round(obtRelExpected*confInt*100, digits=3), "%", sep=""))
    }
    segments(xsAll[,i]-0.05, env[,1], xsAll[,i]+0.05, env[,1], col=2, lwd=1.5)
    segments(xsAll[,i]-0.05, env[,2], xsAll[,i]+0.05, env[,2], col=2, lwd=1.5)
    text(x=xsAll[,1], y=env[,2]+0.2, labels=paste(allLogSumP, "%", sep=""))
}
