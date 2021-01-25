#' A really ugly piece of code to plot variables from model run
#' 
#' @param variableData Numeric array. Data to be plotted
#' @param variableName Character string. The name of the variable
#' @param variableDims Numeric Vector. The dimensions of the \code{variableData}
#' @param otherList List. A list of species names, fleet names, file apth etc to help with labeling and saving of plots
#' 
#' @return a saved plot
#' 
#' This is the start of a generic function able to handle many kinds of data
#' Needs an overhaul

sub_plottingSuite <- function(variableData,variableName,variableDims,otherList) {
  
  speciesNames <- otherList$speciesNames
  nSpecies <- length(speciesNames)
  guildNames <- otherList$guildNames
  nGuilds <- length(guildNames)
  fleetNames <- otherList$fleetNames
  nFleets <- length(fleetNames)
  nRuns <- length(variableData) # should be a list
  howToScale <- otherList$scale
  nStepsyr <- otherList$nStepsyr
  nYrs <- otherList$nYrs
  B0 <- otherList$B0
  B0Guild <- otherList$B0Guild
  B0thresholds <- otherList$threshold
  filePath <- otherList$filePath
  
  isGuild <- as.logical(grep("guild",tolower(variableName)))
  if (length(isGuild) == 0) {isGuild <- as.logical(0)}


  dir.create(file.path(filePath, howToScale), showWarnings = FALSE,recursive = TRUE)
  
  # create directory to write files if not already created
  if (!dir.exists(file.path(filePath,howToScale))) {
    dir.create(file.path(filePath,howToScale))
  }
  
  # global plotting parameters and units
  colVec <- c("black","red","blue","yellow","green") 
  cexAxisVal <-  2.5
  subTitleSize <- 2.5
  titleSize <- 2.5
  xyLabelsize <- 3
  mtextSizey <- 2
  mtextSizeTitle <- 3
  
  ylabel <- ""
  if (as.logical(length(grep("index",variableName))))  {
    scaleData <- 1
    ylabel <- "index"
  } else if ((as.logical(length(grep("biomass",variableName)))) | (as.logical(length(grep("SSB",variableName)))) | (as.logical(length(grep("avByr",variableName))))) {
    # convert mt to kt
    scaleData <- 1000
    ylabel <- "kilotons, (kt)"
  } else if ((as.logical(length(grep("N",variableName)))) | (as.logical(length(grep("recruit",variableName))))) {
    scaleData <- 1
    ylabel <- "millions"
  } else {
    scaleData <- 1
  }
  
    
  png(filename=paste0(filePath,howToScale,"/",variableName,".png"),width=900,height=900,units="px")
  
  if (is.null(variableDims[1])) {
    dataMean <- apply(simplify2array(variableData)/scaleData,1,mean)
    dataMedQuant <- apply(simplify2array(variableData)/scaleData,1,quantile,probs=.5)
    dataLowQuant <- apply(simplify2array(variableData)/scaleData,1,quantile,probs=.05)
    dataUppQuant <- apply(simplify2array(variableData)/scaleData,1,quantile,probs=.95)
    rangeOfData <- c(min(dataLowQuant),max(dataUppQuant))
    
    plot(dataMean,type="n",ylim=rangeOfData,ylab="",xlab = "time, t",main=variableName,cex.main=titleSize,cex.axis=cexAxisVal,cex.lab=xyLabelsize)
    lines(dataMedQuant,lty=1)
    lines(dataLowQuant,lty=2)
    lines(dataUppQuant,lty=2)
    print(variableName)


    
  } else if (variableDims[1] == nSpecies) {
    # if there are 10 rows then data are species specific
    par(mfrow=c(floor(sqrt(nSpecies)),ceiling(sqrt(nSpecies))))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,4,4,0))
    
    data <- simplify2array(variableData)
    data[is.na(data)] <- 0  # replace nans with zero. This wont happen in future versions. C++ wont produce nans
    # take the mean over all data runs
    dataMean <- apply(data/scaleData,1:2,mean)
    dataMedQuant <- apply(data/scaleData,1:2,quantile,probs=.5)
    dataLowQuant <- apply(data/scaleData,1:2,quantile,probs=.05)
    dataUppQuant <- apply(data/scaleData,1:2,quantile,probs=.95)
    rangeOfData <- c(min(dataLowQuant),max(dataUppQuant))
    
    for (ispec in 1:nSpecies) {
      if(howToScale == "global") {
        # global range already calculated        
      } else if (howToScale == "local") {
        rangeOfData <- c(min(dataLowQuant[ispec,]),max(dataUppQuant[ispec,]))
      } else {
        stop("no such Scaling Exists - sub_plottingSuite.R")
      }
      plot(dataMean[ispec,] ,type="n",ylim=rangeOfData,main=speciesNames[ispec],cex.main=titleSize,cex.axis=cexAxisVal,cex.lab=xyLabelsize)
       lines(dataMedQuant[ispec,],lty=1)
       lines(dataLowQuant[ispec,],lty=2)
       lines(dataUppQuant[ispec,],lty=2)
       if (variableName == "avByr") {
         # add bimass threshold to plot
         lines(c(0,nYrs),c(B0[ispec]*(.2+B0thresholds[ispec]),B0[ispec]*(.2+B0thresholds[ispec])),lty=3,col="red")
       }
       
    }


    
  } else if ((variableDims[1] == nGuilds) & (isGuild) ){
    
    # take the mean over all data runs
    dataMean <- apply(simplify2array(variableData)/scaleData,1:2,mean)
    dataMedQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.5)
    dataLowQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.05)
    dataUppQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.95)
    rangeOfData <- c(min(dataLowQuant),max(dataUppQuant))
    # plots at guild level
    par(mfrow=c(nGuilds,1))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,4,4,0))
    for (iguild in 1:nGuilds) {
      if(howToScale == "global") {
        # global range already calculated        
      } else if (howToScale == "local") {
        rangeOfData <- c(min(dataLowQuant[iguild,]),max(dataUppQuant[iguild,]))
      } else {
        stop("no such Scaling Exists - sub_plottingSuite.R")
      }
      if (as.logical(length(grep("guild",tolower(variableName))))) {
        plotNames <-  guildNames[iguild]
      } else {
        plotNames <-  fleetNames[iguild]
      }
      plot(dataMean[iguild,] ,type="n",ylim=rangeOfData,main=plotNames,cex.main=titleSize,cex.axis=cexAxisVal,cex.lab=xyLabelsize)
      lines(dataMedQuant[iguild,],lty=1)
      lines(dataLowQuant[iguild,],lty=2)
      lines(dataUppQuant[iguild,],lty=2)
      if (variableName == "est_survey_guild_biomass") {
        # add guild biomass threshold
        lines(c(0,nYrs),c(.2*B0Guild[iguild],.2*B0Guild[iguild]),lty=3,col = "red")
      }
    }
    
  } else if (variableDims[1] == nFleets)  {
    
    # take the mean over all data runs
    dataMean <- apply(simplify2array(variableData)/scaleData,1:2,mean)
    dataMedQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.5)
    dataLowQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.05)
    dataUppQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.95)
    rangeOfData <- c(min(dataLowQuant),max(dataUppQuant))
    # plots at guild level
    par(mfrow=c(nFleets,1))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,4,4,0))
    for (ifleet in 1:nFleets) {
      if(howToScale == "global") {
        # global range already calculated        
      } else if (howToScale == "local") {
        rangeOfData <- c(min(dataLowQuant[ifleet,]),max(dataUppQuant[ifleet,]))
      } else {
        stop("no such Scaling Exists - sub_plottingSuite.R")
      }
      plotNames <-  fleetNames[ifleet]

      plot(dataMean[ifleet,] ,type="n",ylim=rangeOfData,main=plotNames,cex.main=titleSize,cex.axis=cexAxisVal,cex.lab=xyLabelsize)
      lines(dataMedQuant[ifleet,],lty=1)
      lines(dataLowQuant[ifleet,],lty=2)
      lines(dataUppQuant[ifleet,],lty=2)

    }
    
  } else if (variableDims[1] == (nFleets*nGuilds)) {
    # plots at guild level for each fleet
    # take the mean over all data runs
    dataMean <- apply(simplify2array(variableData)/scaleData,1:2,mean)
    dataMedQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.5)
    dataLowQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.05)
    dataUppQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.95)
    rangeOfData <- c(min(dataLowQuant),max(dataUppQuant))
    
    par(mfrow=c(nGuilds,nFleets))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,4,4,0))
    
    for (iGuild in 1:nGuilds) {
      for (iFleet in 1:nFleets) {
        index <- iFleet + ((iGuild-1)*nFleets)
        if(howToScale == "global") {
          # global range already calculated        
        } else if (howToScale == "local") {
          rangeOfData <- c(min(dataLowQuant[index,]),max(dataUppQuant[index,]))
        } else {
          stop("no such Scaling Exists - sub_plottingSuite.R")
        }
        plot(dataMean[index,] ,type="n",ylim=rangeOfData,main=guildNames[iGuild],cex.main=titleSize,cex.axis=cexAxisVal,cex.lab=xyLabelsize)
        legend(0,rangeOfData[2],fleetNames[iFleet])
        lines(dataMedQuant[index,],lty=1)
        lines(dataLowQuant[index,],lty=2)
        lines(dataUppQuant[index,],lty=2)
        
      }
      
    }
    
   
  } else if (variableDims[1] == (nFleets*nSpecies)) {
    # 2 plots are produced. One shows each fleet per species.
    # the other aggregated over species. sp1 # fleet1, #fleet2, #fleet3
    index <- seq(0,nSpecies*nFleets-1,nFleets) 
    variableDataSum <- vector(mode="list",length=nRuns)
    for (iRun in 1:nRuns) {
      varSum <- matrix(0,nFleets,nYrs) # preallocate
      for(j in 1:nFleets) {
        varSum[j,] <- colSums(variableData[[iRun]][index+j,]) # sum over species for fleet total
      }
      variableDataSum[[iRun]] <- varSum
    }
    
    # take the mean over all data runs
    dataMean <- apply(simplify2array(variableDataSum)/scaleData,1:2,mean)
    dataMedQuant <- apply(simplify2array(variableDataSum)/scaleData,1:2,quantile,probs=.5)
    dataLowQuant <- apply(simplify2array(variableDataSum)/scaleData,1:2,quantile,probs=.05)
    dataUppQuant <- apply(simplify2array(variableDataSum)/scaleData,1:2,quantile,probs=.95)
    rangeOfData <- c(min(dataLowQuant),max(dataUppQuant))
    # plots at guild level
    par(mfrow=c(nFleets,1))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,4,4,0))
    for (ifleet in 1:nFleets) {
      if(howToScale == "global") {
        # global range already calculated        
      } else if (howToScale == "local") {
        rangeOfData <- c(min(dataLowQuant[ifleet,]),max(dataUppQuant[ifleet,]))
      } else {
        stop("no such Scaling Exists - sub_plottingSuite.R")
      }
      plot(dataMean[ifleet,] ,type="n",ylim=rangeOfData,main=fleetNames[ifleet],cex.main=titleSize,cex.axis=cexAxisVal,cex.lab=xyLabelsize)
      lines(dataMedQuant[ifleet,],lty=1)
      lines(dataLowQuant[ifleet,],lty=2)
      lines(dataUppQuant[ifleet,],lty=2)
    }
  
    
  } else if (variableDims[1] == (nStepsyr*nYrs*nSpecies)) { # M2,F, N
    par(mfrow=c(floor(sqrt(nSpecies)),ceiling(sqrt(nSpecies))))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,4,4,0))
    nSize <- variableDims[2]
    
    # first we need to aggregate over nStepyr to obtain annual time series
    if (as.logical(length(grep("N",variableName)))) {
      operation <- "mean"
    } else {
      operation <- "sum"
    }
    variableDataAgg <- lapply(variableData,aggregate_t_to_yr,operation)
    variableData <- variableDataAgg
    
    # abundance data over size class for each timestep within a year
    #dataMean <- apply(simplify2array(variableData),1:2,mean)
    dataMedQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.5)
    #dataLowQuant <- apply(simplify2array(variableData),1:2,quantile,probs=.05)
    #dataUppQuant <- apply(simplify2array(variableData),1:2,quantile,probs=.95)
    rangeOfData <- c(min(dataMedQuant),max(dataMedQuant))
    for (ispec in 1:nSpecies) {
#      index <- 1 + ((ispec-1)*nStepsyr*nYrs)
      index <- 1 + ((ispec-1)*nYrs)
      
      if(howToScale == "global") {
        # global range already calculated        
      } else if (howToScale == "local") {
      #  rangeOfData <- c(min(dataMedQuant[(index:(ispec*nStepsyr*nYrs)),]),max(dataMedQuant[(index:(ispec*nStepsyr*nYrs)),]))
        rangeOfData <- c(min(dataMedQuant[(index:(ispec*nYrs)),]),max(dataMedQuant[(index:(ispec*nYrs)),]))
      } else {
        stop("no such Scaling Exists - sub_plottingSuite.R")
      }

      for (isize in 1:nSize) {
        if (isize == 1) {
        #  plot(dataMedQuant[(index:(ispec*nStepsyr*nYrs)),isize] ,type="l",ylim=rangeOfData,main=speciesNames[ispec],col=colVec[isize],cex.main=titleSize,cex.axis=cexAxisVal,cex.lab=xyLabelsize)
          plot(dataMedQuant[(index:(ispec*nYrs)),isize] ,type="l",ylim=rangeOfData,main=speciesNames[ispec],col=colVec[isize],cex.main=titleSize,cex.axis=cexAxisVal,cex.lab=xyLabelsize)
        } else {
        #  lines(dataMedQuant[(index:(ispec*nStepsyr*nYrs)),isize] ,type="l",ylim=rangeOfData,col=colVec[isize])
          lines(dataMedQuant[(index:(ispec*nYrs)),isize] ,type="l",ylim=rangeOfData,col=colVec[isize])
        }
      
      #lines(dataMedQuant[ispec,],lty=1)
      #lines(dataLowQuant[ispec,],lty=2)
      #lines(dataUppQuant[ispec,],lty=2)
      }
      legend(0,rangeOfData[2],c("1","2","3","4","5"),lty=c(1,1,1),col=colVec)
    }

  } else if (variableDims[1] == (nYrs*nSpecies)) { #predation_mortality_size 10*53*5
    # data are by size class for each species over time series
    par(mfrow=c(floor(sqrt(nSpecies)),ceiling(sqrt(nSpecies))))
    par(mar=c(2,2,3,2)+0.1)
    par(oma=c(2,4,4,0))
    nSize <- variableDims[2]
    
    dataMedQuant <- apply(simplify2array(variableData)/scaleData,1:2,quantile,probs=.5)
    rangeOfData <- c(min(dataMedQuant),max(dataMedQuant))
    for (ispec in 1:nSpecies) {
      #      index <- 1 + ((ispec-1)*nStepsyr*nYrs)
      index <- 1 + ((ispec-1)*nYrs)
      
      if(howToScale == "global") {
        # global range already calculated        
      } else if (howToScale == "local") {
        #  rangeOfData <- c(min(dataMedQuant[(index:(ispec*nStepsyr*nYrs)),]),max(dataMedQuant[(index:(ispec*nStepsyr*nYrs)),]))
        rangeOfData <- c(min(dataMedQuant[(index:(ispec*nYrs)),]),max(dataMedQuant[(index:(ispec*nYrs)),]))
      } else {
        stop("no such Scaling Exists - sub_plottingSuite.R")
      }
      
      for (isize in 1:nSize) {
        if (isize == 1) {
          #  plot(dataMedQuant[(index:(ispec*nStepsyr*nYrs)),isize] ,type="l",ylim=rangeOfData,main=speciesNames[ispec],col=colVec[isize],cex.main=titleSize,cex.axis=cexAxisVal,cex.lab=xyLabelsize)
          plot(dataMedQuant[(index:(ispec*nYrs)),isize] ,type="l",ylim=rangeOfData,main=speciesNames[ispec],col=colVec[isize],cex.main=titleSize,cex.axis=cexAxisVal,cex.lab=xyLabelsize)
        } else {
          #  lines(dataMedQuant[(index:(ispec*nStepsyr*nYrs)),isize] ,type="l",ylim=rangeOfData,col=colVec[isize])
          lines(dataMedQuant[(index:(ispec*nYrs)),isize] ,type="l",ylim=rangeOfData,col=colVec[isize])
        }
        
        #lines(dataMedQuant[ispec,],lty=1)
        #lines(dataLowQuant[ispec,],lty=2)
        #lines(dataUppQuant[ispec,],lty=2)
      }
      legend(0,rangeOfData[2],c("1","2","3","4","5"),lty=c(1,1,1),col=colVec)
    }
    
    
  } else {
    dev.off()
    stop(paste("You are trying to plot",variableName,"with dims",variableDims[1],".This hasnt been coded for ..."))
  }
  mtext(ylabel, outer=T, side=2, line=2, cex=mtextSizey)
  mtext(variableName, line = 0,outer=T, side=3,cex = mtextSizeTitle)
  dev.off()
  
}

aggregate_t_to_yr <- function(x,operation) {
  y <- matrix(data=NA,dim(x)[1]/dim(x)[2],dim(x)[2])
  #for each sizeclass
  for (isize in 1:dim(x)[2]) {
    y[,isize] <- tapply(x[,isize],(seq_along(1:dim(x)[1])-1) %/% dim(x)[2],operation)
  }
  
  return(y)
}