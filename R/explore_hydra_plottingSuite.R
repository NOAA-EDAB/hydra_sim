#' Creates plots of output variables
#' 
#' @param scale Character string. How to scale the subplots. "local" or "global".
#' All subplots on the same global y axis = "global". All sup plots on separate axis = "local. Default ="local"
#' @param plottingLevel Character string. Subset of plots to create, Options are "all" (Default), "species", "catch",
#' "survey","guild",indices"
#' @param wd Character string. Full path to folder where model run data is stored
#' @param dataFolder Character string. Folder in which output files are stored. Default = "diagnostics"
#' 
#' @return Plots are created and saved withing the \code{dataFolder} folder.
#' 


if (!exists("sub_reptoRlist", mode="function")){source(here::here("R","sub_reptoRlist.R"))}
if (!exists("sub_plottingSuite", mode="function")){source(here::here("R","sub_plottingSuite.R"))}

explore_hydra_plottingSuite <- function(scale="local",plottingLevel="all",wd=here::here(),dataFolder="diagnostics"){
  
  if (!dir.exists(wd)) {
   stop("You can not run this script unless you pass a correct argument of wd. wd = \"path to the root of your diagnostic folder\"")
  }
  ################################ READ IN THE DATA FILES ################################
  # path to where the simulated data files are located *.out
  filePath <- paste0(wd,"/",dataFolder,"/")
  
  # list of all of the .out files. 
  fname<-list.files(path=filePath,pattern="\\.out$")
  if (length(fname)== 0) {
    stop(paste0("There are no output files to process. Make sure they have been copied to the folder -",filePath))
  }
  # create a list of names. full path
  f <- list()
  for (ifile in fname) {
    f[ifile] <- paste0(filePath,"/",ifile)
  }
  nRuns <- length(f)
  # process all files using function sub_reptoRlist
  print("Please wait ... processing files")
  print(paste0("in folder ... ",filePath))
  
  A<-lapply(f, sub_reptoRlist) #all in one huge list, no loop

 #########################################################################################
 # pick out some constants + error checking
  nStepsyr <- unique(simplify2array(lapply(A,'[[',"Nstepsyr"))) # number of time steps in a year
  if (length(nStepsyr) != 1) { # check to make sure same over all runs.
    stop("The number of time increments within a year is inconsistent between runs")
  }
  nYrs <- unique(simplify2array(lapply(A,'[[',"Nyrs"))) # number of years
  if (length(nYrs) != 1) { # check to make sure same over all runs.
    stop("The number of years is inconsistent between runs")
  }

  #used for plotting
  speciesNames <- c("spiny_dogfish", "winter_skate", "Aherring", "Acod", "haddock", "ytail_fl", "wint_fl", "Amackerel", "silverhake", "goosefish")
  guildNames  <- c("Piscivores","Planktivores","Benthivores","Elasmobranchs")
  fleetNames <- c("Benthic","Pelagic","Longline","smallMesh","gillnet")

  # selects groups of variables to plot
  variablesToPlot <- plotLevel(plottingLevel)
  nVariables <- length(variablesToPlot)
  # following used to plot thresholds on relevant plots
  B0 <- unlist(A[[1]]["B0"])/1000 # convert to kilotons for plotting
  B0Guild <- unlist(A[[1]]["B0_guilds"])/1000
  threshold_species <- unlist(A[[1]]["threshold_species"])
  
  print("sorting and plotting data ...")
  # create a list of stuff to pass as common to plotting function
  otherList <- list(speciesNames=speciesNames,guildNames=guildNames,fleetNames=fleetNames,
                    scale=scale,nStepsyr=nStepsyr,nYrs=nYrs,filePath=filePath,B0=B0,B0Guild=B0Guild,threshold=threshold_species)


  # loops through all variable to create plots of data
  for (iplot in 1:nVariables) {
    print(variablesToPlot[iplot])
    # creates variable to pass to plotting function
    variableName <- variablesToPlot[iplot]
    variableData <- lapply(A,'[[',variablesToPlot[iplot])
    variableDims <- dim(variableData[[1]])

    print(variableDims)

    if (length(variableData[[1]])==0) {
        print(paste(variablesToPlot[iplot], "do not exist in this output file"))
        next
    } 
    
    sub_plottingSuite(variableData,variableName,variableDims,otherList)
    
  }
  

}

plotLevel <- function(plottingLevel) {
  individualLevel <- c("N","recruitment")
  speciesLevel <- c("recruitment","avByr","SSB","M2","F","Z","TrueF","coastWideF","FishingAndDiscardMortality","obs_survey_biomass","est_survey_biomass",
                    "catch_biomass","obs_catch_biomass","est_catch_biomass","otherDead_Discardbiomass","otherDead_Catchbiomass",
                    "obs_survey_biomass","est_survey_biomass","eaten_biomass","discard_biomass","M1_biomass","otherDead_biomass","total_biomass",
                    "predation_mortality","predation_mortality_size","fishing_mortality","fishing_mortality_size")
  surveyLevel <- c("obs_survey_biomass","est_survey_biomass","est_survey_guild_biomass","est_survey_guild_biomass_assessment")
  catchLevel <- c("catch_biomass","obs_catch_biomass","est_catch_biomass","fleet_catch_biomass","est_catch_guild_biomass",
                  "est_fleet_catch_biomass","est_fleet_catch_guild_biomass")
  guildLevel <-  c("est_catch_guild_biomass","est_survey_guild_biomass","est_fleet_catch_guild_biomass")
  indicesLevel <- c("index_LFI_Biomass","index_LFI_Catch","index_LFI_N","index_Simpsons_Nrecip","index_Simpsons_Crecip",
                    "index_predToPreyRatio","index_plankToPiscRatio","index_stdev_catch","index_stdev_biomass","index_ExploitationRate",
                    "index_SystemExploitationRate")
  dataLevel <- c("obs_temp","obs_effort","obs_effort_total","obs_effort_proportion","obs_effortAssess","exploitation_update")
  allLevel <- unique(c(individualLevel,speciesLevel,guildLevel,surveyLevel,catchLevel,indicesLevel,dataLevel))
  
  switch (plottingLevel,all=allLevel,catch=catchLevel,survey=surveyLevel,species=speciesLevel,
          guild=guildLevel,individual=individualLevel,indices=indicesLevel,data=dataLevel)
  

}