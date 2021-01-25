#' Run hydra and plot some output
#' 
#' Run Hydra many times, move output into project subfolder and create plots
#' 
#' @param pathToTPL Character string. Path to where the compiled code resides
#' @param rootFolder Character string. The name of the folder to be created to contain all of the output the name of th 
#' 
#' @return model run files and plots of some key variables
#' 
#' @section Note:
#' 
#' You must keep the hydra_sim.dat and .pin files in the root of the directory. they will be copied into \code{rootFolder}
#' 
#' Two folders will be created inside of \code{rootFolder}. A "diagnostics" and an "indices" folder 
#' Model run output will be saved in each. The diagnostics folder contains full output. 
#' All plots will be saved within this folder.
#' 
#' @section Dependencies: 
#' 
#' You will need to install following packages
#' 
#' here


source(here::here("R","explore_hydra_plottingSuite.R"))


run_model <- function(pathToTPL="C:/Users/Andrew.Beet/Documents/MyWork/Hydra/Beet-Hydra/ADMB",rootFolder) {

  hydraVersion <- "hydra_sim"
  datpinPath <- here::here(rootFolder)
  
    
  if(!dir.exists(here::here(rootFolder))) {
    dir.create(here::here(rootFolder))
    dir.create(here::here(rootFolder,"indices"))
    dir.create(here::here(rootFolder,"diagnostics"))
  }

  a <- file.copy(from = here::here("hydra_sim.pin"),to = here::here(rootFolder,"hydra_sim.pin"))
  a <- file.copy(from = here::here("hydra_sim.dat"),to = here::here(rootFolder,"hydra_sim.dat"))
  
  message("running Hydra")
  for(i in 1:100) {
    isamp <- sample(1E6,1)
    message(paste0("run ",i," of 100"))
    hydramse::run_single_hydra(isamp,exePath = paste0(pathToTPL,"/",hydraVersion,".exe"),pinPath = paste0(datpinPath,"/",hydraVersion,".pin"),datPath = paste0(datpinPath,"/",hydraVersion,".dat"))
  }
  
  message("copy output files from root")
  if (Sys.info()['sysname']=="Windows") {
    filesToMove <- list.files(here::here(),pattern="\\.txt$")
    a <- file.rename(from = here::here(filesToMove),to = here::here(rootFolder,"indices",filesToMove))
    filesToMove <- list.files(here::here(),pattern="\\.out$")
    a <- file.rename(from = here::here(filesToMove),to = here::here(rootFolder,"diagnostics",filesToMove))
  } else {
    message("Not coded copying output files fr Linux or Mac os")
    
  }
  
  message("creating historic plot")
  
  explore_hydra_plottingSuite(wd=here::here(rootFolder),dataFolder="diagnostics")
  
}
