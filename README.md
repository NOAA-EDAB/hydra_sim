# hydra_sim
Code development for the hydra multispecies model

- Code used in: [Gaichas et al. (2017)](https://github.com/NOAA-EDAB/hydra_sim/releases/tag/0.1.0). Combining stock, multispecies, and ecosystem level fishery objectives within an operational management procedure: simulation to start the conversation. ICES 74:2, 552-565. doi: [10.1093/icesjms/fsw119](https://doi.org/10.1093/icesjms/fsw119)

Hydra is a length-based multispecies, multifleet, spatial model designed to provide simulated data for performance 
testing of simpler (non-size structured) multispecies assessment models and managemnet procedures for the Northeast US Continental Shelf.

For more information please visit the [wiki](https://github.com/NOAA-EDAB/hydra_sim/wiki). Assumptions made in the model are detailed there.

## Usage

* Clone the repo as an R project
* Move the *.tpl file to its own folder and compile it
* Source the `run_model.r` file
* Run the `run_model.r` file

```r
run_model(pathToTPL = "full path to where you compiled the model",rootFolder="name of folder to store model output")
```

A sample set of parameter files `hydra_sim.dat` and `hydra_sim.pin` are included and should remain in the projects root after cloning. They will be copies to the `rootFolder`.
The parameterizarion reflects what is termed a "historic" run where (somewhat) realistic fishing effort is used to drive the model.


## Results

The model runs 100 times. The output files (*.out and *.txt) are temporarily stored in the projects root folder. Once all runs have completed they will be moved to the `rootFolder` that you specified in the `run_model`call. These files will then be processed and a suite of plots will be made. They will be saved in the "diagnostics" folder  



*This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.*
