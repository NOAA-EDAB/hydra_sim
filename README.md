# hydra_sim

[![gitleaks](https://github.com/NOAA-EDAB/hydra_sim/actions/workflows/secretScan.yml/badge.svg)](https://github.com/NOAA-EDAB/hydra_sim/actions/workflows/secretScan.yml)

Code development for the hydra multispecies model

- Code used in: [Gaichas et al. (2017)](https://github.com/NOAA-EDAB/hydra_sim/releases/tag/0.1.0). Combining stock, multispecies, and ecosystem level fishery objectives within an operational management procedure: simulation to start the conversation. ICES 74:2, 552-565. doi: [10.1093/icesjms/fsw119](https://doi.org/10.1093/icesjms/fsw119)

Hydra is a length-based multispecies, multifleet, spatial model designed to provide simulated data for performance 
testing of simpler (non-size structured) multispecies assessment models and managemnet procedures for the Northeast US Continental Shelf.

For more information please visit the [wiki](https://github.com/NOAA-EDAB/hydra_sim/wiki). Assumptions made in the model are detailed there.

## Usage

The Hydra model is coded in Automatic Differentiation Model Builder (ADMB). You will either need to install ADMB locally or use Docker/podman and install ADMB inside the container. A Dockerfile is provided to achieve this.

### On PC

* Clone the repo as an R project
* Move the *.tpl file to its own folder and compile it
* Install dependencies (`remotes::install_deps(".")`)
* Source the `run_model.r` file
* Run the `run_model.r` file

```r
run_model(pathToTPL = "full path to where you compiled the model",rootFolder="name of folder to store model output")
```
For example, if you move the tpl file into a folder called `ADMB` and want output saved into a folder called `out`

```r
run_model(pathToTPL = here::here("ADMB"), rootFolder = "out")
```

A sample set of parameter files `hydra_sim.dat` and `hydra_sim.pin` are included and should remain in the projects root after cloning. They will be copied to the `rootFolder`.
The parameterizarion reflects what is termed a "historic" run where (somewhat) realistic fishing effort is used to drive the model.

## Results

The model runs 100 times. The output files (*.out and *.txt) are temporarily stored in the projects root folder. Once all runs have completed they will be moved to the `rootFolder` that you specified in the `run_model`call. These files will then be processed and a suite of plots will be made. They will be saved in the "diagnostics" folder 

### Container

A Dockerfile is provided to run Hydra inside a container.

* Clone the repo
* Build the Docker image (here it is named `hydra`)

  `docker build -t hydra .`

* Create a directory in which to write model output, for example `output`
* Copy `hydra_sim.dat`,`hydra_sim.pin` to this folder
* Run Hydra 

  `docker run --rm --mount "type=bind,src=/path_to_output_folder/,dst=/HYDRA/mount" hydra`

The output from the model should appear in the `output` folder. Running Hydra in this manner uses a random number seed of 1, always. This is a good first step to ensure your image built correctly. The preferred method to run Hydra is to pass a value of the seed as an argument to the docker run command: 

* `docker run --rm --mount "type=bind,src=/path_to_output_folder/,dst=/HYDRA/mount" hydra hydra_sim.dat hydra_sim.pin 88`

This runs Hydra using dat and pin files with specified names (these must reside in the `path_to_output_folder`) and a seed of 88

* Process output

## Contact

| [andybeet](https://github.com/andybeet)        
| ----------------------------------------------------------------------------------------------- 
| [![](https://avatars1.githubusercontent.com/u/22455149?s=100&v=4)](https://github.com/andybeet) | 



#### Legal disclaimer

*This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.*
