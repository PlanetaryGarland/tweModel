# tweModel.py

:ringed_planet: 

tweModel.py is a package used to create models of giant planet atmospheres by extrapolating cloud-top temperature and wind measurements to the deep troposphere using the geostrophic-balance form of the thermal wind equation. 

:ringed_planet:

### Dependancies
- python 3.8 (should work on older versions)
- numpy
- scipy (.signal)
- matplotlib

### Example Run

To perform a sample run and verify that everything has been downloaded correctly, run the following example:

1.) Navigate to the [examples/jupiterExample/](examples/jupiterExample/) directory

2.) Run [jupiterExample.py](examples/jupiterExample/jupiterExample.py) in your favorite IDE or with `python jupiterExample.py` from the [examples/jupiterExample/](examples/jupiterExample/) directory. 

If you see INSTALL SUCCESSFUL, everything is working properly. Example .csv outputs and plots will be in the examples/jupiterExample/results/jupiter/runResults/ directory.

### How to Generate a New Model

1.) Place your input cloud-top zonal wind profile and temperatures in the [data/](data/) directory under a \<planet>/ directory. For example, for a Jupiter model you would have two .csv files in the directory data/jupiter/. Both .csv files should be in a two column format: the 1st column containing the latitude in degrees and the 2nd the wind/temperature value at that latitude. 

2.) Create a configuration file named config.ini in the same directory as [tweModel.py](tweModel.py). A blank template ([configTemplate.ini](examples/configTemplate.ini)) and defaults containing some physical parameters of the 4 outer planets may be found in the [examples/](examples/) directory. Each user input is described in the template config.ini. This file will be copied and saved with your outputs.

3.) Run [tweModel.py](tweModel.py) in your favorite IDE or with `python tweModel.py` from the main directory.

## Outputs

Output pressure, latitudes, temperature, density, gravity, and wind values as well as plots will be saved in the [results/](results/) directory under \<planet>/\<date>. 

The temperature, density, and wind outputs are 2D arrays with rows of pressures from p_0 to pFinal and columns of latitudes from phiStart to phiStop. The gravity output is a 1D array from phiStart to phiStop.

### Validation

To ensure the output model is a physical, stable atmosphere run the validation routine with these steps:

1.)  Create a validation configuration file named validationConfig.ini in the same directory as [validation.py](validation.py). A blank template ([validationConfigTemplate.ini](examples/validationConfigTemplate.ini)) may be found in the [examples/](examples/) directory. Each user input is described in the template validationConfig.ini. The code will check the model located in the results/\<planet>/\<date> directory.

2.) Run [validation.py](validation.py) in your favorite IDE or with `python validation.py` from the main directory.

If you see VALIDATE SUCCESS! for both the Brunt-Väisälä and potential temperature validations, the model is a stable atmosphere. 

[tweModel.py](tweModel.py) will also save .csv files containing the square of the Brunt-Väisälä frequency, potential temperature, and change in the potential temperature with pressure in a similar structure to the 2D outputs from [tweModel.py](tweModel.py) described above.

### Documentation

The physics behind each step of model generation is described in [theory.md](docs/theory.md) in the [docs/](docs/) directory.

### Author
Justin Garland

Ph.D. Candidate, Hampton University, Department of Atmospheric and Planetary Sciences

planetary@garland.run
