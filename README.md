# climate_state_thesis
repository containing files pertinent to Brandon Smith's M.Sc. Thesis from the University of Victoria, and to the end of creating a publication.

This repository contains python scripts that call upon CESM model output data stored on ([FRDR](https://doi.org/10.20383/103.0652)). If you wish to utilize the scripts in this repository to re-create the figures contained in the thesis, please download the model output from their database and adjust the file paths called in each script to where you store the copy of said model output. For example, from the file "map_TS_response.py" line 27:

<inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME #location of parent data set>

"inpath" should be rewritten to be wherever you have decided to store the model output on your machine. 

Please make the necessary adjustments to directory calls for where you wish to store your temporary output files and figures in each script. Scripts should run as-is, with the exception of the file paths. 

These figures were created using model output from the Community Earth System Model (CESM) version 1.2.2. , configured with a slab ocean model and pre-industrial solar insolation and CO2 concentrations. The component set used for these experiments was E_1850_CAM5. Solar insolation and CO2 concentrations were altered for each background climate simulation to maintain near-equivalent global mean surface temperatures (GMST). The available model output contains the 30 year climatology runs from the 60 year simulations that were created for each control, 2xCO2 and 4xCO2 simulation (the prior 30 years of each simulation being treated as spin-up runs), along with their respective restart files (the restart files for the 4xCO2 run at 90% solar luminosity were regrettably lost). If one wishes to use these runs as a starting point for future CESM simulations of the same configuration, only the restart files will be needed. 

