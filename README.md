# climate_state_thesis
repository containing files pertinent to Brandon Smith's M.Sc. Thesis from the University of Victoria, and to the end of creating a publication.

This repository contains python scripts that call upon CESM model output data stored on FDRD (insert link). If you wish to utilize the scripts in this repository to re-create the figures contained in the thesis, please download the model output from their database and adjust the file paths called in each script to where you store the copy of said model output. For example, from the file "map_TS_response.py" line 27:

<inpath = '/home/brandonsmith/modeloutput/'+run+'/'+CASENAME #location of parent data set>

inpath should be rewritten to be wherever you have decided to store the model output on your machine. 

Please make the necessary adjustments to directory calls for where you wish to store your temporary output files and figures in each script. Scripts should run as-is, with the exception of the file paths. 

The following lines contain instructions for repeating the experiment using the Community Earth System Model (v.1.2.2):

