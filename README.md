# EDJO-identification
Python code for the Eddy-Driven Jet Object (EDJO) identification methodology described in the paper:

A new characterization of the North Atlantic eddy-driven jet using 2-dimensional moment analysis. 

## Folder contents 
In the folder EDJO_finding_scripts, you will find 5 files: 

1. edjo_diagnostic_finding_script.py 
2. functions.py
3. inputs.py
4. my_regions.py
5. jet_env.yml

If you only want to apply the method then you need to install the conda environment using the jet_env.yml, then update the inputs.py with the relevant files (wind and grid area data) and processing (low pass filtering, domain slicing, choice of region size, choice of Ucrit value). Once these are updated you can run the edjo_diagnostic_finding_script.py which will save two files, one containing all of the EDJOs found and another only containg the EDJOs with the largest mass. 

The other files contain the code that finds the objects (find_blobs.py) and the EDJO class (my_regions.py) which is used to calculate the various diagnostics. 
