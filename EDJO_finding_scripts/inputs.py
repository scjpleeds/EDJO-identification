### INPUTS FOR EDJO FINDING SCRIPT ### 

PATH = ' ' # path for where to save the files 
output_fname ='ERA5_19592020_diagnostics_Ucrit8' # two output files are made, <output_fname>_full.npz contains all data with no distinction between regions 
                                                 # <output_fname>_lm.npz contains the data of the regions with the largest umass on days with multiple objects

var_cube_fname =  'U.nc' # Zonal wind cube file path 
grid_cube_fname = 'gridarea.nc' # grid area cube file path 


cube_constraints = [(-60,0),(15,75),(850,850)] # Cube constraints, if already done with input cube enter as None 

filtering = True # Apply a low pass Lanczos filter to the wind data 
window, length = 61,10


### Constraints for the EDJO finding algorithm 
flood_val = 8
min_length = 1661
min_zonal_length = 20



                  
