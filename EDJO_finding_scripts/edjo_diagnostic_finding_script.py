import numpy as np 
import iris
from functions import blob_finder,constrain_data,low_pass_filter
from my_regions import EDJregion
import warnings
# from inputs import * 
import tqdm as tqdm
warnings.filterwarnings('ignore')


PATH = '/Users/admin/OneDrive - University of Leeds/ARC4/scripts/' # path for where to save the files 
# two output files are made, <output_fname>_full.npz contains all data with no distinction between regions 
                  # <output_fname>_lm.npz contains the data of the regions with the largest umass on days with multiple objects

var_cube_fname =  '/Users/admin/OneDrive - University of Leeds/ARC4/ERA/U19592020/U.nc' # Zonal wind cube file path 
grid_cube_fname = '/Users/admin/OneDrive - University of Leeds/ARC4/scripts/gridarea.nc' # grid area cube file path 


cube_constraints = [(-60,0),(15,75),(850,850)] # Cube constraints, if already done with input cube enter as None 

filtering = True # Apply a low pass Lanczos filter to the wind data 
window, length = 61,10


# Importing data
var_cube =  iris.load_cube(var_cube_fname)
grid_cube = iris.load_cube(grid_cube_fname)


# Applying constraints to the data if provided 
if cube_constraints!=None:
    lon_constraint,lat_constraint,pressure_constraint=cube_constraints
    grid_cube = constrain_data(grid_cube,None,[lat_constraint[0],lat_constraint[1]],[lon_constraint[0],lon_constraint[1]])
    var_cube = constrain_data(var_cube,[pressure_constraint[0],pressure_constraint[1]],[lat_constraint[0],lat_constraint[1]],[lon_constraint[0],lon_constraint[1]])

# Applying low-pass Lanczos filter if True
if filtering == True: 
    var_cube = low_pass_filter(var_cube,window,length)

#collapse pressure coordinate if more than one level is in the cube
if np.diff(pressure_constraint)!=0: 
    var_cube.collapsed('air_pressure',iris.analysis.MEAN)


### Constraints for the EDJO finding algorithm 

min_length = 1661
min_zonal_length = 20
flood_val = 8
# Running EDJO finding algorithm 
    



output_fname ='ERA5_19592020_diagnostics_filtered' 
regions_store,flood_store,region_maxima_coords= blob_finder(var_cube,flood_val,grid_cube,min_length,min_zonal_length)


#----------- Extracting Diagnotics ----------- # 

# Number of objects defined on each day
num_of_labels = []
N = len(var_cube.data)
for i in range(0,N):
    num_of_labels.append(len(regions_store[i]))

lon,lat = var_cube.coord('longitude').points,var_cube.coord('latitude').points


largest_mass_regions = np.zeros(N,dtype=object)
largest_area_regions = np.zeros(N,dtype=object)
largest_mass_flood = np.zeros(N,dtype=object)
largest_mass_orientation = np.zeros(N)
largest_mass_maxima = np.zeros(N,dtype=object)

max_index_store = []
com_y_full,com_x_full = [],[]
orientation_values_full = []
mean_intensity_full = [] 
major_axis_length_full = []
minor_axis_length_full = []
area_full = []
mass_full = []

for i in np.arange(0,N):
    mass_store = []
    area_store = []
    com_x_store,com_y_store=[],[]
    orientation_values_store = []
    mean_intensity_store = []
    major_axis_length_store,minor_axis_length_store = [],[]

    if regions_store[i] == []:
        com_y_full.append([0])
        com_x_full.append([0])
        mass_full.append([0])
        orientation_values_full.append([0])
        area_full.append([0])
        mean_intensity_full.append([0])
        major_axis_length_full.append([0])
        minor_axis_length_full.append([0])
    else:
        for k,regions in enumerate(regions_store[i]):
            region_edj = EDJregion(regions[0],flood_store[i][k],var_cube.data[i],grid_cube.data,lon,lat)
            com_y_store.append(region_edj.phibar())
            com_x_store.append(region_edj.lambdabar())
            mean_intensity_store.append(region_edj.mean_intensity())
            mass_store.append(region_edj.mass())
            area_store.append(region_edj.region_area())
            major_axis_length_store.append(region_edj.get_axis_length())
            minor_axis_length_store.append(region_edj.get_axis_width())
            orientation_values_store.append(region_edj.alpha())

        com_y_full.append(com_y_store)
        com_x_full.append(com_x_store)
        orientation_values_full.append(orientation_values_store)
        mean_intensity_full.append(mean_intensity_store)
        area_full.append(area_store)
        mass_full.append(mass_store)
        max_mass_index = np.argmax(mass_store)
        major_axis_length_full.append(major_axis_length_store)
        minor_axis_length_full.append(minor_axis_length_store)
        max_index_store.append(max_mass_index)
    
        largest_mass_regions[i] = regions_store[i][max_mass_index]
        largest_mass_flood[i] = flood_store[i][max_mass_index]


com_y_full = np.array(com_y_full,dtype=object)
com_x_full = np.array(com_x_full,dtype=object)     
orientation_values_full = np.array(orientation_values_full,dtype=object)
mean_intensity_full = np.array(mean_intensity_full,dtype=object)
area_full = np.array(area_full,dtype=object)
mass_full = np.array(mass_full,dtype=object)
major_axis_length_full = np.array(major_axis_length_full,dtype=object)
minor_axis_length_full = np.array(minor_axis_length_full,dtype=object)



# Extract out regions based on the largest mass object

mass_values= []
area_values = []
orientation_values = []
com_y = []
com_x = []
mean_intensity = []
major_axis_length_values,minor_axis_length_values = [],[]


for i,region in enumerate(largest_mass_regions):
    if region == 0:
        mass_values.append(0) 
        area_values.append(0)
        orientation_values.append(0)
        com_y.append(0)
        com_x.append(0)
        mean_intensity.append(0)
        major_axis_length_values.append(0)
        minor_axis_length_values.append(0)
    else:
    
        region_edj = EDJregion(region[0],largest_mass_flood[i],var_cube.data[i],grid_cube.data,lon,lat)
        mass_values.append(region_edj.mass())
    
        com_y.append(region_edj.phibar())
        com_x.append(region_edj.lambdabar())
        orientation_values.append(region_edj.alpha())
    
        area_values.append(region_edj.region_area())
        mean_intensity.append(region_edj.mean_intensity())
        major_axis_length_values.append(region_edj.get_axis_length())
        minor_axis_length_values.append(region_edj.get_axis_width())
    

mass_values = np.array(mass_values,dtype=object)
area_values = np.array(area_values,dtype=object)
orientation_values = np.array(orientation_values,dtype=object)
com_y = np.array(com_y,dtype=object)
com_x = np.array(com_x,dtype=object)
mean_intensity = np.array(mean_intensity,dtype=object)
major_axis_length = np.array(major_axis_length_values,dtype=object)
minor_axis_length = np.array(minor_axis_length_values,dtype=object)

np.savez(output_fname+'_lm.npz',phibar=com_y,lambdabar=com_x,alpha=orientation_values,umean=mean_intensity,umass=mass_values,area=area_values,major_axis_length=major_axis_length,minor_axis_length=minor_axis_length,labels=num_of_labels,regions=largest_mass_regions,flood=largest_mass_flood)
np.savez(output_fname+'_full.npz',phibar=com_y_full,lambdabar=com_x_full,alpha=orientation_values_full,umean=mean_intensity_full,area=area_full,umass=mass_full,major_axis_length=major_axis_length_full,minor_axis_length=minor_axis_length_full,labels=num_of_labels,regions=regions_store,flood=flood_store)
