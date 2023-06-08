import numpy as np 
from skimage.measure import regionprops,label
from skimage.morphology import flood
from skimage.feature import peak_local_max
import iris
from my_regions import EDJregion

def constrain_data(cube,pressure,lat,lon):
    """
    Constrains the cube
    Args:
        cube (iris.cube): cube of data
        pressure (array): array of pressure range in ascending order
        lat (array): array of latitude in ascending order
        lon (array): array of longitude in ascending order

    Returns:
        cube : same cube but now constrained
    """
    if pressure == None: 
        cube = cube.extract(iris.Constraint(latitude=lambda x: lat[0]<= x <= lat[1]))
        cube = cube.extract(iris.Constraint(longitude=lambda x: lon[0]<= x <= lon[1]))
    elif lat == None: 
        cube = cube.extract(iris.Constraint(air_pressure=lambda x: pressure[0]<= x <= pressure[1]))
        cube = cube.extract(iris.Constraint(longitude=lambda x: lon[0]<= x <= lon[1]))
    elif lon == None: 
        cube = cube.extract(iris.Constraint(air_pressure=lambda x: pressure[0]<= x <= pressure[1]))
        cube = cube.extract(iris.Constraint(latitude=lambda x: lat[0]<= x <= lat[1]))
    else: 
        cube = cube.extract(iris.Constraint(air_pressure=lambda x: pressure[0]<= x <= pressure[1]))
        cube = cube.extract(iris.Constraint(latitude=lambda x: lat[0]<= x <= lat[1]))
        cube = cube.extract(iris.Constraint(longitude=lambda x: lon[0]<= x <= lon[1]))
    return cube


def blob_finder(wind,Ucrit,gridarea,region_length,zonal_length):
    """
    Args:
        wind (iris.Cube): Iris cube of your zonal wind field 
        Ucrit (int or float): Sets the maxima finding and the flooding criteria
        gridarea (iris.Cube): Area of each grid cell **IMPORTANT** dimensions have to match that of the wind cube
        region_length (int or float): Minimum length in km of each region
        zonal_length (int or float): Minimum zonal extent in longitude of each region 

    Returns:
        3-tuple : (regions found on each day, binary mask of the region, maxima for each region)
    """
    
    wind_copy = wind.copy()
    days = np.arange(0,len(wind.data))
    regions = np.zeros(len(days),dtype=object)
    floods = np.zeros(len(days),dtype=object)
    maximas = np.zeros(len(days),dtype=object)

    lon,lat = wind.coord('longitude').points,wind.coord('latitude').points
    if np.max(lon)!=180: 
        print('converting longitudes')
        lon = (lon + 180) % 360 - 180

    for day in days: 
        # find all maxima 
        wind_day = wind_copy.data[day]
        wind_intense = wind.data[day]

        region_maxima = peak_local_max(wind_day,threshold_abs=Ucrit,exclude_border=False)
        
        region_maxima_value = []
        for m in region_maxima: 
            region_maxima_value.append(wind_day[m[0],m[1]])
        
        #sort into descending order starting with the largest
        try:
            sorted_index = np.argsort(np.array(region_maxima_value)) 
            region_maxima_coords = (np.array(region_maxima)[sorted_index])[::-1]
        except IndexError: 
            print(sorted_index)
       
        #begin flood
        regions_store = []
        flood_store = []
        maximas_store = []


        for maxima in region_maxima_coords: 
        
            xmax,ymax = maxima[0],maxima[1] #maxima used for seed point
            max_vel = wind_day[xmax,ymax] #max velocity used for flooding tolerance
             
            if max_vel == -1.0:
                continue
            else:
                
                tol =  max_vel - Ucrit
                flood_max = flood(wind_day, (xmax,ymax),tolerance=tol,connectivity=2)
                label_flood_max = label(flood_max,connectivity=2)
                region =  regionprops(label_flood_max,intensity_image=wind_intense)
                wind_day[flood_max] = -1 #remove region


                regionEDJ= EDJregion(region[0],flood_max,wind_intense,gridarea,lon,lat)
                length = regionEDJ.length()
                
                if region == []: #empty label check, change this to appending a 0 
                    continue
                if  (regionEDJ.phibar()<np.min(lat) or regionEDJ.lambdabar()< np.min(lon)) or (regionEDJ.phibar()> np.max(lat) or regionEDJ.lambdabar()> np.max(lon)): #weighted centroid out of domain (weird bug)    
                    continue
                if  length <= region_length: #length of region check
                    continue
                if  regionEDJ.zonal_length() <= zonal_length: 
                    continue
                else:
                    flood_store.append(flood_max)
                    regions_store.append(region)
                    maximas_store.append(maxima)
        


                
        maximas[day] = maximas_store
        regions[day] = regions_store
        floods[day] = flood_store

    return regions,floods,maximas
    
