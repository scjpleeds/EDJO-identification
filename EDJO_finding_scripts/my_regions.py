import numpy as np 
from scipy.interpolate import interp1d
from shapely.geometry import Polygon,Point
from geopy.distance import geodesic
from scipy.ndimage import binary_fill_holes
from skimage.measure import find_contours

def pix_to_lat(pixel,lat):

    dlat,latmax = abs(lat[0]-lat[-1])/len(lat),np.max(abs(lat))
    lat = latmax-dlat*pixel    
    return lat

def pix_to_lon(pixel,lon):

    dlon,lonmax = abs(lon[0]-lon[-1])/len(lon),np.max(abs(lon))
    lon = dlon*pixel-lonmax

    return lon

class EDJregion: 
    def __init__(self,blob,flood,field,area,lon,lat):
        global R; 
        R = 6371e3
        self.blob = blob
        self.field = field
        self.area = area
        self.lon = lon
        self.lat = lat
        self.flood = flood

    def region_area(self): 
        A=self.area 
        BLOB_COORDS=self.blob.coords
        return np.sum(np.array([A[y,x] for y,x in BLOB_COORDS]))
     
    def mass(self): 
        F=self.field
        A=self.area 
        BLOB_COORDS=self.blob.coords
        return np.sum([F[y,x]*A[y,x] for y,x in BLOB_COORDS])
    
    def mean_intensity(self):

        return self.mass()/self.region_area()

    def phibar(self): 
        F = self.field
        A =  self.area
        LAT = self.lat
        BLOB_COORDS = self.blob.coords
        MASS = self.mass()

        mu01= np.sum([F[y,x]*A[y,x]*LAT[y] for y,x in BLOB_COORDS])
        return mu01/MASS

    def lambdabar(self): 
        F = self.field
        A =  self.area
        LON = self.lon
        BLOB_COORDS = self.blob.coords
        MASS = self.mass()
        mu10 = np.sum([F[y,x]*A[y,x]*LON[x] for y,x in BLOB_COORDS])
        
        return mu10/MASS

    def alpha(self): 
    # Calculation here is taken from skimage.regionprops 
        a,b,b,c = self.get_inertia_tensor().flat

        if a-c==0:
            if b<0:
                return -np.pi/4
            else:
                return np.pi/4
        else: 
            return np.degrees(0.5*np.arctan2(-2*b,a-c))%180-90

    def get_inertia_tensor(self):
        F = self.field
        A = self.area
        BLOB_COORDS = self.blob.coords
        LAT,LON = R*np.radians(self.lat),R*np.radians(self.lon)
        PHIBAR,LAMBDABAR = R*np.radians(self.phibar()),R*np.radians(self.lambdabar())

        mu11 = np.sum([F[y,x]*A[y,x]*(LAT[y]-PHIBAR)*(LON[x]-LAMBDABAR) for y,x in BLOB_COORDS])
        mu20 = np.sum([F[y,x]*A[y,x]*(LAT[y]-PHIBAR)**2 for y,x in BLOB_COORDS])
        mu02 = np.sum([F[y,x]*A[y,x]*(LON[x]-LAMBDABAR)**2 for y,x in BLOB_COORDS]) 

        return np.array([[mu20,mu11],[mu11,mu02]])
    
    def get_axis_length(self):
        return (3/1e3)*np.sqrt(np.max(np.linalg.eigvals(self.get_inertia_tensor()))/self.mass())
    
    def get_axis_width(self):
        return (3/1e3)*np.sqrt(np.min(np.linalg.eigvals(self.get_inertia_tensor()))/self.mass())

    def length(self):
        #Use Rtree or SRTree to speed this up ? - https://stackoverflow.com/questions/14697442/faster-way-of-polygon-intersection-with-shapely
        
        BLOB = self.blob
        LON = self.lon
        LAT = self.lat
        FLOOD = self.flood

        _,min_col,_,max_col = BLOB.bbox
        points = np.linspace(min_col-100,max_col+100,1000) #+- 100 is done to ensure that the line extrapolation goes outside of the region, guaranteeing that it cuts the edges of the blob
        y0, x0 = BLOB.centroid
        x0,y0=x0+1,y0+1 #shift centre by 1 to account for padding
        orientation = BLOB.orientation
        x2 = x0 - np.sin(orientation) * 0.5 * BLOB.major_axis_length 
        y2 = y0 - np.cos(orientation) * 0.5 * BLOB.major_axis_length 

        f = interp1d((x0,x2),(y0,y2),fill_value='extrapolate')  # extrapolate a line through the major axis 
        line_points = f(points)
        
        xy_points = np.array(list(zip(points,line_points)))

        flood = np.pad(FLOOD,1,mode='constant',constant_values=0)

        co = np.concatenate(find_contours(binary_fill_holes(flood),level=0)[:])
        co = [(t[1],t[0]) for t in co]
        poly  = Polygon(co)
        
        points_inside = []
        split = np.split(xy_points,2)
        xy_top=split[0]
        xy_bottom=split[1]
        
        split_line = np.r_[xy_top,xy_bottom]

        for p in split_line: 
            p = Point(p)
            if p.intersects(poly):
                points_inside.append((p.x,p.y))

        if points_inside == []: 
            return 0.0
        
        points_inside.sort(key=lambda tup: tup[0])

        yp0,xp0,yp1,xp1 = points_inside[0][1],points_inside[0][0],points_inside[-1][1],points_inside[-1][0]
        lon0,lat0 = pix_to_lon(xp0,LON),pix_to_lat(yp0,LAT)
        lon1,lat1 = pix_to_lon(xp1,LON),pix_to_lat(yp1,LAT)
        length = geodesic((lat0,lon0),(lat1,lon1)).km

        return length
            
    def zonal_length(self):
        BLOB_BOX = self.blob.bbox
        LON = self.lon

        _,min_row,_,max_row = BLOB_BOX

        return abs(abs(LON[min_row])-abs(LON[max_row-1]))

    def plot_region(): 

        return 0
