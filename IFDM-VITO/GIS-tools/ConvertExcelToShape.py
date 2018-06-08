
#load packages
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString


#define coordinate system
crs='epsg:4326'

# read input
input=pd.read_csv('Lijnbroninvoer.txt',sep= '\t',skiprows=1,header=None)
input.columns=['lat1','lon1','lat2','lon2', 'width','hight','spread','NOx',' Wind','Time','NO/NOx','type']

#start conversion
input['geometry']=""
for i in input.index:
    geometry=LineString([(input.loc[i,'lon1'],input.loc[i,'lat1']),(input.loc[i,'lon2'],input.loc[i,'lat2'])]) #note: check column headers!
    input.loc[i,'geometry']=geometry
    
#save
geo_df = gpd.GeoDataFrame(input, crs=crs, geometry=input['geometry'])
geo_df.to_file('LineSources.shp')
print 'Saved as Shapefile'