import geopandas as gpd
import pandas as pd
from shapely.geometry import LineString
import os

file_in='Emissions2.xlsx'

data=pd.read_excel(file_in)
data['geometry']=''

for j in data.index:
    data.loc[j,'geometry']= LineString([(data.loc[j,'Lat1']*1e3,data.loc[j,'Lon1']*1e3),\
            (data.loc[j,'Lat2']*1e3,data.loc[j,'Lon2']*1e3)])
            
crs='epsg:31370'
geo_df = gpd.GeoDataFrame(data, crs=crs,geometry=data['geometry'])
geo_df.to_file(file_in.replace('xlsx','shp'))
