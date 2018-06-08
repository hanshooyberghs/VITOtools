import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

data_in=pd.read_csv('CompleteGrid_GIS.txt',sep=' ',header=None)
data_in=data_in[[3,10]]
data_in.columns=[['x','y']]
data_in['geo']=''

for i in range(len(data_in)):
    data_in.loc[i,'geo']=Point(data_in.loc[i,'x'],data_in.loc[i,'y'])

    
crs='epsg:31370'
geo_df = gpd.GeoDataFrame(data_in, crs=crs, geometry=data_in['geo'])
geo_df.to_file('Grid.shp')