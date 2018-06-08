import geopandas as gpd
import pandas as pd
from shapely.geometry import LineString
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import mplleaflet


tunnels=pd.read_excel('TunnelsBridges.xlsx',sheet='Tunnels')
tunnels['geometry']=''
tunnels.geometry=tunnels.apply(lambda p: LineString([(p.XA,p.YA),(p.XB,p.YB)]),axis=1)
tunnels_gpd=gpd.GeoDataFrame(tunnels,geometry=tunnels.geometry,crs={'init':'epsg:31370'})
tunnels_gpd.to_file('Tunnels.shp')


colorslist=['green','red','blue']
CustomCmap = ListedColormap(colorslist)
tunnels_gpd.plot(column='Type',categorical='True',colormap=CustomCmap,legend=True,linewidth=5)
mplleaflet.show(path='tunnels.html', crs={'init':'epsg:31370'})