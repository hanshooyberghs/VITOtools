import geopandas as gpd
import pandas as pd
import sys

#############
##  files  ##
#############

file_LEZ_region='LEZ_Antwerpen.shp'
file_emissions_LEZ='EnkelLEZ2020.shp'
file_emissions_NOLEZ='NoLEZ2020.shp'
outputfile='Emissions_processed.shp'

##########################
## read data and checks ##
##########################

# read emissions
emissions_LEZ=gpd.read_file(file_emissions_LEZ)
emissions_NOLEZ=gpd.read_file(file_emissions_NOLEZ)

# check geo
if ~emissions_LEZ.geometry.equals(emissions_NOLEZ.geometry):
    print 'ERROR in input: geometry of both emission files differ'
    sys.exit()


# read polygon
region=gpd.read_file(file_LEZ_region)


##################
##  processing  ##
##################

# inside geometry
intersects=emissions_LEZ.geometry.intersects(region.iloc[0].geometry)
emissions_inside=pd.DataFrame(emissions_LEZ[intersects])

# outside geometry
emissions_outside=pd.DataFrame(emissions_NOLEZ[~intersects])

#combine
emissions=pd.concat([emissions_inside,emissions_outside],ignore_index=True)


################
##   saving   ##
################

gdf = gpd.GeoDataFrame(emissions, crs='epsg:31370', geometry = emissions.geometry)
gdf.to_file(outputfile)







