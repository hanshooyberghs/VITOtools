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
if ~emissions_LEZ.geom_equals(emissions_NOLEZ).all():
    print 'ERROR in input: geometry of both emission files differ'
    sys.exit()

# process geo data of emissions
emissions_geo=emissions_LEZ.geometry
emissions_LEZ.drop('geometry',inplace=True,axis=1)
emissions_NOLEZ.drop('geometry',inplace=True,axis=1)

# read polygon
region=gpd.read_file(file_LEZ_region)


##################
##  processing  ##
##################

# inside geometry
intersects=emissions_geo.geometry.intersection(region.iloc[0].geometry)
intersects2=pd.DataFrame(intersects[intersects.length>0],columns=['geometry'])

# outside geometry
difference=emissions_geo.geometry.difference(region.iloc[0].geometry)
difference2=pd.DataFrame(difference[difference.length>0],columns=['geometry'])

#inside emissions
emissions_inside=intersects2.merge(emissions_LEZ,left_index=True,right_index=True,how='left')

#outside emissions
emissions_outside=difference2.merge(emissions_NOLEZ,left_index=True,right_index=True,how='left')

#combine
emissions=pd.concat([emissions_inside,emissions_outside],ignore_index=True)

# recalculate X and Y coordinates
for i in range(len(emissions)):
    coord_segment=emissions.loc[i,'geometry'].coords[:]
    emissions.loc[i,'X1']=coord_segment[0][0]
    emissions.loc[i,'Y1']=coord_segment[0][1]
    emissions.loc[i,'X2']=coord_segment[1][0]
    emissions.loc[i,'Y2']=coord_segment[1][1]

################
##   saving   ##
################

gdf = gpd.GeoDataFrame(emissions, crs='epsg:31370', geometry = emissions.geometry)
gdf.to_file(outputfile)







