import geopandas as gpd
import pandas as pd
from shapely.geometry import LineString
import os


##########################################
## convert VVC network to correct input ##
##########################################
data=pd.read_csv('VVC_netwerk.txt',sep='\t')
data['geometry']=''

for j in data.index:
    data.loc[j,'geometry']= LineString([(data.loc[j,'Xa'],data.loc[j,'Ya']),\
            (data.loc[j,'Xb'],data.loc[j,'Yb'])])


kolom_select=['A','B', u'DISTANCE', u'LINK',u'LINKTYPE', u'LT2SOORT',\
		u'URBANISATI',   u'VTGKM_PW',   u'VTGKM_VL',\
	   u'VTGKM_VR',   u'VTGKM_VZ',   u'V_FREEFLOW',
	  u'dist_eucl',   u'geometry']
data_out=data[kolom_select]

#PAE => bevat snelheid (anders niet mee te geven)
kolom_correct=['A','B', u'DISTANCE', u'LINK',u'LINKTYPE', u'LT2SOORT',\
		u'URBANISATI',   u'VTGKM_PW',   u'VTGKM_VL',\
	   u'VTGKM_VZ',   u'VTGKM_vrac',   u'PAE',
	  u'length',   u'geometry']
data_out.columns=kolom_correct



crs='epsg:31370'
geo_df = gpd.GeoDataFrame(data_out, crs=crs,geometry=data_out['geometry'])
geo_df.to_file('rekenfolder\\VVC.shp')


##########################################
##              kromtrekker             ##
##########################################

os.system('rekenfolder\\run.bat')


##########################################
##              Verwerking              ##
##########################################

# verander wat kolomnamen voor de duidelijkheid
data_verwerkt = gpd.read_file('rekenfolder\\outnet.shp')
kolommen_select=[u'LINK',   u'LINKTYPE',   u'LT2SOORT',\
             u'PAE',    u'URBANISATI','length',\
        u'VTGKM_PW',   u'VTGKM_VL',   u'VTGKM_VZ', u'VTGKM_vrac','geometry']
data_verwerkt=data_verwerkt[kolommen_select]
kolommen_nieuw=[u'LINK',   u'LINKTYPE',   u'LT2SOORT',\
             u'V_FREEFLOW',    u'URBANISATI','length',\
        u'VTGKM_PW',   u'VTGKM_VL',   u'VTGKM_VZ', u'VTGKM_vrac','geometry']
data_verwerkt.columns=kolommen_nieuw

crs='epsg:31370'
geo_df = gpd.GeoDataFrame(data_verwerkt, crs=crs,geometry=data_verwerkt['geometry'])
geo_df.to_file('OutputStap1.shp')
