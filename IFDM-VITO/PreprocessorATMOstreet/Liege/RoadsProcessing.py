import geopandas as gpd
import pandas as pd
import os
import numpy as np
from shapely.geometry import LineString


file_in='PICC_voirie_axe_sz.shp'
file_buildings='PICC_constr_batiemprise_sz.shp'


############### 
## read data ##
###############

data_buildings=gpd.read_file(file_buildings)
roads_orig=gpd.read_file(file_in)
crs={'init' :'epsg:31370'}
if not os.path.exists('tmp'):
            os.makedirs('tmp')
if not os.path.exists('output'):
            os.makedirs('output')
data_buildings=data_buildings[['geometry','NATUR_CODE']]
            
########################## 
## first iteration: 10m ##
##########################
print 'simplification 10m'
# simplification
simplified_10=gpd.GeoDataFrame(roads_orig.simplify(10))
simplified_10=simplified_10.join(roads_orig.drop('geometry',axis=1))
simplified_10=gpd.GeoDataFrame(simplified_10,geometry=simplified_10[0],crs=crs)

# Line Polygon intersection
joined=gpd.sjoin(simplified_10,data_buildings,op='intersects',lsuffix='',rsuffix='right')
list_intersection=joined.OBJECTID.unique()

data_ok10m=simplified_10[~simplified_10.OBJECTID.isin(list_intersection)]  # already OK
data_tofurthersimplify_10=roads_orig[simplified_10.OBJECTID.isin(list_intersection)] # still to simplify
data_ok10m[['OBJECTID','geometry']].to_file('tmp/10m_ok.shp')
data_tofurthersimplify_10[['OBJECTID','geometry']].to_file('tmp/10m_nok.shp')


########################## 
## second iteration: 3m ##
##########################
print 'simplification 3m'

# simplification
simplified_3=gpd.GeoDataFrame(data_tofurthersimplify_10.simplify(3))
simplified_3=simplified_3.join(data_tofurthersimplify_10.drop('geometry',axis=1))
simplified_3=gpd.GeoDataFrame(simplified_3,geometry=simplified_3[0],crs=crs)

# Line Polygon intersection
joined=gpd.sjoin(simplified_3,data_buildings,op='intersects')
list_intersection=joined.OBJECTID.unique()

data_ok3m=simplified_3[~simplified_3.OBJECTID.isin(list_intersection)]  # already OK
data_tofurthersimplify_3=data_tofurthersimplify_10[simplified_3.OBJECTID.isin(list_intersection)] # still to simplify
data_ok3m[['OBJECTID','geometry']].to_file('tmp/3m_ok.shp')
data_tofurthersimplify_3[['OBJECTID','geometry']].to_file('tmp/3m_nok.shp')



########################## 
## Third iteration: 1m ##
##########################
print 'simplification 1m'

# simplifications
simplified_1=gpd.GeoDataFrame(data_tofurthersimplify_3.simplify(1))
simplified_1=simplified_1.join(data_tofurthersimplify_3.drop('geometry',axis=1))
simplified_1=gpd.GeoDataFrame(simplified_1,geometry=simplified_1[0],crs=crs)

# Line Polygon intersection
joined=gpd.sjoin(simplified_1,data_buildings,op='intersects')
list_intersection=joined.OBJECTID.unique()

data_ok1m=simplified_1[~simplified_1.OBJECTID.isin(list_intersection)]  # already OK
data_nok=data_tofurthersimplify_3[simplified_1.OBJECTID.isin(list_intersection)] # still to simplify
data_ok1m[['OBJECTID','geometry']].to_file('tmp/1m_ok.shp')
data_nok[['OBJECTID','geometry']].to_file('tmp/1m_nok.shp')


########################## 
##    Merge all data    ##
##########################

data_ok10m['simplification']='simpl_10m'
data_ok3m['simplification']='simpl_3m'
data_ok1m['simplification']='simpl_1m'
data_nok['simplification']='original'
entire_network=pd.concat([data_ok10m,data_ok3m,data_ok1m,data_nok],ignore_index=True)
entire_network.drop([0],axis=1,inplace=True)

columns=['OBJECTID','RUE_ID1','RUE_ID2','RUE_NOM1','RUE_NOM2','simplification','geometry']
entire_network=entire_network[columns]

geo_df = gpd.GeoDataFrame(entire_network, crs=crs)
geo_df.to_file('output\\finalresult.shp')

################################# 
##    Split in single lines    ##
#################################
dfs=[]
kolommen=entire_network.columns.tolist()
for i in entire_network.index:
    if i%1000 == 0:
        print i
    if entire_network.loc[i,'geometry'].type == 'MultiLineString':
        coord_segment=entire_network.loc[i,'geometry'].geoms[0].coords[:]
        #for j in entire_network.loc[i,'geometry'].geoms:
        #    coord_segment.extend(j.coords[:])        
        #    print j.coords[:]
    else:
        coord_segment=entire_network.loc[i,'geometry'].coords[:]
    lengte_segment=len(coord_segment)
    splitted=pd.DataFrame(columns=['X1','Y1','X2','Y2']+kolommen,\
        index=range(lengte_segment-1))
    for j in range(lengte_segment-1):
        splitted.loc[j,kolommen]=entire_network.loc[i,kolommen]
        splitted.loc[j,'X1']=coord_segment[j][0]
        splitted.loc[j,'Y1']=coord_segment[j][1]
        splitted.loc[j,'X2']=coord_segment[j+1][0]
        splitted.loc[j,'Y2']=coord_segment[j+1][1]
    dfs.append(splitted)   
    
entire_network_splitted=pd.concat(dfs,ignore_index=True)
entire_network_splitted.drop('geometry',axis=1,inplace=True)
entire_network_splitted['geometry']=entire_network_splitted.apply(lambda p: 
    LineString([(p['X1'],p['Y1']),(p['X2'],p['Y2'])]),axis=1)
    
    

###################################  
##    Saving results as shape    ##
################################### 

#shape file
geo_df = gpd.GeoDataFrame(entire_network_splitted, crs=crs)
geo_df.to_file('output/finalresult_splitted.shp')
