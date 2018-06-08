# PP voor 1 domein
# input variable: nummer van domein

####################
## Load modules   ##
####################
import sys
import os
import shutil
import geopandas as gpd
import pandas as pd
from shapely.geometry import LineString
pd.options.mode.chained_assignment = None

saga_cmd='saga_cmd'
ogr='ogr2ogr'
file_gebouwen='../Data/Gebouwen/Gebouwen_GDI_final_version.shp'
file_lijnbronnen='../Stap3_FASTRACE/VerwerkteEmissies.shp'
folder_domeinen= '../Data/Domeinen/'


#####################
## Input verwerken ##
#####################

input=sys.argv
nummer=input[1]



##########################
##  FOLDER KLAARZETTEN  ##
##########################
folder_run='PP_runs/'+str(nummer).zfill(2)+'/'
# remove existing folder
if os.path.exists(folder_run):
   shutil.rmtree(folder_run) 

# copy files of PP_basic to new dir
shutil.copytree('PP_Basic',folder_run)

# temporary folder
folder_tmp=folder_run+'/tmp/'
os.makedirs(folder_tmp)


#####################
##   LIJNBRONNEN   ##
#####################
file_lijnbronnen_zone=folder_tmp+'/lijnbronnen_zone.shp'
file_lijnbronnen_zone_tmp=folder_tmp+'/lijnbronnen_zone_tmp.shp'


# clippen
file_zone=folder_domeinen+'domein'+str(nummer)+'.shp'
cmd=saga_cmd+' shapes_polygons 11 -CLIP '+file_zone+\
    ' -MULTIPLE 0 -S_INPUT '+file_lijnbronnen+\
    ' -S_OUTPUT '+file_lijnbronnen_zone_tmp
print cmd
os.system(cmd)

# inlezen en correcte coordinaten eraan hangen
lijnbronnen=gpd.read_file(file_lijnbronnen_zone_tmp)
for i in lijnbronnen.index:
    coord_segment=lijnbronnen.loc[i,'geometry'].coords[:]
    lijnbronnen.loc[i,'X1']=coord_segment[0][0]
    lijnbronnen.loc[i,'Y1']=coord_segment[0][1]
    lijnbronnen.loc[i,'X2']=coord_segment[1][0]
    lijnbronnen.loc[i,'Y2']=coord_segment[1][1]

crs='epsg:31370'
geo_df = gpd.GeoDataFrame(lijnbronnen, crs=crs, geometry=lijnbronnen['geometry'])
geo_df.to_file(file_lijnbronnen_zone)


#####################
##    GEBOUWEN     ##
#####################

# buffer maken voor gebouwen
file_buffer=folder_tmp+'buffer.shp'
cmd=saga_cmd+' shapes_tools 18 -SHAPES '+file_lijnbronnen_zone+' -BUFFER '+ file_buffer+ \
    ' -DIST_FIELD_DEFAULT 100 -DISSOLVE 1'
print cmd
os.system(cmd)

# polygon intersect
file_gebouwen_zone=folder_tmp+'gebouwen.shp'
cmd=saga_cmd+' shapes_polygons 14 -A '+file_gebouwen+' -B '+ file_buffer+ ' -RESULT '+file_gebouwen_zone
print cmd
os.system(cmd)

# ogr bewerkingen
file_gebouwen_zone_2=folder_tmp+'gebouwen_2.shp'
cmd=ogr+' '+file_gebouwen_zone_2+' '+file_gebouwen_zone
print cmd
os.system(cmd)

file_gebouwen_correct=folder_tmp+'gebouwen_OK.shp'
cmd=ogr+' -nlt POLYGON '+file_gebouwen_correct+' '+file_gebouwen_zone_2
print cmd
os.system(cmd)

file_lijnbronnen_correct=folder_tmp+'lijnbronnen_OK.shp'
cmd=ogr+' '+file_lijnbronnen_correct+' '+file_lijnbronnen_zone
print cmd
os.system(cmd)


##########################
##  OSPM PP KLAARZETTEN ##
##########################

ext_list=['shp','dbf','shx']
for ext in ext_list:
    shutil.copy(file_lijnbronnen_correct.replace('shp',ext),folder_run+'/data/')
    shutil.copy(file_gebouwen_correct.replace('shp',ext),folder_run+'/data/')

    
##########################
##   OSPM PP UITVOEREN  ##
##########################
runcommando='./ospmpp data/lijnbronnen_OK.shp data/gebouwen_OK.shp '+\
    ' --out OSPM_punten.txt --fmt 0 --hattrib altitude --ds 20 --plen 50 --hmode 2 > out_ospmpp.txt'
cmd='cd '+folder_run+'; ulimit -s unlimited; chmod +x ospmpp;'+runcommando 
print cmd
os.system(cmd)
