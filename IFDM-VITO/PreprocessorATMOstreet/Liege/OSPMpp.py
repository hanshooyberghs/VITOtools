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
file_gebouwen='GebouwenHeightCentimer.shp'
file_lijnbronnen='../RoadsProcessing/output/finalresult_splitted.shp'




##########################
##  FOLDER KLAARZETTEN  ##
##########################
folder_run='runfolder_centimeter/'


# remove existing folder
if os.path.exists(folder_run):
    shutil.rmtree(folder_run) 
    os.makedirs(folder_run)
else:
    os.makedirs(folder_run)

shutil.copy('ospmpp',folder_run)

#####################
##    GEBOUWEN     ##
#####################


# ogr bewerkingen
file_gebouwen_zone_2=folder_run+'gebouwen_2.shp'
cmd=ogr+' '+file_gebouwen_zone_2+' '+file_gebouwen
print cmd
os.system(cmd)

file_gebouwen_correct=folder_run+'gebouwen_OK.shp'
cmd=ogr+' -nlt POLYGON '+file_gebouwen_correct+' '+file_gebouwen_zone_2
print cmd
os.system(cmd)

file_lijnbronnen_correct=folder_run+'lijnbronnen_OK.shp'
cmd=ogr+' '+file_lijnbronnen_correct+' '+file_lijnbronnen
print cmd
os.system(cmd)


    
##########################
##   OSPM PP UITVOEREN  ##
##########################
runcommando='./ospmpp lijnbronnen_OK.shp gebouwen_OK.shp'+\
    ' --out OSPM_punten.txt --fmt 0 --hattrib MEDIAN --ds 20 --plen 50 --hmode 2 > out_ospmpp.txt'
cmd='cd '+folder_run+'; ulimit -s unlimited; chmod +x ospmpp;'+runcommando 
print cmd
os.system(cmd)
