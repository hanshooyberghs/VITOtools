""""
usage: ETDM.py country_code

ATTENTION! The script will use 32 cores!

Welcome to the European Traffic Data Model

The model will perform the following actions:
    1. Process the OTM data. The following steps will be performed:
    * Remove roads without traffic counts
    * Simply roads with 0.0005 degree (approx. 50m) resolution 
    * Split lines at vertices
    * Other columns (height, speed...) added  
    
    2. Spread the COPERT-vehicle kms according to the OTM data, using capacity

positional arguments:
  country_code          Country Code for which data will be read

Usage: conda activate /projects/concawe/ETDM/anaconda/
  
"""

# Note: issue with input files exists:
# * Empty files => IO error mentioned. Can be neglected
# * Incorrect format => still to be checked

import geopandas as gpd
import glob, argparse, multiprocessing,os,pyproj,gdal
import pandas as pd
from shapely.geometry import LineString
from joblib import Parallel, delayed
import numpy as np
pd.options.mode.chained_assignment = None 

column_list=[ u'inspireid',   u'sourceid', u'beginlifes', 
        u'direction',u'nationalro',u'functional',  u'formofway', u'roadsurfac', u'speedlimit',
         u'capacity', u'trafficvol',u'geometry'] 

"""
Aux Routines
"""

#############################
## split lines aux routine ##
#############################

def split_line(index,row):
    coord_segment=row['geometry'].coords[:]
    lengte_segment=len(coord_segment)
    splitted=pd.DataFrame(columns=['id','lat1','lon1','lat2','lon2'],index=range(lengte_segment-1))
    for j in range(lengte_segment-1):
        splitted.loc[j,'id']=row['id']
        splitted.loc[j,'lat1']=coord_segment[j][1]
        splitted.loc[j,'lon1']=coord_segment[j][0]
        splitted.loc[j,'lat2']=coord_segment[j+1][1]
        splitted.loc[j,'lon2']=coord_segment[j+1][0]
    return(splitted)

###########################
## read file aux routine ##
###########################
    
def read_file(file):
    global column_list
    try:
        data=gpd.read_file(file)
        
        if not ('trafficvol' in data):
            print('Different format file: '+file)
            print('Try to combine csv and shape file')
            return pd.DataFrame(columns=column_list)
        else: 
            temp= data.loc[data.trafficvol>0,column_list]
            return temp
    except:
        print('Empty file: '+file)
        return pd.DataFrame(columns=column_list)

"""
Main routine
"""
parser = argparse.ArgumentParser(description='Processing Open Transport Maps')
parser.add_argument('country', metavar='country_code', type=str, help='Country code')
parser.add_argument('-u','--urban',dest='file_urban',metavar='string', type=str, 
    help='File with JRC Global Settlement. Default: JRC_GlobalSettlement/GHS_WGS84.tif',
    default='JRC_GlobalSettlement/GHS_WGS84.tif')
parser.add_argument('-r','--roads',dest='otm',metavar='string', type=str, 
    help='Folder with raw OTM input/ Default: otminput/',
    default='otminput/')
parser.add_argument('-c','--copert',dest='copert',metavar='string', type=str, 
    help='Folder with COPERT input. Default: COPERTinput/',
    default='COPERTinput/')
parser.add_argument('-o','--output',dest='output',metavar='string', type=str, 
    help='OutputFolder. Default: Output/',
    default='Output/')

args = parser.parse_args()
    

country=args.country
print('Welcome to the European Traffic Data Model. \n'+
    'We hope you enjoy your time using this VITO Software.')
print('----------------------------------------------------\n')

print('ETDM is now running for country: '+country)
print('Processing may take some time, so grab a coffee!')
    
###################
## Reading files ##
###################   

# read and combine
file_list=glob.glob(args.otm+'roadlinks_'+country+'*.shp')

dfs=[]
print('Start reading files')
print('Number of files: '+str(len(file_list)))
pool = multiprocessing.Pool(processes=32)  
results=pool.map(read_file, file_list) 
pool.close()
pool.join()
combined=pd.concat(results)


##################################
## Selecting roads with traffic ##
##################################  
print('Combine files')
combined=combined[combined.trafficvol>0]
combined.drop_duplicates('inspireid',inplace=True)

###########################
## Simplifying the roads ##
###########################  
print('Simplify roads')
simplified=combined.simplify(50*0.00001)
combined['geometry']=simplified
combined['id']=combined.inspireid
keep_columns=['id','trafficvol','geometry','functional','speedlimit']
combined=combined[keep_columns]

#########################
## Splitting the roads ##
#########################  
print('Split segments')
dfs=[]
results = Parallel(n_jobs=32)(delayed(split_line)(index,row) for index,row in combined.iterrows())
for i in range(len(results)):
    dfs.append(results[i])
    
splitted_concat=pd.concat(dfs)     
splitted_concat=splitted_concat.astype(float)
splitted_data=pd.merge(combined,splitted_concat,left_on='id',right_on='id')
splitted_data['geometry']=splitted_data.apply(lambda p:LineString([(p['lon1'],p['lat1']),\
    (p['lon2'],p['lat2'])]),axis=1)
geod = pyproj.Geod(ellps='WGS84')
splitted_data['length_m']=splitted_data.apply(lambda p: geod.inv(p['lon1'],p['lat1'],p['lon2'],p['lat2'])[2],axis=1)


###################
## Add road type ##
###################
# Load settlement data
print('Add road types')
data=gdal.Open(args.file_urban)
gt=data.GetGeoTransform()
settlement=data.ReadAsArray()
x = ((splitted_data.geometry.centroid.x - gt[0])/gt[1]).astype(int)
y = ((splitted_data.geometry.centroid.y - gt[3])/gt[5]).astype(int)
splitted_data['settlement']=settlement[y,x] 

# Add the values 
loc_highway=(splitted_data.functional=='mainRoad')
loc_urban=(splitted_data.settlement.isin([2,3]))
splitted_data['roadtype']='R'
splitted_data.loc[loc_urban,'roadtype']='U'
splitted_data.loc[loc_highway,'roadtype']='H'


#############################
## Read COPERT vehicle kms ##
#############################
print('Read COPERT vehicle km')
COPERT_counts=pd.read_csv(args.copert+'/'+country+'.txt',sep=';')

# rename the columns for final output
COPERT_counts.Fuel=COPERT_counts.Fuel.str[0]+'_'
old_names=['Passenger Cars','Light Commercial Vehicles','Heavy Duty Trucks','Buses','L-Category']
new_names=['CAR_','LDV_','HDV_','BUS_','TWH_']
COPERT_counts.Category.replace(old_names,new_names,inplace=True)
COPERT_counts.Euro=COPERT_counts.Euro.str.replace('Euro ','')
COPERT_counts.Euro=COPERT_counts.Euro.str.replace('0-','<')
COPERT_counts['OutputCategory']=COPERT_counts[['Fuel','Category','Euro']].sum(axis=1).str.lower()
print('Total vehicle kms year based on COPERT:')
print('\tHighway:'+str(COPERT_counts[COPERT_counts.RoadType=='H'].km.sum()))
print('\tRural:  '+str(COPERT_counts[COPERT_counts.RoadType=='R'].km.sum()))
print('\tUrban:  '+str(COPERT_counts[COPERT_counts.RoadType=='U'].km.sum()))
print('\tU+R:    '+str(COPERT_counts[COPERT_counts.RoadType!='H'].km.sum()))



###########################
## Spread traffic counts ##
###########################
print('Calculate spreading factors based on OTM')
splitted_data['veh_kms_year']=splitted_data['trafficvol']*splitted_data['length_m']/1e3*365
total_veh_kms_year_H=splitted_data.loc[splitted_data.roadtype=='H','veh_kms_year'].sum()
total_veh_kms_year_R=splitted_data.loc[splitted_data.roadtype=='R','veh_kms_year'].sum()
total_veh_kms_year_U=splitted_data.loc[splitted_data.roadtype=='U','veh_kms_year'].sum()
total_veh_kms_year_other=splitted_data.loc[splitted_data.roadtype!='H','veh_kms_year'].sum()

print('Total vehicle kms year based on OTM:')
print('\tHighway:'+str(total_veh_kms_year_H))
print('\tRural:  '+str(total_veh_kms_year_R))
print('\tUrban:  '+str(total_veh_kms_year_U))
print('\tOther:   '+str(total_veh_kms_year_other))

print('Ratios OTM / COPERT')
print('\tHighway: '+str(total_veh_kms_year_H/COPERT_counts[COPERT_counts.RoadType=='H'].km.sum()))
print('\tOther:   '+str(total_veh_kms_year_other/COPERT_counts[COPERT_counts.RoadType!='H'].km.sum()))


print('Actual spreading for all parameters')
vehicle_types=COPERT_counts['OutputCategory'].unique()    
loc_nonh_otm=splitted_data.roadtype!='H'
loc_h_otm=splitted_data.roadtype=='H'

for param in vehicle_types:
    # highway
    splitted_data[param+'_vehkms']=np.nan
    loc_h_copert=((COPERT_counts.OutputCategory==param)&(COPERT_counts.RoadType=='H'))
    splitted_data.loc[loc_h_otm,param+'_vehkms']=float(COPERT_counts.loc[loc_h_copert,'km'])* \
        splitted_data.loc[loc_h_otm,'veh_kms_year']/total_veh_kms_year_H

    # non highway
    loc_nonh_copert=((COPERT_counts.OutputCategory==param)&(COPERT_counts.RoadType!='H'))
    splitted_data.loc[loc_nonh_otm,param+'_vehkms']=COPERT_counts.loc[loc_nonh_copert,'km'].sum()* \
        splitted_data.loc[loc_nonh_otm,'veh_kms_year']/total_veh_kms_year_other
        
    # check => larger than 1
    check1=splitted_data[param+'_vehkms'].sum()
    check2=COPERT_counts[COPERT_counts.OutputCategory==param].sum()['km']
    if check1*check2>0:
        if ((check1-check2)/check1)>1e-4:
            print('ERROR: Check program')
            stop      
    
    # convert to number of vehicles
    splitted_data[param]=splitted_data[param+'_vehkms']/(splitted_data['length_m']/1e3)
    

            
#################
## Road speeds ##
#################
print('Add speeds')
splitted_data['speed']=splitted_data[ 'speedlimit']
# add copert values for missing 
speed_copert=pd.read_csv(args.copert+'/roadtype_speed_COPERT-EU28_'+country+'.txt',sep=';')
road_types=['H','U','R']
for road_type in road_types:
    missing=(splitted_data.speed==0)&(splitted_data.roadtype==road_type)
    select=(speed_copert.Category=='Passenger Cars')&(speed_copert.RoadType==road_type)
    splitted_data.loc[missing,'speed']=int(speed_copert.loc[select,'Speed'])
 
        
#######################
## Add other columns ##
#######################
print('Add other columns')
splitted_data['height']=0
splitted_data['country']=country


###########################
## Save results: CONCAWE ##
###########################
print('Save final result')
OutputFolder=args.output+'/'+country+'/'
if not os.path.exists(OutputFolder):
    os.makedirs(OutputFolder)

columns_correct=['id','lat1','lon1','lat2','lon2']+list(vehicle_types)\
    +['height','roadtype','speed','length_m','geometry','country']
tosave=splitted_data[columns_correct]
tosave.loc[:,'id']=range(len(tosave))
tosave.to_file(OutputFolder+country+'_LineSegments_CONCAWE.shp')



############################
## Save results: ATMOplan ##
############################
print('Save final result in ATMOplan version')
# calculate sums for ATMOplan vehicle types
vehicle_types_ATM=['ldv','hdv','car', 'bus']
splitted_data = splitted_data[splitted_data.columns.drop(list(splitted_data.filter(regex='vehkms')))]

for type in vehicle_types_ATM:
    splitted_data[type]=splitted_data.loc[:,splitted_data.columns.str.contains(type)].sum(axis=1)

columns_correct=['id','lat1','lon1','lat2','lon2']+list(vehicle_types_ATM)\
    +['height','roadtype','speed','length_m','geometry','country']
tosave=splitted_data[columns_correct]
tosave.loc[:,'id']=range(len(tosave))
tosave.to_file(OutputFolder+country+'_LineSegments_ATMOplan.shp')


