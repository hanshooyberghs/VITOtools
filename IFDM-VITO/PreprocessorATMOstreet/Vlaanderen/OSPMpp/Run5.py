import pandas as pd
import os
from shapely.geometry import Point
import geopandas as gpd
from joblib import Parallel, delayed
import shutil


#########################
## verwerking per zone ##
#########################
def calc_per_zone(nummer):
    crs='epsg:31370'
    kolom_save=['X','Y','HEIGTH','WIDTH','LENGTH1','LENGTH2','THETA','FID']

    nr_str=str(nummer+1).zfill(2) 
    print 'Processing zone: ' +nr_str
    folder_run='PP_runs/'+nr_str+'/'
    
    # outputfiles per zone
    nosc_zone=folder_run+'nosc_z'+nr_str+'.csv'
    sc_zone=folder_run+'sc_z'+nr_str+'.csv'
    
    # read file
    input=pd.read_csv(folder_run+'OSPM_punten.txt',sep='\t')
  
    
    # add geometry
    input['geometry']=''
    for i in input.index:
        input.loc[i,'geometry']=Point(input.loc[i,'X'],input.loc[i,'Y'])

    # find street canyon locations
    loc_sc=(input['WIDTH']>0)&(input['HEIGTH']>200)

    # save street canyon
    sc_data=input[loc_sc]   
    sc_data.to_csv(sc_zone)
    geo_df = gpd.GeoDataFrame(sc_data[kolom_save], crs=crs,geometry=sc_data['geometry'])
    geo_df.to_file(sc_zone.replace('.csv','.shp'))

    
    # save non-street canyon + conversion to 
    nosc_data=input[~loc_sc]
    nosc_data.to_csv(nosc_zone)
    geo_df = gpd.GeoDataFrame(nosc_data[kolom_save], crs=crs,geometry=nosc_data['geometry'])
    geo_df.to_file(nosc_zone.replace('.csv','.shp'))
    
    print 'Finished zone: ' +nr_str
    return (input,sc_data,nosc_data)

########################
##   Combine results  ##
########################

# set-up
nummers=range(32)
dfs_sg=[]
dfs_sc=[]
dfs_nosc=[]
kolom_save=['X','Y','HEIGTH','WIDTH','LENGTH1','LENGTH2','THETA','FID']

# par processing
results = Parallel(n_jobs=32)(delayed(calc_per_zone)(nummer) for nummer in nummers)
for i in range(len(results)):
    dfs_sg.append(results[i][0])
    dfs_sc.append(results[i][1])
    dfs_nosc.append(results[i][2])


# concat
nosc_combined=pd.concat(dfs_nosc,ignore_index=True)
sc_combined=pd.concat(dfs_sc,ignore_index=True)
sg_combined=pd.concat(dfs_sg,ignore_index=True)


# outputfiles
outputfile_sg='straatgegegevens.txt'
outputfile_nosc='nosc.csv'
outputfile_sc='sc.csv'

# conversion for straatgegevens file
output=sc_combined.copy()
output.drop('geometry', 1,inplace=True)
output[['X','Y']]/=1000.0   # rescale coordinates
output['HEIGTH']/=100.0          # rescale height column
kolommen_hoogte=[x for x in sg_combined.columns if '_H' in x]
output[kolommen_hoogte]/=100.0  # rescale height exceptions
output[output <0]=-1            # correction for rescaling height exceptions
output.loc[output.LENGTH1<30,'LENGTH1']=30  # minimal length
output.loc[output.LENGTH2<30,'LENGTH2']=30  # minimal length
file='straatgegevens.txt'
with open(file, 'w') as f:
    f.write(str(len(output))+'\n')
    output.to_csv(f,index=False,header=None,sep='\t')


# save combined files sc / nosc (slow step => sometimes off)
crs='epsg:31370'
geo_df = gpd.GeoDataFrame(sc_combined[kolom_save], crs=crs,geometry=sc_combined['geometry'])
geo_df.to_file('Output/sc.shp')
geo_df = gpd.GeoDataFrame(nosc_combined[kolom_save], crs=crs,geometry=nosc_combined['geometry'])
geo_df.to_file('Output/nosc.shp')

#############################
##   config_exc uitvoeren  ##
#############################

# lijnbronnen klaarzetten
lijnbronfile='../Stap3_FASTRACE/VerwerkteEmissies.xlsx'
lijnbronnen=pd.read_excel(lijnbronfile)
file='Lijnbronnen.txt'
kolommen=['XA','YA','XB','YB','width','height']
with open(file, 'w') as f:
    f.write(str(len(lijnbronnen))+'\n')
    lijnbronnen[kolommen].to_csv(f,index=False,header=None,sep='\t')
    

# runnen fortran script
os.system('./config_exc.out')

# move output to correct directory
shutil.move('straatgegevens.txt','Output/')
shutil.move('straatgegevens_config_exc.txt','Output/')
shutil.move('Lijnbronnen.txt','Output/')

