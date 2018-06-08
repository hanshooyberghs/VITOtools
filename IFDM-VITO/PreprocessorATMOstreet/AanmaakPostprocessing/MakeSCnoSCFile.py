import pandas as pd
import os
from shapely.geometry import Point
import geopandas as gpd
from joblib import Parallel, delayed
import shutil
import numpy as np

xmin=18700
xmax=259000
ymin=147100
ymax=246100


#########################
## verwerking per zone ##
#########################
def calc_per_zone(nummer,folder_input):
    crs='epsg:31370'
    kolom_save=['X','Y','HEIGTH','WIDTH','LENGTH1','LENGTH2','THETA','FID']

    nr_str=str(nummer+1).zfill(2) 
    print 'Processing zone: ' +nr_str
    folder_run=folder_input+nr_str+'/'
    
    # outputfiles per zone
    nosc_zone=folder_run+'nosc_z'+nr_str+'.csv'
    sc_zone=folder_run+'sc_z'+nr_str+'.csv'
    
    # read file
    input=pd.read_csv(folder_run+'OSPM_punten.txt',sep='\t')
 
    
    # add geometry
    input['geometry']=''
    input.geometry=input.apply(lambda p: Point(p['X'],p['Y']),axis=1)

    # find street canyon locations
    loc_sc=DetermineStreetCanyonLocations(input)

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

def DetermineStreetCanyonLocations(input):

    # where some width is found
    input['check']=((input.WIDTH>0)&(input.HEIGTH>200))
    
    # loop over streets
    streets=input.FID.unique()
    loc=[]
    for street in streets:
        
        data_street=input[input.FID==street]
        
        if len(data_street)>2:
            # find locations with minimal three points in line sc
            controle=data_street['check'].rolling(window=3,center=True).mean()
            sc_loc=controle==1
                    
            # add points before and after for these points as well
            locations=sc_loc[sc_loc].index
            before=locations[locations>sc_loc.index.min()]-1
            sc_loc[before]=True
            after=locations[locations<sc_loc.index.max()]+1
            sc_loc[after]=True
            
            # find locations with a lonely no-sc, make them sc
            controle=sc_loc.rolling(window=3,center=True).mean()
            lonely_sc=(controle>0.5)&(sc_loc==False)
            sc_loc[lonely_sc]=True
            
        else:
            controle=data_street['check'].mean()
            sc_loc=pd.Series(index=data_street.index)
            sc_loc[:]=controle==1

        loc.append(sc_loc)
    
    
    loc=pd.concat(loc)
    return loc
    
########################
##   Combine results  ##
########################

def ProcessingSCnoSC(folder_input,folder_output):
    print 'Determine street canyon locations'
    print '----------------------------------'

    global x_min
    global x_max
    global y_min
    global y_max

    
    if os.path.exists(folder_output):
       shutil.rmtree(folder_output) 
    os.makedirs(folder_output)

    # set-up
    nummers=range(32)
    dfs_sg=[]
    dfs_sc=[]
    dfs_nosc=[]
    kolom_save=['X','Y','HEIGTH','WIDTH','LENGTH1','LENGTH2','THETA','FID']

    # par processing
    results = Parallel(n_jobs=32)(delayed(calc_per_zone)(nummer,folder_input) for nummer in nummers)
    for i in range(len(results)):
        dfs_sg.append(results[i][0])
        dfs_sc.append(results[i][1])
        dfs_nosc.append(results[i][2])

    print 'Start concatenating'
    # concat
    nosc_combined=pd.concat(dfs_nosc,ignore_index=True)
    sc_combined=pd.concat(dfs_sc,ignore_index=True)
    sg_combined=pd.concat(dfs_sg,ignore_index=True)
    


    # combine sc and no sc
    crs='epsg:31370'
    geo_df_sc = gpd.GeoDataFrame(sc_combined[kolom_save], crs=crs,geometry=sc_combined['geometry'])
    geo_df_sc.to_file(folder_output+'sc.shp')
    geo_df_nosc = gpd.GeoDataFrame(nosc_combined[kolom_save], crs=crs,geometry=nosc_combined['geometry'])
    geo_df_sc['LocationsSC']=1
    geo_df_nosc['LocationsSC']=0
    combined=pd.concat([geo_df_sc,geo_df_nosc],ignore_index=True)
    
    # streetcanyon file
    output=sc_combined.copy()
    output=output[(output.WIDTH>0)&(output.HEIGTH>200)]
    output.drop('geometry', 1,inplace=True)
    output[['X','Y']]/=1000.0   # rescale coordinates
    output['HEIGTH']/=100.0          # rescale height column
    kolommen_hoogte=[x for x in sg_combined.columns if '_H' in x]
    output[kolommen_hoogte]/=100.0  # rescale height exceptions
    output[output <0]=-1            # correction for rescaling height exceptions
    output.loc[output.LENGTH1<30,'LENGTH1']=30  # minimal length
    output.loc[output.LENGTH2<30,'LENGTH2']=30  # minimal length
    file=folder_output+'straatgegevens.txt'
    with open(file, 'w') as f:
        f.write(str(len(output))+'\n')
        output.to_csv(f,index=False,header=None,sep='\t')


    # change points close to each other (within 10m distance from each other)
    buffer_points=combined.geometry.buffer(10)
    buffer_points=gpd.GeoDataFrame(buffer_points,columns=['geometry'],crs=crs)
    buffer_points.index=buffer_points.index.rename('left_index')
    combined.index=combined.index.rename('right_index')
    joined=gpd.sjoin(buffer_points,combined,how='left',op='intersects')
    type_mean=joined.groupby(by=joined.index).mean()['LocationsSC']
    type_mean[type_mean>0.5]=1
    type_mean[type_mean<=0.5]=0

    combined.drop(['LocationsSC'],axis=1,inplace=True)
    combined=combined.merge(pd.DataFrame(type_mean),left_index=True,right_index=True)

    combined.to_csv(folder_output+'LocationsSC.csv')
    #combined.to_file(folder_output+'LocationsSC.shp')
    
    vrt_file='<OGRVRTDataSource>\n'+\
         '<OGRVRTLayer name="LocationsSC">\n'+\
         '<SrcDataSource>'+folder_output+'LocationsSC.csv</SrcDataSource> \n'+\
         '<GeometryType>wkbPoint</GeometryType> \n'+\
         '<LayerSRS>EPSG:31370</LayerSRS> <GeometryField encoding="PointFromColumns" x="X" y="Y" z="LocationsSC"/> \n'+\
         '</OGRVRTLayer> \n'+\
         '</OGRVRTDataSource>'
    f=open('gridding.vrt','w')
    f.write(vrt_file)
    f.close()    

    outputfile=folder_output+'LocationsSC.tif'
    cmd='gdal_grid  -a NEAREST:radius1=100:radius2=100:nodata=-9999 -ot Byte '+\
        '-txe '+str(xmin)+' '+str(xmax)+' -tye '+str(ymax)+' '+str(ymin)+' -outsize '+\
        '24030 9900 -zfield LocationsSC '+\
        '-co COMPRESS=LZW gridding.vrt '+outputfile
    os.system(cmd)
    
    print 'Determined street canyon locations'
    print '----------------------------------'

