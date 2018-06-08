# Buffer file aanmaken voor 1 domein
# input variable: nummer van domein

####################
## Load modules   ##
####################
import sys
import os
import shutil
import geopandas as gpd
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
import gdal



x_min=18700
x_max=259000
y_min=147100
y_max=246100

####################################
##    Start parallel processing   ##
####################################    


def runAllBufferRoads(folder_input,folder_output):
    print 'Make buffers around roads'
    print '-------------------------'
    
    if os.path.exists(folder_output):
       shutil.rmtree(folder_output) 
    os.makedirs(folder_output)

    nummers=range(1,33)
    results = Parallel(n_jobs=32)(delayed(runFolder)(nummer,folder_input,folder_output) for nummer in nummers)
    
    print 'Start merging'
    
    # read size of data
    gdata=gdal.Open(results[0])
    data=gdata.ReadAsArray().astype(float)
    union=np.zeros_like(data)
    gdata=None
    
    for file in results:
        print 'Read '+file
        gdata=gdal.Open(file)
        data=gdata.ReadAsArray().astype(float)
        union=union+data
        
    # read geo information
    gt = gdata.GetGeoTransform()
    proj=gdata.GetProjection()

    # some conversion
    print 'Save'
    union[union>0]=1
    union[union==0]=-9999
    print 'Merged'

    # saving
    file_buffer=folder_output+'/Buffer.tif'
    print 'Save results in GeoTiff'
    [cols, rows] = union.shape
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(file_buffer, rows, cols, 1, gdal.GDT_Byte,options = [ 'COMPRESS=LZW' ])
    outdata.SetGeoTransform(gt)
    outdata.SetProjection(proj)
    outdata.GetRasterBand(1).WriteArray(union)
    outdata.GetRasterBand(1).SetNoDataValue(-9999)
    outdata.FlushCache() ##saves to disk!!
    outdata = None
    
    print 'Made buffers around roads'
    print '-------------------------'





def runFolder(nummer,folder_input,folder_output):
    global x_min
    global x_max
    global y_min
    global y_max
    grid=str(x_min)+' '+str(y_min)+' '+str(x_max)+' '+str(y_max)
    crs='epsg:31370'
    
    print 'Start processing of folder '+str(nummer)
    #####################
    ## Input verwerken ##
    #####################

    folder_ospmpp=folder_input+str(nummer).zfill(2)+'/'

    file_gebouwen=folder_ospmpp+'tmp/gebouwen_OK.shp'
    file_lijnbronnen=folder_ospmpp+'tmp/lijnbronnen_OK.shp'
    file_possc=folder_ospmpp+'sc_z'+str(nummer).zfill(2)+'.shp'

    folder_run=folder_output+'/'+str(nummer).zfill(2)+'/'
    # remove existing folder
    if os.path.exists(folder_run):
       shutil.rmtree(folder_run) 
    os.makedirs(folder_run)
       
       
    ##############################################
    ## Koppelen van sc breedte aan lijnbronnen  ##
    ##############################################

    # Buffer maken van 25m rond de wegen
    file_buffer_25=folder_run+'buffer25.shp'
    lijnbronnen=gpd.read_file(file_lijnbronnen)
    lijnbronnen.crs=crs
    buffer_25=gpd.GeoDataFrame(lijnbronnen.buffer(25),columns=['geometry'],crs=crs)
     
    # gemiddelde wegbreedte per buffer
    pos_sc=gpd.read_file(file_possc)
    pos_sc.crs=crs
    joined=gpd.sjoin(buffer_25,pos_sc,how='left',op='contains')
    joined['lijnbron']=joined.index
    avgwidth_per_linesource_joined=joined.groupby('lijnbron',as_index=False).mean()[['WIDTH','lijnbron']]
    lijnbronnen_metbuffer=avgwidth_per_linesource_joined=joined.\
        groupby('lijnbron',as_index=False).mean()[['WIDTH','lijnbron']]
    lijnbronnen_metbuffer=lijnbronnen_metbuffer.join(lijnbronnen['geometry'])

    # in geval er geen sc is => neem breedte 25 als default (kan misschien weg later)
    lijnbronnen_metbuffer['WIDTH'].fillna(25,inplace=True)

    # bepaal lengte
    lijnbronnen_metbuffer['WIDTH']=lijnbronnen_metbuffer['WIDTH']/2.0

    ##################################
    ## Definitieve buffer aanmaken  ##
    ##################################

    # buffer van avg width / 2 (including dissolve)
    geo=lijnbronnen_metbuffer.apply(lambda p: p.geometry.buffer(p.WIDTH),axis=1)
    buffer=gpd.GeoDataFrame(geo,columns=['geometry'],crs=crs)
    buffer['dissolve_column']=1
    buffer=buffer.dissolve('dissolve_column')

    buffer['ID']=1
    buffer_out=buffer.dissolve('ID')
    buffer_out['ID']=1.0

    # save
    file_buffer_OK=folder_run+'buffer_OK.shp'
    buffer_out.to_file(file_buffer_OK)


    ###############
    ## Gridding  ##
    ###############
    file_out=folder_run+'Buffer.tif'
    cmd='gdal_rasterize -a ID -co "COMPRESS=LZW" -ot Byte -te '+grid+' -tr 10 10 '+file_buffer_OK+' '+file_out +' > /dev/null'
    os.system(cmd)
    
    print 'End processing of folder '+str(nummer)
    
    return file_out
    

    
