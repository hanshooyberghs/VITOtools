"""
Postprocessing routine for IFDM, adaptation for LiboVITO-system
===============================================================

Basic description: visualize one IFDM-run using Python and GDAL


Python-Dependencies:
* Python routines:
    -os
    -pandas
* External software:
    - gdal

* To add: read time series / select location of OPAQ-points + coupling with OPAQ
"""

############################
# Input parameters
############################

# folder with IFDM results
Runfolder='/projects/uitstel-no2/IFDM_runs/2020_BAU/AllSectors/01'

# folder in which GeoTIFF is saved (should exist)
Outputfolder='/projects/uitstel-no2/IFDM_runs/2020_BAU/AllSectors/01'

# list of pollutants
Pollutants=["NO2","PM25","PM10"]

# gdal_grid command (should not be modified if GDAL is installed correctly)
gdal_grid='gdal_grid'

# output resolution in meter
resolution=100  

#################################      
# Import modules
#################################

import os
import pandas

Pollutants = [w.ljust(4,'a') for w in Pollutants]    


for pollutant in Pollutants:
    # define file names
    runfile=Runfolder+'/output'+pollutant+'_saga.txt'
    result=Outputfolder+'/'+pollutant+'.tif'   
    csv_file=Runfolder+'/results'+pollutant+'.csv'
    vrtfile=Runfolder+'/outputformat'+pollutant+'.vrt'
    
    # read data and store in format required for gdal_grid
    data=pandas.read_csv(runfile,sep=';')
    data.columns=['x','y',pollutant]
    data.to_csv(csv_file,index=False)
    
    # deduce number of grids in x and y direction
    x_min=min(data['x'])
    x_max=max(data['x'])
    y_min=min(data['y'])
    y_max=max(data['y'])
    print str(x_min)+' '+str(x_max)+' '+str(y_min)+' '+str(y_max)
    x_diff=(x_max-x_min)/resolution
    y_diff=(y_max-y_min)/resolution

    # make vrt file (determines format for gdal_grid)    
    f= open(vrtfile, 'w')
    content='<OGRVRTDataSource> \n <OGRVRTLayer name="results'+pollutant+'"> \n <SrcDataSource>'+csv_file \
        +'</SrcDataSource> \n <GeometryType>wkbPoint</GeometryType> \n ' \
        +'<LayerSRS>EPSG:31370</LayerSRS> <GeometryField encoding="PointFromColumns" x="X" y="Y" z="' \
        +pollutant+'"/> \n </OGRVRTLayer> \n </OGRVRTDataSource>'
    f.write(content)
    f.close()

    # grid data
    commando=gdal_grid+' -a LINEAR:radius=0:nodata=-9999 -outsize '+str(x_diff)+\
        ' '+str(y_diff)+' -zfield '+pollutant+' '+vrtfile+' '+result
    print commando
    os.system(commando)  


            
            
