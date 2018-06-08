import gdal
import numpy as np
import pandas as pd
import os
from shapely.geometry import Point
import geopandas as gpd


####################
##       input    ##
####################

# inputfiles: output of IFDM and OSPM in csv-file format, with headers ['X','Y','values']
result_ifdm='OutputRuns/IFDM_res.csv'
result_OSPM='OutputRuns/OSPM_res.csv'

# folder where output will be stored
Outputfolder='Gridded_Maps/'

# name of final output file
out_combined=Outputfolder+'RIO-IFDM-OSPM_NO2.tif'


######################
##   fixed input    ##
######################
resolution = 10

# default extent of grid (should not be changed unless buffer file is also changed)
x_min=18700
x_max=259000
y_min=147100
y_max=246100

# positions of street canyons GeoTiff
posities_streetcanyons='/projects/IFDM/PostprocessingOSPM/FinalRoadMap/Output_2017/OSPM.tif'

#crs
crs={'init':'epsg:31370'}

# gdal programs
gdal_grid='gdal_grid'
gdalwarp='gdalwarp'

# max coords => if the numer of coordinates is larger, outputfiles will be read in different parts
maxcoords=100000

if not os.path.exists(Outputfolder):
    os.makedirs(Outputfolder)
else:
    os.system('rm '+Outputfolder+'*')

# convert filenames
result_ifdm_file=result_ifdm.split('/')[-1].split('.')[0]
result_ospm_file=result_OSPM.split('/')[-1].split('.')[0]
    
# read canyon grid file
gdata_canyon = gdal.Open(posities_streetcanyons)
gt_canyon = list(gdata_canyon.GetGeoTransform())
data_canyon = gdata_canyon.ReadAsArray().astype(float)
gdata_canyon = None   

# temporary folder
tmpfolder='tmp/'
if not os.path.exists(tmpfolder):
    os.makedirs(tmpfolder)
else:
    os.system('rm '+tmpfolder+'*')

#############################
##  gridding ifdm results  ##
#############################


print 'Grid IFDM resultaten'
    
    
print str(x_min)+' '+str(y_min)+' '+str(x_max)+' '+str(y_max)
x_diff=(x_max-x_min)/resolution
y_diff=(y_max-y_min)/resolution
print str(x_diff)+' '+str(y_diff)

# make vrt file (determines format for gdal_grid)
vrtfile=tmpfolder+'outputformatIFDM.vrt'
f= open(vrtfile, 'w')
content='<OGRVRTDataSource> \n <OGRVRTLayer name="'+result_ifdm_file+'"> \n <SrcDataSource>'+result_ifdm+ ''\
    +'</SrcDataSource> \n <GeometryType>wkbPoint</GeometryType> \n ' \
    +'<LayerSRS>EPSG:31370</LayerSRS> <GeometryField encoding="PointFromColumns" x="X" y="Y" z="values"' \
    +'/> \n </OGRVRTLayer> \n </OGRVRTDataSource>'
f.write(content)
f.close()

# grid data
gridded_ifdm=Outputfolder+'/IFDM.tif'   
gridded=gdal.Grid(gridded_ifdm,vrtfile, algorithm = 'linear:radius=0:nodata=-9999',\
    outputBounds = [x_min, y_max, x_max , y_min],  creationOptions = ['BIGTIFF=yes', 'COMPRESS=deflate'],\
    outputType = gdal.GDT_Float32,width = x_diff, height = y_diff, \
    zfield = pollutant+'_'+indicator)
del gridded



############################
##  gridding ospm results ##
############################
print 'Grid OSPM resultaten'
# take mean over 10m
ospm=pd.read_csv(result_OSPM)
ospm['geometry']=ospm.apply(lambda p:Point(p.X,p.Y),axis=1)
ospm=gpd.GeoDataFrame(ospm,crs=crs)

buffer_points=ospm.geometry.buffer(10)
buffer_points=gpd.GeoDataFrame(buffer_points,columns=['geometry'],crs=crs)
buffer_points.index=buffer_points.index.rename('left_index')
ospm.index=ospm.index.rename('right_index')
joined=gpd.sjoin(buffer_points,ospm,how='left',op='intersects')
mean=pd.DataFrame(joined.groupby(by=joined.index).mean()['values'],columns=['values'])
ospm_save=ospm.drop('values',axis=1).merge(mean,left_index=True,right_index=True)
ospm_save.to_csv(result_OSPM.replace('.csv','_smoothed.csv'))

del mean
del joined
del buffer_points
del ospm_save
del ospm

# make vrt file (determines format for gdal_grid)
vrtfile=tmpfolder+'/outputformatOSPM.vrt'
f= open(vrtfile, 'w')
content='<OGRVRTDataSource> \n <OGRVRTLayer name="'+result_ospm_file+'_smoothed"> \n <SrcDataSource>'+result_OSPM.replace('.csv','_smoothed.csv')+ ''\
    +'</SrcDataSource> \n <GeometryType>wkbPoint</GeometryType> \n ' \
    +'<LayerSRS>EPSG:31370</LayerSRS> <GeometryField encoding="PointFromColumns" x="X" y="Y" z="values"' \
    +'/> \n </OGRVRTLayer> \n </OGRVRTDataSource>'
f.write(content)
f.close()

# grid data
gridded_ospm=tmpfolder+'OSPM.tif'   
gridded=gdal.Grid(gridded_ospm,vrtfile, algorithm = 'NEAREST:radius1=100:radius2=100:nodata=-9999',\
    outputBounds = [x_min, y_max, x_max , y_min], \
    outputType = gdal.GDT_Float32,width = x_diff, height = y_diff, \
    zfield = pollutant+'_'+indicator, creationOptions = ['COMPRESS=deflate','BIGTIFF=YES'])
del gridded



##############################
##  combining IFDM and OSPM ##
##############################
print 'Combine IFDM and OSPM'

# read ospm data
gdata_ospm = gdal.Open(gridded_ospm)
gt_ospm = gdata_ospm.GetGeoTransform()
data_ospm = gdata_ospm.ReadAsArray().astype(float)
proj=gdata_ospm.GetProjection()
gdata_ospm = None

# read ifdm data
gdata_ifdm = gdal.Open(gridded_ifdm)
gt_ifdm = gdata_ifdm.GetGeoTransform()
data_ifdm = gdata_ifdm.ReadAsArray().astype(float)
gdata_ifdm = None

# read ifdm data
gdata_ifdm = gdal.Open(gridded_ifdm)
gt_ifdm = gdata_ifdm.GetGeoTransform()
data_ifdm = gdata_ifdm.ReadAsArray().astype(float)
gdata_ifdm = None

# checks on grid
if np.shape(data_ifdm)!=np.shape(data_canyon):
    print 'Error: grid with canyons has wrong resolution'
    sys.exit() 

# error in the orientation of the canyon file
if gt_ospm[5]*gt_canyon[5]<0:
    tmp=data_canyon.copy()
    del data_canyon
    data_canyon=np.flipud(tmp)
    del tmp
    gt_canyon[5]*=(-1)
    
# combination
final_result=data_ifdm
loc=((data_canyon==1)&(~np.isnan(data_ospm)))&(data_ospm>data_ifdm)
final_result[loc]=data_ospm[loc]


# saving results
print 'Save results in GeoTiff'
[cols, rows] = final_result.shape
driver = gdal.GetDriverByName("GTiff")
outdata = driver.Create(out_combined, rows, cols, 1, gdal.GDT_Float32,options = [ 'COMPRESS=LZW' ])
outdata.SetGeoTransform(gt_ospm)
outdata.SetProjection(proj)
outdata.GetRasterBand(1).WriteArray(final_result)
outdata.GetRasterBand(1).SetNoDataValue(-9999)
outdata.FlushCache() ##saves to disk!!
outdata = None


del final_result
del data_ifdm
del data_ospm
del loc


print '----------------------'
print 'Processing done'
print '----------------------'
