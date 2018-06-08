from MakeBufferRoads import runAllBufferRoads
from MakeSCnoSCFile import ProcessingSCnoSC
import gdal
import numpy as np
import os,shutil

inputfolder='OSPM_PP_runs_2017/' 
outputfolder='Output_2017/'

# make sc / no-sc division
ProcessingSCnoSC(inputfolder,outputfolder+'/SC_noSC/')
# make buffer around roads
runAllBufferRoads(inputfolder,outputfolder+'/RoadsBuffer/')

print 'Combine all datasets'
print '--------------------'


# combination of both
file_buffer=outputfolder+'/RoadsBuffer/Buffer.tif'
print file_buffer
gdata = gdal.Open(file_buffer)
streets = gdata.ReadAsArray().astype(float)

file_sc=outputfolder+'/SC_noSC/LocationsSC.tif'
print file_sc
gdata = gdal.Open(file_sc)
gt = gdata.GetGeoTransform()
sc = gdata.ReadAsArray().astype(float)
proj=gdata.GetProjection()

ospm=np.full_like(sc,np.nan)
loc=(streets==1) & (sc==1)
ospm[loc]=1

# save  
[cols, rows] =np.shape(ospm)
driver = gdal.GetDriverByName("GTiff")
outdata = driver.Create(outputfolder+'OSPM_locations.tif', rows, cols, 1, gdal.GDT_Byte,options = [ 'COMPRESS=LZW' ])
outdata.SetGeoTransform(gt)
outdata.SetProjection(proj)
outdata.GetRasterBand(1).WriteArray(ospm)
outdata.GetRasterBand(1).SetNoDataValue(-9999)
outdata.FlushCache() ##saves to disk!!
outdata = None
