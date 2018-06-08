from netCDF4 import Dataset
import gdal
import numpy as np

file='DECSO-CHINA.NOx_2016_OMI.v5.1.nc'
input=Dataset(file)
data=input.variables['NOx'][:]
lat=input.variables['lat'][:]
lon=input.variables['lon'][:]

wkt='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'
data_jan=data[1,:,:]

dlon=lon[1]-lon[0]
dlat=lat[0]-lat[1]

driver = gdal.GetDriverByName('GTiff')
dataset = driver.Create(file.replace('nc','tif'),len(lon),len(lat),1, gdal.GDT_Float32)
dataset.SetGeoTransform((min(lon)-dlon/2.0,dlon,0,max(lat)+dlat/2.0,0,dlat) )
dataset.SetProjection(wkt)
dataset.GetRasterBand(1).WriteArray(np.flipud(data_jan))
dataset.FlushCache()  
