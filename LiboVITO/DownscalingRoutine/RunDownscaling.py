# -*- coding: utf-8 -*-
"""
Run the preprocessing of IFDM. 
==============================

In this file, the input parameters for the IFDM PreProcessor should be presented.
More details on the input parameters can be found in the Preprocessor.py file.

@author: hooyberh (hans.hooyberghs@vito.be)
"""

# first initialize the run
from DownscalingClass import DownscalingClass
import datetime

############################
#Compulsory input paramters. 
############################




# outputfiles
fileEmissions='D:/ifdm/IFDMwerkfolder/Ping/Emissions_downscaling.xlsx'
fileMeteo='D:/ifdm/IFDMwerkfolder/Ping/Meteo_downscaling.xlsx'
file_NewRaster='D:/ifdm/IFDMwerkfolder/Ping/New_InputRaster.txt'
file_CompleteGrid='D:/ifdm/IFDMwerkfolder/Ping/CompleteGrid.txt'

#input variables
pollutants=["NOX","BC","PM25","PM10"]
IFDMdata='D:/ifdm/IFDMwerkfolder/Ping/IFDMinput/'
ProxyData='D:/ifdm/IFDMwerkfolder/DataChina/'
Domains=['D:/ifdm/IFDMwerkfolder/Ping/pingdingshan_domain.shp']



saga='D:/saga_3.0.0_x64/saga_cmd.exe'
gdal='D:/gdal/gdal/apps/gdalwarp.exe'
ogr='D:/gdal/gdal/apps/ogr2ogr.exe'

startdate=datetime.datetime(2017,01,01,00)
enddate=datetime.datetime(2017,12,31,23)
meteoname='meteo2017_ERAint.nc'


############################
# Actual routines. 
############################
Run = DownscalingClass(Domains,ProxyData,fileEmissions,fileMeteo,file_NewRaster,
    file_CompleteGrid,IFDMdata,pollutants,saga,gdal,ogr,startdate,enddate,meteoname)
Run.GISprocessing()    
Run.convertEmissionsToExcel()
Run.convertMeteoToExcel()
Run.MakeRaster()
