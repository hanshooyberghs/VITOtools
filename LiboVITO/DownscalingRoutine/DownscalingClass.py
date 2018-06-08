# -*- coding: utf-8 -*-
"""
Preprocessing routine for the use of IFDM in China.
===================================================

@author: hooyberh (hans.hooyberghs@vito.be)

Basic description:
Starting from (MarcoPolo-)emissions and general proxy data for China, the module 
constructs the input emissions for IFDM. The meteorological input required by IFDM
is composed starting from NetCDF data. 


Input parameters for constructor:
* domain: link to shape file with the domain that will be studied
* inputfolder_proxy: link to the folder with Emissions and proxy data for China
    IMPORTANT: this folder should have the correct structure, otherwith you have
        to modify the input parameters in the constructor
* outputEmissions: file in which Emissions are stored
* outputMeteo: file in which Meteo is stored
* tempfolder: Temporary folder used by the program
* pollutants: list of pollutants that should be treated
    available are (case-sensitive!!):
        * NOX
        * SO2
        * PM25
        * PM10
        * BC

Routines:
* Constructor
 *convertEmissionsToExcel
    Make Excel Emissions (uses the GISProcessing routine)
* convertMeteoToExcel
    Make Excel Meteo
* GISProcessing
    GIS processing of the Emissions and proxies
* Auxiliary routines:
    * getCenterCoordinates: aux routine to get the center coordinates
    * CalculateAreaCell: aux routine to calculate area of a grid cell    
    
Dependencies:
* Python routines:
    -os
    -glob
    -pandas
    -netCDF4
    -simpledbf
    -datetime
    -calendar
    -math
* External software:
    - SAGA-GIS
    - GDAl
        
Output files: 2 Excel files, containing
* Meteo input for IFDM preprocessor
* Emissions input for IFDM preprocessor
These files can readily be processed by the Python IFDM Preprocessor.

Other input parameters: see constructor. Only change these if you know what you
are doing (don't tell me that I didn't warn you ;) )

===

* to do:
    * correct implementation of domain size
    * treatment of center coordinates

"""
import os
import glob
import pandas
from netCDF4 import Dataset
from simpledbf import Dbf5
import datetime
import calendar
import math
import datetime
import shutil
import sys
pandas.options.mode.chained_assignment = None

class DownscalingClass:
    def __init__(self,domainfiles,inputfolder_proxy,outputEmissions,outputMeteo,outputnewraster,
        outputCompleteGrid,IFDMdatafolder,pollutants,saga,gdal,ogr,startdate,enddate,meteo_name):        
		# constants
        self.resolution_population = 0.01 #in degrees, the output resolution for residential proxy
        self.resolution_globcover = 0.01 #in degrees, the output resolution for globcover proxy
        self.meteoStart=startdate
        self.meteoEnd=enddate
        self.year=self.meteoStart.year
        self.heightBridge=8 # height of a bridge   
        
        # output files
        self.output_emissions=outputEmissions
        self.output_meteo=outputMeteo
        self.output_NewRaster=outputnewraster
        self.output_CompleteGrid=outputCompleteGrid
        
        # input files
        self.input_domainfiles=domainfiles
        self.input_roads=inputfolder_proxy+'OpenStreetMapRoads/OpenStreetMap.shp'
        self.input_roads_constants=inputfolder_proxy+'OpenStreetMapRoads/RoadValues.xlsx'
        self.input_globcover=inputfolder_proxy+'GlobCover/Globcover.tif'
        self.input_population=inputfolder_proxy+'PopulationData/WorldPopulation.tif'
        self.input_powerplants=inputfolder_proxy+'Powerplants/Powerplants_Enipedia.shp'
        self.input_emissions=inputfolder_proxy+'Emissions_MarcoPolo/'   #use MarcoPolo or MEIC?
        self.input_defaultfile=inputfolder_proxy+'Emissions/DefaultValuesEmissions.xlsx'
        self.input_meteo=inputfolder_proxy+'Meteo/'+meteo_name
        self.input_pollutants=pollutants

        
        # temporary folder
        self.IFDMfolder=IFDMdatafolder
        print self.IFDMfolder
        print IFDMdatafolder
        print  os.path.exists(self.IFDMfolder)
        if not os.path.exists(self.IFDMfolder):
            print 'jan'
            os.makedirs(self.IFDMfolder)
        self.tempfolder=self.IFDMfolder+'PP_china/'
        if os.path.exists(self.tempfolder):
            shutil.rmtree(self.tempfolder)
        os.makedirs(self.tempfolder)
        self.input_domain=self.tempfolder+'Domain.shp'

        self.gdal=gdal
        self.saga=saga
        self.ogr=ogr
        
        #temporary files
        self.outputfile_emissions=self.tempfolder+'/Emissions.shp'
        
        # find coordinates of centroid
        self.getCenterCoordinates()
        
        # dete ine size of the grid cells (in m2)
        self.cell_area_input=self.CalculateAreaCell(0.25)
        
        # check files
        for domain_i in self.input_domainfiles:
            self.checkPresence(domain_i)
        self.checkPresence(inputfolder_proxy)

    
    def convertEmissionsToExcel(self):

        
        print "==========================================="
        print 'Processing of the emissions, part II'
        print 'Convert to Excel'
        print "==========================================="
        ####
        # Define basic input
        ####
        DefaultValuesLineSources=pandas.read_excel(self.input_defaultfile,sheetname='LineSources')
        DefaultValuesPointSources=pandas.read_excel(self.input_defaultfile,sheetname='PointSources')
        DefaultValuesSurfaceSources=pandas.read_excel(self.input_defaultfile,sheetname='SurfaceSources')
        #define order
        temp=['Lat1','Lon1','Lat2','Lon2','Width (m)','Height (m)','Spread (m)']
        temp.extend(self.input_pollutants)
        temp.extend(['Wind/Temp dependence','Time factor','NO/NOx ratio','Type','Sector'])
        ExpectedColumnsLineSources=temp 
           
        temp=['Lat','Lon','Height (m)','Diameter (m)','Temperature (K)','Volume stream (Nm3/s)']
        temp.extend(self.input_pollutants)
        temp.extend(['Wind/Temp dependence','Time factor','NO/NOx ratio','Type','Sector'])
        ExpectedColumnsPointSources=temp
        
        temp=['Lat1','Lon1','Lat2','Lon2']
        temp.extend(self.input_pollutants)
        temp.extend(['Height (m)','Spread (m)','Wind/Temp dependence','Time factor','NO/NOx ratio','Type','Sector'])
        ExpectedColumnsSurfaceSources=temp
        
        temp=['Lat','Lon','Area (m2)']
        temp.extend(self.input_pollutants)
        temp.extend(['Height (m)','Spread (m)','Wind/Temp dependence','Time factor','NO/NOx ratio','Type','Sector'])
        ExpectedColumnsAltSurfaceSources=temp
        
        
        ####
        # Make dataframes
        ####
        # convert emissions to database object
        dbf=Dbf5(self.outputfile_emissions.replace('shp','dbf'))
        emissions=dbf.to_dataframe()
        emissions.columns=[w.split(' ')[0] for w in emissions.columns]    
        
        # convert residential sources to database object (and delete rows without population)
        residential=pandas.DataFrame()
        for i in range(len(self.input_domainfiles)):
			file_in=self.file_population_points.replace('.shp',str(i)+'.dbf')
			if os.path.exists(file_in):
				dbf=Dbf5(file_in)
				residential_in=dbf.to_dataframe()
				residential_in=residential_in[residential_in[residential_in.columns[-1]]>0]
				residential_in['Sector']='Residential'+str(i)
				residential=pandas.concat([residential,residential_in],ignore_index=True)
        

        temp=residential[residential.columns[-2]]/sum(residential[residential.columns[-2]])*len(residential.index) 
        
        residential['fraction']=temp
        
        # convert industrial sources to database object (and delete rows without built-up areas)
        industrial=pandas.DataFrame()
        for i in range(len(self.input_domainfiles)):
			file_in=self.file_Globcover_points.replace('.shp',str(i)+'.dbf')
			if os.path.exists(file_in):
				dbf=Dbf5(file_in)
				industrial_in=dbf.to_dataframe()
			   # industrial_in=industrial_in[industrial_in[industrial_in.columns[-1]]>0]
				industrial_in['Sector']='Industry'+str(i)
				industrial=pandas.concat([industrial,industrial_in],ignore_index=True)
				
        
        # convert road map to database 
        ####object (and delete rows without built-up areas)
        roads=pandas.DataFrame()
        for i in range(len(self.input_domainfiles)):
			file_in=self.outputfile_roads.replace('.shp',str(i)+'.dbf')
			if os.path.exists(file_in):
				dbf=Dbf5(file_in)
				roads_in=dbf.to_dataframe()
				roads_in['Sector']='Traffic'+str(i)
				roads=pandas.concat([roads,roads_in],ignore_index=True)
        
        
        # convert powerplant sources to database object
        powerplants=pandas.DataFrame()
        for i in range(len(self.input_domainfiles)):
            file_in=self.outputfile_powerplants.replace('.shp',str(i)+'.dbf')
            if os.path.exists(file_in):
                dbf=Dbf5(file_in)
                powerplants_in=dbf.to_dataframe()
                temp=powerplants_in['OutputMWh']/sum(powerplants_in['OutputMWh'])
                powerplants_in['fraction']=temp
                powerplants_in['Sector']='Energy'+str(i)
                powerplants=pandas.concat([powerplants,powerplants_in],ignore_index=True)
        
        ####
        # Add all data per sector
        ####
        # make output database
        surfaceSources_industry=pandas.DataFrame()
        surfaceSources_residential=pandas.DataFrame()
        surfaceSources=pandas.DataFrame()
        lineSources=pandas.DataFrame()
        pointSources=pandas.DataFrame()       
       
        
        #residential
        if residential.empty==False:
            surfaceSources_residential['Lat1']=residential['Y']-self.resolution_population/2
            surfaceSources_residential['Lat2']=residential['Y']+self.resolution_population/2
            surfaceSources_residential['Lon1']=residential['X']-self.resolution_population/2
            surfaceSources_residential['Lon2']=residential['X']+self.resolution_population/2
            surfaceSources_residential['Sector']=residential['Sector']
            surfaceSources_residential['Height (m)']=DefaultValuesSurfaceSources['Height (m)']['Residential']
            surfaceSources_residential['Spread (m)']=DefaultValuesSurfaceSources['Spread (m)']['Residential']
            surfaceSources_residential['Wind/Temp dependence']=DefaultValuesSurfaceSources['Wind/Temp Dependence']['Residential']
            surfaceSources_residential['Time factor']=DefaultValuesSurfaceSources['Time factor']['Residential']
            surfaceSources_residential['NO/NOx ratio']=DefaultValuesSurfaceSources['NO/Nox ratio']['Residential']
            surfaceSources_residential['Type']='Latlon'
            for pol in self.input_pollutants:
                surfaceSources_residential[pol]=0          
                em_res=emissions['Res_'+pol][0]*1000000/(365*24)/self.cell_area_input
                surfaceSources_residential[pol]=residential['fraction']*em_res
        else:
            surfaceSources_residential=pandas.DataFrame(columns=ExpectedColumnsSurfaceSources)
        
        
        #industry
        if industrial.empty==False:
            temp=industrial[industrial.columns[-2]]/sum(industrial[industrial.columns[-2]])
            industrial['fraction']=temp*len(industrial.index)
            
            surfaceSources_industry['Lat1']=industrial['Y']-self.resolution_globcover/2
            surfaceSources_industry['Lat2']=industrial['Y']+self.resolution_globcover/2
            surfaceSources_industry['Lon1']=industrial['X']-self.resolution_globcover/2
            surfaceSources_industry['Lon2']=industrial['X']+self.resolution_globcover/2
            surfaceSources_industry['Sector']=industrial['Sector']
            surfaceSources_industry['Height (m)']=DefaultValuesSurfaceSources['Height (m)']['Industry']
            surfaceSources_industry['Spread (m)']=DefaultValuesSurfaceSources['Spread (m)']['Industry']
            surfaceSources_industry['Wind/Temp dependence']=DefaultValuesSurfaceSources['Wind/Temp Dependence']['Industry']
            surfaceSources_industry['Time factor']=DefaultValuesSurfaceSources['Time factor']['Industry']
            surfaceSources_industry['NO/NOx ratio']=DefaultValuesSurfaceSources['NO/Nox ratio']['Industry']
            surfaceSources_industry['Type']='Latlon'
            for pol in self.input_pollutants:
                surfaceSources_industry[pol]=0
                em_pol=emissions['Ind_'+pol][0]*1000000/(365*24)/self.cell_area_input/2.0
                #industry: spread 50% uniform, and 50% according to proxy
                surfaceSources_industry[pol]=industrial['fraction']*em_pol+em_pol
                
        else:
            surfaceSources_industry=pandas.DataFrame(columns=ExpectedColumnsSurfaceSources)

        
        
        
        #powerplants
        if powerplants.empty==False:
            pointSources['Lat']=powerplants['lat']
            pointSources['Lon']=powerplants['lon']
            pointSources['Sector']=powerplants['Sector']
            pointSources['Height (m)']=DefaultValuesPointSources['Height (m)']['Energy']
            pointSources['Diameter (m)']=DefaultValuesPointSources['Diameter (m)']['Energy']
            pointSources['Temperature (K)']=DefaultValuesPointSources['Temperature (K)']['Energy']
            pointSources['Volume stream (Nm3/s)']=DefaultValuesPointSources['Volume stream (Nm3/s)']['Energy']
            pointSources['Wind/Temp dependence']=DefaultValuesPointSources['Wind/Temp Dependence']['Energy']
            pointSources['Time factor']=DefaultValuesPointSources['Time factor']['Energy']
            pointSources['NO/NOx ratio']=DefaultValuesPointSources['NO/Nox ratio']['Energy']
            pointSources['Type']='Latlon'      
            for pol in self.input_pollutants:
                pointSources[pol]=0
                pointSources[pol]=powerplants['fraction']*emissions['Ene_'+pol][0]/(365*24)*1000
        else:
            pointSources=pandas.DataFrame(columns=ExpectedColumnsPointSources)
  
        
        
        #roads       
        if roads.empty==False:
            lineSources['Lat1']=roads['y1']
            lineSources['Lat2']=roads['y2']
            lineSources['Lon1']=roads['x1']
            lineSources['Lon2']=roads['x2']
            lineSources['Sector']=roads['Sector']
            lineSources['Height (m)']=DefaultValuesLineSources['Height (m)']['Traffic']
            lineSources['Height (m)'][roads.bridge==1]=lineSources['Height (m)'][roads.bridge==1]+self.heightBridge
            lineSources['Wind/Temp dependence']=DefaultValuesLineSources['Wind/Temp Dependence']['Traffic']
            lineSources['Spread (m)']=DefaultValuesLineSources['Spread (m)']['Traffic']
            lineSources['Time factor']=DefaultValuesLineSources['Time factor']['Traffic']
            lineSources['NO/NOx ratio']=DefaultValuesLineSources['NO/Nox ratio']['Traffic']
            lineSources['Type']='Latlon'
            
            DefaultValuesRoads=pandas.read_excel(self.input_roads_constants)
            lineSources['Width (m)']=0
            roads['weighted_length']=0
            road_types=DefaultValuesRoads['factor'].index
            for i in road_types:
                lineSources['Width (m)'][roads.type==i]=DefaultValuesRoads['width'][i]
                roads['weighted_length'][roads.type==i]=roads['length_km'][roads.type==i]*DefaultValuesRoads['factor'][i]
            sum_weighted_length_total=sum(roads['weighted_length'])  
            fraction_road_type=pandas.DataFrame(columns=road_types,index=['0'])
            sum_length_type=pandas.DataFrame(columns=road_types,index=['0'])
            for i in road_types:
                fraction_road_type[i]=sum(roads['weighted_length'][roads['type']==i])/sum_weighted_length_total
                sum_length_type[i]=sum(roads['length_km'][roads['type']==i])
                
            for pol in self.input_pollutants:
                lineSources[pol]=0
                for i in road_types:
                    lineSources[pol][roads.type==i]=fraction_road_type[i][0]*emissions['Tra_'+pol][0]*1000/(365*24)/sum_length_type[i][0]

            lineSources=lineSources[lineSources['NOX']>0]
        
        else:
            lineSources=pandas.DataFrame(columns=ExpectedColumnsLineSources)	
        
        ####
        # Save data
        ####        
        # merge surface sources
        surfaceSources=pandas.concat([surfaceSources_industry,surfaceSources_residential],ignore_index=True)      
        
        #set in correct order
        surfaceSources=surfaceSources[ExpectedColumnsSurfaceSources]
        lineSources=lineSources[ExpectedColumnsLineSources]
        pointSources=pointSources[ExpectedColumnsPointSources]
        
        #make empty alternative surface sources
        altSurfaceSources=pandas.DataFrame(columns=ExpectedColumnsAltSurfaceSources)
        
        #save data
        if os.path.exists(self.output_emissions):
            os.remove(self.output_emissions)
        writer = pandas.ExcelWriter(self.output_emissions)
        surfaceSources.to_excel(writer,'SurfaceSources',index=False)
        lineSources.to_excel(writer,'LineSources',index=False)
        pointSources.to_excel(writer,'PointSources',index=False)
        altSurfaceSources.to_excel(writer,'AltSurfaceSources',index=False)
        writer.save()
        writer.close()
        
    
    def convertMeteoToExcel(self):
        ####
        #Processing of the ERA-interim meteo data
        ####
        print "==========================================="
        print 'Processing of Meteo data'
        print "==========================================="
                
        meteo_file=Dataset(self.input_meteo)
        temperature_input=pandas.Panel(meteo_file.variables['t2m'][:])
        uwind_input=pandas.Panel(meteo_file.variables['u10'][:])
        vwind_input=pandas.Panel(meteo_file.variables['v10'][:])
        
        lat=pandas.Series(meteo_file.variables['latitude'][:])  
        loc_lat=lat[self.lat_center_domain>lat].index[0] #meenemen: loc_lat en loc_lat-1
        lon=pandas.Series(meteo_file.variables['longitude'][:])
        loc_lon=lon[self.lon_center_domain>lon].index[-1] #meenemen: loc_lon en loc_lon+1
        
        ic=(lon[loc_lon+1]-self.lon_center_domain)/(lon[loc_lon+1]-lon[loc_lon])
        jc=(lat[loc_lat-1]-self.lat_center_domain)/(lat[loc_lat-1]-lat[loc_lat])
        
        output_temperature=(1-ic)*(1-jc)*temperature_input[:,loc_lat-1,loc_lon+1]+ \
            ic*(1-jc)*temperature_input[:,loc_lat-1,loc_lon]+ \
            (1-ic)*jc*temperature_input[:,loc_lat+1,loc_lon+1]+ \
            ic*jc*temperature_input[:,loc_lat,loc_lon]
        
        uwind_output=(1-ic)*(1-jc)*uwind_input[:,loc_lat-1,loc_lon+1]+ \
            ic*(1-jc)*uwind_input[:,loc_lat-1,loc_lon]+ \
            (1-ic)*jc*uwind_input[:,loc_lat+1,loc_lon+1]+ \
            ic*jc*uwind_input[:,loc_lat,loc_lon]
            
        vwind_output=(1-ic)*(1-jc)*vwind_input[:,loc_lat-1,loc_lon+1]+ \
            ic*(1-jc)*vwind_input[:,loc_lat-1,loc_lon]+ \
            (1-ic)*jc*vwind_input[:,loc_lat+1,loc_lon+1]+ \
            ic*jc*vwind_input[:,loc_lat,loc_lon]
        
        output=pandas.DataFrame()
        output['WE wind speed (m/s)']=uwind_output
        output['NS wind speed (m/s)']=vwind_output
        output['Temperature (K)']=output_temperature;
        output['Lat']=self.lat_center_domain
        output['Lon']=self.lon_center_domain
        output['Height (m)']=10
        output['Time parameter']=0
        print 'Assembled MeteoData'
        
        # set date
        output.index=pandas.date_range(start=datetime.datetime(self.year,01,01,0),\
			periods=len(output),freq='H')
        output['Year']=output.index.year
        output['Month']=output.index.month
        output['Day']=output.index.day
        output['Hour']=output.index.hour
        print 'Added time data'
        
        # select correct period (default: complete year)
        output=output[self.meteoStart:self.meteoEnd]

        
        # put in correct order and save
        output=output[['WE wind speed (m/s)','NS wind speed (m/s)','Height (m)','Temperature (K)',\
            'Year','Month','Day','Hour','Lat','Lon','Time parameter']]
        writer = pandas.ExcelWriter(self.output_meteo)
        output.to_excel(writer,'Meteo',index=False)
        writer.save()
        writer.close()
        meteo_file.close()
        print 'Finished'
        
        
    
    def GISprocessing(self):
        print "==========================================="
        print 'Processing of the emissions, part I'
        print 'Start geo-preprocessing using gdal and saga'
        print "==========================================="
        
        
                
        print "Emissions for complete domain"
        inputfiles_emissions=glob.glob(self.input_emissions+'*.sgrd')
        filelist='"'
        for file in inputfiles_emissions:
            filelist+=file
            filelist+=';'
        filelist=filelist[:-1]+'"'
        commando=self.saga+' shapes_grid 2 -POLYGONS '+self.input_domain+' -GRIDS '+filelist+' -RESULT '+ \
            self.outputfile_emissions +' -SUM 0 -COUNT 0 -MIN 0 -MAX 1 -RANGE 0 -MEAN 0 -VAR 0 -STDDEV 0 -METHOD 3'
        os.system(commando)
        print "GIS processing for entire domain"
        print "Industrial proxy: GlobCover"
        #Globcover (industrial proxy)
        commandogdal_basis=self.gdal+'  -crop_to_cutline -t_srs "EPSG:4326" -r near -cutline '+self.input_domain+' ' 
        outputfile_Globcover=self.tempfolder+'/Globcover.tif'
        commando_Globcover=commandogdal_basis+self.input_globcover+' '+outputfile_Globcover
        print commando_Globcover
        os.system(commando_Globcover) #clipping
        outputfile_Globcover_resampled=outputfile_Globcover.replace('.tif','_resampled.tif')
        commandogdal_basis=self.gdal+' -r average -tr '
        commando= commandogdal_basis+ str(self.resolution_globcover)+' '+str(self.resolution_globcover) \
            + ' ' + outputfile_Globcover + ' ' +outputfile_Globcover_resampled
        os.system(commando) #regridding
        file_Globcover_saga=outputfile_Globcover_resampled.replace('.tif','.sgrd')
        commando=self.saga+'  io_gdal 0 -FILES '+outputfile_Globcover_resampled+' -GRIDS '+file_Globcover_saga
        os.system(commando) #convert to saga
        self.file_Globcover_points=file_Globcover_saga.replace('.sgrd','.shp')
        commando=self.saga+'  shapes_grid 3 -GRIDS '+file_Globcover_saga+' -SHAPES '+self.file_Globcover_points
        os.system(commando) #convert to points shape
        
        print "Population proxy"
        #Population
        outputfile_population=self.tempfolder+'/Population.tif'
        commandogdal_basis=self.gdal+'  -crop_to_cutline -t_srs "EPSG:4326" -r near -cutline '+self.input_domain+' ' 
        commando_population=commandogdal_basis+self.input_population+' '+outputfile_population
        os.system(commando_population)  #clipping
        commandogdal_basis=self.gdal+' -r average -tr '
        outputfile_population_resampled=outputfile_population.replace('.tif','_resampled.tif')
        commandogdal_basis=self.gdal+' -r average -tr '
        commando= commandogdal_basis+ str(self.resolution_population)+' '+str(self.resolution_population) \
            + ' ' + outputfile_population + ' ' +outputfile_population_resampled
        os.system(commando) #regridding
        file_population_saga=outputfile_population_resampled.replace('.tif','.sgrd')
        commando=self.saga+'  io_gdal 0 -FILES '+outputfile_population_resampled+' -GRIDS '+file_population_saga
        os.system(commando) #convert to saga
        self.file_population_points=file_population_saga.replace('.sgrd','.shp')
        commando=self.saga+'  shapes_grid 3 -GRIDS '+file_population_saga+' -SHAPES '+self.file_population_points
        os.system(commando) #convert to points shape
        
        print "Road map: Open Street Map"
        print "Powerplants: Enipedia"      
        self.outputfile_powerplants=self.tempfolder+'/Powerplants.shp'
        self.outputfile_roads=self.tempfolder+'/Roads.shp'        
        commando= self.ogr+' -clipsrc '+self.input_domain+' '+ self.outputfile_powerplants+' '+ self.input_powerplants
        os.system(commando)
        
        commando= self.ogr+' -clipsrc '+self.input_domain+' '+ self.outputfile_roads+' '+ self.input_roads
        os.system(commando)           
            
            
        print "Clipping on subdomains"
        for i in range(len(self.input_domainfiles)):
            #globcover
            commando=self.saga+'  shapes_points 8 -POINTS '+self.file_Globcover_points+' -POLYGONS '\
                + self.input_domainfiles[i]+' -CLIPS '+self.file_Globcover_points.replace('.shp',str(i)+'.shp')
            os.system(commando)
            
            #population
            commando=self.saga+'  shapes_points 8 -POINTS '+self.file_population_points+' -POLYGONS '\
                + self.input_domainfiles[i]+' -CLIPS '+self.file_population_points.replace('.shp',str(i)+'.shp')
            os.system(commando)
            
            #powerplants
            commando=self.saga+'  shapes_points 8 -POINTS '+self.outputfile_powerplants+' -POLYGONS '\
                + self.input_domainfiles[i]+' -CLIPS '+self.outputfile_powerplants.replace('.shp',str(i)+'.shp')
            os.system(commando)
            
            #roads
            commando=self.saga+'  shapes_polygons 11 -M_INPUT "" -S_INPUT '+self.outputfile_roads+' -CLIP '\
                + self.input_domainfiles[i]+' -S_OUTPUT '+self.outputfile_roads.replace('.shp',str(i)+'.shp')\
                +' -MULTIPLE 0'
            print commando
            os.system(commando)
    
    def MakeRaster(self):
        print "Make Raster"
        
        # location of default IFDM input files
        self.inputfileTimeSeries=self.IFDMfolder+'Input_Timeseries.txt'   
        self.inputfileExecutable=self.IFDMfolder+'ifdmmodel_windows.exe'
        
        # location of default input config files
        self.inputfileConfig=self.IFDMfolder+'ifdm_add.conf'

        # location of default input raster parameter files
        self.inputfileRasterParameters=self.IFDMfolder+'InputRaster.txt'
        self.inputfileRasterParameters2=self.IFDMfolder+'InputRaster2_point.txt'
        self.inputfileRasterParameters3=self.IFDMfolder+'InputRaster3.txt'

        # make raster folder
        self.RasterFolder=self.IFDMfolder+'Rasterfolder/'
        if not os.path.exists(self.RasterFolder):
            os.makedirs(self.RasterFolder)


        # check presence of all input files        
        self.checkPresence(self.inputfileConfig)
        self.checkPresence(self.inputfileTimeSeries)
        self.checkPresence(self.inputfileExecutable)    
        self.checkPresence(self.inputfileRasterParameters)
        self.checkPresence(self.inputfileRasterParameters2)
        self.checkPresence(self.inputfileRasterParameters3)
        
        # make rasterfolder versions
        self.tempfileConfig=self.RasterFolder+'ifdm.conf'
        self.tempfileRasterParameters=self.RasterFolder+'InputRaster.txt'
        
        """
        Changes to config and RasterParameters file
        """
        self.pollutantsIFDM = [w.ljust(4,'a') for w in self.input_pollutants] 
        # make list for config file
        self.pollutantList=''        
        for i in range(len(self.pollutantsIFDM)):
            self.pollutantList=self.pollutantList+self.pollutantsIFDM[i]+','
        self.pollutantList=self.pollutantList.rstrip(',')
        self.modifyFile(self.inputfileConfig,self.tempfileConfig,'pollutantslist',self.pollutantList)
        print 'Modified config file to include pollutants'
        
        filereader=open(self.tempfileConfig,'a')
        filereader.write('\n IFDM_RASTER_Create_Grid \t = \t 1')
        filereader.close()
        print 'Modified config file to start RasterMakingRoutine'
        
        centercoordinates=str(self.lat_center_domain)+'\t'+str(self.lon_center_domain)
        self.modifyFile(self.inputfileRasterParameters,self.tempfileRasterParameters,
                'center_coordinates',centercoordinates)  
        print 'Modified InputRasterParameter to include center coordinates'
              
            
            
        """
        Make emissions input 
        """
        xls = pandas.ExcelFile(self.output_emissions)
        sheets = xls.sheet_names
        
        if 'LineSources' in sheets:
            LineSources=pandas.read_excel(self.output_emissions,sheetname='LineSources')
            tempfileEmissions=self.RasterFolder+'InputLinesources.txt'
            number_rows=len(LineSources.index)
            f=open(tempfileEmissions,'w')
            f.write('%d' % number_rows)
            f.write('\n')
            f.close()  
            LineSources.to_csv(tempfileEmissions,sep='\t',header=False,index=False,mode='a')           
            
        else:
            print 'Problem: no line sources'
            print 'I will stop'
            sys.exit()
            
        
        if 'PointSources' in sheets:
            PointSources=pandas.read_excel(self.output_emissions,sheetname='PointSources')
            tempfileEmissions=self.RasterFolder+'InputPointsources.txt'
            number_rows=len(PointSources.index)
            f=open(tempfileEmissions,'w')
            f.write('%d' % number_rows)
            f.write('\n')
            f.close()  
            PointSources.to_csv(tempfileEmissions,sep='\t',header=False,index=False,mode='a')           
            
        else:
            print 'No point sources'
            print 'I will only build a line source following grid'
            filereader=open(self.tempfileConfig,'a')
            filereader.write('\n IFDM_RASTER_Raster_type \t = \t 10')
            filereader.close()
            
        print 'Added emissions to RasterFolder'
        
        """
        Copy files to raster folder
        """
        shutil.copy(self.inputfileRasterParameters3,self.RasterFolder)
        shutil.copy(self.inputfileRasterParameters2,self.RasterFolder)
        shutil.copy(self.inputfileTimeSeries,self.RasterFolder)
        shutil.copy(self.inputfileExecutable,self.RasterFolder+'ifdm.exe')
        print 'Copied files to RasterFolder'
        print 'Raster folder is located at '+self.RasterFolder
        
        
        print 'Raster will be made'
        print 'This may take some time'
        commando='cd '+self.RasterFolder+' & ifdm.exe'
        print commando
        os.system(commando)
        
        print '\n \n Made Raster'
        shutil.copy(self.RasterFolder+'CompleteGrid.txt',self.output_CompleteGrid)
        shutil.copy(self.RasterFolder+'New_InputRaster.txt',self.output_NewRaster)
        
        

    def getCenterCoordinates(self):
        #routine to find the center coordinates
        
        print "Merge all subdomains"
        
        # Merge and dissolve if necessary
        if len(self.input_domainfiles)>1:    

            # buffer (no gaps)
            for i in range(len(self.input_domainfiles)):
                commando=self.saga+'  shapes_tools 18 -SHAPES '+self.input_domainfiles[i]+' -BUFFER ' +\
                    self.input_domainfiles[i].replace('.shp','_buffered.shp')+' -DIST_FIELD_DEFAULT 0.01'
                os.system(commando)

        
            commando=self.saga+'  shapes_tools 2 -MERGED '+self.input_domain.replace('.shp','_merged.shp')+' -INPUT "' \
                +self.input_domainfiles[0].replace('.shp','_buffered.shp')+ ';'            
            for i in range(1,len(self.input_domainfiles)):
                commando+=self.input_domainfiles[i].replace('.shp','_buffered.shp')
                commando+=';'
            
            commando = commando[:-1]
            commando = commando+'"'
            print commando
            os.system(commando)        

            #Dissolve (first self-interaction to get rid of overlaps)
            commando=self.saga+'  shapes_polygons 12 -POLYGONS '+self.input_domain.replace('.shp','_merged.shp')+ \
                ' -INTERSECT '+self.input_domain.replace('.shp','_selfint.shp')
            print commando
            os.system(commando)
            
            
            commando=self.saga+'  shapes_polygons 5 -POLYGONS '+self.input_domain.replace('.shp','_selfint.shp')+' -DISSOLVED '+self.input_domain 
            print commando
            os.system(commando)
            
        else:
            self.input_domain=self.input_domainfiles[0]    

        print "Get center coordinates"
        centroids_shape=self.tempfolder+'centroids.shp'
        print self.input_domain
        commando=self.saga+'  shapes_polygons 1 -POLYGONS '+self.input_domain+' -CENTROIDS '+centroids_shape
        os.system(commando)
        
        commando=self.saga+'  shapes_points 6 -INPUT ' +centroids_shape
        os.system(commando)
        
        
        dbf=Dbf5(centroids_shape.replace('shp','dbf'))
        print 'Reading '+centroids_shape
        centroids_data=dbf.to_dataframe()

        self.lat_center_domain=(centroids_data['Y'][0]).item()
        self.lon_center_domain=(centroids_data['X'][0]).item()
        
        print 'lon-coordinate of centroid: '
        print self.lon_center_domain
        print 'lat-coordinate of centroid: '
        print self.lat_center_domain

    def CalculateAreaCell(self,resolution):
        #calculate size of grid cell in m2
        Re=6371000 #radius of earth in m
        
        phi_2=(self.lat_center_domain+resolution/2)*math.pi/180
        phi_1=(self.lat_center_domain-resolution/2)*math.pi/180
        
        area=Re**2.0*math.acos(math.sin(phi_1)*math.sin(phi_2)+math.cos(phi_1)*math.cos(phi_2)) * \
            math.acos(math.sin(phi_1)**2+math.cos(phi_1)**2*math.cos(resolution*math.pi/180))
        return area

        
    def checkPresence(self,file_in):        
        if not os.path.exists(file_in):
            print 'Computer says "no".'
            print 'Input file not found error'
            print 'File '+file_in+' does not exist'            
            print 'Exiting program...'
            sys.exit()

            
    def modifyFile(self,inputFile,outputFile,text_in,text_out):        
        file_in=open(inputFile)        
        inputtext=file_in.read()
        file_in.close()        
        outputtext=inputtext.replace(text_in,text_out)
        filesave=open(outputFile,'w')
        filesave.write(outputtext)
        filesave.close()
