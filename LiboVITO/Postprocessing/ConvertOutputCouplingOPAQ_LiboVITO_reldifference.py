"""
Postprocessing routine for IFDM, adaptation for LiboVITO-system
===============================================================

Basic description: postprocessing of IFDM using Python and GDAL


Python-Dependencies:
* Python routines:
    -os
    -glob
    -pandas
    -shutil
    -simpledbf
    -numpy
    -matplotlib
* External software:
    - gdal

"""

############################
# Input parameters
############################



Runfolder='/usr/ifdm/czb/pingdingshan/run_base/run/'
Outputfolder='/usr/ifdm/czb/pingdingshan/run_base/results/'
Outputfolder0='/usr/ifdm/czb/pingdingshan/run_base/results/'
Pollutants=["NOX","BC","PM25","PM10","SO2"]
folder_OPAQ_results = '/usr/local/pingdingshan_opaq/pingdingshan_7CST/fc/output_RMSE/'    # folder where OPAQ results can be found
#date='20161222'    # date of the OPAQ run

gdal_grid='/usr/local/bin/gdal_grid'
gdal_translate='gdal_translate'
gdal_warp='gdalwarp'
gdal_calc='gdal_calc.py'

resolution=0.0005
x_min=113.06
x_max=113.47
y_min=33.65
y_max=33.84


#################################      
# Import modules
#################################

import os
import shutil
import gdal
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas
import numpy
import glob
import sys


######################################      
# Determine number of days and sectors
######################################

date=sys.argv[1]

sectorList=os.walk(Runfolder).next()[1]
print "Sectorlist: "
print sectorList    

list_in=os.walk(Runfolder+sectorList[0]).next()[1]
days_list=[i for i in list_in if i.isdigit()]
days_list2=['day'+i for i in list_in if i.isdigit()]
print "Dayslist: "
print days_list

Pollutants_orig=Pollutants
Pollutants = [w.ljust(4,'a') for w in Pollutants]    

if not os.path.exists(Outputfolder):
    os.makedirs(Outputfolder)

#sectorList=["AllSectors"]
    
#################################      
# Loops over all folders
#################################
for sector in sectorList:

    # empty variables for locations of measurement stations
    results_timeseries={}
    for pollutant in Pollutants:
        results_timeseries[pollutant]=pandas.DataFrame(columns=days_list2)


    for day in days_list:
        folder=Runfolder+sector+'/'+day+'/'
        print '\n'
        print day+' '+sector
        for pollutant in Pollutants:
            runfile=folder+'output'+pollutant+'_saga.txt'
            result=Outputfolder+pollutant+'_'+sector+'_day'+day+'.tif'   
            csv_file_end='results'+pollutant+'.csv'
            csv_file=folder+csv_file_end
            shutil.copy(runfile,csv_file)
            os.system("sed -i 's/;/,/g' "+ csv_file)
            os.system("sed -i 's/ //g' "+ csv_file)
            print str(x_min)+' '+str(x_max)+' '+str(y_min)+' '+str(y_max)
            x_diff=(x_max-x_min)/resolution
            y_diff=(y_max-y_min)/resolution

            # make vrt file (determines format for gdal_grid)
            vrtfile='outputformat'+pollutant+'.vrt'
            f= open(vrtfile, 'w')
            content='<OGRVRTDataSource> \n <OGRVRTLayer name="'+csv_file_end.replace('.csv','')+'"> \n <SrcDataSource>'+csv_file \
                +'</SrcDataSource> \n <GeometryType>wkbPoint</GeometryType> \n ' \
                +'<LayerSRS>EPSG:4326</LayerSRS> <GeometryField encoding="PointFromColumns" x="X" y="Y" z="' \
                +pollutant+'"/> \n </OGRVRTLayer> \n </OGRVRTDataSource>'
            f.write(content)
            f.close()

            # grid data
            if os.path.isfile(result):
                print 'file removed'
                os.remove(result)
            commando=gdal_grid+' -a LINEAR:radius=0:nodata=-9999 -txe '+str(x_min)+' '+str(x_max)+ \
                 ' -tye '+str(y_min)+' '+str(y_max)+\
                 ' -outsize '+str(x_diff)+' '+str(y_diff)+' -zfield values '+vrtfile+' '+result
            print commando
            os.system(commando)  
            
            # conversion to asc file
            # upside-down error         
            tmp_file=result.replace('.tif','_tmp.tif')  
            if os.path.isfile(tmp_file):
                print 'file removed'
                os.remove(tmp_file)            
            
            commando=gdal_warp+' '+result+' '+tmp_file
            print commando
            os.system(commando)
            
            
                        
            # read time series
            timeseries_file=folder+'output_time_series'+pollutant+'.txt'
            input_data_series=pandas.read_fwf(timeseries_file)
            input_data_series.columns=['value','longitude','latitude','dummy3','location']    
            mean_series_days=input_data_series.groupby('location').mean()
            results_timeseries[pollutant]['day'+day]=mean_series_days['value']*1e9
            
    # save time series data
    for pollutant in Pollutants:
        outputfile=Outputfolder+pollutant+'_'+sector+'_MeasStations.csv'  
        results_timeseries[pollutant].loc['mean']=results_timeseries[pollutant].mean(axis=0)
        results_timeseries[pollutant].to_csv(outputfile)
        
#################################      
# Coupling to OPAQ
#################################
background_total=pandas.DataFrame(columns=days_list2,index=Pollutants)
for i in range(len(Pollutants)):
    pollutant_OPAQ=Pollutants_orig[i].lower()
    if pollutant_OPAQ != 'bc':    
        file_OPAQ=folder_OPAQ_results+date+'/ovlfc_'+pollutant_OPAQ+'_'+date+'.txt'
        print file_OPAQ
        OPAQ=pandas.read_csv(file_OPAQ,sep='\t',comment='#',header=None)    
        OPAQ.columns=['location','date']+days_list2
    else:
        OPAQ=pandas.DataFrame(index=['AVG_URB'],columns=['location','date']+days_list2)
        OPAQ.location='AVG_URB'
        OPAQ.fillna(0,inplace=True)
    
    pollutant_IFDM=Pollutants[i]
    file_IFDM=Outputfolder0+pollutant_IFDM+'_AllSectors_MeasStations.csv' 
    IFDM=pandas.read_csv(file_IFDM)
    
    for day in days_list2:
        OPAQ_day=OPAQ[OPAQ.location=='AVG_URB'][day]
        IFDM_day=IFDM[IFDM.location=='mean'][day]
        
        # calculate background
        background=max(float(OPAQ_day)-float(IFDM_day),0)
        background_total[day][pollutant_IFDM]=background
        
        # gdal calc to add background to all sector files
        total_without_bg=Outputfolder+pollutant_IFDM+'_AllSectors_'+day+'.tif'   
        total_with_bg=Outputfolder+pollutant_IFDM+'_AllSectorsWithBG_'+day+'.tif'   
        
        # make new pictures
        if os.path.isfile(total_with_bg):
            print 'file removed'
            os.remove(total_with_bg)
        commando=gdal_calc+' -A '+total_without_bg +' --outfile='+total_with_bg \
               + ' --calc="A+'+str(background)+'"'+' --NoDataValue=-9999'
        print commando
        os.system(commando)
        
        
        
        # calculate relative difference
        file_basecase=Outputfolder0+pollutant_IFDM+'_AllSectorsWithBG_'+day+'.tif'   
        file_scenario=Outputfolder+pollutant_IFDM+'_AllSectorsWithBG_'+day+'.tif'  
        file_relativedifference=Outputfolder+pollutant_IFDM+'_AllSectorsWithBG_'+day+'_RelativeDifference.tif'  
        if os.path.isfile(file_relativedifference):
            print 'file removed'
            os.remove(file_relativedifference)
        commando=gdal_calc+' -A '+file_basecase +' -B  '+file_scenario\
            +' --outfile='+file_relativedifference + ' --calc="(B-A)/A*100" '\
            +'--NoDataValue=-9999'
        os.system(commando)
        
        # upside-down error            
        tmp_file=file_relativedifference.replace('.tif','_tmp.tif')
        if os.path.isfile(tmp_file):
            print 'file removed'
            os.remove(tmp_file)
        commando=gdal_warp+' '+file_relativedifference+' '+tmp_file
        print commando
        os.system(commando)
        

        
        
        
        

    
    
    
            
            

