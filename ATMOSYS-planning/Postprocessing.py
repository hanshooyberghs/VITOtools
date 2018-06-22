#####################################################
### POSTPROCESSING TOOL ATMOSYS SCENARIO ANALYSER ###
#####################################################


## USAGE ##
###########
#
# Input arguments:
#   * folder with IFDM runs (string). In the folder, the subfolder NO2, PM10 and PM25 with
#           IFDM-runs should be present
#   * Number of runs per pollutant (integer)
#   * Outputfolder. (string) Note: folder should exist 
#   * Pollutant: either NO2, PM10 or PM2.5 (string)
# example: python Postproccessing.py DemoOutput/ 6 outputfolder/ pollutant


## REQUIREMENTS ##
##################
#   python installation with following packages:
#       pandas, numpy, os, sys
#       Note: the installation in /tools/python/2.7.12/ is ok.
#
#   gdal, minimal version 2.1.1
#       Note: the version in /tools/gdal/2.1.1/ is ok.


##############
##############
## CONSTANT ##
##############
##############
#thresholds for exceedances
hour_limit_NO2=200
day_limit_pm10=50
day_limit_pm25=25

# programs
gdal_grid='gdal_grid'
#gdal_translate='gdal_translate'
#gdal_warp='gdalwarp'
#gdal_calc='gdal_calc.py'


# list of statistics to plot
statistics_output=['MeanConcentration','P9040']
# options: MeanConcentration','P50','P66','P80','P90','P9040','P95','P98','P99','P9979','exceedances']

# output resolution (in degrees)
resolution=0.0002


##############
##############
## PACKAGES ##
##############
##############
import pandas as pd
import os
import sys
import numpy as np
import gdal


############
############
## SET UP ##
############
############

# read input parameters
fp = open("postproc.config","r")
config = fp.read().splitlines()
fp.close()
folder_runs = config[0]
nr_runs = int(config[1])
folder_output = config[2]
pol = config[3]
print folder_runs
print nr_runs
print folder_output
print pol

#make output folder
if not os.path.exists(folder_output):
    os.makedirs(folder_output)

# read time stamps
timestamps=pd.Series()
for i in range(nr_runs):
    # read meteo file and select datatime
    file=folder_runs+'/NO2/'+str(i)+'/Meteotype1invoer.txt'
    data_in=pd.read_csv(file,sep=' ',header=None,skiprows=1)
    data_in=data_in[range(4,8)]     
    
    # convert to datetime
    data_in.columns=['Year','Month','Day','Hour']
    dt_in=pd.to_datetime(data_in)
    
    # add to timestamps
    timestamps=timestamps.append(dt_in,ignore_index=True)
ntime=len(timestamps)
ntime_perrun=int(ntime/nr_runs)


# read coordinates
file=folder_runs+'/NO2/0/uitvoerNO2a_atmosys.txt'    
data_in=pd.read_csv(file,sep=';',header=None)
data_in.columns=[['Lon','Lat','tmp']]
coordinates=data_in[['Lon','Lat']]
ncoords=len(coordinates)


# read locations of time series
# read locations of time series
TimeSeries=True
file=folder_runs+'/NO2/0/Tijdsreeksuitvoer.txt'    
with open(file, 'r') as f:
    first_line = int(f.readline())
if first_line == 0:
    TimeSeries=False
    print 'No time series found. No output for time series will be generated'
if TimeSeries:
    file=folder_runs+'/NO2/0/Tijdsreeksuitvoer.txt'    
    data_in=pd.read_fwf(file,skiprows=1,header=None)
    data_in.columns=[['Lon','Lat']]
    coordinates_timeseries=data_in[['Lon','Lat']]
    nseries=len(coordinates_timeseries)

    conc_timeseries=pd.DataFrame(index=timestamps,columns=range(nseries))


# read coordinate type
type='Latlon' # currently fixed, possible to read it from file in the future


###############################
###############################
## PROCESSING PER POLLLUTANT ##
###############################
###############################


pol_ifdm=pol.ljust(4,'a')   

## READ RAW IFDM OUTPUT ##
##########################

# loop over runs
dfs=[] # list of dataframes
cnt=0
print 'Read data: '+pol
for i in range(nr_runs):
    print 'Reading: '+pol+' '+str(i)
    folder=folder_runs+'/'+pol+'/'+str(i)+'/'
    # read concentrations: loop over situations
    for j in range(ntime_perrun):
        # reading
        file=folder+'uitvoer'+pol_ifdm+str(j+1).zfill(4)+'_atmosys.txt'
        data_in=pd.read_csv(file,header=None)
        data_in=data_in.transpose()
        data_in.index=pd.Series(timestamps[cnt])
        
        # add dataframe to list with all dataframes
        dfs.append(data_in)        
        cnt = cnt+1
    
    if TimeSeries:
        # read time series
        file=folder+'uitvoerTijdsreeksen'+pol_ifdm+'.txt'
        data_in=pd.read_fwf(file,skiprows=1,header=None)
        data_in.columns=['value','longitude','latitude','dummy3','location']    
        
        
        # select individual time series
        times=conc_timeseries.index[range(i*ntime_perrun,(i+1)*ntime_perrun)]
        for j in range(nseries):        
            where=data_in['location']==(j+1)
            values=data_in.loc[where,'value']           
            values.index=times     
            value_frame=values.rename(j).to_frame()
            conc_timeseries.loc[times,j]=values*1e9
        
# merge all concentrations
concentrations=pd.concat(dfs)

## CALCULATE PERCENTILES ##
###########################
print '\nProcessing: '+pol+'\n' 
# convert to numpy    
concentrations_np=concentrations.values 

#mean concentrations
mean_conc=concentrations_np.mean(axis=0)

# sort
sorted=np.sort(concentrations_np,axis=0)

# concentrations
p50=sorted[int(ntime*0.5),:]
p66=sorted[int(ntime*0.66),:]
p80=sorted[int(ntime*0.80),:]
p90=sorted[int(ntime*0.90),:]
p9040=sorted[int(ntime*0.9040),:]
p95=sorted[int(ntime*0.95),:]
p98=sorted[int(ntime*0.98),:]
p99=sorted[int(ntime*0.99),:]
p9979=sorted[int(ntime*0.9979),:]

## CALCULATE EXCEEDANCES ##
###########################

if pol == 'NO2':
    exceedance=(concentrations>hour_limit_NO2).sum()
if pol == 'PM10':
    conc_daily=concentrations.resample('D').sum()/24
    exceedance=(conc_daily>day_limit_pm10).sum()
if pol == 'PM25':
    conc_daily=concentrations.resample('D').sum()/24
    exceedance=(conc_daily>day_limit_pm25).sum()    

## CALCULATE DAILY SERIES ##
############################
if TimeSeries:
    conc_daily_timeseries=conc_timeseries.resample('D').sum()/24


## SAVE OUTPUT ##
#################
# for all grid cells
print '\nSaving: '+pol+'\n' 
column_names=['Lat','Lon','MeanConcentration','P50','P66','P80','P90','P9040','P95','P98','P99','P9979','exceedances']
output=pd.DataFrame(columns=column_names,index=coordinates.index)
output[['Lat','Lon']]=coordinates[['Lat','Lon']]
output.MeanConcentration=mean_conc
output.P50=p50
output.P66=p66
output.P80=p80
output.P90=p90
output.P9040=p9040
output.P95=p95
output.P98=p98
output.P99=p99
output.P9979=p9979
output.exceedances=exceedance
output=output[output.Lat!=-10]

# save concentrations
file=folder_output+'/outputIndicators'+pol+'.csv'
with open(file, 'w') as f:
    f.write('Unit: ug/m3 \nYear: '+str(timestamps[0].year)+'\n\n')
    output.to_csv(f,index=False)
    
if TimeSeries:
    # for time series
    name_series=['Location'+str(i+1) for i in range(nseries)]
    output_ts=conc_timeseries
    output_ts['Timestamp']=output_ts.index.strftime('%Y%m%d%H')    
    columns=name_series+['Timestamp']
    output_ts.columns=columns

    # for daily time series
    output_ts_daily=conc_daily_timeseries
    output_ts_daily['Timestamp']=output_ts_daily.index.strftime('%Y%m%d')    
    columns=name_series+['Timestamp']
    output_ts_daily.columns=columns

    # save time series
    coordinates_timeseries.index=name_series
    file=folder_output+'/outputTimeSeries'+pol+'.csv'
    with open(file, 'w') as f:
        coordinates_timeseries[['Lat','Lon']].to_csv(f)
        f.write('Unit: ug/m3\n\n')
        output_ts[['Timestamp']+name_series].to_csv(f)        

    # save daily time series
    file=folder_output+'/outputTimeSeriesDaily'+pol+'.csv'
    with open(file, 'w') as f:
        coordinates_timeseries[['Lat','Lon']].to_csv(f)
        f.write('Unit: ug/m3\n\n')
        output_ts_daily[['Timestamp']+name_series].to_csv(f)  
else:
    file=folder_output+'/outputTimeSeries'+pol+'.csv'
    with open(file, 'w') as f:
        f.write('No POI provided.')
    
    file=folder_output+'/outputTimeSeriesDaily'+pol+'.csv'    
    with open(file, 'w') as f:
        f.write('No POI provided.')


## VISUALISATION ##
###################
print '\nConversion to GeoTIFF: '+pol+'\n' 
# for the moment: only visualise yearly mean concentration
# temporary folder
folder_tmp=folder_output+'/tmp'+pol+'/'
if not os.path.exists(folder_tmp):
    os.makedirs(folder_tmp)

# save data in format readible for gdal
mean_conc_df=output[['Lon','Lat']+statistics_output]
tmp=(np.random.rand(len(mean_conc_df),1)-0.5)*1e-7
mean_conc_df['rand']=tmp
mean_conc_df['Lon']=mean_conc_df['Lon']+mean_conc_df['rand']

tmpfile=folder_tmp+'results'+pol+'.csv'
mean_conc_df.to_csv(tmpfile,index=False)

# loop over statistics
for stat in statistics_output:        
    
    # find extreme coordinates
    x_min=min(mean_conc_df['Lon'])
    x_max=max(mean_conc_df['Lon'])
    y_min=min(mean_conc_df['Lat'])
    y_max=max(mean_conc_df['Lat'])
    print str(x_min)+' '+str(x_max)+' '+str(y_min)+' '+str(y_max)
    x_diff=(x_max-x_min)/resolution
    y_diff=(y_max-y_min)/resolution
    print str(x_diff)+' '+str(y_diff)

    # make vrt file (determines format for gdal_grid)
    vrtfile=folder_tmp+'outputformat'+pol+stat+'.vrt'
    f= open(vrtfile, 'w')
    content='<OGRVRTDataSource> \n <OGRVRTLayer name="results'+pol+'"> \n <SrcDataSource>'+tmpfile+ ''\
        +'</SrcDataSource> \n <GeometryType>wkbPoint</GeometryType> \n ' \
        +'<LayerSRS>EPSG:4326</LayerSRS> <GeometryField encoding="PointFromColumns" x="Lon" y="Lat" z="'+ \
        stat+'"/> \n </OGRVRTLayer> \n </OGRVRTDataSource>'
    f.write(content)
    f.close()

    # grid data
    result=folder_output+'/'+stat+pol+'.tif'   
    #commando=gdal_grid+' -a LINEAR:radius=0:nodata=-9999 -outsize '+str(int(x_diff))+' '+str(int(y_diff))+' -zfield '+stat+' #'+vrtfile+' '+result
    #print commando
    #os.system(commando)  
    gridded=gdal.Grid(result,vrtfile, algorithm = 'linear:radius=0:nodata=-9999',\
        outputBounds = [x_min, y_max, x_max , y_min],  creationOptions = ['BIGTIFF=yes', 'COMPRESS=deflate'],\
        outputType = gdal.GDT_Float32,width = x_diff, height = y_diff, \
        zfield = stat)
    del gridded
            
# remove temporary folder
#for i in os.listdir(folder_tmp):
#    os.remove(folder_tmp+i)
#os.rmdir(folder_tmp) 

print ' OK'   
    
