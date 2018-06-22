"""
####################################################
##          GetMeteoBackground.py                 ##
####################################################

Script to acquire Meteo and CAMS Background for IFDM Atmosys tool

usage: GetMeteo-Background_ERA5.py [-h] [-c string] [-lat float] [-lon float]

GetMeteoBackground_ERA5.py

optional arguments:
  -h, --help            show this help message and exit
  -c string, --city string
                        Name of city. Default: Mol
  -lat float, --latitude float
                        Latitude. Default: 51.22
  -lon float, --longitude float
                        Longitude. Default: 5.09
"""

#Load packages and set-up
###############################
from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import pandas as pd
import argparse

####################################################
##         READ COMMAND LINE INPUT                ##
####################################################


parser = argparse.ArgumentParser(description='GetMeteoBackground_ERA5.py')
parser.add_argument('-c','--city',dest='city',metavar='string', type=str, 
    help='Name of city. Default: Mol',default='Mol')
parser.add_argument('-lat','--latitude',dest='latitude',metavar='float', \
    type=float, help='Latitude between 30 and 70 degrees. Default: 51.22',default=51.22)
parser.add_argument('-lon','--longitude',dest='longitude',metavar='float',\
    type=float,     help='Longitude between -25 and 35 degrees. Default: 5.09',default=5.09)
args = parser.parse_args()
city=args.city       # name of city
latitude=args.latitude         # longitude of city
longitude=args.longitude        # latitude of city

#check input
if (latitude>70)|(latitude<30)|(longitude>35)|(longitude<-25):
    print 'Illegal input'
    raise ValueError('Incorrect longitude or latitude.')


print 'Set-up:'
print 'City: '+city
print 'Lon:  '+str(longitude)
print 'Lat:  '+str(latitude)


####################################################
##         FIXED INPUT PARAMETERS                 ##
####################################################

inputfolder_CAMS='/projects/cams95f/CAMS_data/Input/' # folder with CAMS hourly data in NetCDF format 
inputfile_ERA5='/projects/cams95f/DownloadMeteo/MeteoEurope/2015/ERA5_Europe_2015.nc'
startdate='20150101'    # startdate (format YYYYMMDD)
enddate='20160101'      # enddate (format YYYYMMDD)
pollutants=['NO2','PM10','PM25','O3']  
year=pd.to_datetime(startdate).year

###############################
###############################
##          METEO            ##
###############################
###############################

# read coordinates and find city
###############################
# read coordinates
print '\nProcessing Meteo'
print 'Reading ERA5 data'
print inputfile_ERA5
meteo_file=Dataset(inputfile_ERA5)
lon=meteo_file.variables['longitude'][:]
lat=meteo_file.variables['latitude'][:]
xx,yy=np.meshgrid(lon,lat)

#find closest point to city
distance=(xx-longitude)**2+(yy-latitude)**2
loc_x,loc_y=np.unravel_index(distance.argmin(),distance.shape)

# select 3x3 around city
x_min=max(0,loc_x-1)
x_max=min(loc_x+2,xx.shape[0])
y_min=max(0,loc_y-1)
y_max=min(loc_y+2,xx.shape[1])
xx=xx[x_min:x_max,y_min:y_max]
yy=yy[x_min:x_max,y_min:y_max]
print 'Domain meteo:'
print xx
print yy
print '\n'
        
# read meteorological data
###############################
# load data

time=meteo_file.variables['time'][:]
T=meteo_file.variables['t2m'][:,x_min:x_max,y_min:y_max]
U=meteo_file.variables['u10'][:,x_min:x_max,y_min:y_max]
V=meteo_file.variables['v10'][:,x_min:x_max,y_min:y_max]
Rs=meteo_file.variables['ssrd'][:,x_min:x_max,y_min:y_max]
meteo_file.close()

# convert time to Python format
time_python=pd.datetime(1900,1,1)+pd.to_timedelta(time[:],'h')

# radiation per second
Rs=Rs/10800.

 # Spatial interpolation
###############################
print 'Spatial interpolation meteo'
# Interpolation to city location
T_city= np.zeros(shape=(len(time_python),))
U_city= np.zeros(shape=(len(time_python),))
V_city= np.zeros(shape=(len(time_python),))
Rs_city= np.zeros(shape=(len(time_python),))


for i in range(len(time_python)):
    f = interpolate.interp2d(xx, yy, T[i,:,:], kind='linear')
    T_city[i]=f(longitude,latitude)
    f = interpolate.interp2d(xx, yy, U[i,:,:], kind='linear')
    U_city[i]=f(longitude,latitude)
    f = interpolate.interp2d(xx, yy, V[i,:,:], kind='linear')
    V_city[i]=f(longitude,latitude)
    f = interpolate.interp2d(xx, yy, Rs[i,:,:], kind='linear')
    Rs_city[i]=f(longitude,latitude)

    
# Correction
Rs_city[Rs_city<0]=0


# Interpolation to hourly data
###############################
print 'Temporal interpolation meteo'

# store input in dataframe
input=pd.DataFrame(columns=['Temperature (K)','WE wind speed (m/s)','NS wind speed (m/s)','Solar radiation (J /m**2)'])
input['Temperature (K)']=T_city
input['WE wind speed (m/s)']=U_city
input['NS wind speed (m/s)']=V_city
input['Solar radiation (W / m**2)']=Rs_city
input.index=time_python

# make output dataframe for hourly data
start=pd.to_datetime(startdate)
end=pd.to_datetime(enddate)
output_meteo=pd.DataFrame(index=pd.date_range(start=start,end=end,freq='1H'))

# store 3h data in 1h output frame
output_meteo=pd.merge(output_meteo,input,how='left',left_index=True,right_index=True)

# interpolation to hourly results
output_meteo=output_meteo.interpolate()   

# neglect last hour (is next year)
output_meteo=output_meteo.ix[:-1]

# fill first and last hours if empty
output_meteo=output_meteo.fillna(method='backfill') 

# Finalize IFDM fromat
###############################

# add parameters for IFDM
output_meteo['Year']=output_meteo.index.year
output_meteo['Month']=output_meteo.index.month
output_meteo['Day']=output_meteo.index.day
output_meteo['Hour']=output_meteo.index.hour
output_meteo['Height (m)']=10
output_meteo['Lat']=latitude
output_meteo['Lon']=longitude
output_meteo['Time parameter']=0
output_meteo=output_meteo[['WE wind speed (m/s)','NS wind speed (m/s)','Height (m)','Temperature (K)',\
     'Year','Month','Day','Hour','Lat','Lon','Time parameter','Solar radiation (W / m**2)']]

print 'Meteo done\n\n'



###############################
###############################
##        BACKGROUND         ##
###############################
###############################

# read coordinates and find city
###############################
# read coordinates
print 'Processing background'
print 'Reading CAMS data'
CAMS_files=Dataset(inputfolder_CAMS+'ENSa.'+str(year)+'.NO2.yearlyrea.nc')
lon=CAMS_files.variables['lon'][:]
lat=CAMS_files.variables['lat'][:]
xx,yy=np.meshgrid(lon,lat)

#find closest point to city
distance=(xx-longitude)**2+(yy-latitude)**2
loc_x,loc_y=np.unravel_index(distance.argmin(),distance.shape)

# select 3x3 around city
x_min=max(0,loc_x-1)
x_max=min(loc_x+2,xx.shape[0])
y_min=max(0,loc_y-1)
y_max=min(loc_y+2,xx.shape[1])
xx=xx[x_min:x_max,y_min:y_max]
yy=yy[x_min:x_max,y_min:y_max]

print 'Domain background:'
print xx
print yy
print '\n'

#read and process pollutants 
###############################
output_cams=pd.DataFrame(columns=pollutants+['NO'])
for pollutant in pollutants:
    print 'Processing '+pollutant
    print 'Reading data'
    CAMS_files=Dataset(inputfolder_CAMS+'ENSa.'+str(year)+'.'+pollutant+'.yearlyrea.nc')
    data=CAMS_files.variables[pollutant][:,x_min:x_max,y_min:y_max]
            
    # interpolation
    print 'Interpolation'
    data_interpolated= np.zeros(shape=(len(data),))
    for i in range(len(data)):
        f = interpolate.interp2d(xx, yy, data[i,:,:], kind='linear')
        data_interpolated[i]=f(longitude,latitude)
    
    output_cams.loc[:,pollutant]=data_interpolated
    
# Add time format
output_cams.index=output_meteo.index
    
## NO calculation
###############################
print 'Calculating NO'
# corrections for small ozone values
output_cams.O3[output_cams.O3<1]=1

# molar masses
mNO2 = 46.0056 
mNO  = 30.0061 
mO3  = 47.9983 
  
# Avogadro
Na   = 6.02214e23 # avogadro's constant

# Conversion to molecules/cm3
nNO2=output_cams['NO2']*1e-12*Na/mNO2
nO3=output_cams['O3']*1e-12*Na/mO3

# Reaction constants
J=0.8e-3*np.exp(-10.0/output_meteo['Solar radiation (W / m**2)'])+7.4e-6 * output_meteo['Solar radiation (W / m**2)']
K= 2.2e-12 * np.exp( - 1430.0 /output_meteo['Temperature (K)']  )
nNO = pd.DataFrame(0, index=output_meteo.index, columns=['NO'])
nNO[output_meteo['Solar radiation (W / m**2)']>0]=J/K*nNO2/nO3
output_cams.NO=nNO  / 1.0e-12 * mNO  / Na
output_cams['NOx']=output_cams.NO+output_cams.NO2



###############################
###############################
##          SAVING           ##
###############################
###############################
# save in ATMOSYS tool background - meteo format
print 'Saving file for ATMOSYS tool'
# merge dataframes
output=pd.merge(output_meteo,output_cams,how='left',left_index=True,right_index=True)
output['Winddirection']=180+np.arctan2(output['WE wind speed (m/s)'],output['NS wind speed (m/s)'])*180/np.pi
output['Windspeed']=np.sqrt(output['NS wind speed (m/s)']**2+output['WE wind speed (m/s)']**2)
output['Temperature (C)']=output['Temperature (K)']-273.15
output[[u'Year', u'Month', u'Day', u'Hour']]=output[[u'Year', u'Month', u'Day', u'Hour']].astype(int)
output['Winddirection']=output['Winddirection'].astype(int)
output['BC']=0

columns=[u'Year', u'Month', u'Day', u'Hour',
       u'Windspeed', u'Winddirection', 
       u'Temperature (C)', u'Height (m)', u'Lat', u'Lon',
       u'Time parameter']
output_save=output[columns]
columns_save=[u'year', u'month', u'day', u'hour',
       u'wind_speed (m/s)', u'wind_direction (degree)', 
       u'temperature (C)', u'height (m)', u'lat', u'lon',
       u'time_zone (h)']
output_save.columns=columns_save
output_save.to_csv('MeteoFile_ATMOSYS_'+city+'.txt',sep='\t',index=False)


columns=[u'Year', u'Month', u'Day', u'Hour',
       u'O3',u'NOx',u'NO2',u'PM25',u'PM10','BC']
output_save=output[columns]
columns_save=[u'year', u'month', u'day', u'hour',
       u'O3 (ug/m3)',u'NOX (ug/m3)',u'NO2 (ug/m3)',u'PM25 (ug/m3)',u'PM10 (ug/m3)','BC (ug/m3)']
output_save.columns=columns_save
output_save.to_csv('BackgroundFile_ATMOSYS_'+city+'.txt',sep='\t',index=False)
       
print 'Done'
