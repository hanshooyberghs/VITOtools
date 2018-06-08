
######################################
## Usage
######################################
# example: python ConvertMeteo.py 2016 12 19

######################################
## INPUT
######################################
import pandas
import datetime
import sys

## input
input=sys.argv
startyear=int(float(input[1]))
startmonth=int(float(sys.argv[2]))
startday=int(float(sys.argv[3]))
inputfile='meteo.csv'

## input: changei for different location and different set-up
nrdays=5					# nr days
latitude=38.051771			# latitude
longitude=114.496545		# longitude			

######################################
## Convert meteo.csv to correct format
######################################


# read file
input=pandas.read_csv(inputfile)
input.index=pandas.to_datetime(input['me_date'])

# make converted variable
startdate=datetime.datetime(startyear,startmonth,startday,0,0)	# startdate
print startdate
enddate=startdate+datetime.timedelta(days=nrdays-1,hours=21)  #end-time: 21h
data=pandas.DataFrame(columns=['temperature','WE_wind_speed','NS_wind_speed'])

list_values=['value'+str(w) for w in range(0,22,3)]

for i in range(nrdays):
    day=startdate+datetime.timedelta(days=i)
    temp=input[input.index==day]
    temp.index=temp.meteo_type
    input_day=temp[list_values].transpose()
    input_day.index=pandas.date_range(day, day+datetime.timedelta(hours=21), freq='3h')
    data=data.append(input_day) 
    
######################################
## INTERPOLATION
######################################


data_out=data.resample('h').interpolate()
data_out['Height']=10
data_out['Year']=data_out.index.year
data_out['Month']=data_out.index.month
data_out['Day']=data_out.index.day
data_out['Hour']=data_out.index.hour
data_out['Lat']=latitude
data_out['Lon']=longitude
data_out['TimeParameter']=0

data_out=data_out[['WE_wind_speed','NS_wind_speed','Height','temperature','Year','Month','Day','Hour','Lat','Lon','TimeParameter']]
data_out.columns=['WE wind speed (m/s)','NS wind speed (m/s)','Height (m)','Temperature (K)',\
            'Year','Month','Day','Hour','Lat','Lon','Time parameter']


######################################
## Save data
######################################

writer = pandas.ExcelWriter('Meteo.xlsx')		# Change name of output file, if you want it in another directory
data_out.to_excel(writer,'Meteo',index=False)
writer.save()
writer.close()
