""""
Download meteo using the ECMWF API

usage: DownloadMeteo.py [-h] str str float float

Download ERA-Interim data


optional arguments:
  -h, --help  show this help message and exit
"""

from ecmwfapi import ECMWFDataServer
import os, argparse, calendar, sys, fileinput

print 'Start meteo download'

parser = argparse.ArgumentParser(description='Download ERA-Interim data')
parser.add_argument('configfiles', metavar='str', type=str, help='List of configfiles, comma separated')
parser.add_argument('matlabscript', metavar='str', type=str, help='Location of matlabscript for further processing')
args = parser.parse_args()

configfiles=args.configfiles.split(',')

for configfile in configfiles:
    print configfile
    # read config file and get input    
    params = {}
    required = ['clat', 'clon', 'fname_ecmwf','fname_atmo2D','fname_atmo3D',
        'year','startmonth','endmonth']
    with open(configfile, 'r') as file:
        for line in file:
            if len(line.strip()) <> 0 and line.strip()[0] != '#':
                splitted = line.rsplit('=')
                var = splitted[0].strip()
                val = splitted[1].strip()
                params[var] = val
  
	for param in required:
		if param not in params:
			print param, ' not found in config file!'
			raise Exception('Not all necessary parameters could be found, check your configuration file!')

    latitude=float(params['clat'])
    longitude=float(params['clon'])
    surfacefile=params['fname_atmo2D']
    profilesfile=params['fname_atmo3D']
    meteofile=params['fname_ecmwf']
    startmonth=int(params['startmonth'])
    endmonth=int(params['endmonth'])
    year=int(params['year'])
    
    # process start- and enddate
    startdate=str(year)+'-'+str(startmonth).zfill(2)+'-01'
    enddate=str(year)+'-'+str(endmonth).zfill(2)+'-'+str(calendar.monthrange(year,endmonth)[1])

    # process input location
    area=str(latitude+0.6)+'/'+str(longitude-0.6)+'/'+\
        str(latitude-0.6)+'/'+str(longitude+0.6)        
    

    # download surface file
    server = ECMWFDataServer()
    mars_surface={
        'class':"ea",
        'stream':"oper",
        'date':startdate+'/to/'+enddate,
        'dataset':"era5",
        'expver':"1",
        'repres':"ll",
        'type':"fc",
        'levtype':"sfc",
        'time':"06/18",
        'step':"1/2/3/4/5/6/7/8/9/10/11/12",
        'param':"129.128/134.128/142.128/143.128/146.128/147.128/172.128/139.128/170.128/183.128/236.128/169.128/175.128/39.128/40.128/41.128/42.128/180.128/181.128/34.128",
        'grid':"0.3/0.3",
        'area':area,		      
        'target':surfacefile,
        'format':"netcdf"	
        }
    server.retrieve(mars_surface)


    # download vertical profiles
    mars_profiles={
        'class':"ea",
        'stream':"oper",
        'date':startdate+'/to/'+enddate,
        'dataset':"era5",
        'expver':"1",
        'repres':"ll",
        'type':"fc",
        'levtype':"ml",
        'time':"06/18",
        'step':"1/2/3/4/5/6/7/8/9/10/11/12",
        'levelist':"100/to/137",
        'param':"130.128/131.128/132.128/133.128",
        'grid':"0.3/0.3",
        'area':area,		      
        'target':profilesfile,
        'format':"netcdf"	
        }

    server.retrieve(mars_profiles)

    
    # start processing in matlab
    commando=args.matlabscript+" /tools/matlab/2016a/x86_64/ "+configfile+' > '+meteofile.replace('.nc','.out')  
    os.system(commando)

        
