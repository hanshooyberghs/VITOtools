""""
usage: PostprocessingIFDMOSPM_IRCEL.py [-h] [-r INT] [-f folder]
                                       folder_run folder_output

Postprocessing IFDM-OSPM results

positional arguments:
  folder_run            Folder with results of IFDM and OSPM runs
  folder_output         Folder where processed files will be stored

optional arguments:
  -h, --help            show this help message and exit
  -r INT, --resolution INT
                        Resolution. Default:10
  -f file, --file_streetcanyon_data file
                        Resolution. Default: ./Inputfiles/

"""

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import os,gdal,sys,datetime,argparse,Indicators
import numpy as np
gdal.UseExceptions()


######################
## input parameters ##
######################
# list of pollutants
pol_list=['NO2','PM10','PM25','BC','O3']

# list of indicators
indicators_NO2 = ["YearlyMean","ExcHour_200","P9978","MeanDayMax1Hour"]
indicators_PM10 = ["YearlyMean","ExcDay_505","P9014_daily"]
indicators_BC = ["YearlyMean"]
indicators_PM25 = ["YearlyMean","P9014_daily"]
indicators_O3 =  ["YearlyMean","NET60","WHO","AOT40veg", "AOT40for","AOT60","meanmax8h"]
indicators_lists = {'NO2': indicators_NO2, 'PM10': indicators_PM10, 'PM25': indicators_PM25,\
    'BC': indicators_BC,'O3': indicators_O3}

# timeframe
startDate=datetime.datetime(2016,01,01,00)  
endDate=datetime.datetime(2016,12,31,23)  

######################
##   fixed input    ##
######################
# do not change the following input unless you know what you are doing!
# read command line input
parser = argparse.ArgumentParser(description='Postprocessing IFDM-OSPM results')
parser.add_argument('runfolder', metavar='folder_run', type=str, help='Folder with results of IFDM and OSPM runs')
parser.add_argument('outputfolder', metavar='folder_output', type=str, help='Folder where processed files will be stored')
parser.add_argument('-r','--resolution',dest='resolution',metavar='INT', type=int, 
    help='Resolution. Default:10',default=10)
parser.add_argument('-f','--file_streetcanyon_locations',dest='file_streetcanyon_locations',metavar='folder', type=str, 
    help='Folder with streetcanyon data. Default: ./Inputfiles/LocationsSC.tif',default='Inputfiles/')
args = parser.parse_args()
    

Runfolder=args.runfolder
Outputfolder=args.outputfolder+'/'
resolution=args.resolution
posities_streetcanyons=args.file_streetcanyon_locations


# default extent of grid (should not be changed unless buffer file is also changed)
x_min=20000
x_max=300000
y_min=21000
y_max=245000

# coordinate
crs={'init':'epsg:31370'}


# max coords => if the numer of coordinates is larger, outputfiles will be read in different parts
maxcoords=100000


#############################
##  read large input files ##
#############################
print 'Load street canyon locations'
# read canyon grid file
gdata_canyon = gdal.Open(posities_streetcanyons)
gt_canyon = list(gdata_canyon.GetGeoTransform())
data_canyon = gdata_canyon.ReadAsArray().astype(float)
gdata_canyon = None         

#make output folder
if not os.path.exists(Outputfolder):
    os.makedirs(Outputfolder)
else:
    os.system('rm '+Outputfolder+'*')
    
# temporary folder
tmpfolder='tmp/'
if not os.path.exists(tmpfolder):
    os.makedirs(tmpfolder)
else:
    os.system('rm '+tmpfolder+'*')
    
# daterange
datetimes=pd.date_range(start=startDate,end=endDate,freq='H')

    
for pollutant in pol_list:
    print '----------------------'
    print 'Start processing '+pollutant
    
    print '----------------------'
    
    # selection of indicators
    indicators= indicators_lists[pollutant]
    
    ##############################
    ## convert data IFDM to csv ##
    ##############################

    output_array=[]  # array to combine different parts
    read_start=1
    kolom_list=['X','Y']+[pollutant+'_'+w for w in indicators]
    
    # determine maximal size by reading one file
    file=Runfolder+'/output'+pollutant.ljust(4,'a')+'_saga_'+datetimes[0].strftime('%Y%m%d%H')+'.txt'
    data_in=pd.read_csv(file,sep=';',header=None)
    ncoords=len(data_in)
    fileindex=data_in.index
    
    # read part by part 
    print 'Load IFDM output. This might take some time...'
    while read_start<ncoords:
        # define pratical end of reading
        read_end=min(read_start+maxcoords,ncoords)
        print 'Read from '+str(read_start)+' to '+str(read_end)
        toskip=np.argwhere(~fileindex.isin(range(read_start,read_end)))
        
        # reading
        dfs=[]
        for dt in datetimes: 
            if (dt.day==1)&(dt.hour==0):
                print 'Reading: '+dt.isoformat()
            # read file
            
            file=Runfolder+'/output'+pollutant.ljust(4,'a')+'_saga_'+dt.strftime('%Y%m%d%H')+'.txt'
            data_in=pd.read_csv(file,sep=';',header=None,skiprows=toskip)
            data_in.columns=['x','y',dt.strftime('%Y%m%d%H')]
            #data_in=data_in.iloc[range(read_start,read_end)]
            dfs.append(data_in[dt.strftime('%Y%m%d%H')])
            
            # get coordinates
            if not ('coords' in locals()):
                coords=data_in[['x','y']]
            del data_in
            
        # concatenate values
        result=pd.concat(dfs,axis=1)       
        result.columns=pd.to_datetime(result.columns+'00')

        # free memory
        while len(dfs)>0:
            del dfs[0]
        del dfs
        
            
        
        ## indicators IFDM ##
        #####################
        print 'Calculate Indicators'
        ind_IFDM=[]
        for indicator in indicators:        
            function=getattr(Indicators,indicator) 
            ind_IFDM.append(function(result))
            
        # add coordinates
        final_IFDM=pd.concat(ind_IFDM,axis=1)
        final_IFDM=pd.concat([coords,final_IFDM],axis=1)
        final_IFDM.columns=kolom_list
        
        # add result to array 
        output_array.append(final_IFDM)
        del final_IFDM
        del result
        del coords
    
        while len(ind_IFDM)>0:
            del ind_IFDM[0]
        del ind_IFDM

        
        # update loop 
        read_start=read_start+maxcoords
            
    # combine subloops
    final_IFDM=pd.concat(output_array)
    
    # wegschrijven IFDM resultaten
    result_ifdm=Outputfolder+pollutant+'_IFDM.csv'

    print 'Save results as CSV-file'
    
    final_IFDM.columns=kolom_list
    final_IFDM.to_csv(result_ifdm,index=False)
    
    
    # free memory
    while len(output_array)>0:
        del output_array[0]
    del output_array
    del final_IFDM
    

    ##############################
    ## convert data OSPM to csv ##
    ##############################
    if pollutant != 'O3':
        output_array=[]  # array to combine different parts
        read_start=1
        kolom_list=['X','Y']+[pollutant+'_'+w for w in indicators]
        
        # determine maximal size by reading one file
        file=Runfolder+'/OSPMfinal'+pollutant.ljust(4,'a')+'_saga_'+datetimes[0].strftime('%Y%m%d%H')+'.txt'
        data_in=pd.read_csv(file,sep=';',header=None)
        ncoords=len(data_in)
        fileindex=data_in.index
        
        # read part by part 
        print 'Load OSPM output. This might take some time...'
        while read_start<ncoords:        
            # define pratical end of reading
            read_end=min(read_start+maxcoords,ncoords)
            print 'Read from '+str(read_start)+' to '+str(read_end)
            toskip=np.argwhere(~fileindex.isin(range(read_start,read_end)))
            
            dfs=[]
            for dt in datetimes:  
                if (dt.day==1)&(dt.hour==0):
                    print 'Reading: '+dt.isoformat()
                # read file
                file=Runfolder+'/OSPMfinal'+pollutant.ljust(4,'a')+'_saga_'+dt.strftime('%Y%m%d%H')+'.txt'
                
                data_in=pd.read_csv(file,sep=';',header=None,skiprows=toskip)
                #data_in=data_in.iloc[range(read_start,read_end)]
                data_in.columns=['x','y',dt.strftime('%Y%m%d%H'),'tmp1']
                dfs.append(data_in[dt.strftime('%Y%m%d%H')])
                # get coordinates
                if not ('coords' in locals()):
                    coords_OSPM=data_in[['x','y']]
                del data_in


            # concatenate values
            result_OSPM=pd.concat(dfs,axis=1)
            result_OSPM.columns=pd.to_datetime(result_OSPM.columns+'00')

            # free memory
            while len(dfs)>0:
                del dfs[0]
            del dfs
            
            
            ## indicators OSPM ##
            #####################
            
            ind_OSPM=[]
            for indicator in indicators:        
                function=getattr(Indicators,indicator) 
                ind_OSPM.append(function(result_OSPM))    
            
            
            # add coordinates
            final_OSPM=pd.concat(ind_OSPM,axis=1)
            final_OSPM=pd.concat([coords_OSPM,final_OSPM],axis=1)
            final_OSPM.columns=kolom_list
            
            # add result to array 
            output_array.append(final_OSPM)
            del final_OSPM
            del coords_OSPM
            del result_OSPM
            while len(ind_OSPM)>0:
                del ind_OSPM[0]
            del ind_OSPM

            
            # update loop 
            read_start=read_start+maxcoords
                
        # combine subloops
        
        final_OSPM=pd.concat(output_array)
        
        # output van het script
        result_ospm=Outputfolder+pollutant+'_OSPM.csv'        
        final_OSPM.to_csv(result_ospm,index=False)
        
        # remove points that are close to each other
        print 'Remove points that are close to each other'
        final_OSPM['geometry']=final_OSPM.apply(lambda p:Point(p.X,p.Y),axis=1)
        final_OSPM=gpd.GeoDataFrame(final_OSPM,crs=crs)

        buffer_points=final_OSPM.geometry.buffer(10)
        buffer_points=gpd.GeoDataFrame(buffer_points,columns=['geometry'],crs=crs)
        buffer_points.index=buffer_points.index.rename('left_index')
        final_OSPM.index=final_OSPM.index.rename('right_index')
        columns=[pollutant +'_'+w for w in indicators]
        joined=gpd.sjoin(buffer_points,final_OSPM,how='left',op='intersects')
        mean=pd.DataFrame(joined.groupby(by=joined.index).mean()[columns],columns=columns)
        final_OSPM_save=final_OSPM.drop(columns,axis=1).merge(mean,left_index=True,right_index=True)
        final_OSPM_save.to_csv(result_ospm.replace('.csv','_smoothed.csv'))
        print 'Save results as CSV-file and move on to gridding'
            
        # free memory
        while len(output_array)>0:
            del output_array[0]
        del output_array
        del mean
        del joined
        del buffer_points
        del final_OSPM
        del final_OSPM_save
        

    
    #############################
    #############################
    ##       Geoprocessing     ##
    #############################
    #############################

    for indicator in indicators:    
        
        output=Outputfolder+'/'+pollutant+'_'+indicator+'_IFDM-OSPM.tif'
        gridded_ifdm=Outputfolder+'/'+pollutant+'_'+indicator+'_IFDM.tif'
        
        ############################
        ##  gridding ifdm results ##
        ############################
        print 'Grid IFDM resultaten'           
            
        print str(x_min)+' '+str(y_min)+' '+str(x_max)+' '+str(y_max)
        x_diff=(x_max-x_min)/resolution
        y_diff=(y_max-y_min)/resolution
        print str(x_diff)+' '+str(y_diff)

        # make vrt file (determines format for gdal_grid)
        vrtfile=tmpfolder+'outputformatIFDM.vrt'
        f= open(vrtfile, 'w')
        content='<OGRVRTDataSource> \n <OGRVRTLayer name="'+pollutant+'_IFDM"> \n <SrcDataSource>'+result_ifdm+ ''\
            +'</SrcDataSource> \n <GeometryType>wkbPoint</GeometryType> \n ' \
            +'<LayerSRS>EPSG:31370</LayerSRS> <GeometryField encoding="PointFromColumns" x="X" y="Y" z="'+pollutant+'_'+indicator+'"' \
            +'/> \n </OGRVRTLayer> \n </OGRVRTDataSource>'
        f.write(content)
        f.close()

        # grid data

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
        if pollutant != 'O3':
            
            # make vrt file (determines format for gdal_grid)
            vrtfile=tmpfolder+'/outputformatOSPM.vrt'
            f= open(vrtfile, 'w')
            content='<OGRVRTDataSource> \n <OGRVRTLayer name="'+pollutant+'_OSPM_smoothed"> \n <SrcDataSource>'+result_ospm.replace('.csv','_smoothed.csv')+ ''\
                +'</SrcDataSource> \n <GeometryType>wkbPoint</GeometryType> \n ' \
                +'<LayerSRS>EPSG:31370</LayerSRS> <GeometryField encoding="PointFromColumns" x="X" y="Y" z="'+pollutant+'_'+indicator+'"' \
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

            # checks on grid
            if np.shape(data_ifdm)!=np.shape(data_canyon):
                print 'Error: grid with canyons has wrong resolution'
                sys.exit() 

            # error in the orientation of the canyon file
            if gt_ospm[5]*gt_canyon[5]<0:
                data_canyon=np.flipud(data_canyon)
                gt_canyon[5]*=(-1)
                
            # combination
            print 'actual combination' 
            final_result=data_ifdm
            loc=((data_canyon==1)&(~np.isnan(data_ospm)))&(data_ospm>data_ifdm)
            final_result[loc]=data_ospm[loc]


            # saving results
            print 'Save results in GeoTiff'
            [cols, rows] = final_result.shape
            driver = gdal.GetDriverByName("GTiff")
            outdata = driver.Create(output, rows, cols, 1, gdal.GDT_Float32,options = [ 'COMPRESS=LZW' ])
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
    print 'Processing '+pollutant+' done'
    print '----------------------'

