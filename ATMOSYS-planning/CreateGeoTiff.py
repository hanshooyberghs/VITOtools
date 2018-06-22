import gdal,argparse,os
import pandas as pd


## USAGE ##
###########
#
# Converting of csv to GeoTiff file using gdal.Grid

# Gridding GeoTiff

# optional arguments:
  # -h, --help            show this help message and exit
  # -i string, --inputFile string
                        # Inputfile. Default: indicators.csv
  # -o string, --outputFolder string
                        # Outputfolder. Default: ./
  # -r float, --resolution float
                        # Resolution of output file. Default: 0.0002



## REQUIREMENTS ##
##################
#   python installation with following packages:
#       pandas
#
#   gdal, minimal version 2.1.1


## Parse input ##
#################

parser = argparse.ArgumentParser(description='Gridding GeoTiff')
parser.add_argument('-i','--inputFile',dest='inputfile',metavar='string', type=str, 
    help='Inputfile. Default: indicators.csv',
    default='indicators.csv')
parser.add_argument('-o','--outputFolder',dest='outputfolder',metavar='string', type=str, 
    help='Outputfolder. Default: ./',
    default='./')
parser.add_argument('-r','--resolution',dest='resolution',metavar='float', type=float, 
    help='Resolution of output file. Default: 0.0002',
    default='0.0002')
args = parser.parse_args()

    

folder_output=args.outputfolder
if not os.path.exists(folder_output):
    os.makedirs(folder_output)
folder_tmp=folder_output+'/tmp/'
if not os.path.exists(folder_tmp):
    os.makedirs(folder_tmp)

resolution=args.resolution  
    
## Determine column headers and boundaries ##
#############################################
# read header of inputfile
inputfile=args.inputfile
data=pd.read_csv(inputfile)
indicators=data.columns.drop(['Lat','Lon'])

x_min=min(data['Lon'])
x_max=max(data['Lon'])
y_min=min(data['Lat'])
y_max=max(data['Lat'])
print str(x_min)+' '+str(x_max)+' '+str(y_min)+' '+str(y_max)
x_diff=(x_max-x_min)/resolution
y_diff=(y_max-y_min)/resolution
print str(x_diff)+' '+str(y_diff)

for indicator in indicators:
    
    ## Make VRT file ##
    ###################
    filename=inputfile.split('/')[-1].split('.')[0]

    # make vrt file (determines format for gdal_grid)
    vrtfile=folder_tmp+'outputformat'+indicator+'.vrt'
    f= open(vrtfile, 'w')
    content='<OGRVRTDataSource> \n <OGRVRTLayer name="'+filename+'"> \n <SrcDataSource>'+inputfile+ ''\
        +'</SrcDataSource> \n <GeometryType>wkbPoint</GeometryType> \n ' \
        +'<LayerSRS>EPSG:4326</LayerSRS> <GeometryField encoding="PointFromColumns" x="Lon" y="Lat" z="'+ \
        indicator+'"/> \n </OGRVRTLayer> \n </OGRVRTDataSource>'
    f.write(content)
    f.close()


    ## Gridding ##
    ##############
    outputfile=folder_output+'/'+indicator+'.tif'   
    gridded=gdal.Grid(outputfile,vrtfile, algorithm = 'linear:radius=0:nodata=-9999',\
        outputBounds = [x_min, y_max, x_max , y_min],  creationOptions = ['BIGTIFF=yes', 'COMPRESS=deflate'],\
        outputType = gdal.GDT_Float32,width = x_diff, height = y_diff, \
        zfield = indicator)