import sys
import shutil
import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
from PIL import Image
import matplotlib.pyplot as plt
plt.ioff()
from netCDF4 import Dataset
from matplotlib import colors
import gdal

from collections import Counter
import tifffile

def terrain2urbclim(configfile_input):
    # ******************************************************************
    # ******************************************************************
    #   READ CONFIG FILE (KEY VALUE PAIRS)
    # ******************************************************************
    # ******************************************************************

    params = {}
    required = ['nx', 'ny', 'clat', 'clon', 'dxy', 'EPSG', 'pname_modis', 'veg_option',             'fname_dem', 'fname_soil',
                'fname_land', 'fname_land_parms', 'fname_ahf', 'fname_surface',
                'use_buildings_z0m','fname_buildings','building_height','urban_tree_fraction']
    with open(configfile_input, 'r') as config_file:
        for line in config_file:
            if len(line.strip()) <> 0 and line.strip()[0] != '#':
                splitted = line.rsplit('=')
                var = splitted[0].strip()
                val = splitted[1].strip()
                params[var] = val

    for param in required:
        if param not in params:
            print param, ' not found in config file!'
            raise Exception('Not all necessary parameters could be found, check your configuration file!')

    nx = int(params['nx'])
    ny = int(params['ny'])
    dxy = long(params['dxy'])
    logfolder=os.path.dirname(params['fname_surface'])+'/'
    if os.path.exists(logfolder):
        shutil.rmtree(logfolder)
    os.makedirs(logfolder)

    # *****************************************************************
    # *****************************************************************
    #    GENERATE GRID
    # *****************************************************************
    # *****************************************************************

    # The grid is constructed with clon and clat used to define the nearest round number in the choosen coordinate system,
    # and then taking (nx-1)/2 and (ny-1)/2 grid cells in every direction. nx and ny should be uneven to ensure that the coordinates denote the centers of the grid cells.
    # The resolution is set to the choosen value.

    if nx % 2 == 0 or ny % 2 == 0:
        raise Exception('nx or ny is an even number, only uneven numbers are allowed')

    # DEFINE CENTER COORDINATES
    str1 = 'echo ' + params['clon'] + ' ' + params['clat'] + ' | gdaltransform -s_srs EPSG:4326 -t_srs EPSG:' + params[
        'EPSG']
    result = execute(str1)
    result = result.split(' ')
    ccorx = round(float(result[0]))
    ccory = round(float(result[1]))

    print 'Center coordinates: ', ccorx, ccory

    # CALCULATE GRID

    xmin = ccorx - ((nx - 1) / 2) * dxy
    ymin = ccory - ((ny - 1) / 2) * dxy
    xmax = ccorx + ((nx - 1) / 2) * dxy
    ymax = ccory + ((ny - 1) / 2) * dxy

    xcor = xmin + np.arange(nx) * dxy
    ycor = ymin + np.arange(ny) * dxy

    xcortemp = []
    for i in range(int(ny)):
        xcortemp.append(xcor)
    ycortemp = []
    for i in range(int(nx)):
        ycortemp.append(ycor)
    xcor = xcortemp
    ycor = ycortemp

    ycor = np.transpose(ycor)

    #xcor=np.transpose(np.tile(xcortmp,(ny,1)))
    #ycor = np.fliplr(np.tile(ycortmp,(nx,1)))

    xmins = xmin - 0.5 * dxy
    xmaxs = xmax + 0.5 * dxy
    ymins = ymin - 0.5 * dxy
    ymaxs = ymax + 0.5 * dxy

    print xmins, xmaxs, ymins, ymaxs
    # #CREATE LAT/LON GRID

    with open(logfolder+'/temp_coord.txt', 'w') as coordfile:
        for j in range(len(xcor[0])):
            for i in range(len(xcor)):
                coordfile.write("" + "{:.5f}".format(xcor[i][j]) + " " + "{:.5f}".format(ycor[i][j]) + "\n")

    str2 = 'gdaltransform -s_srs EPSG:' + params['EPSG'] + ' -t_srs EPSG:4326 < '+logfolder+'/temp_coord.txt > '+logfolder+'/temp_lonlat.txt'
    execute(str2)

    lat = np.zeros((ny, nx), float)
    lon = np.zeros((ny, nx), float)

    with open(logfolder+'/temp_lonlat.txt', 'r') as lonlatfile:
        for j in range(nx):
            for i in range(ny):
                line = lonlatfile.readline()
                result = line.rsplit(" ")

                lon[i][j] = float(result[0])
                lat[i][j] = float(result[1])

    lnmn = np.min(lon)
    lnmx = np.max(lon)
    ltmn = np.min(lat)
    ltmx = np.max(lat)

    # ******************************************************************
    # ******************************************************************
    # ******************************************************************
    #                          LAND USE
    # ******************************************************************
    # ******************************************************************
    # ******************************************************************

    str3 = 'gdalwarp -s_srs EPSG:3035 -t_srs EPSG:' + params['EPSG'] + ' -te ' + str(xmins) + ' ' + str(ymins) + ' ' + str(
        xmaxs) + ' ' + str(ymaxs) + ' -tr ' + str(dxy) + ' ' + str(dxy) + ' -r mode ' + params[
               'fname_land'] +' '+ logfolder+'/temp_cor.tif'
    execute(str3)

    vtyp = read_tiff(logfolder+'/temp_cor.tif')

    # Read land parameters from external file
    # ---------------------------------------
    albv1 = []
    albs1 = []
    epsv1 = []
    epss1 = []
    z0m1 = []
    rst1 = []
    rp1 = []
    rootb1 = []
    corine = []
    nroftypes = 15

    count = 0

    with open(params['fname_land_parms'], 'r') as params_file:
        for line in params_file:
            if len(line.strip()) <> 0 and line.strip()[0] != '#':
                splitted = filter(None, line.split(" "))
                count += 1
                albv1.append(float(splitted[1]))
                albs1.append(float(splitted[2]))
                epsv1.append(float(splitted[3]))
                epss1.append(float(splitted[4]))
                z0m1.append(float(splitted[5]))
                rst1.append(float(splitted[6]))
                rp1.append(float(splitted[7]))
                rootb1.append(float(splitted[8]))
                corine.append(splitted[9])

    if count <> nroftypes:
        raise Exception('Could not find input data for all URBCLIM land use classes, check your land_parameters.dat file!')

    # associate CORINE classes with URBCLIM land use types
    # -----------------------------------------------------

    typval = []
    for i in range(50):
        typval.append(15)

    for i in range(15):
        types = corine[i].split(',')
        for x in types:
            typval[int(x)] = i + 1

    #
    # for x in range(0,len(vtyp)):
    #     for y in range(0, len(vtyp[x])):
    #         vtyp[x][y] = typval[vtyp[x][y]]

    sfctyp = replace_values(vtyp, typval, int)
    sfctyp[0][0] = 14
    # set up colour scheme for visualization
    dbord = 245
    lbord = 235
    purple = 20
    orange = 220
    lorange = 190
    lblue = 70
    dblue = 40
    medgrn = 135
    lgreen = 150
    dgreen = 85
    yellow = 165

    # color for each land use type (1-15)
    color = [0, dbord, lbord, purple, dgreen, orange, medgrn, medgrn, lgreen, yellow, dgreen, dgreen, dgreen, lblue, lblue,
             dblue]

    # Check if every grid point belongs to a class + make land and sea mask
    # ----------------------------------------------------------------------
    print np.unique(sfctyp)
    landmask = np.full_like(sfctyp, 0)
    landmask[np.where(sfctyp < 13)] = 1

    seemask = np.full_like(sfctyp, 0)
    seemask[np.where(sfctyp >= 13)] = 1

    lan = sfctyp
    if np.min(lan) == 0:
        raise Exception('Problem in the land use image')

    albv = replace_values(lan-1, albv1, float)
    albs = replace_values(lan-1, albs1, float)
    epsv = replace_values(lan-1, epsv1, float)
    epss = replace_values(lan-1, epss1, float)
    z0m = replace_values(lan-1, z0m1, float)
    rst = replace_values(lan-1, rst1, float)
    rp = replace_values(lan-1, rp1, float)
    rootb = replace_values(lan-1, rootb1, float)

    # ******************************************************************
    # ******************************************************************
    # ******************************************************************
    #                      VEGETATION NDVI
    # ******************************************************************
    # ******************************************************************
    # ******************************************************************
    # MODIS 250m NDVI (MOD13Q1)

    cover = []
    # nn = 0

    for ii in range(1, 13):
        month = str(ii).zfill(2)
        vgtfile = params['pname_modis'] + '/mm62' + month + 'k.tif'
        str4 = 'gdalwarp -s_srs EPSG:3035 -t_srs EPSG:' + params['EPSG'] + ' -te ' + str(xmins) + ' ' + str(
            ymins) + ' ' + str(xmaxs) + ' ' + str(ymaxs) + ' -tr ' + str(dxy) + ' ' + str(
            dxy) + ' ' + vgtfile + ' '+logfolder+'/temp_veg_' + month + '.tif'
        execute(str4)
        dvi = read_tiff(logfolder+'/temp_veg_' + month + '.tif')
        if np.min(dvi) < 10 or np.max(dvi) > 250:  # missing pixels
            dvi[dvi < 10] = 1
            dvi[dvi > 250] = 1
        dvi = dvi * landmask
        cnt = sum(row.count(1) for row in dvi.tolist())
        if cnt > 0:
            problem = np.where(dvi == 1)
            for i in range(len(problem[0])):
                val = get_smoothed_value(dvi, 5, problem[0][i], problem[1][i], 1)
                dvi[problem[0][i], problem[1][i]] = val
        dvi = 0.004 * dvi - 0.08
        dvi[dvi < 0] = 0
        dvimax = 0.85  # np.max(dvi)
        dvimin = np.ma.array(dvi, mask=np.logical_not(landmask)).min()
        veg = (dvi - dvimin) / (dvimax - dvimin)
        veg[veg < 0] = 0
        veg[veg > 1] = 1
        if params['veg_option'] == 'MEAN':
            for kk in range(1, 16):
                sites = np.where(lan == kk)
                veg[sites] = np.median(veg[sites])
        
        veg = veg * landmask
        cover.append(veg)
        
        
    # ******************************************************************
    # ******************************************************************
    #                        BUILDINGS ROUGHNESS LENGTH
    # ******************************************************************
    # ******************************************************************
    if int(params['use_buildings_z0m']) == 1:
        print('Use building height for calculation of the roughness length')
        str3 = 'gdalwarp -s_srs EPSG:3035 -t_srs EPSG:' + params['EPSG'] + ' -te ' + str(xmins) + ' ' + str(ymins) + ' ' + str(
            xmaxs) + ' ' + str(ymaxs) + ' -tr ' + str(dxy) + ' ' + str(dxy) + ' -r average -ot Float32 ' +\
            params['fname_buildings'] +' '+ logfolder+'/temp_buildings.tif'
       
        execute(str3)
        buildings = read_tiff(logfolder+'/temp_buildings.tif')
        green=cover[6]
        # combination of buildings and tree 
        a = float(params['building_height'])
        b = float(params['urban_tree_fraction'])
#        z0m_buildings=(buildings*params['building_height']/10.)+\
#            (params['urban_tree_fraction']*green*2.) 
        z0m_buildings=(buildings*a/10.)+\
            (b*green*2.) 
        z0m_buildings[z0m_buildings >= 2.0]=2.0
        z0m_buildings[z0m_buildings <= 0.1]=0.1
 #       z0m[(sfctyp <= 3.) & (z0m_buildings >= 0.1)]=\
 #           z0m_buildings[(sfctyp <= 3.) & (z0m_buildings >= 0.1)]
        z0m[(sfctyp <= 3.)]=\
            z0m_buildings[(sfctyp <= 3.)]

    # ******************************************************************
    # ******************************************************************
    #                        DEM TERRAIN HEIGHT
    # ******************************************************************
    # ******************************************************************


    str5 = 'gdalwarp -s_srs EPSG:4326 -t_srs EPSG:' + params['EPSG'] + ' -te ' + str(xmins) + ' ' + str(ymins) + ' ' + str(
        xmaxs) + ' ' + str(ymaxs) + ' -tr ' + str(dxy) + ' ' + str(dxy) + ' -r bilinear ' + params[
               'fname_dem'] +' '+ logfolder+'/temp_dem.tif'
    execute(str5)

    with tifffile.TiffFile(logfolder+'/temp_dem.tif') as tif:
        dem = tif.asarray()

    outofdem = np.where(dem < 0)
    dem[outofdem] = 0
    dem = dem.astype(float)

    # ******************************************************************
    # ******************************************************************
    #                       SOIL SEALING
    # ******************************************************************
    # ******************************************************************

    str7 = 'gdalwarp -s_srs EPSG:3035 -t_srs EPSG:' + params['EPSG'] + ' -te ' + str(xmins) + ' ' + str(ymins) + ' ' + str(
        xmaxs) + ' ' + str(ymaxs) + ' -tr ' + str(dxy) + ' ' + str(dxy) + ' ' + params['fname_soil'] + ' '+logfolder+'/temp_soil.tif'

    execute(str7)

    soil = read_tiff(logfolder+'/temp_soil.tif')

    # rescaling
    soil=soil/100.0
    soil = soil * landmask

    problem = np.where(soil > 1.0)
    soil[problem] = 10
    if len(problem[0]) > 0:
        for i in range(len(problem[0])):
            val = get_smoothed_value(soil, 10,problem[0][i], problem[1][i], 10)
            soil[problem[0][i], problem[1][i]] = val

    soil = soil * landmask

    # ******************************************************************
    # ******************************************************************
    #                       ANTROPOGENIC HEAT FLUX
    # ******************************************************************
    # ******************************************************************

    str8 = 'gdalwarp -s_srs EPSG:4326 -t_srs EPSG:' + params['EPSG'] + ' -te ' + str(xmins) + ' ' + str(ymins) + ' ' + str(
        xmaxs) + ' ' + str(ymaxs) + ' -tr ' + str(dxy) + ' ' + str(dxy) + ' -r cubicspline ' + params[
               'fname_ahf'] + ' '+logfolder+'/temp_ahf.tif'

    execute(str8)

    ahf = read_tiff(logfolder+'/temp_ahf.tif')

    # ******************************************************************
    # ******************************************************************
    #                           VISUALISATION
    # ******************************************************************
    # ******************************************************************

    f, axarr = plt.subplots(2, 3, figsize=(20, 20))
    im = axarr[1, 1].imshow(soil, extent=[lnmn, lnmx, ltmn, ltmx], aspect='auto')
    axarr[1, 1].set_title('Soil sealing (%)')
    im.set_cmap('rainbow')
    f.colorbar(im, ax=axarr[1, 1], orientation='horizontal')

    im = axarr[1, 0].imshow(dem, extent=[lnmn, lnmx, ltmn, ltmx], aspect='auto')
    axarr[1, 0].set_title('Terrain height (m)')
    im.set_cmap('rainbow')
    f.colorbar(im, ax=axarr[1, 0], orientation='horizontal')
    
    im = axarr[0, 2].imshow(z0m, extent=[lnmn, lnmx, ltmn, ltmx], aspect='auto')
    axarr[0, 2].set_title('Roughness Length (m)')
    im.set_cmap('rainbow')
    f.colorbar(im, ax=axarr[0, 2], orientation='horizontal')

    im = axarr[1, 2].imshow(soil * ahf, extent=[lnmn, lnmx, ltmn, ltmx], aspect='auto')
    axarr[1, 2].set_title('AHF (W/m)')
    im.set_cmap('rainbow')
    f.colorbar(im, ax=axarr[1, 2], orientation='horizontal')

    im = axarr[0, 0].imshow(cover[7], extent=[lnmn, lnmx, ltmn, ltmx], aspect='auto')
    axarr[0, 0].set_title('Vegetation cover (August)')
    im.set_cmap('rainbow')
    f.colorbar(im, ax=axarr[0, 0], orientation='horizontal')

    cmap = colors.ListedColormap(
        ['#AF324B', '#C83264', '#C800C8', '#008C00', '#FF7D00', '#96C800', '#96C800', '#C8FF78', "#FFFF00", "#008C00",
         "#008C00", '#008C00', '#6464FF', '#6464FF', '#646496'])
    bounds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    im = axarr[0, 1].imshow(lan, cmap=cmap, norm=norm, extent=[lnmn, lnmx, ltmn, ltmx], aspect='auto')
    axarr[0, 1].set_title('Land Use Type')
    f.colorbar(im, ax=axarr[0, 1], orientation='horizontal')


    # Fine-tune figure; make subplots farther from each other.
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    f.savefig(logfolder+'/surface_visualisation.png')

    # ******************************************************************
    # ******************************************************************
    #                      WRITE TO OUTPUT FILE
    # ******************************************************************
    # ******************************************************************


    # Define NetCDF file
    w_nc_fid = Dataset(params['fname_surface'], 'w', format='NETCDF4')

    # Create dimensions
    w_nc_fid.createDimension('nx', nx)
    w_nc_fid.createDimension('ny', ny)
    w_nc_fid.createDimension('months', 12)


    # make variable
    vid0 = w_nc_fid.createVariable('X', 'f8', ('ny', 'nx'))
    vid1 = w_nc_fid.createVariable('Y', 'f8', ('ny', 'nx'))
    vid2 = w_nc_fid.createVariable('LAT', 'f8', ('ny', 'nx'))
    vid3 = w_nc_fid.createVariable('LON', 'f8', ('ny', 'nx'))
    vid4 = w_nc_fid.createVariable('DEM', 'f8', ('ny', 'nx'))
    vid5 = w_nc_fid.createVariable('LAN', 'f8', ('ny', 'nx'))
    vid6 = w_nc_fid.createVariable('SOIL', 'f8', ('ny', 'nx'))
    vid7 = w_nc_fid.createVariable('ALBV', 'f8', ('ny', 'nx'))
    vid8 = w_nc_fid.createVariable('ALBS', 'f8', ('ny', 'nx'))
    vid9 = w_nc_fid.createVariable('EPSV', 'f8', ('ny', 'nx'))
    vid10 = w_nc_fid.createVariable('EPSS', 'f8', ('ny', 'nx'))
    vid11 = w_nc_fid.createVariable('Z0M', 'f8', ('ny', 'nx'))
    vid12 = w_nc_fid.createVariable('RST', 'f8', ('ny', 'nx'))
    vid13 = w_nc_fid.createVariable('RP', 'f8', ('ny', 'nx'))
    vid14 = w_nc_fid.createVariable('ROOTB', 'f8', ('ny', 'nx'))
    vid15 = w_nc_fid.createVariable('LANDMASK', 'f8', ('ny', 'nx'))
    vid16 = w_nc_fid.createVariable('SEEMASK', 'f8', ('ny', 'nx'))
    vid17 = w_nc_fid.createVariable('AHF', 'f8', ('ny', 'nx'))
    vid18 = w_nc_fid.createVariable('COV', 'f8', ('months', 'ny', 'nx'))
    vid19 = w_nc_fid.createVariable('SST', 'f8', ('ny', 'nx'))


    vid0.setncatts({'long_name': 'X-coordinates of the urbclim grid (User-defined projection)', 'units': 'm'})
    vid1.setncatts({'long_name': 'Y-coordinates of the urbclim grid (User-defined projection)', 'units': 'm'})
    vid2.setncatts({'long_name': 'Latitude of urbclim grid', 'units': 'degrees'})
    vid3.setncatts({'long_name': 'Longitude of urbclim grid', 'units': 'degrees'})
    vid4.setncatts({'long_name': 'Surface height of urbclim grid', 'units': 'meters'})
    vid5.setncatts({'long_name': 'Land use value', 'units': 'scale of 1-15'})
    vid6.setncatts({'long_name': 'Soil sealing fraction of grid cell', 'units': '%'})
    vid7.setncatts({'long_name': 'Vegetation albedo', 'units': '-'})
    vid8.setncatts({'long_name': 'Ground surface albedo', 'units': '-'})
    vid9.setncatts({'long_name': 'Vegetation emissivity', 'units': '-'})
    vid10.setncatts({'long_name': 'Ground surface emissivity', 'units': '-'})
    vid11.setncatts({'long_name': 'Roughness length for momentum', 'units': 'm'})
    vid12.setncatts({'long_name': 'Stomatal resistance', 'units': 's/m'})
    vid13.setncatts({'long_name': 'Total plant resistance', 'units': 'm'})
    vid14.setncatts({'long_name': 'Root distribution coefficient', 'units': '-'})
    vid15.setncatts({'long_name': 'Land mask', 'units': '-'})
    vid16.setncatts({'long_name': 'Sea mask', 'units': '-'})
    vid17.setncatts({'long_name': 'Antropogenic Heat Flux', 'units': 'W/m'})
    vid18.setncatts({'long_name': 'Land cover of grid cell', 'units': '%'})
    vid19.setncatts({'long_name': 'Sea surface temperature', 'units': '%'})

    # add variable data
    w_nc_fid.variables['X'][:,:] = xcor
    w_nc_fid.variables['Y'][:,:] = ycor
    w_nc_fid.variables['LAT'][:,:] = lat
    w_nc_fid.variables['LON'][:,:] = lon
    w_nc_fid.variables['DEM'][:,:] =np.flipud(dem)
    w_nc_fid.variables['LAN'][:,:] = np.flipud(lan)
    w_nc_fid.variables['SOIL'][:,:] = np.flipud(soil)
    w_nc_fid.variables['ALBV'][:,:] = np.flipud(albv)
    w_nc_fid.variables['ALBS'][:,:] = np.flipud(albs)
    w_nc_fid.variables['EPSV'][:,:] = np.flipud(epsv)
    w_nc_fid.variables['EPSS'][:,:] = np.flipud(epss)
    w_nc_fid.variables['Z0M'][:,:] = np.flipud(z0m)
    w_nc_fid.variables['RST'][:,:] = np.flipud(rst)
    w_nc_fid.variables['RP'][:,:] = np.flipud(rp)
    w_nc_fid.variables['ROOTB'][:,:] = np.flipud(rootb)
    w_nc_fid.variables['LANDMASK'][:,:] = np.flipud(landmask)
    w_nc_fid.variables['SEEMASK'][:,:] = np.flipud(seemask)
    w_nc_fid.variables['AHF'][:,:] = np.flipud(ahf)
    cover_save=np.zeros((12,ny,nx))
    for i in range(12):
        cover_save[i,:,:]=np.flipud(cover[i])
    w_nc_fid.variables['COV'][:,:,:] = cover_save

    # close the file
    w_nc_fid.close()



def read_tiff(filename):
    #im = Image.open(filename)
    #tiff = np.array(im)
    data=gdal.Open(filename)
    tiff= data.ReadAsArray()
    return tiff


def get_smoothed_value(arr, window, x, y, nan):
    offset = window / 2

    xmin = x - offset if x - offset > 0 else 0
    ymin = y - offset if y - offset > 0 else 0
    xmax = x + offset if x + offset < len(arr) else len(arr) - 1
    ymax = y + offset if y + offset < len(arr[0]) else len(arr[0]) - 1

    conv_arr = arr[xmin:xmax + 1, ymin: ymax + 1]
    sum = 0.0
    count = 0

    for (i, row) in enumerate(conv_arr):
        for (j, value) in enumerate(row):
            if (i == x and y == j) or value == nan:
                continue;
            sum += value
            count += 1

    if count == 0:
        print 'Problem with the map, area with to many undefined values'
        raise Exception('Problem with the map, area with to many undefined values')
        
    return sum / count


def replace_values(array, values, type):

    array = np.copy(array).astype(type)
    for x in range(0, len(array)):
        for y in range(0, len(array[x])):
            array[x][y] = type(values[int(array[x][y])])
    return array


def execute(command):
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ""
    for line in p.stdout.readlines():
        output += line + "\n"
    retval = p.wait()
    #print command, "executed with return value ", retval
    return output


