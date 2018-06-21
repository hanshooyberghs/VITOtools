from netCDF4 import Dataset
import argparse, sys, fileinput,os, shutil

# script to add the number of time steps to the config file

parser = argparse.ArgumentParser(description='Add number of time steps to config file')
parser.add_argument('config', metavar='file', type=str, help='Config_file')
args = parser.parse_args()

config_file_input=args.config
# read config file for meteo-file
with open(config_file_input, 'r') as config_file:
    for line in config_file:
        if len(line.strip()) <> 0 and line.strip()[0] != '#':
            splitted = line.rsplit('=')
            var = splitted[0].strip()
            val = splitted[1].strip()
            if var == 'fname_ecmwf':
                meteofile=val      

# read meteo file
file = Dataset(meteofile)
nrsteps=len(file.variables['YEAR'][:])               

# save changes
for line in fileinput.input(config_file_input, inplace=True):
    if line.strip().startswith('nt'):
        line = 'nt = '+str(nrsteps)+'\n'
    sys.stdout.write(line)