import argparse, sys, fileinput,os, shutil
import pandas as pd
sys.path.append('Scripts/Terrain')
from terrain2urbclim import terrain2urbclim 


"""
usage: Masterscript.py [-h] folder

Get UrbClim data

positional arguments:
  folder      Runfolder

optional arguments:
  -h, --help  show this help message and exit

  
to do 
-----
* versie maken waarbij data naar VSC te verhuizen is (relatieve paden gebruiken tov. startpunt)
"""


########################
### Global variables ###
########################

executable=os.path.abspath('Inputfiles/urbclim')
input_excel='InputSimulations.xlsx'
config_input=os.path.abspath('Inputfiles/generic_config.cfg')
script_meteo_matlab=os.path.abspath('Scripts/Meteo/run_era52urbclim_bash.sh')
script_meteo_matlab_aux=os.path.abspath('Scripts/Meteo/era52urbclim_bash')
script_meteo_python=os.path.abspath('Scripts/Meteo/DownloadMeteo.py')
script_meteo_add_nt=os.path.abspath('Scripts/Meteo/Add_nt.py')

####################
### Main Routine ###
####################

def main(args):
    # compose config files and process input
    print 'Gather information and compose config files\n\n'
    runfolder=os.path.abspath(args.runfolder)
    simulations = SetUp(runfolder)
    
    # Compose Terrain data
    print 'Make terrain files'
    ComposeTerrainInput(simulations)
    
    # get meteorological data (only once for every city)
    print 'Compose download of meteorological data \n\n'
    ComposeMeteoInput(simulations)   
    
    
    print 'Make Launchscripts \n\n'
    MakeRuns(simulations,runfolder)
    

##########################
### Startup processing ###
########################## 
def SetUp(runfolder):
    global input_excel
    # read Excel file
    input_simulations=pd.read_excel(input_excel)   
        
    # run- and meteofolders
    input_simulations['Runfolder']=runfolder+'/'+input_simulations['City']+'_'+ \
        input_simulations['Resolution'].astype(str)+'/'
    input_simulations['Meteofolder']=runfolder+'/Meteo/'+input_simulations['City']+'/'
    input_simulations['Configfolder']=input_simulations['Runfolder']+'/configfiles/'
    input_simulations['Terrainfolder']=input_simulations['Runfolder']+'/terrain/'
    input_simulations['Results']=input_simulations['Runfolder']+'results/'
    input_simulations['Scriptsfolder']=input_simulations['Runfolder']+'/scripts/'
    
    # years and months
    input_simulations['Years']=input_simulations.apply(lambda p: range(p.Startyear,p.Endyear+1),axis=1)
    input_simulations['Months']=input_simulations.apply(lambda p: range(p.Startmonth,p.Endmonth+1),axis=1)
    
    # configuration files
    input_simulations['Configlist']=''
    for index,row in input_simulations.iterrows():
        configlist= MakeConfigfiles(row)
        configlist_string=",".join(configlist)
        input_simulations.loc[index,'Configlist']= configlist_string

        
    return input_simulations
         
######################
###  Terain Input  ###
###################### 
def ComposeTerrainInput(simulations):
    for index,simulation in simulations.iterrows():
        print simulation['City']+' '+str(simulation['Resolution'])
        # make folder and copy files
        terrain_folder=simulation['Terrainfolder']
        makeFolder(terrain_folder)
        terrain_configfile=terrain_folder+'config_general.cfg'
        shutil.copy(simulation.Configfolder+'config_general.cfg',terrain_configfile)
        terrain2urbclim(terrain_configfile)
        
 
######################
###  Meteo input   ###
###################### 
 
def ComposeMeteoInput(simulations):
    # gather links to scripts
    global script_meteo_matlab
    global script_meteo_python
    global script_meteo_matlab_aux
    
    # chmodding
    os.system('chmod +x '+script_meteo_matlab) 
    os.system('chmod +x '+script_meteo_matlab_aux)
    
    # get unique cities
    cities=simulations.drop_duplicates('City')
    for index,city in cities.iterrows():        
        # make folder for meteo
        makeFolder(city.Meteofolder)
        
        # output of script
        out_file=city.Meteofolder+'out.txt'
        
        # make bash script content
        content=MakeMeteoScript(city.Configlist,script_meteo_python,script_meteo_matlab,out_file)
        
        # bash script locations
        location_script=city.Meteofolder+'DownloadScript.sh'
        
        # save script
        file=open(location_script,'w')
        file.write(content)
        file.close()
        
        # launch donwload of meteo. If on=> automatically download meteo.
        #os.system('qsub '+location_script)

######################
###   Start runs   ###
######################
def MakeRuns(simulations,runfolder):
    # get executable
    global executable
    global script_meteo_add_nt
    
    # chmod executable
    os.system('chmod +x '+executable) 
    
    # launchscript for all simulations
    launchscript_combined='' 
    
    for index,simulation in simulations.iterrows():
        # make resultsfolder
        makeFolder(simulation.Results)
    
        # folder for scripts
        scriptsfolder=simulation.Scriptsfolder
        makeFolder(scriptsfolder)
                
        # combined launchscript
        launchscripts=''       
        
        for year in simulation.Years:
            # script for year
            config_file=simulation.Configfolder+'config_'+str(year)+'.cfg'
            out_file=scriptsfolder+'out_'+str(simulation.Resolution)+'_'+str(year)+'.txt'
            script_content=MakeScript(config_file,out_file,executable,script_meteo_add_nt)
            script_file=scriptsfolder+'script_'+str(year)
            
            # add to combined launchscript
            launchscripts=launchscripts+'qsub '+script_file+'\n'
            launchscript_combined=launchscript_combined+'qsub '+script_file+'\n'            
            
            # save script
            file=open(script_file,'w')
            file.write(script_content)
            file.close()
        
        # save combined Launchscript
        launchscript_file=scriptsfolder+'LaunchAllRuns.sh'
        file=open(launchscript_file,'w')
        file.write(launchscripts)
        file.close()
        
    # save combined Launchscript
    launchscript_file=runfolder+'/LaunchAllRuns.sh'
    file=open(launchscript_file,'w')
    file.write(launchscript_combined)
    file.close()

"""
MAJOR AUXILIARY ROUTINES#
"""      
   
#########################
### ComposeConfigFile ###
#########################
def MakeConfigfiles(rundata):
    global config_input
    # make folder and generic file
    makeFolder(rundata.Configfolder)
    configfile_general=rundata.Configfolder+'config_general.cfg'
    shutil.copy(config_input,configfile_general)
    
    # add parameters that are the same for all years
    AddKeyWord(configfile_general,'city',rundata.City)
    AddKeyWord(configfile_general,'clon',rundata.Longitude)
    AddKeyWord(configfile_general,'clat',rundata.Latitude)
    AddKeyWord(configfile_general,'nx',rundata.nx)
    AddKeyWord(configfile_general,'ny',rundata.ny)
    AddKeyWord(configfile_general,'dxy',rundata.Resolution)
    AddKeyWord(configfile_general,'startmonth',rundata.Startmonth)
    AddKeyWord(configfile_general,'endmonth',rundata.Endmonth)
    AddKeyWord(configfile_general,'building_height',rundata.BuildingHeight)
    AddKeyWord(configfile_general,'urban_tree_fraction',rundata.UrbanTreeFraction)
    AddKeyWord(configfile_general,'timecor',rundata.TimeFactor)
    fname_surface=rundata.Terrainfolder+'urbclim_surface_'+rundata.City+'_'+str(rundata.Resolution)+'.nc'
    AddKeyWord(configfile_general,'fname_surface',fname_surface)
    configlist=[]
    for year in rundata.Years:
        # add parameters that differ per year
        configfile_year=rundata.Configfolder+'config_'+str(year)+'.cfg'
        configlist.append(configfile_year)
        shutil.copy(configfile_general,configfile_year)
        AddKeyWord(configfile_year,'year',year)
        
        fname_output=rundata.Results+'urbclim_output_'+rundata.City+'_'+str(rundata.Resolution)+'_'+str(year)+'.nc'
        AddKeyWord(configfile_year,'fname_output',fname_output)
        fname_ecmwf=rundata.Meteofolder+'urbclim_meteo_'+rundata.City+'_'+str(year)+'.nc'
        AddKeyWord(configfile_year,'fname_ecmwf',fname_ecmwf)
        
        fname_ecmwf_surface=rundata.Meteofolder+'urbclim_meteo_'+rundata.City+'_'+str(year)+'_2d.nc'
        AddKeyWord(configfile_year,'fname_atmo2D',fname_ecmwf_surface)
        fname_ecmwf_profiles=rundata.Meteofolder+'urbclim_meteo_'+rundata.City+'_'+str(year)+'_3d.nc'
        AddKeyWord(configfile_year,'fname_atmo3D',fname_ecmwf_profiles)
        
    return configlist

    
##########################
### ComposeMeteoscript ###
##########################

def MakeMeteoScript(configlist,script_meteo_python,script_meteo_matlab,out_file):
    content='#!/bin/bash \n'+\
                '#PBS -N meteodownload \n'+\
                '#PBS -m e \n'+\
                'python '+script_meteo_python+' '+configlist+' '+script_meteo_matlab+' &> '+out_file
    
    return content

###########################
### ComposeLaunchscript ###
###########################
        
def MakeScript(config_file,out_file,executable,script_meteo_add_nt):
    content= '#!/bin/bash \n'+\
                '#PBS -N '+config_file+'\n'+\
                '#PBS -m e \n'+\
                'export LD_LIBRARY_PATH=/tools/maphdf5/1.0.1/x86_64/lib:/tools/hdf5/1.8.8/x86_64/lib:/opt/intel/compiler/lib/intel64:/tools/netcdf/4.1.3/x86_64/lib \n'+\
                'ulimit -s unlimited \n'+\
                'python '+script_meteo_add_nt+' '+config_file+'\n'+\
                executable+' '+config_file+' > '+out_file
                
    return content

        

"""
MINOR AUXILIARY ROUTINES#
"""
 
def AddKeyWord(filename,keyword,keydata):
    for line in fileinput.input(filename, inplace=True):
        if line.strip().startswith(keyword):
            line = keyword+' = '+str(keydata)+'\n'
        sys.stdout.write(line)

def makeFolder(folder):        
    if not os.path.exists(folder):
        os.makedirs(folder)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get OSM data')
    parser.add_argument('runfolder', metavar='folder', type=str, help='Runfolder')
    args = parser.parse_args()

    main(args)
