# -*- coding: utf-8 -*-
"""
Preprocessing routine for IFDM
==============================

@author: hooyberh (hans.hooyberghs@vito.be)

Basic description: Preprocessing routine for IFDM-OSPM. This preprocessor sets up
the IFDM-OSPM runs and starts running them if executed on a Linux torque system. 


Input parameters passed to constructor:
    None. Parameters must be set by hand (see example RunPreprocesser script
    
Parameters linked to class: see constructor. There are different types of input parameters:
- Compulsory that be provided by hand in the RunPreprocessor.py script 
- Optional parameters. They describe the IFDM-run that will be started. Some of these
    become compulsory in certain cases (see below for details) 
- Other input parameters: see constructor. Only change these if you know what you
    are doing (don't tell me that I didn't warn you ;) )

Main routines:
* Constructor
* setUp: check and process input parameters that must be said by hand
* makeMeteo: makes meteo files using meteo class
* makeEmissions: make emissions using emissions class
* makeFileTree: makes the file tree for the IFDM runs
* RunIFDM: starts the IFDM runs on the cluster
* runAll: do all of the above (except the constructor)

Major auxiliary routines: (description: see below)
* makeFileTreeSector
* makeFolderForRun
* runIFDMSector
* makeSymbolicLinksCTM
* makeSymbolicLinksIFDM

Minor auxiliary routines: see below
    
Dependencies:
* Python routines:
    -os
    -glob
    -pandas
    -shutil
    -math
    -Emissions
    -Meteo
* External software:
    - torque cluster (for starting the runs)
        
Output: 
* Emissions input for IFDM


===

TO DO LIST
----------
If we really have time
* solution for empty point sources file in building the raster
* allow meteo input
"""

# load Python libraries
import shutil
from Emissions import Emissions
from Meteo import Meteo
import os
import sys
import math
import pandas as pd
import datetime


class Preprocessor:
    
#################################      
# Constructor and main parameters
#################################

    def __init__(self):
# Always set the following parameters in RunPreprocessor
 #######################################################
        # variables that must be set by hand in the RunPostprocessor script
        self.runfolder=False         # folder where IFDM runs will take place
        self.inputfolder=False          # folder where inputdata can be found
        self.outputfolder=None         # folder where outputfiles are copied to
        self.pollutants=False           # list of pollutants
        self.lengthOfPeriod=744      # length of meteorological periods in hours

# Optionally set the following set of parameters
 #######################################################    
        
        self.allSectorsRun=True         # runs for all sectors combined (default is True)
        
        # variables for OSPM
        self.OSPM=False                 # OSPM? (default is False)
        self.generateWarning=True
        
        # variables for rugosity
        self.Rugosity=False
                
        # variables for coupling to background
        self.processBackground=False    # Processing of background needed?
        self.CTMfolder=False             # folder with CTM details
        self.FixedBackground=None
        
        
        # variables for start and enddate
        self.startDate=None
        self.endDate=None
        
        # variables for MeteoPreprocessing
        self.meteoPreprocess=False
        self.meteoCopyAssimil=False
        
        #shipping?
        self.Shipping=False

        
        
# Only changes these parameters if you know what you're doing
 ############################################################
        # location of important input files
        self.inputfileEmissions='Emissions.xlsx'      
        self.inputfileMeteo='Meteo.xlsx'
        self.inputfileRaster='CompleteGrid.txt'
        self.inputfileNewRaster='New_InputRaster.txt'
        self.inputfileRasterOSPM='InputRasterOSPM.txt'
        self.inputfileTimeSeries='Input_Timeseries.txt'  
        
        # location of default IFDM input files
        self.inputfileTimeFactors='default/Timefactors.txt'
        self.inputfileBuilding='default/InputBuilding.txt'
        self.inputfileExecutable='default/ifdm'
        self.inputfileLimitValues='default/InputLimitValues.txt'
        self.inputfileStraatgegevens='default/straatgegevens_config_exc.txt'
        self.inputfileRugosity='default/Rugosity250m.txt'
        
        # location of default input config files
        self.inputfileConfig='default/ifdm_add.conf'
        self.inputfileConfigCoupling='default/ifdm_coupling.conf'
        self.inputfileConfigOSPM='default/ifdm_ospm.conf'
        
        # location of default input cluster launchscript file
        self.inputfileLaunchscript='default/Launchscript.txt'
        self.inputfileLaunchscript_ospm='default/Launchscript_OSPM.txt'
        self.inputfileLaunchscript_ospm_fixedBackground='default/Launchscript_OSPM_fixedBackground.txt'
        self.inputfileLaunchscript_VSC='default/Launchscript_VSC.txt'
        self.inputfileLaunchscript_VSC_ospm='default/Launchscript_VSC_OSPM.txt'
              
        
    ##############################        
    # Run All scripts
    ##############################
       
    
    def runAll(self):
        # routine to run all the important subfunctions
        self.setUp()
        self.makeMeteo()
        self.makeEmissions()
        self.makeFileTree()     
        self.RunIFDM()
        

    ##############################        
    # Set-up of the preprocessor
    ##############################

    def setUp(self):
        """
        Process input files and check whether all the files are present
        """
        # make empty sectorlist
        self.sectorList=pd.DataFrame()        
        
        # Process inputfiles to take inputfolder into account
        self.inputfileEmissions=self.inputfolder+self.inputfileEmissions
        self.inputfileMeteo=self.inputfolder+self.inputfileMeteo
        self.inputfileTimeFactors=self.inputfolder+self.inputfileTimeFactors
        self.inputfileBuilding=self.inputfolder+self.inputfileBuilding
        self.inputfileConfig=self.inputfolder+self.inputfileConfig
        self.inputfileConfigCoupling=self.inputfolder+self.inputfileConfigCoupling
        self.inputfileConfigOSPM=self.inputfolder+self.inputfileConfigOSPM
        self.inputfileTimeSeries=self.inputfolder+self.inputfileTimeSeries
        self.inputfileExecutable=self.inputfolder+self.inputfileExecutable
        self.inputfileRaster=self.inputfolder+self.inputfileRaster
        self.inputfileRasterOSPM=self.inputfolder+self.inputfileRasterOSPM    
        self.inputfileNewRaster=self.inputfolder+self.inputfileNewRaster
        self.inputfileLaunchscript=self.inputfolder+self.inputfileLaunchscript
        if self.FixedBackground is None:            
            self.inputfileLaunchscript_ospm=self.inputfolder+self.inputfileLaunchscript_ospm
        else:
            self.inputfileLaunchscript_ospm=self.inputfolder+self.inputfileLaunchscript_ospm_fixedBackground  
        self.inputfileLaunchscript_VSC=self.inputfolder+self.inputfileLaunchscript_VSC
        self.inputfileLaunchscript_VSC_ospm=self.inputfolder+self.inputfileLaunchscript_VSC_ospm
        self.inputfileLimitValues=self.inputfolder+self.inputfileLimitValues
        self.inputfileStraatgegevens=self.inputfolder+self.inputfileStraatgegevens    
        self.inputfileRugosity=self.inputfolder+self.inputfileRugosity
        
        # Make temporary folders
        self.tempfolderbase=self.inputfolder+'temp/'
        self.tempoutputfolderMeteo=self.tempfolderbase+'/Meteo/'
        self.tempoutputfolderEmissions=self.tempfolderbase+'/Emissions/'
        self.tempoutputfolderRIO=self.tempfolderbase+'/RIO/'	
        self.tempfileConfig=self.tempfolderbase+'ifdm_add.conf'
        self.tempfileConfigCoupling=self.tempfolderbase+'ifdm_coupling.conf'
        self.tempfileConfigOSPM=self.tempfolderbase+'ifdm_ospm.conf'
        
        # make the output folders
        
        self.makeFolder(self.tempoutputfolderMeteo)
        self.makeFolder(self.tempoutputfolderEmissions)
        self.makeFolder(self.tempoutputfolderRIO) 	
        
        # check presence of all input files        
        self.checkPresence(self.inputfileTimeFactors)
        self.checkPresence(self.inputfileBuilding)
        self.checkPresence(self.inputfileConfig)
        self.checkPresence(self.inputfileTimeSeries)
        self.checkPresence(self.inputfileExecutable)
        
        # make outputfolder if needed
        if self.outputfolder is not None:
            self.makeFolder(self.outputfolder)
        
        #check rugosity
        if os.path.exists(self.inputfileRugosity):
            self.Rugosity=True
        
                               
        
        """
        Add pollutants to config file
        """
        # make run pollutant names
        self.pollutantsIFDM = [w.ljust(4,'a') for w in self.pollutants]        
        # make list for config file
        self.pollutantList=''        
        for i in range(len(self.pollutantsIFDM)):
            self.pollutantList=self.pollutantList+self.pollutantsIFDM[i]+','
        self.pollutantList=self.pollutantList.rstrip(',')
        self.modifyFile(self.inputfileConfig,self.tempfileConfig,'pollutantslist',self.pollutantList)
        print 'Modified config file to include pollutants'
        
        
        """ 
        Make changes to config files for coupling or scenario run
        """
        # changes that are only needed for coupling run
        self.modifyFile(self.inputfileConfigCoupling,self.tempfileConfigCoupling,'pollutantslist',self.pollutantList)
        if self.FixedBackground is None:    
            self.checkPresence(self.CTMfolder+'concentraties/gridRIO.grid_def')

        """
        Add pollutants to config files for OSPM step
        """
        if self.OSPM:
            self.modifyFile(self.inputfileConfigOSPM,self.tempfileConfigOSPM,'pollutantslist',self.pollutantList)           
            print 'Modified config files for OSPM routine'
      
        """
        InputRaster files
        """     
        # check whether input files are present
        self.checkPresence(self.inputfileRaster)
        self.checkPresence(self.inputfileNewRaster)
         
        if self.OSPM:
            self.checkPresence(self.inputfileRasterOSPM)
                
            # two options for OSPM files
            if not os.path.exists(self.inputfileStraatgegevens):
                self.inputfileStraatgegevens=self.inputfileStraatgegevens.replace('_exc','')
                print 'Did not find Straatgegevens_config_exc.txt'
                print 'Will now try Straatgegevens_config.txt'
                self.checkPresence(self.inputfileStraatgegevens)          

    ##############################        
    # Call meteo class
    ##############################
    
    def makeMeteo(self): 
        # initialize: check if start and enddate are provided
        met=Meteo(self.inputfileMeteo,self.tempoutputfolderMeteo,self.tempoutputfolderRIO,self.lengthOfPeriod,self.meteoCopyAssimil,self.FixedBackground,self.startDate,self.endDate)
              
        # write meteo input to tempfolder
        met.writeMeteo()
        
        
        if self.meteoPreprocess:
            met.preprocessAssimilMeteo(self.meteoLocation)
            
            
        if self.processBackground:
            met.preprocessRIObackground(self.CTMfolder)
        
        # get some important values for the remainder of the script
        self.hoursPerPeriod=self.lengthOfPeriod
        self.hoursLastPeriod=met.hoursLastPeriod
        self.numberOfPeriods=met.NumberOfPeriods
        self.periodStart=met.periodStart
        
        print 'Meteo has been made'
        print '----------------------------------------------'
        print ' '
        
    ##############################        
    # Call emissions class
    ##############################
        
    def makeEmissions(self):
        # initialize
        print self.inputfileEmissions
        em=Emissions(self.inputfileEmissions,self.tempoutputfolderEmissions,self.pollutants,self.Shipping,self.OSPM)
        
        # write emission input to tempfolder
        # emissions input for run with all sectors
        em.makeEmissionsAll()
            
        
        print 'Emissions have been made'
        print '----------------------------------------------'
        print ' '



    ##############################        
    # Make file tree
    ##############################
    def makeFileTree(self):  
               
        """
        Make file tree and runs scripts (use major auxiliary routines)
        """
        print "Make file trees"
        self.makeFileTreeSector('AllSectors')
                            
        print "File trees have been made"


    ####################################
    # Submit the runs to a cluster queue   
    ####################################
    def RunIFDM(self):
    
        print 'Starting runs'
        print 'Note: this will start the actual simulations, do only proceed if you want to run on a cluster!'
        print 'Do you want to proceed?'
        var = raw_input('Enter "yes" to continue: ')
        print "Your choice:", var
        if var == 'yes':  
            sectorlabel='AllSectors'
            self.runIFDMSector(sectorlabel)
            
            print 'All runs finished'
                
        else:
            print 'Exiting program.'
            exit()
            
            
    ####################################
    # Submit the runs to the VSC cluster   
    ####################################

    def RunIFDM_VSC(self):
    
        var='yes'
        if var == 'yes':  
            sectorlabel='AllSectors'
            self.runIFDMSector_VSC(sectorlabel)                
                            
        else:
            print 'Exiting program.'
            exit()
            
    
################################       
# Major auxiliary routines
################################       

    
        
    def makeFileTreeSector(self,sectorlabel_run):
        # make the file tree for a specific sector    
        savefolder_runscripts=self.runfolder+'/'+sectorlabel_run+'/Scripts/'
        self.makeFolder(savefolder_runscripts)
        launchscript_allruns_data=''
        if self.OSPM:        
            launchscript_allruns_data_ospm=''

        for i in range(self.numberOfPeriods):
            # make folder
            nrString=self.periodStart[i]
            meteostart_hours=i*self.hoursPerPeriod+1 
            numberMeteo=self.hoursPerPeriod
            if i == self.numberOfPeriods-1:
                numberMeteo = self.hoursLastPeriod
            # copy all files to runfolder (most processing in subroutine)
            savefolder_run=self.runfolder+sectorlabel_run+'/'+nrString+'/'
            self.CopyFilesForRun(sectorlabel_run,nrString,savefolder_run)
                   
            # Changes for coupling to config files
            meteostart_ctm=str(meteostart_hours*3600)
            config_add=savefolder_run+'ifdm_add.conf'
            self.modifyFile(self.tempfileConfig,config_add,'ctm_start',str(meteostart_ctm))
            config_coupling=savefolder_run+'ifdm_coupling.conf'
            self.modifyFile(self.tempfileConfigCoupling,config_coupling,'meteostart',str(meteostart_hours))
            
            # copy files for coupling to assimil meteo.
            if self.meteoCopyAssimil:
                self.CopyMeteoFiles(savefolder_run,self.periodStart[i],numberMeteo)
                meteofile_run=self.tempoutputfolderMeteo+'Meteofile'+self.periodStart[i]+'.txt' 
                shutil.copy(meteofile_run,savefolder_run+'InputMeteoType1.txt')
            else:
                self.modifyFile(config_add,config_add,'meteostart',str(meteostart_hours))
                meteofile_run=self.tempoutputfolderMeteo+'Meteofile'+self.periodStart[i]+'.txt' 
                shutil.copy(meteofile_run,savefolder_run+'InputMeteoType1.txt')
            
            # background coupling
            if self.FixedBackground is not None:
                backgroundfile_run=self.tempoutputfolderMeteo+'Background'+self.periodStart[i]+'.txt' 
                shutil.copy(backgroundfile_run,savefolder_run+'InputBackground.txt')
            else:
                self.makeSymbolicLinksCTM(savefolder_run,nrString,numberMeteo)
            
            # runscript             

            savefile_launchscript=savefolder_runscripts+'Launch_'+sectorlabel_run+'_'+nrString                   
            self.modifyFile(self.inputfileLaunchscript,savefile_launchscript,'folder_base',self.runfolder+'/'+sectorlabel_run)
            self.modifyFile(savefile_launchscript,savefile_launchscript,'outputlocation',self.outputfolder)
            self.modifyFile(savefile_launchscript,savefile_launchscript,'date',nrString)
            launchscript_allruns_data=launchscript_allruns_data+'echo '+nrString+'\n'
            launchscript_allruns_data=launchscript_allruns_data+'bash '+savefile_launchscript+'\n'
            
            
            # in case of OSPM run    
            if self.OSPM: 
                # make folder                
                savefolder_runospm=self.runfolder+sectorlabel_run+'/'+nrString+'ospm/'                 
                self.CopyFilesForRun(sectorlabel_run,nrString,savefolder_runospm)
                self.makeSymlink(self.inputfileStraatgegevens,savefolder_runospm+'straatgegevens_config_exc.txt')
                shutil.copy(self.inputfileRasterOSPM,savefolder_runospm+'InputRaster.txt')
                
                # changes to config file
                config_OSPM=savefolder_runospm+'ifdm_ospm.conf'
                self.modifyFile(self.tempfileConfigOSPM,config_OSPM,'meteostart',str(meteostart_hours))    
                self.modifyFile(config_OSPM,config_OSPM,'ctm_start',str(meteostart_ctm))                
                
                # symbolic links between IFDM and OSPM folder
                if self.FixedBackground is not None:
                    backgroundfile_run=self.tempoutputfolderMeteo+'Background'+self.periodStart[i]+'.txt' 
                    shutil.copy(backgroundfile_run,savefolder_runospm+'InputBackground.txt')
                else:
                    self.makeSymbolicLinksIFDM(savefolder_runospm,savefolder_run,numberMeteo)  

                
                # copy files for coupling to assimil meteo.
                if self.meteoCopyAssimil:
                    self.CopyMeteoFiles(savefolder_runospm,self.periodStart[i],numberMeteo)
                    meteofile_run=self.tempoutputfolderMeteo+'Meteofile'+self.periodStart[i]+'.txt' 
                    shutil.copy(meteofile_run,savefolder_runospm+'InputMeteoType1.txt')
                else:
                    self.modifyFile(config_add,config_add,'meteostart',str(meteostart_hours))
                    meteofile_run=self.tempoutputfolderMeteo+'Meteofile'+self.periodStart[i]+'.txt' 
                    shutil.copy(meteofile_run,savefolder_runospm+'InputMeteoType1.txt')               
                
                # runscript
                savefile_launchscript_ospm=savefolder_runscripts+'Launch_'+sectorlabel_run+'_'+nrString+'_ospm'  
                self.modifyFile(self.inputfileLaunchscript_ospm,savefile_launchscript_ospm,'folder_base',self.runfolder+'/'+sectorlabel_run)
                self.modifyFile(savefile_launchscript_ospm,savefile_launchscript_ospm,'outputlocation',self.outputfolder)
                self.modifyFile(savefile_launchscript_ospm,savefile_launchscript_ospm,'date',nrString)
                launchscript_allruns_data_ospm=launchscript_allruns_data_ospm+'echo '+nrString+'\n'
                launchscript_allruns_data_ospm=launchscript_allruns_data_ospm+'bash '+savefile_launchscript_ospm+'\n'
       
       # save runscripts
        launchscript_allruns= open(self.runfolder+'Launch_'+sectorlabel_run+'_ifdm','w')
        launchscript_allruns.write(launchscript_allruns_data)
        launchscript_allruns.close()        
        
        if self.OSPM:    
            launchscript_allruns= open(self.runfolder+'Launch_'+sectorlabel_run+'_ifdm_ospm','w')
            launchscript_allruns.write(launchscript_allruns_data_ospm)
            launchscript_allruns.close()
        
    
    def CopyFilesForRun(self,sectorLabel,meteoLabel,savefolder):  
        # make a single run folder
        self.makeFolder(savefolder)    
        
        #files to copy
        shutil.copy(self.inputfileNewRaster,savefolder+'InputRaster.txt')
        shutil.copy(self.inputfileTimeFactors,savefolder+'Timefactors.txt')
        shutil.copy(self.inputfileLimitValues,savefolder+'InputLimitValues.txt')
        shutil.copy(self.inputfileBuilding,savefolder+'InputBuilding.txt')     
        shutil.copy(self.inputfileTimeSeries,savefolder+'Input_Timeseries.txt')
        shutil.copy(self.inputfileExecutable,savefolder+'ifdm')

        # files to symlink (large files)   
        surfacesourcesT1_run=self.tempoutputfolderEmissions+'InputSurfaceSourcesT1_'+sectorLabel+'.txt'
        self.makeSymlink(surfacesourcesT1_run,savefolder+'InputSurfaceSourcesT1.txt')
        surfacesourcesT2_run=self.tempoutputfolderEmissions+'InputSurfaceSourcesT2_'+sectorLabel+'.txt'
        self.makeSymlink(surfacesourcesT2_run,savefolder+'InputSurfaceSourcesT2.txt')
        self.makeSymlink(self.inputfileRaster,savefolder+'CompleteGrid.txt')        
        linesources_run=self.tempoutputfolderEmissions+'InputLineSources_'+sectorLabel+'.txt'
        self.makeSymlink(linesources_run,savefolder+'InputLineSources.txt')        
        pointsources_run=self.tempoutputfolderEmissions+'InputPointSources_'+sectorLabel+'.txt'
        self.makeSymlink(pointsources_run,savefolder+'InputPointSources.txt')
        
        if self.Rugosity:
            self.makeSymlink(self.inputfileRugosity,savefolder+'Rugosity250m.txt')

    # run IFDM for a specific sector        
    def runIFDMSector(self,sectorlabel):        
        if self.OSPM:
            runfile=self.runfolder+'Launch_'+sectorlabel+'_ifdm_ospm'
            os.system('chmod +x '+runfile)
            os.system(runfile)     
        else:
            runfile=self.runfolder+'Launch_'+sectorlabel+'_ifdm'
            os.system('chmod +x '+runfile)
            os.system(runfile)         
            
            
    # run IFDM for a specific sector on the VSC-cluster        
    def runIFDMSector_VSC(self,sectorlabel):
        folder_sector=self.runfolder+'/'+sectorlabel+'/'
        launchscript_VSC=folder_sector+'Launchscript_VSC.txt'
        if self.OSPM:
            self.modifyFile(self.inputfileLaunchscript_VSC_ospm,launchscript_VSC,'folder_base',folder_sector)
        else:
            self.modifyFile(self.inputfileLaunchscript_VSC,launchscript_VSC,'folder_base',folder_sector)
        self.modifyFile(launchscript_VSC,launchscript_VSC,'outputlocation',self.outputfolder)
        self.periodStart.name='date'
        self.periodStart.to_csv(folder_sector+'dates.txt',index=False,header=True)
        launch_command='wsub -A lp_h_ircel  -batch '+launchscript_VSC+' -data '+folder_sector+'dates.txt'
        print 'Submit runs command: '+launch_command
        os.system(launch_command)

            
    # make symbolic links to CTM-output
    def makeSymbolicLinksCTM(self,savefolder_run,nrString,numberMeteo):
        startdate=datetime.datetime.strptime(nrString,'%Y%m%d%H')
        enddate=startdate+datetime.timedelta(hours=numberMeteo-1)
        cnt=3600
        for day in pd.date_range(startdate,enddate,freq='H'):
            uitvoerfile=savefolder_run+"Aurora_conc"+str(cnt).zfill(6)
            filename=day.strftime('%Y%m%d%H')+'Aurora_conc'
            self.checkPresence(self.tempoutputfolderRIO+filename)	
            shutil.move(self.tempoutputfolderRIO+filename,uitvoerfile)	
            cnt=cnt+3600
            
        folder_griddef=savefolder_run+'concentraties/'
        self.makeFolder(folder_griddef)
        uitvoerfile=folder_griddef+"gridRIO.grid_def"
        self.makeSymlink(self.tempoutputfolderRIO+'concentraties/gridRIO.grid_def',uitvoerfile)
        
        
    # copy files for assimil meteo
    def CopyMeteoFiles(self,savefolder_run,nrString,numberMeteo):
        file_in=open(savefolder_run+'lijst','w')        
        file_in.write('27 \n')
        self.makeSymlink(self.meteoLocation+'GridOlivier.txt',savefolder_run+'GridOlivier.txt')
        startdate=datetime.datetime.strptime(nrString,'%Y%m%d%H')
        enddate=startdate+datetime.timedelta(hours=numberMeteo-1)
        for day in pd.date_range(startdate,enddate,freq='H'):
            filename='MeteoAssimil_'+day.strftime('%Y%m%d%H')+'.txt'
            self.makeSymlink(self.tempoutputfolderMeteo+filename,savefolder_run+filename)
            file_in.write("'"+filename+"' \n")
        file_in.close()


    # make symbolic links to parent folder for OSPM or scenario runs
    def makeSymbolicLinksIFDM(self,savefolder,savefolder_coupling,numberMeteo):       
        # symbolic links for normal pollutants
        for pol in self.pollutantsIFDM:
            for i in range(numberMeteo):
                invoerfile=savefolder_coupling+"uitvoer_AuroPod_"+pol+"_"+str(i+1).zfill(4)+".txt"
                uitvoerfile=savefolder+"uitvoer_AuroPod_"+pol+"_"+str(i+1).zfill(4)+".txt"
                self.checkPresenceSymlink(uitvoerfile)
                self.makeSymlink(invoerfile,uitvoerfile)
        
        # symbolic links for chemistry (ozone, NO2 and coupling constants)
        for i in range(numberMeteo):
            invoerfile=savefolder_coupling+"uitvoer_AuroPod_NO2"+str(i+1).zfill(4)+".txt"
            uitvoerfile=savefolder+"uitvoer_AuroPod_NO2"+str(i+1).zfill(4)+".txt"

            self.makeSymlink(invoerfile,uitvoerfile)
            
            invoerfile=savefolder_coupling+"uitvoer_AuroPod_O3"+str(i+1).zfill(4)+".txt"
            uitvoerfile=savefolder+"uitvoer_AuroPod_O3"+str(i+1).zfill(4)+".txt"

            self.makeSymlink(invoerfile,uitvoerfile)
            
            invoerfile=savefolder_coupling+"uitvoer_AuroPod_rjk"+str(i+1).zfill(4)+".txt"
            uitvoerfile=savefolder+"uitvoer_AuroPod_rjk"+str(i+1).zfill(4)+".txt"

            self.makeSymlink(invoerfile,uitvoerfile)
   
   
################################       
# Minor auxiliary routines
################################       
        
    # make a folder (check first for existence to reduce warnings)            
    def makeFolder(self,folder):        
        if not os.path.exists(folder):
            os.makedirs(folder)
            
    # check presence of files (and exits the program if the file can't be found)        
    def checkPresence(self,file_in):        
        if not os.path.exists(file_in):
            print 'Computer says "no".'
            print 'Input file not found error'
            print 'File '+file_in+' does not exist'            
            print 'Exiting program...'
            sys.exit()
            
    # replace a word in a file    
    def modifyFile(self,inputFile,outputFile,text_in,text_out):        
        file_in=open(inputFile)        
        inputtext=file_in.read()
        file_in.close()        
        outputtext=inputtext.replace(text_in,text_out)
        filesave=open(outputFile,'w')
        filesave.write(outputtext)
        filesave.close()
        
    def checkPresenceSymlink(self,file_in):
        if os.path.islink(file_in):
            os.remove(file_in)
            if self.generateWarning:
                print 'WARNING: symlink exists, overwriting...'
                print 'This warning will only appear once. \n \n'
                self.generateWarning=False
                
                
    def makeSymlink(self,file_in,file_out):
        self.checkPresenceSymlink(file_out)
        os.symlink(file_in,file_out)
                
            

