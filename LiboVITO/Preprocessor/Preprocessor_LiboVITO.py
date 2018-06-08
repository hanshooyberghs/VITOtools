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
    -Meteofall
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
        self.outputfolder=False         # folder where IFDM runs will take place
        self.inputfolder=False          # folder where inputdata can be found
        self.pollutants=False           # list of pollutants
        self.lengthOfPeriod=744      # length of meteorological periods in hours

# Optionally set the following set of parameters
 #######################################################    
        
        # variables for run that could be changed but have default values
        self.individualSectorRuns=True  # runs per sector (default is True)
        self.withoutIndividualSectorRuns=False  # runs without sector (default if False)
        self.allSectorsRun=True         # runs for all sectors combined (default is True)
        
        
# Only changes these parameters if you know what you're doing
 ############################################################
        # location of important input files
        self.inputfileEmissions='Emissions.xlsx'      
        self.inputfileMeteo='Meteo.xlsx'
        self.inputfileRaster='CompleteGrid.txt'
        self.inputfileNewRaster='New_InputRaster.txt'
        self.inputfileRasterOSPM='InputRasterOSPM.txt'
        self.inputfileBuilding='default/InputBuilding.txt'
        
        # location of default IFDM input files
        self.inputfileTimeFactors='default/Timefactors.txt'
        self.inputfileTimeSeries='default/Input_Timeseries.txt'   
        self.inputfileExecutable='default/ifdm'
        self.inputfileLimitValues='default/InputLimitValues.txt'
        
        # location of default input config files
        self.inputfileConfig='default/ifdm_add.conf'
        
        # location of default input cluster launchscript file
        self.inputfileLaunchscript='default/Launchscript.txt'
              
        
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
        self.inputfileBuilding=self.inputfolder+self.inputfileBuilding
        self.inputfileTimeFactors=self.inputfolder+self.inputfileTimeFactors
        self.inputfileConfig=self.inputfolder+self.inputfileConfig
        self.inputfileTimeSeries=self.inputfolder+self.inputfileTimeSeries
        self.inputfileExecutable=self.inputfolder+self.inputfileExecutable
        self.inputfileRaster=self.inputfolder+self.inputfileRaster
        self.inputfileNewRaster=self.inputfolder+self.inputfileNewRaster
        self.inputfileLaunchscript=self.inputfolder+self.inputfileLaunchscript
        self.inputfileLimitValues=self.inputfolder+self.inputfileLimitValues
        
        # Make temporary folders
        self.tempfolderbase=self.inputfolder+'temp/'
        self.tempoutputfolderMeteo=self.tempfolderbase+'/Meteo/'
        self.tempoutputfolderEmissions=self.tempfolderbase+'/Emissions/'
        self.tempfileConfig=self.tempfolderbase+'ifdm_add.conf'
        
        # make the output folders
        self.makeFolder(self.tempoutputfolderMeteo)
        self.makeFolder(self.tempoutputfolderEmissions)
        
        # check presence of all input files        
        self.checkPresence(self.inputfileTimeFactors)
        self.checkPresence(self.inputfileBuilding)
        self.checkPresence(self.inputfileConfig)
        self.checkPresence(self.inputfileTimeSeries)
        self.checkPresence(self.inputfileExecutable)
        
        
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
        InputRaster files
        """     
        # check whether input files are present
        self.checkPresence(self.inputfileRaster)
        self.checkPresence(self.inputfileNewRaster)   

    ##############################        
    # Call meteo class
    ##############################
    
    def makeMeteo(self):
        # initialize
        met=Meteo(self.inputfileMeteo,self.tempoutputfolderMeteo,self.lengthOfPeriod)
        
        # write meteo input to tempfolder
        met.writeMeteo()
        
        # get some important values for the remainder of the script
        self.hoursPerPeriod=met.hoursPerPeriod
        self.hoursLastPeriod=met.hoursLastPeriod
        self.numberOfPeriods=met.NumberOfPeriods
        
        print 'Meteo has been made'
        print '----------------------------------------------'
        print ' '
        
    ##############################        
    # Call emissions class
    ##############################
        
    def makeEmissions(self):
        # initialize
        print self.inputfileEmissions
        em=Emissions(self.inputfileEmissions,self.tempoutputfolderEmissions,self.pollutants)
        self.sectorList=em.SectorList
        
        # write emission input to tempfolder
        # emissions input for run with all sectors
        em.makeEmissionsAll()
            
        # emission input for runs per sector
        if self.individualSectorRuns:    
            em.makeEmissionsRunsPerSector()
        
        # emission input for run leaving one sectro out
        if self.withoutIndividualSectorRuns:
            em.makeEmissionsRunsWithoutSector()     
        
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
        if self.allSectorsRun:    
            print "Make file tree for run with all sectors"
            self.makeFileTreeSector('AllSectors')
            
                
        if self.individualSectorRuns:
            print "Make file tree for sector runs"
            for sector in self.sectorList['Sector']:
                sectorlabel_run=sector
                self.makeFileTreeSector(sectorlabel_run)
            
        if self.withoutIndividualSectorRuns:
            print "Make file tree for runs leaving one sector out"
            for sector in self.sectorList['Sector']:
                sectorlabel_run='without'+sector
                self.makeFileTreeSector(sectorlabel_run)
                
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
            if self.allSectorsRun:   
                sectorlabel='AllSectors'
                self.runIFDMSector(sectorlabel)
                
            if self.individualSectorRuns:
                for sector in self.sectorList['Sector']:
                    sectorlabel_run=sector 
                    self.runIFDMSector(sectorlabel_run)
            
            if self.withoutIndividualSectorRuns:
                for sector in self.sectorList['Sector']:
                    sectorlabel_run='without'+sector 
                    self.runIFDMSector(sectorlabel_run)
            
            print 'All runs have started'
            print 'Good luck with IFDM!'
                
        else:
            print 'Exiting program.'
            exit()
    
    
################################       
# Major auxiliary routines
################################       

    
        
    def makeFileTreeSector(self,sectorlabel_run):
        # make the file tree for a specific sector    
        savefolder_runscripts=self.outputfolder+'/'+sectorlabel_run+'/Scripts/'
        self.makeFolder(savefolder_runscripts)
        launchscript_allruns_data=''

        for i in range(self.numberOfPeriods):
            # make folder
            nrString=str(i+1).zfill(2)
            meteolabel_run=nrString            
            meteostart_hours=(int(meteolabel_run)-1)*self.hoursPerPeriod+1
            
            # copy all files to runfolder (most processing in subroutine)
            savefolder_run=self.outputfolder+sectorlabel_run+'/'+nrString+'/'
            self.CopyFilesForRun(sectorlabel_run,meteolabel_run,savefolder_run)
            self.modifyFile(self.tempfileConfig,savefolder_run+'ifdm_add.conf','meteostart',str(meteostart_hours)) 
            shutil.copy(self.inputfileNewRaster,savefolder_run+'InputRaster.txt')
            shutil.copy(self.inputfileRaster,savefolder_run+'CompleteGrid.txt')

            # runscript
            savefile_launchscript=savefolder_runscripts+'Launch_'+sectorlabel_run+'_'+nrString            
            self.modifyFile(self.inputfileLaunchscript,savefile_launchscript,'folder',savefolder_run)
            runname=savefolder_run.replace('/','_')+'add'
            self.modifyFile(savefile_launchscript,savefile_launchscript,'runname',runname)
            launchscript_allruns_data=launchscript_allruns_data+'qsub '+savefile_launchscript+'\n'

            
       # save runscripts
        launchscript_allruns= open(self.outputfolder+'Launch_'+sectorlabel_run+'_add','w')
        launchscript_allruns.write(launchscript_allruns_data)
        launchscript_allruns.close()        
    
    def CopyFilesForRun(self,sectorLabel,meteoLabel,savefolder):  
        # make a single run folder
        self.makeFolder(savefolder)    
        
        shutil.copy(self.inputfileTimeFactors,savefolder+'Timefactors.txt')
        shutil.copy(self.inputfileLimitValues,savefolder+'InputLimitValues.txt')
        shutil.copy(self.inputfileBuilding,savefolder+'InputBuilding.txt')        
        shutil.copy(self.inputfileTimeSeries,savefolder+'Input_Timeseries.txt')
        shutil.copy(self.inputfileExecutable,savefolder+'ifdm')
     
        meteofile_run=self.tempoutputfolderMeteo+'Meteofile'+meteoLabel+'.txt' 
        shutil.copy(meteofile_run,savefolder+'InputMeteoType1.txt')
        
        linesources_run=self.tempoutputfolderEmissions+'InputLineSources_'+sectorLabel+'.txt'
        shutil.copy(linesources_run,savefolder+'InputLineSources.txt')
        surfacesourcesT1_run=self.tempoutputfolderEmissions+'InputSurfaceSourcesT1_'+sectorLabel+'.txt'
        shutil.copy(surfacesourcesT1_run,savefolder+'InputSurfaceSourcesT1.txt')
        surfacesourcesT2_run=self.tempoutputfolderEmissions+'InputSurfaceSourcesT2_'+sectorLabel+'.txt'
        shutil.copy(surfacesourcesT2_run,savefolder+'InputSurfaceSourcesT2.txt')
        pointsources_run=self.tempoutputfolderEmissions+'InputPointSources_'+sectorLabel+'.txt'
        shutil.copy(pointsources_run,savefolder+'InputPointSources.txt')
  

    # run IFDM for a specific sector        
    def runIFDMSector(self,sectorlabel):        
            
        runfile=self.outputfolder+'Launch_'+sectorlabel+'_add'
        os.system('chmod +x '+runfile)
        os.system(runfile)
            
   
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
            

