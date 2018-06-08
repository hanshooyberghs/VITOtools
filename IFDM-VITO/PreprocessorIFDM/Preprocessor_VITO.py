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
        
        # variables for run that could be changed but have default values
        self.individualSectorRuns=True  # runs per sector (default is True)
        self.withoutIndividualSectorRuns=False  # runs without sector (default if False)
        self.allSectorsRun=True         # runs for all sectors combined (default is True)
        
        # variables for OSPM
        self.OSPM=False                 # OSPM? (default is False)
        self.generateWarning=True
        
        # variables for rugosity
        self.Rugosity=False
                
        # variables for coupling to background
        self.couplingBackground=False   # coupling to CTM (default is False)
        self.linkToCTM=False            # location of CTM files (has to be provided by hand in case of coupling with CTM)
        self.CTMfile=False
        self.CTMgriddef=False
        
        # variables for scenario run        
        self.scenario=False             # scenario run (default is False)
        self.scenarioFolder=False       # location of base case for scenario run(has to be provided by hand in case of scenario run)
        
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
        
        # location of default IFDM input files
        self.inputfileTimeFactors='default/Timefactors.txt'
        self.inputfileBuilding='default/InputBuilding.txt'
        self.inputfileTimeSeries='default/Input_Timeseries.txt'   
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
        self.inputfileLaunchscript_VSC='default/Launchscript_VSC.txt'
              
        
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
        self.inputfileLaunchscript_VSC=self.inputfolder+self.inputfileLaunchscript_VSC
        self.inputfileLimitValues=self.inputfolder+self.inputfileLimitValues
        self.inputfileStraatgegevens=self.inputfolder+self.inputfileStraatgegevens    
        self.inputfileRugosity=self.inputfolder+self.inputfileRugosity
        
        # Make temporary folders
        self.tempfolderbase=self.inputfolder+'temp/'
        self.tempoutputfolderMeteo=self.tempfolderbase+'/Meteo/'
        self.tempoutputfolderEmissions=self.tempfolderbase+'/Emissions/'
        self.tempfileConfig=self.tempfolderbase+'ifdm_add.conf'
        self.tempfileConfigCoupling=self.tempfolderbase+'ifdm_coupling.conf'
        self.tempfileConfigOSPM=self.tempfolderbase+'ifdm_ospm.conf'
        
        # make the output folders
        
        self.makeFolder(self.tempoutputfolderMeteo)
        self.makeFolder(self.tempoutputfolderEmissions)
        
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
        
           
        if not(self.CTMfile) and self.couplingBackground:
            print 'You have demanded linkage with a CTM, but did not provide the file of the CTM'
            print 'Unresolvable error'
            print 'Check your input (probably providing the vairable "linkToCTM" solves the issue)'
            print 'Exiting'
            exit()
        
        if not(self.CTMgriddef) and self.couplingBackground:
            print 'You have demanded linkage with a CTM, but did not provide the grid definition of the CTM'
            print 'Unresolvable error'
            print 'Check your input (probably providing the vairable "linkToCTM" solves the issue)'
            print 'Exiting'
            exit()
                    
        # check consistency of input values related to scenario run
        if self.scenario:
            self.checkPresence(self.scenarioFolder)        
        if self.scenario and self.couplingBackground:
            print 'Your input parameters do not make sense:'
            print 'Coupling to background and scenario'
            print 'You have to make two runs: the one with the coupling, and then one with the scenario'
            print 'Exiting'
            exit()
        
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
        if self.couplingBackground:
            self.modifyFile(self.inputfileConfigCoupling,self.tempfileConfigCoupling,'pollutantslist',self.pollutantList)
            self.modifyFile(self.tempfileConfigCoupling,self.tempfileConfigCoupling,'folder_of_ctm',"./")
            self.modifyFile(self.tempfileConfigCoupling,self.tempfileConfigCoupling,'name_of_ctm','background.')

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
                print 'Did not find straatgegevens_config_exc.txt'
                print 'Will now try straatgegevens_config.txt'
                self.checkPresence(self.inputfileStraatgegevens)          

    ##############################        
    # Call meteo class
    ##############################
    
    def makeMeteo(self): 
        # initialize: check if start and enddate are provided
        met=Meteo(self.inputfileMeteo,self.tempoutputfolderMeteo,self.lengthOfPeriod,self.meteoCopyAssimil,self.startDate,self.endDate)
              
        # write meteo input to tempfolder
        met.writeMeteo()
        
        
        if self.meteoPreprocess:
            met.preprocessAssimilMeteo(self.meteoLocation)
        
        # get some important values for the remainder of the script
        self.hoursPerPeriod=met.hoursPerPeriod
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
            
            
    ####################################
    # Submit the runs to the VSC cluster   
    ####################################

    def RunIFDM_VSC(self):
    
        var='yes'
        if var == 'yes':  
            if self.allSectorsRun:   
                sectorlabel='AllSectors'
                self.runIFDMSector_VSC(sectorlabel)
                
            if self.individualSectorRuns:
                for sector in self.sectorList['Sector']:
                    sectorlabel_run=sector 
                    self.runIFDMSector_VSC(sectorlabel_run)
            
            if self.withoutIndividualSectorRuns:
                for sector in self.sectorList['Sector']:
                    sectorlabel_run='without'+sector 
                    self.runIFDMSector_VSC(sectorlabel_run)
            
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
        savefolder_runscripts=self.runfolder+'/'+sectorlabel_run+'/Scripts/'
        self.makeFolder(savefolder_runscripts)
        launchscript_allruns_data=''
        if self.couplingBackground:        
            launchscript_allruns_data_coupling=''
            launchscript_allruns_data_all=''
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
                   
            # Changes for coupling to conf file
            meteostart_ctm=str(meteostart_hours*3600)
            config_add=savefolder_run+'ifdm_add.conf'
            self.modifyFile(self.tempfileConfig,config_add,'ctm_start',str(meteostart_ctm))
                        
            # copy files for coupling to assimil meteo.
            if self.meteoCopyAssimil:
                self.CopyMeteoFiles(savefolder_run,nrString,numberMeteo)
            else:
                self.modifyFile(config_add,config_add,'meteostart',str(meteostart_hours))

            # runscript
            savefile_launchscript=savefolder_runscripts+'Launch_'+sectorlabel_run+'_'+nrString            
            self.modifyFile(self.inputfileLaunchscript,savefile_launchscript,'folder',savefolder_run)
            if self.outputfolder is not None:
                self.modifyFile(savefile_launchscript,savefile_launchscript,'outputlocation',self.outputfolder)
                self.modifyFile(savefile_launchscript,savefile_launchscript,'date',nrString)

            runname=savefolder_run.replace('/','_')+'add'  
            self.modifyFile(savefile_launchscript,savefile_launchscript,'runname',runname)
            launchscript_allruns_data=launchscript_allruns_data+'qsub '+savefile_launchscript+'\n'
            
            
            #in case of coupling to background
            if self.couplingBackground:
                # coupling in config file      
                config_coupling=savefolder_run+'ifdm_coupling.conf'
                self.modifyFile(self.tempfileConfigCoupling,config_coupling,'ctm_start',str(meteostart_ctm))  
                self.modifyFile(config_coupling,config_coupling,'meteostart',str(meteostart_hours))
                
                #Symbolic links
                self.makeSymbolicLinksCTM(savefolder_run)
                
                # runscript
                savefile_launchscript_coupling=savefolder_runscripts+'Launch_'+sectorlabel_run+'_'+nrString+'_coupling'
                self.modifyFile(savefile_launchscript,savefile_launchscript_coupling,'add','coupling')                
                self.modifyFile(savefile_launchscript_coupling,savefile_launchscript_coupling,'true','false')
                launchscript_allruns_data_coupling=launchscript_allruns_data_coupling+'qsub '+savefile_launchscript_coupling+'\n'
                launchscript_allruns_data_all=launchscript_allruns_data_all+'echo '+nrString+'\n'
                launchscript_allruns_data_all=launchscript_allruns_data_all+'qsub '+savefile_launchscript_coupling+'\n'
                launchscript_allruns_data_all=launchscript_allruns_data_all+' sleep 2s \n'
                launchscript_allruns_data_all=launchscript_allruns_data_all+'qsub '+savefile_launchscript+'\n'
            
            # in case of scenario run                    
            if self.scenario:
                # make symbolic links to scenario folder
                folder_scenario=self.scenarioFolder+'/'+nrString+'/'
                self.checkPresence(folder_scenario)
                self.makeSymbolicLinksIFDM(savefolder_run,folder_scenario) 
                
            # in case of OSPM run    
            if self.OSPM: 
                # make folder                
                savefolder_runospm=self.runfolder+sectorlabel_run+'/'+nrString+'ospm/'                 
                self.CopyFilesForRun(sectorlabel_run,nrString,savefolder_runospm)
                self.makeSymlink(self.inputfileStraatgegevens,savefolder_runospm+'straatgegevens_config_exc.txt')
                shutil.copy(self.inputfileRasterOSPM,savefolder_runospm+'InputRaster.txt')
                
                # changes to config file
                config_OSPM=savefolder_runospm+'ifdm_ospm.conf'
                self.modifyFile(self.tempfileConfigOSPM,config_OSPM,'ctm_start',str(meteostart_ctm))
                if self.meteoCopyAssimil:
                    self.CopyMeteoFiles(savefolder_runospm,nrString,numberMeteo)
                else:
                    self.modifyFile(config_OSPM,config_OSPM,'meteostart',str(meteostart_hours))          
                
                # symbolic links between IFDM and OSPM folder
                self.makeSymbolicLinksIFDM(savefolder_runospm,savefolder_run)                
                
                # runscript
                savefile_launchscript_ospm=savefolder_runscripts+'Launch_'+sectorlabel_run+'_'+nrString+'_ospm'  
                self.modifyFile(self.inputfileLaunchscript,savefile_launchscript_ospm,'folder',savefolder_runospm)
                runname=savefolder_runospm.replace('/','_')
                self.modifyFile(savefile_launchscript_ospm,savefile_launchscript_ospm,'runname',runname)
                self.modifyFile(savefile_launchscript_ospm,savefile_launchscript_ospm,'add','ospm')
                if self.outputfolder is not None:
                    self.modifyFile(savefile_launchscript,savefile_launchscript,'outputlocation',self.outputfolder)
                    self.modifyFile(savefile_launchscript,savefile_launchscript,'date',nrString)
                    self.modifyFile(savefile_launchscript,savefile_launchscript,'saga','ospm')

                launchscript_allruns_data_ospm=launchscript_allruns_data_ospm+'qsub '+savefile_launchscript_ospm+'\n'
       
       # save runscripts
        launchscript_allruns= open(self.runfolder+'Launch_'+sectorlabel_run+'_add','w')
        launchscript_allruns.write(launchscript_allruns_data)
        launchscript_allruns.close()
        
        if self.couplingBackground:
            launchscript_allruns= open(self.runfolder+'Launch_'+sectorlabel_run+'_coupling','w')
            launchscript_allruns.write(launchscript_allruns_data_coupling)
            launchscript_allruns.close()
            
            launchscript_allruns= open(self.runfolder+'Launch_'+sectorlabel_run+'_couplingadd','w')
            launchscript_allruns.write(launchscript_allruns_data_all)
            launchscript_allruns.close()

        
        if self.OSPM:    
            launchscript_allruns= open(self.runfolder+'Launch_'+sectorlabel_run+'_ospm','w')
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
        meteofile_run=self.tempoutputfolderMeteo+'Meteofile'+meteoLabel+'.txt' 
        shutil.copy(meteofile_run,savefolder+'InputMeteoType1.txt')
        
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
        if self.couplingBackground:
            runfile=self.runfolder+'Launch_'+sectorlabel+'_couplingadd'
            os.system('chmod +x '+runfile)
            os.system(runfile)         
        else:
            runfile=self.runfolder+'Launch_'+sectorlabel+'_add'
            os.system('chmod +x '+runfile)
            os.system(runfile)
        
        if self.OSPM:
            runfile=self.runfolder+'Launch_'+sectorlabel+'_ospm'
            os.system('chmod +x '+runfile)
            os.system(runfile)
            
            
    # run IFDM for a specific sector on the VSC-cluster        
    def runIFDMSector_VSC(self,sectorlabel):
        folder_sector=self.runfolder+'/'+sectorlabel+'/'
        launchscript_VSC=folder_sector+'Launchscript_VSC.txt'
        self.modifyFile(self.inputfileLaunchscript_VSC,launchscript_VSC,'folder_base',folder_sector)
        self.periodStart.name='date'
        self.periodStart.to_csv(folder_sector+'dates.txt',index=False,header=True)
        launch_command='wsub  -batch '+launchscript_VSC+' -data '+folder_sector+'dates.txt'
        print 'Submit runs command: '+launch_command
        os.system(launch_command)

            
    # make symbolic links to CTM-output
    def makeSymbolicLinksCTM(self,savefolder):
        uitvoerfile=savefolder+"Aurora_conc"
        self.makeSymlink(self.CTMfile,uitvoerfile)
        
        folder_griddef=savefolder+'concentraties/'
        self.makeFolder(folder_griddef)
        uitvoerfile=folder_griddef+"background.grid_def"
        self.makeSymlink(self.CTMgriddef,uitvoerfile)
        
        
    # copy files for assimil meteo
    def CopyMeteoFiles(self,savefolder_run,nrString,numberMeteo):
        file_in=open(savefolder_run+'lijst','w')        
        file_in.write('27 \n')
        self.makeSymlink(self.meteoLocation+'GridOlivier.txt',savefolder_run+'GridOlivier.txt')
        startdate=datetime.datetime.strptime(nrString,'%Y%m%d%H')
        enddate=startdate+datetime.timedelta(hours=numberMeteo-1)
        for day in pd.date_range(startdate,enddate,freq='H'):
            filename='MeteoAssimil_'+day.strftime('%Y%m%d%H')+'.txt'
            shutil.copy(self.meteoLocation+filename,savefolder_run+filename)
            file_in.write("'"+filename+"' \n")
        file_in.close()


    # make symbolic links to parent folder for OSPM or scenario runs
    def makeSymbolicLinksIFDM(self,savefolder,savefolder_coupling):       
        # symbolic links for normal pollutants
        for pol in self.pollutantsIFDM:
            invoerfile=savefolder_coupling+"uitvoer_AuroPod_"+pol+".txt"
            uitvoerfile=savefolder+"uitvoer_AuroPod_"+pol+".txt"
            self.makeSymlink(invoerfile,uitvoerfile)
        
        # symbolic links for chemistry (ozone, NO2 and coupling constants)
        invoerfile=savefolder_coupling+"uitvoer_AuroPod_NO2.txt"
        uitvoerfile=savefolder+"uitvoer_AuroPod_NO2.txt"
        self.makeSymlink(invoerfile,uitvoerfile)
        
        invoerfile=savefolder_coupling+"uitvoer_AuroPod_O3.txt"
        uitvoerfile=savefolder+"uitvoer_AuroPod_O3.txt"
        self.makeSymlink(invoerfile,uitvoerfile)
        
        invoerfile=savefolder_coupling+"uitvoer_AuroPod_rjk.txt"
        uitvoerfile=savefolder+"uitvoer_AuroPod_rjk.txt"
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
                
            

