# -*- coding: utf-8 -*-
"""
Run the preprocessing of IFDM. 
==============================

In this file, the input parameters for the IFDM PreProcessor should be presented.
More details on the input parameters can be found in the Preprocessor.py file.

@author: hooyberh (hans.hooyberghs@vito.be)
"""

# first initialize the run
import datetime
from Preprocessor_VITO import Preprocessor
Run = Preprocessor()

############################
#Compulsory input paramters. 
############################

Run.inputfolder='/projects/IFDM/IRCEL/InstallationPythonChain/Testfolder/'
Run.runfolder='/projects/IFDM/IRCEL/InstallationPythonChain/Testfolder/Runs/'

Run.pollutants=["NO2","PM10","PM25"]
Run.lengthOfPeriod=120

Run.couplingBackground=True
Run.CTMfile='/projects/N78C0/IFDM/RIO_IFDM_OSPM_2015_Antwerpen/RIO_achtergronden/Aurora_conc_All_RIO2015'
Run.CTMgriddef='/projects/N78C0/IFDM/RIO_IFDM_OSPM_2015_Antwerpen/RIO_achtergronden/rio2015.grid_def'

Run.individualSectorRuns=False
Run.allSectors=True
Run.OSPM=True

#Run.scenario=True
#Run.scenarioFolder='/projects/ifdm-jrc/test/'


#Run.startDate=datetime.datetime(2015,01,01,0)  #optional
#Run.endDate=datetime.datetime(2015,01,05,23)   #optional 

######################
#Preprocessing
######################


Run.setUp()
Run.makeMeteo()
Run.makeEmissions()
Run.makeFileTree()
Run.RunIFDM()