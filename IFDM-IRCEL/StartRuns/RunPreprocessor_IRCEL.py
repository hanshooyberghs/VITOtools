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
from Preprocessor_IRCEL import Preprocessor
Run = Preprocessor()

############################
#Compulsory input paramters. 
############################

Run.inputfolder='/projects/IFDM/IRCEL/InstallationPythonChain/Testfolder/'
Run.runfolder='/projects/IFDM/IRCEL/InstallationPythonChain/Testfolder/Runs/'
Run.outputfolder='/projects/IFDM/IRCEL/InstallationPythonChain/Testfolder/'


Run.pollutants=["NO2","PM10","PM25"]
Run.lengthOfPeriod=1

# coupling to RIO
Run.processBackground=True
Run.CTMfolder='/projects/IFDM/IRCEL/InstallationPythonChain/PythonChain/RIO/'
    # Needed in this folder: ConvertRIO.out, RIO-data,lambert4x4.dat
# fixed backbround
#Run.FixedBackground='Cams.xlsx'
#Run.processBackground=True

Run.meteoPreprocess=True
Run.meteoCopyAssimil=True
Run.meteoLocation='/projects/IFDM/IRCEL/InstallationPythonChain/PythonChain/MeteoAssimil/'
    # Needed in this folder: Assimildata, ASSIMILgrid.txt,readdataOlivier.out,readdataOlivier_grid.out

Run.startDate=datetime.datetime(2015,04,14,9)  #optional
Run.endDate=datetime.datetime(2015,04,14,5)   #optional 

# Shipping yer or no?
Run.Shipping=True

Run.OSPM=True


######################
#Preprocessing
######################


Run.setUp()
Run.makeMeteo()
Run.makeEmissions()
Run.makeFileTree()
Run.RunIFDM()