# -*- coding: utf-8 -*-
"""
Run the preprocessing of IFDM. 
==============================

In this file, the input parameters for the IFDM PreProcessor should be presented.
More details on the input parameters can be found in the Preprocessor.py file.

@author: hooyberh (hans.hooyberghs@vito.be)
"""

# first initialize the run
from Preprocessor_LiboVITO import Preprocessor
Run = Preprocessor()

############################
#Compulsory input paramters. 
############################

Run.inputfolder='/projects/MarcoPolo/ifdm_demo_Shijiazhuang/IFDM_input/Data/'
Run.outputfolder='/projects/MarcoPolo/ifdm_demo_Shijiazhuang/IFDM_runs/'

Run.pollutants=["NOX","BC","PM25","PM10"]
Run.lengthOfPeriod=24


######################
#Preprocessing
######################

Run.setUp()
Run.makeMeteo()
Run.makeEmissions()
Run.makeFileTree()
Run.RunIFDM()
