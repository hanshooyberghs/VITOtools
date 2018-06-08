# -*- coding: utf-8 -*-
"""
Emissions preprocessing routine for IFDM
========================================

@author: hooyberh (hans.hooyberghs@vito.be)

Basic description: Class to preprocess the emissions for IFDM. Using the functions
in this class, the Excel file containing the IFDM input is processed and text files
for the runs with the different pollutants are composed and stored in an output folder.
The routine also includes some checks to infer whether the correct input is provided
in the Excel file (check for column headers, check for missing data).


Input parameters for constructor:
* name: name of the input Excel file. Note: the file must contain the correct work sheets.
* nameOutput: name of the output directory where the text files will be stored
* pollutants: list of pollutants (in order to check the input file)


Routines:
* Constructor: Construct the object, read the Excel file(using the read routine) 
    and perform some checks concerning the input format
* makeEmissionsAll: make emissions for runs with all sectors
* makeEmissionsRunsPerSector: make emissions for runs per sector
* makeEmissionsRunsPerSector: make emissions for runs leaving one sector out
* Automatically invoked routines:
    - readEmissions: reads the emission file (automatically invoked by the constructor)
    - getSectorList: provide list with sectors (invoked by readEmissions)
* Auxiliary routines
    - see below
    
Dependencies:
* Python routines:
    -pandas
    -numpy
* External software:
    nihil
        
Output: 
* Emission input for IFDM, stored in the provided OutputDirectory

Other input parameters: see constructor. Only change these if you know what you
are doing (don't tell me that I didn't warn you ;) )

===
To do list
----------
* remove rows with missing coordinates and emissions
* Check if excel file sheets exists. If not: make empty data array
"""


import pandas as pd

class Emissions:
    # constructor
    def __init__(self,name,nameOutput,pollutants,shipping,ospm):
        self.InputEmissionFile=name
        self.OutputDirectory=nameOutput
        self.LineSources=pd.DataFrame()
        self.PointSources=pd.DataFrame()
        self.SurfaceSources=pd.DataFrame()
        self.AltSurfaceSources=pd.DataFrame()
        self.SectorList=pd.DataFrame()
        self.Pollutants=pollutants
        self.Shipping=shipping
        self.OSPM=ospm
        
        print 'Initialised emissions class'
        
        # read emissions   
        self.readEmissions()
        
    # make emissions for runs including all sectors       
    def makeEmissionsAll(self):    
        self.__OutputfilePoint=self.OutputDirectory+'InputPointSources_AllSectors.txt'
        self.__OutputfileLine=self.OutputDirectory+'InputLineSources_AllSectors.txt'
        self.__OutputfileSurface=self.OutputDirectory+'InputSurfaceSourcesT2_AllSectors.txt'
        self.__OutputfileAltSurface=self.OutputDirectory+'InputSurfaceSourcesT1_AllSectors.txt'
        self.__writeData(self.LineSources,self.__OutputfileLine)
        self.__writeData(self.PointSources,self.__OutputfilePoint)
        self.__writeData(self.SurfaceSources,self.__OutputfileSurface)
        self.__writeData(self.AltSurfaceSources,self.__OutputfileAltSurface)
        
        print 'Outputfiles for all emissions have been made'
    
 
 
################################       
# Automatically invoked routines
################################    
    # read emissions       
    def readEmissions(self):        
        # check the presence of all the sheets. if some are missing, make empty arrays
        xls = pd.ExcelFile(self.InputEmissionFile)
        sheets = xls.sheet_names
        
        if 'LineSources' in sheets:
            self.readLineSources()
        else:
            print 'No Line Sources sheet in Emissions files'
            print 'I will use empty line sources'
            self.LineSources=pd.DataFrame()
            
        if 'PointSources' in sheets:
            self.readPointSources()
        else:
            print 'No Point Sources sheet in Emissions files'
            print 'I will use empty point sources'
            self.PointSources=pd.DataFrame()
            
        if 'SurfaceSources' in sheets:
            self.readSurfaceSources()
        else:
            print 'No Surface Sources sheet in Emissions files'
            print 'I will use empty surface sources'
            self.SurfaceSources=pd.DataFrame()
            
        if 'AltSurfaceSources' in sheets:
            self.readAltSurfaceSources()
        else:
            print 'No Alternative SurfaceSources in Emissions files'
            print 'I will use empty alternative surface sources'
            self.AltSurfaceSources=pd.DataFrame()
                        
        print 'Emissions have been read'
        
    
    # read line sources
    def readLineSources(self):
        self.LineSources=pd.read_excel(self.InputEmissionFile,sheetname='LineSources')
        if self.LineSources.empty==False:
            nrColumns=12+len(self.Pollutants)
            if self.Shipping:
                nrColumns=nrColumns+3
            if self.OSPM:
                nrColumns=nrColumns+4
            self.LineSources=self.LineSources[list(self.LineSources.columns[0:nrColumns])] 
            self.__checkMissingValues(self.LineSources)
    
    #read point sources
    def readPointSources(self):
        self.PointSources=pd.read_excel(self.InputEmissionFile,sheetname='PointSources')
        if self.PointSources.empty==False:
            nrColumns=11+len(self.Pollutants)
            if self.Shipping:
                nrColumns=nrColumns+3
            if self.OSPM:
                nrColumns=nrColumns+4
            self.PointSources=self.PointSources[list(self.PointSources.columns[0:nrColumns])]
            self.__checkMissingValues(self.PointSources)
    
    # read surface sources
    def readSurfaceSources(self):
        self.SurfaceSources=pd.read_excel(self.InputEmissionFile,sheetname='SurfaceSources')
        if self.SurfaceSources.empty==False:
            self.__checkMissingValues(self.SurfaceSources)
    
    # read alternative surface sources
    def readAltSurfaceSources(self):
        self.ALtSurfaceSources=pd.read_excel(self.InputEmissionFile,sheetname='AltSurfaceSources')
        if self.AltSurfaceSources.empty==False:
            self.__checkMissingValues(self.ALtSurfaceSources)
     


################################       
# Auxiliary routines
################################       
        
    # write data to a file
    def __writeData(self,DataToWrite,OutputFile):
        if DataToWrite.empty:
            f=open(OutputFile,'w')
            f.write('%d' % 0)
            f.close()
        else:        
            number_rows=len(DataToWrite.index)
            f=open(OutputFile,'w')
            f.write('%d' % number_rows)
            f.write('\n')
            f.close()  
            DataToWrite.to_csv(OutputFile,sep='\t',header=False,index=False,mode='a')
    
        # check for missing values
    def __checkMissingValues(self,EmInput):
        if EmInput.empty:
            print 'Input OK'
        else:
            self.__cntMissingValues=0
            self.__isnulllist=EmInput.isnull()
            if self.__isnulllist.isnull().sum().sum() > 0:
                print 'Found missing values in the emission input'
                print 'Statistics:'
                print self.__isnulllist.isnull().sum()
                print 'Please fix your emission input'
                print 'Exiting'
                exit()

                
