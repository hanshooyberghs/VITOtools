# -*- coding: utf-8 -*-
"""
Meteopreprocessing routine for IFDM
===================================

@author: hooyberh (hans.hooyberghs@vito.be)

Basic description: Class to preprocess the meteo input for IFDM. Using the functions
in this class, the Excel file containing the IFDM input is processed and text files
for the runs are composed and stored in an output folder. The routine also includes some checks 
to infer whether the correct input is provided in the Excel file (check for column 
headers).

Input parameters for constructor:
* name: name of the input Excel file. Note: the file must contain the correct work sheets.
* nameOutput: name of the output directory where the text files will be stored
* nrPeriods: number of meteorological periods for the simulation

Routines:
* Constructor: Construct the object, read the Excel file and perform some checks 
    concerning the input format
* writeMeteo: make meteo files    
* Automatically invoked routines:
    - readMeteo: reads the meteo file (automatically invoked by the constructor)
    - getSectorList: provide list with sectors (invoked by readEmissions)
* Auxiliary routines:
    - see below
    
Dependencies:
* Python routines:
    -pandas
    -numpy
        
Output: 
* Meteo input for IFDM, stored in the provided OutputDirectory

Other input parameters: see constructor. Only change these if you know what you
are doing (don't tell me that I didn't warn you ;) )

===

To do:
* Add missing values in time series (f.i. by interpolation)
"""

import pandas as pd

class Meteo:
    # constructor
    def __init__(self,name,nameOutput,lengthOfPeriod):
        self.InputMeteoFile=name
        self.OutputDirectory=nameOutput
        self.hoursPerPeriod=lengthOfPeriod         
        self.MeteoData=pd.DataFrame()        
        self.expectedColumns=['WE wind speed (m/s)','NS wind speed (m/s)','Height (m)','Temperature (K)',\
            'Year','Month','Day','Hour','Lat','Lon','Time parameter']
           
        print 'Initialised meteo class'
        
        self.readMeteo()
        self.__CalculateNumberOfHours()
    
    # save the meteo to the output file    
    def writeMeteo(self):
        # loop over all but last file
        for i in range(self.NumberOfPeriods-1):
            FileSplitted=self.MeteoData[i*self.hoursPerPeriod:(i+1)*self.hoursPerPeriod]
            Outputfile=self.OutputDirectory+'Meteofile'+str(i+1).zfill(2)+'.txt'         
            self.__writeIndividualFile(Outputfile,FileSplitted)
        
        # last file
        FileSplitted=self.MeteoData[(self.NumberOfPeriods-1)*self.hoursPerPeriod:self.NumberOfHours]
        Outputfile=self.OutputDirectory+'Meteofile'+str(self.NumberOfPeriods).zfill(2)+'.txt'         
        self.__writeIndividualFile(Outputfile,FileSplitted)                    
            
            
        print 'Splitted meteo file in '+str(self.NumberOfPeriods)+' files'
        print 'Saved all files'
        
        
################################       
# Automatically invoked routines
################################   
        
    # read the meteo file    
    def readMeteo(self):
        self.MeteoData=pd.read_excel(self.InputMeteoFile)
        #perform check on the naming of the columns
        if (self.MeteoData.columns==self.expectedColumns).all():
            print 'Input Meteo is OK'
        else:
            print 'Columns of meteo file are not in correct order'
            print "I will continue, don't blame me if things go wrong"
        
        self.NumberOfHours=len(self.MeteoData['Year'])
        
        
################################       
# Auxiliary routines
################################       

    # write a file      
    def __writeIndividualFile(self,OutputFile,DataToWrite):
        number_rows=len(DataToWrite['Year'])
        f=open(OutputFile,'w')
        f.write('%d' % number_rows)
        f.write('\n')
        f.close()         
        DataToWrite.to_csv(OutputFile,sep='\t',header=False,index=False,mode='a')
    
    # calculate the number of hours in the periods    
    def __CalculateNumberOfHours(self):
        self.NumberOfPeriods=self.NumberOfHours/self.hoursPerPeriod +1 #yields always integer
        self.hoursLastPeriod=self.NumberOfHours-self.hoursPerPeriod*(self.NumberOfPeriods-1)
        if self.hoursLastPeriod == 0:
           self.hoursLastPeriod = self.hoursPerPeriod
           self.NumberOfPeriods=self.NumberOfPeriods-1
        
        print 'Total number of hours: '+str(self.NumberOfHours)
        print 'Hours per period: '+str(self.hoursPerPeriod)
        print 'Hours last period: '+str(self.hoursLastPeriod)
        print 'Number of periods:'+str(self.NumberOfPeriods)
        
    
