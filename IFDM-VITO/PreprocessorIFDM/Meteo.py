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

"""

import pandas as pd
import datetime
import shutil
import os

class Meteo:
    # constructor
    def __init__(self,name,nameOutput,lengthOfPeriod,meteoCopyAssimil,startDate=None,endDate=None):
        self.InputMeteoFile=name
        self.OutputDirectory=nameOutput
        self.hoursPerPeriod=lengthOfPeriod         
        self.MeteoData=pd.DataFrame()        
        self.timeInformation=pd.DataFrame()
        self.meteoCopyAssimil=meteoCopyAssimil
           
        print 'Initialised meteo class'
        
        
        ## read meteo data
        if not(self.meteoCopyAssimil):
            self.MeteoData=pd.read_excel(self.InputMeteoFile)
            self.MeteoData.index=pd.to_datetime(self.MeteoData[['Year','Month','Day','Hour']])
                    
            ## Check on input data
            if startDate is None:
                print 'No startdate provided, use end of Meteofile'
                self.startDate = self.MeteoData.index[0]
            else:
                self.startDate = startDate
                
            if endDate is None:
                print 'No enddate provided, use end of Meteofile'
                self.endDate = self.MeteoData.index[-1]   
            else:
                self.endDate = endDate            
                       
            ## select correct period (default: entire excel file)
            self.MeteoData=self.MeteoData[self.startDate:self.endDate]  
        else:
            if startDate is None:
                print 'No startdate provided, use end of Meteofile'
                self.MeteoData=pd.read_excel(self.InputMeteoFile)
                self.MeteoData.index=pd.to_datetime(self.MeteoData[['Year','Month','Day','Hour']])
                startDate = self.MeteoData.index[0]
                
            if endDate is None:
                print 'No enddate provided, use end of Meteofile'
                self.MeteoData=pd.read_excel(self.InputMeteoFile)
                self.MeteoData.index=pd.to_datetime(self.MeteoData[['Year','Month','Day','Hour']])
                endDate = self.MeteoData.index[-1]   
                
            timesteps=pd.date_range(startDate,endDate,freq='h')
            columns=['WEwind','NSwind','Height','Temperature','Year','Month','Day','Hour','Lat','Lon','TimeParameter']
            self.MeteoData=pd.DataFrame(index=timesteps,columns=columns)
            self.MeteoData.Lat=51.037861 
            self.MeteoData.Lon=4.240528
            self.MeteoData.Year=self.MeteoData.index.year
            self.MeteoData.Month=self.MeteoData.index.month
            self.MeteoData.Day=self.MeteoData.index.day
            self.MeteoData.Hour=self.MeteoData.index.hour
            self.MeteoData.Temperature=300
            self.MeteoData.fillna(value=1,inplace=True)
            print 'Coupling to ASSIMIL meteo. Will use default values for Meteo type 1.'
        
        ## calculate number of hours
        self.NumberOfHours=len(self.MeteoData['Year'])
        self.__CalculateNumberOfHours()
        self.periodStart=pd.Series(index=range(self.NumberOfPeriods),dtype='string')
    
    # save the meteo to the output file    
    def writeMeteo(self):
        # loop over all but last file
        for i in range(self.NumberOfPeriods-1):
            self.periodStart[i]=self.MeteoData.index[i*self.hoursPerPeriod].strftime('%Y%m%d%H')
            FileSplitted=self.MeteoData[i*self.hoursPerPeriod:(i+1)*self.hoursPerPeriod]
            Outputfile=self.OutputDirectory+'Meteofile'+self.periodStart[i]+'.txt'         
            self.__writeIndividualFile(Outputfile,FileSplitted)
        
        # last file
        self.periodStart[self.NumberOfPeriods-1]=self.MeteoData.index[(self.NumberOfPeriods-1)*self.hoursPerPeriod].strftime('%Y%m%d%H')        
        FileSplitted=self.MeteoData[(self.NumberOfPeriods-1)*self.hoursPerPeriod:self.NumberOfHours]
        Outputfile=self.OutputDirectory+'Meteofile'+self.periodStart[self.NumberOfPeriods-1]+'.txt'         
        self.__writeIndividualFile(Outputfile,FileSplitted)                    
        
            
        print 'Splitted meteo file in '+str(self.NumberOfPeriods)+' files'
        print 'Saved all files'
        
        
    def preprocessAssimilMeteo(self,folderMeteo):
        script_grid=folderMeteo+'readdataOlivier_grid.out'
        file_grid_in=folderMeteo+'ASSIMILgrid.txt'
        script_assimil=folderMeteo+'readdataOlivier.out'
        self.checkPresence(script_grid)
        self.checkPresence(script_assimil)
        self.checkPresence(file_grid_in)
        
        print 'Converting ASSIMMIL Meteo to correct format'
        os.system('cd '+folderMeteo+'; chmod +x readdataOlivier_grid.out; ./readdataOlivier_grid.out')
        
        
        for datename in self.MeteoData.index:
            file_in=folderMeteo+'Assimil_'+datename.strftime('%Y.%m.%d_%H')+'.txt'
            if os.path.exists(file_in.replace('txt','txt.gz')): 
                shutil.copy(file_in.replace('txt','txt.gz'),folderMeteo+'Assimil.txt.gz')            
                os.system('gunzip '+folderMeteo+'Assimil.txt.gz  -f')
            else:                
                self.checkPresence(file_in)
                shutil.copy(file_in,folderMeteo+'Assimil.txt')
            os.system('cd '+folderMeteo+'; chmod +x readdataOlivier.out; ./readdataOlivier.out')
            shutil.move(folderMeteo+'klaarvoorgebruik.txt',folderMeteo+'MeteoAssimil_'+datename.strftime('%Y%m%d%H')+'.txt')
            
            
            
    def preprocessRIObackground(self,folderRIO):
        script_conversion=folderRIO+'ConvertRIO.out'
        file_grid_in=folderRIO+'lambert4x4.dat'
        self.checkPresence(script_conversion)
        self.checkPresence(file_grid_in)
        list_days=pd.Series(self.MeteoData.index.date).unique()
        
        print 'Converting RIO background to correct format'
        for datename in list_days:
            y=str(datename.year)
            m=str(datename.month).zfill(2)
            d=str(datename.day).zfill(2)
            os.system('cd '+folderRIO+'; chmod +x ConvertRIO.out; ./ConvertRIO.out '+y+' '+m+' '+d)

            
        
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
        
        
        
    def checkPresence(self,file_in):        
        if not os.path.exists(file_in):
            print 'Computer says "no".'
            print 'Input file not found error'
            print 'File '+file_in+' does not exist'            
            print 'Exiting program...'
            sys.exit()

        
    
