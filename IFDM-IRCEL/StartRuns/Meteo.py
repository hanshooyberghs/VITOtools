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
import sys
import glob

class Meteo:
    # constructor
    def __init__(self,name,nameOutput,nameOutputCTM,lengthOfPeriod,meteoCopyAssimil,fixedbackground,startDate=None,endDate=None):
        self.InputMeteoFile=name
        self.OutputDirectory=nameOutput
        self.OutputCTM=nameOutputCTM	#Charlotte
        self.hoursPerPeriod=lengthOfPeriod         
        self.MeteoData=pd.DataFrame()        
        self.timeInformation=pd.DataFrame()
        self.meteoCopyAssimil=meteoCopyAssimil
        self.inputfileBackground=fixedbackground
           
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
            self.MeteoData.Height=10
            self.MeteoData.fillna(value=1,inplace=True)
            print 'Coupling to ASSIMIL meteo. Will use default values for Meteo type 1.'
        
        ## calculate number of hours
        self.NumberOfHours=len(self.MeteoData['Year'])
        self.__CalculateNumberOfHours()
        self.periodStart=pd.Series(index=range(self.NumberOfPeriods),dtype='string')
        
        ## read background
        if self.inputfileBackground is not None:
            self.BackgroundData=pd.read_excel(self.inputfileBackground)
    
    # save the meteo to the output file    
    def writeMeteo(self):
        # loop over all but last file
        for i in range(self.NumberOfPeriods-1):
            self.periodStart[i]=self.MeteoData.index[i*self.hoursPerPeriod].strftime('%Y%m%d%H')
            FileSplitted=self.MeteoData[i*self.hoursPerPeriod:(i+1)*self.hoursPerPeriod]
            Outputfile=self.OutputDirectory+'Meteofile'+self.periodStart[i]+'.txt'         
            self.__writeIndividualFile(Outputfile,FileSplitted)
            
            if self.inputfileBackground is not None:
                FileSplitted=self.BackgroundData[i*self.hoursPerPeriod:(i+1)*self.hoursPerPeriod]
                Outputfile=self.OutputDirectory+'Background'+self.periodStart[i]+'.txt'        
                self.__writeIndividualFileBackground(Outputfile,FileSplitted)

        
        # last file
        self.periodStart[self.NumberOfPeriods-1]=self.MeteoData.index[(self.NumberOfPeriods-1)*self.hoursPerPeriod].strftime('%Y%m%d%H')        
        FileSplitted=self.MeteoData[(self.NumberOfPeriods-1)*self.hoursPerPeriod:self.NumberOfHours]
        Outputfile=self.OutputDirectory+'Meteofile'+self.periodStart[self.NumberOfPeriods-1]+'.txt'         
        self.__writeIndividualFile(Outputfile,FileSplitted)                    
        
        if self.inputfileBackground is not None:
            FileSplitted=self.BackgroundData[(self.NumberOfPeriods-1)*self.hoursPerPeriod:self.NumberOfHours]
            Outputfile=self.OutputDirectory+'Background'+self.periodStart[self.NumberOfPeriods-1]+'.txt'         
            self.__writeIndividualFileBackground(Outputfile,FileSplitted)          

            
        print 'Splitted meteo file in '+str(self.NumberOfPeriods)+' files'
        print 'Saved all files'
        
        
    def preprocessAssimilMeteo(self,folderMeteo):
        script_grid=folderMeteo+'readdataOlivier_grid.out'
        file_grid_in=folderMeteo+'ASSIMILgrid.txt'
        script_assimil=folderMeteo+'readdataOlivier.out'
        self.checkPresence(script_grid)
        self.checkPresence(script_assimil)
        
        # unzip
        
        print 'Converting ASSIMMIL Meteo to correct format'

        # First copy files from folderMeteo to self.OutputDirectory, then execute scripts
        shutil.copy(folderMeteo+'ASSIMILgrid.txt',self.OutputDirectory+'ASSIMILgrid.txt')
        shutil.copy(folderMeteo+'readdataOlivier_grid.out',self.OutputDirectory+'readdataOlivier_grid.out')
        os.system('cd '+self.OutputDirectory+'; ./readdataOlivier_grid.out')
        
        
        for datename in self.MeteoData.index:
            file_in=folderMeteo+'Assimil_'+datename.strftime('%Y.%m.%d_%H')+'.txt'
            if os.path.exists(file_in.replace('txt','txt.gz')): 
                shutil.copy(file_in.replace('txt','txt.gz'),self.OutputDirectory+'Assimil.txt.gz')            
                os.system('gunzip '+self.OutputDirectory+'Assimil.txt.gz  -f')
            else:               
                self.checkPresence(file_in)
                shutil.copy(file_in,self.OutputDirectory+'Assimil.txt')

            shutil.copy(folderMeteo+'readdataOlivier.out',self.OutputDirectory+'readdataOlivier.out')
            os.system('cd '+self.OutputDirectory+'; ./readdataOlivier.out')
            shutil.move(self.OutputDirectory+'klaarvoorgebruik.txt',self.OutputDirectory+'MeteoAssimil_'+datename.strftime('%Y%m%d%H')+'.txt')
            
            
            
    def preprocessRIObackground(self,folderRIO):
        script_conversion=folderRIO+'ConvertRIO.out'
        file_grid_in=folderRIO+'lambert4x4.dat'
        self.checkPresence(script_conversion)
        self.checkPresence(file_grid_in)
        list_days=pd.Series(self.MeteoData.index.date).unique()
        print list_days
        if len(list_days) > 1:
            list_days=pd.Series(self.MeteoData.index.date).unique()
        print list_days
        
        shutil.copy(folderRIO+'ConvertRIO.out',self.OutputCTM+'ConvertRIO.out')
        shutil.copy(folderRIO+'lambert4x4.dat',self.OutputCTM+'lambert4x4.dat')
        if not os.path.exists(self.OutputCTM+'concentraties/'):
        	shutil.copytree(folderRIO+'concentraties/',self.OutputCTM+'concentraties/')

        print 'Converting RIO background to correct format'
        for datename in list_days:
            y=str(datename.year)
            m=str(datename.month).zfill(2)
            d=str(datename.day).zfill(2)
            RIOfiletmp='*'+'_1h_'+y+m+d+'.dat'
            RIOfile=glob.glob(folderRIO+RIOfiletmp)
            for filename in RIOfile:
            	shutil.copy(filename,self.OutputCTM)
            os.system('cd '+self.OutputCTM+'; ./ConvertRIO.out '+y+' '+m+' '+d)
            os.system('cd '+self.OutputCTM+'; rm *_1h_'+y+m+d+'.dat')
            
        
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
        
        
    def __writeIndividualFileBackground(self,OutputFile,DataToWrite):
        with open(OutputFile, 'w') as f:
            DataToWrite[['O3','NOx','NO2']].to_csv(f,index=False,header=False,sep='\t')
            DataToWrite[['PM10']].to_csv(f,index=False,header=False)
            DataToWrite[['PM25']].to_csv(f,index=False,header=False)
            DataToWrite[['BC']].to_csv(f,index=False,header=False)

    
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

        
    
