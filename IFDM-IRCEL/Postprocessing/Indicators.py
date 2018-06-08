"""

List of indicators

"""
import pandas as pd
import datetime
import numpy as np

# YearlyMean
def YearlyMean(input_hourly):
    indicator=input_hourly.mean(axis=1)
    return indicator    
        
# Exceedance of 200 ug/m3 of hourly concentrations        
def ExcHour_200(input_hourly):
    indicator=(input_hourly>200).sum(axis=1)
    return indicator

# 99.78 percentile of hourly concentration    
def P9978(input_hourly):
    indicator=input_hourly.quantile(0.9978,axis=1)
    return indicator

# Mean of daily maximal values     
def MeanDayMax1Hour(input_hourly):
    daymax=input_hourly.resample('d',axis=1).max()
    indicator=daymax.mean(axis=1)
    return indicator    

# Exceedance of 50.5 ug/m3 of daily mean concentrations            
def ExcDay_505(input_hourly):
    input_daily=input_hourly.resample('d',axis=1).mean()
    indicator=(input_daily>50.5).sum(axis=1)
    return indicator    

# 99.78 percentile of daily concentration        
def P9014_daily(input_hourly):
    input_daily=input_hourly.resample('d',axis=1).mean()
    indicator=input_daily.quantile(0.9014,axis=1)
    return indicator

# Exceedance of 25.5 ug/m3 of daily mean concentrations            
def ExcDay_255(input_hourly):
    input_daily=input_hourly.resample('d',axis=1).mean()
    indicator=(input_daily>25.5).sum(axis=1)
    return indicator    

# NET60 for ozone: number of days with daily max 8 hour above 120.5
def NET60(input_hourly):
    input_8hmax=input_hourly.rolling(8,axis=1).mean().resample('d',axis=1).max()
    indicator=(input_8hmax>120.5).sum(axis=1)
    return indicator    

# WHO for ozone: number of days with daily max 8 hour above 100.5
def WHO(input_hourly):
    input_8hmax=input_hourly.rolling(8,axis=1).mean().resample('d',axis=1).max()
    indicator=(input_8hmax>100.5).sum(axis=1)
    return indicator    

# AOT40veg for ozone    
def AOT40veg(input_hourly):
    data_months=input_hourly.iloc[:,np.array([item in range(5,8) for item in input_hourly.columns.month])]
    data_selected=data_months.iloc[:,np.array([item in range(8,20) for item in data_months.columns.hour])]
    indicator=np.maximum(data_selected-80,0).sum(axis=1)    
    return indicator    

# for for ozone      
def AOT40for(input_hourly):
    data_months=input_hourly.iloc[:,np.array([item in range(4,10) for item in input_hourly.columns.month])]
    data_selected=data_months.iloc[:,np.array([item in range(8,20) for item in data_months.columns.hour])]
    indicator=np.maximum(data_selected-80,0).sum(axis=1)
    return indicator    

# AOT60 for ozone    
def AOT60(input_hourly):
    input_8hmax=input_hourly.rolling(8,axis=1).mean().resample('d',axis=1).max()
    indicator=np.maximum(input_8hmax-120,0).sum(axis=1)*8.0
    return indicator    

# yearly mean of the daily 8h max exceedances    
def meanmax8h(input_hourly):
    input_8hmax=input_hourly.rolling(8,axis=1).mean().resample('d',axis=1).max()
    indicator=input_8hmax.mean(axis=1)
    return indicator
    
