"""
Conversion of FastRace output to IFDM input
===============================================================

Operations
* see manual

Execution: Usage: python Fastrace2IFDM.py inputshape outputbase [cut_off use_mult_factor tunnels_file conversionFactors_file]

Actions performed:
* Read FastRace output
* Remove roads with emissions below a threshold (default: 0)
* Correct for roads that appear twice
* Remove roads of length zero
* Conversion from double to single direction
* Save output as shape and Excel
* Provide some statistics

TO ADD
------
* Removal based on emissions in other file (f.i. threshold based on 2012 emissions)
    -> for instance using pre-compiled use / use not shape file (merge)
* LEZ function

"""



"""
=======
Input
=======
"""
import sys
input=sys.argv


if len(input)<3:
   print "Not enough input arguments"         
   print "Usage: python Fastrace2IFDM.py inputshape outputbase [cut_off use_mult_factor tunnels_file conversionFactors_file]"
   print 'Exiting program...'
   sys.exit()
   
input_shape=input[1]
outputbase=input[2]
if len(input)>3:
    cut_off_nox=float(input[4])
    print 'Cut-off NOx: '+str(cut_off_nox)
else:
    print 'No cut-off  supplied. Will use default: 0'
    cut_off_nox=0

if len(input)>4:
    mult_factor=int(float(input[4]))
    if mult_factor == 1:
        print 'Will ask for multiplication factors.'
    else:
        print 'No multiplication factors'
        
else:
    print 'No multiplication factor. Will use default: 1.0'
    mult_factor=0    
    
if len(input)>5:
    tunnels_bridges=input[5]
else:
    print 'NoTunnelsBridges File supplied. Will use default: TunnelsBridges.xlsx'
    tunnels_bridges='TunnelsBridges.xlsx'

if len(input)>6:
    file_factors=input[6]
else:
    print 'No ConversionFactors File supplied. Will use default: ConversionFactors.xlsx'
    file_factors='ConversionFactors.xlsx'

output_excel=outputbase+'.xlsx'
output_shape=outputbase+'.shp'
pollutants=['NOx','PM10','PM25','BC']


import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString
pd.options.mode.chained_assignment = None

"""
=======
Routine
======= 
"""

# LEZ Routine
# TO DO
warning_count=0
pollutants_temp=pollutants+['NO2']
pollutants_in=[w+'_in' for w in pollutants_temp]

# #### Read Input

#which type of input file
type=input_shape[-3:]
if type == 'shp':
    print 'Shape input detected'
    FASTRACE_input=gpd.read_file(input_shape)
    FASTRACE_wegendubbel=FASTRACE_input.loc[:,['X1','X2','Y1','Y2','NO2 kg/km','NOX kg/km', \
        'PM10 kg/km','PM25 kg/km','EC kg/km','road_type','length (m)','geometry','Fastrace_I', 'LV','ZV','V_FREEFLOW']]
    

    FASTRACE_wegendubbel.columns=['XA','XB','YA','YB','NO2_in','NOx_in','PM10','PM25', \
        'EC','roadtype','length','geometry','wegID','LV','ZV','V_FREEFLOW']
    FASTRACE_wegendubbel['NOx_in']*=1/(365.0*24)
    FASTRACE_wegendubbel['NO2_in']*=1/(365.0*24)
    FASTRACE_wegendubbel['BC_in']=FASTRACE_wegendubbel['EC']*1/(365.0*24)*1.5
    FASTRACE_wegendubbel['PM10_in']=(FASTRACE_wegendubbel['PM10'])/(365.0*24)
    FASTRACE_wegendubbel['PM25_in']=(FASTRACE_wegendubbel['PM25'])/(365.0*24)
    FASTRACE_wegendubbel[['LV','ZV']]=FASTRACE_wegendubbel[['LV','ZV']].astype(float)/24.0  # tellingen per uur
    FASTRACE_wegendubbel=FASTRACE_wegendubbel[['XA','XB','YA','YB','NOx_in','NO2_in','PM10_in','PM25_in',\
        'BC_in','roadtype','length','geometry','wegID','LV','ZV','V_FREEFLOW']]
else:
    print 'invalid input format'
    sys.exit()

### Corrections because of error in FASTRACE version
FASTRACE_wegendubbel['length']*=1000 
FASTRACE_wegendubbel[pollutants_in]/=1000 

    
    
print 'Check input NOx (ton/year):  '+str((FASTRACE_wegendubbel['NOx_in']*FASTRACE_wegendubbel['length']).sum()*365*24/1e6)
print 'Check input PM10 (ton/year): '+str((FASTRACE_wegendubbel['PM10_in']*FASTRACE_wegendubbel['length']).sum()*365*24/1e6)
print 'Input number of roads: ' + str(len(FASTRACE_wegendubbel))




### Multiplication with factor
if mult_factor == 1:
    print 'Provide multiplication factors.'
    for pol in pollutants_in:     
        var = raw_input('Provide factor for '+pol+': ')
        var=float(var)
        FASTRACE_wegendubbel[pol]=FASTRACE_wegendubbel[pol]*var
    print 'Multiplied with factor'
    print 'Check new NOx (ton/year):  '+str((FASTRACE_wegendubbel['NOx_in']*FASTRACE_wegendubbel['length']).sum()*365*24/1e6)
    print 'Check new PM10 (ton/year): '+str((FASTRACE_wegendubbel['PM10_in']*FASTRACE_wegendubbel['length']).sum()*365*24/1e6)


# #### Find duplicate entries  (roads that appear twice)
print '\n Start processing'
a=FASTRACE_wegendubbel.groupby(['XA','XB','YA','YB'])
FASTRACE_wegendubbel[pollutants_temp+['LV','ZV']]=a[pollutants_in+['LV','ZV']].transform('sum')
print str(len(FASTRACE_wegendubbel)-len(a))+' duplicate line segments removed'
FASTRACE_wegendubbel.drop_duplicates(subset=['XA','XB','YA','YB'],keep='first',inplace=True)


# #### Remove roads with zero length
old_length=len(FASTRACE_wegendubbel)
FASTRACE_wegendubbel = FASTRACE_wegendubbel[FASTRACE_wegendubbel.length != 0]
FASTRACE_wegendubbel=FASTRACE_wegendubbel[['XA','XB','YA','YB']+pollutants_temp+['roadtype','length','geometry','LV','ZV','V_FREEFLOW']]
print str(old_length-len(FASTRACE_wegendubbel))+' zero length segments removed'    
    

# #### Tunnels
# read tunnels
xls = pd.ExcelFile(tunnels_bridges)
sheets = xls.sheet_names
FASTRACE_wegendubbel['ExtraInfo']='none'
if 'Tunnels' in sheets:
    print '\n Processing Tunnels'
    tunnels=pd.read_excel(tunnels_bridges,sheetname='Tunnels')
    tunnel_list=tunnels['Naam'].unique()
    columns=FASTRACE_wegendubbel.columns

    for tunnel in tunnel_list:
        print 'Processing tunnel: '+tunnel
        tunnel_data=tunnels[tunnels['Naam']==tunnel]
        
        
        ## remove tunnel segments
        remove=tunnel_data[tunnel_data['Type']=='Remove']
        remove=FASTRACE_wegendubbel.merge(remove,on=['XA','XB','YA','YB'],how='right')
        length1=len(FASTRACE_wegendubbel)
        
        # check
        if remove['length'].isnull().any():
            print 'WARNING'
            print 'Some road segments are missing'
            print 'Unfound segments:'
            print remove[remove['length'].isnull()][['XA','YA','XB','YB']]
            warning_count+=1

        else:
           print 'Tunnel OK'
            
        # actual removal    
        temp=FASTRACE_wegendubbel.merge(remove[['XA','XB','YA','YB','Number']],how='left')
        FASTRACE_wegendubbel=temp[temp['Number'].isnull()][columns]
        length2=len(FASTRACE_wegendubbel)
        print 'Removed '+str(length1-length2)+' tunnel segments'
            
        
        ## add rest
        rest=tunnel_data[tunnel_data['Type']=='Rest']
        rest2=rest.merge(remove[pollutants_temp+['LV','ZV','V_FREEFLOW','Number','roadtype']],on='Number')
        rest2['geometry']=""
        rest2['LengthRest']=np.sqrt((rest2['XA']-rest2['XB'])**2+(rest2['YA']-rest2['YB'])**2)
        rest2['Length']=rest2['LengthRest']
        for i in rest2.index:
            geometry=LineString([(rest2.loc[i,'XA'],rest2.loc[i,'YA']),(rest2.loc[i,'XB'],rest2.loc[i,'YB'])])
            rest2.loc[i,'geometry']=geometry

        temp=FASTRACE_wegendubbel.append(rest2)
        temp.loc[temp['LengthRest'].notnull(),'ExtraInfo']='TunnelRest'+tunnel
        FASTRACE_wegendubbel=temp[columns]
        length3=len(FASTRACE_wegendubbel)
        print 'Added '+str(length3-length2)+' rest segments'
        
           
        # add tunnel exit
        uitgang=tunnel_data[tunnel_data['Type']=='Exit']
        temp=rest2.groupby('Number').sum()
        temp['Number']=temp.index
        remove=remove.merge(temp[['Number','LengthRest']],on='Number',how='left')
        
        remove.loc[remove['LengthRest'].isnull(),'LengthRest']=0
        remove['LengthRemoved']=remove['length']-remove['LengthRest']
        if not(uitgang.empty):
            total_em=pd.Series(index=pollutants_temp)
            for pol in pollutants_temp:
                total_em[pol]=sum(remove.loc[:,pol]*remove.loc[:,'LengthRemoved'])                  
                uitgang.loc[:,pol]=total_em[pol]*uitgang.loc[:,'Percentage']/(uitgang.loc[:,'Length'])
            
            uitgang=uitgang.merge(remove[['roadtype','Number']],on='Number')
            uitgang['geometry']=""
            for i in uitgang.index:
                if float(uitgang.loc[i,'XB']-uitgang.loc[i,'XA']) == 0:
                    uitgang.loc[i,'XB']=uitgang.loc[i,'XA']
                    teken=np.sign((uitgang.loc[i,'YB']-uitgang.loc[i,'YA']))
                    uitgang.loc[i,'YB']=uitgang.loc[i,'YA']+teken*float(uitgang.loc[i,'Length'])
                
                else: 
                    rico=float(uitgang.loc[i,'YB']-uitgang.loc[i,'YA'])/float(uitgang.loc[i,'XB']-uitgang.loc[i,'XA'])
                    teken=np.sign(uitgang.loc[i,'XB']-uitgang.loc[i,'XA'])
                    uitgang.loc[i,'XB']=uitgang.loc[i,'XA']+teken*np.sqrt(uitgang.loc[i,'Length']**2/(1+rico**2))
                    teken=np.sign((uitgang.loc[i,'YB']-uitgang.loc[i,'YA'])*rico)
                    uitgang.loc[i,'YB']=uitgang.loc[i,'YA']+teken*rico*np.sqrt(uitgang.loc[i,'Length']**2/(1+rico**2))     
               
                geometry=LineString([(uitgang.loc[i,'XA'],uitgang.loc[i,'YA']),(uitgang.loc[i,'XB'],uitgang.loc[i,'YB'])])
                uitgang.loc[i,'geometry']=geometry
            
            uitgang['length']=uitgang['Length']
            temp=FASTRACE_wegendubbel.append(uitgang)
            temp.loc[temp['Length'].notnull(),'ExtraInfo']='TunnelExit'+tunnel
            FASTRACE_wegendubbel=temp[columns]
            
            length4=len(FASTRACE_wegendubbel)
            print 'Added '+str(length4-length3)+' exit segments'
        
# #### Go over on single roads
print '\n Convert from double to single direction'
pollutants_x=[w+'_x' for w in pollutants_temp+['LV','ZV']]
pollutants_y=[w+'_y' for w in pollutants_temp+['LV','ZV']]
pollutants_total=[w+'_total' for w in pollutants_temp+['LV','ZV']]
FASTRACE_temp=pd.merge(FASTRACE_wegendubbel,FASTRACE_wegendubbel,left_on=['XA','XB','YA','YB'],right_on=['XB','XA','YB','YA'],how='left')
values_x=FASTRACE_temp[pollutants_x].fillna(0)
values_y=FASTRACE_temp[pollutants_y].fillna(0)
values_x.columns=pollutants_total
values_y.columns=pollutants_total
values=values_x+values_y
FASTRACE_temp[pollutants_total]=values
FASTRACE_temp=FASTRACE_temp[['XA_x','XB_x','YA_x','YB_x','roadtype_x','length_x','XA_y','geometry_x','ExtraInfo_x','V_FREEFLOW_x']+pollutants_total]
meenemen0=FASTRACE_temp.XA_y.isnull()
meenemen1=FASTRACE_temp.XA_x>FASTRACE_temp.XB_x
meenemen2= (FASTRACE_temp.XA_x == FASTRACE_temp.XB_x) & (FASTRACE_temp.YA_x >= FASTRACE_temp.YB_x)
meenemen=(meenemen0 | meenemen1) | meenemen2
FASTRACE=FASTRACE_temp.loc[meenemen,['XA_x','XB_x','YA_x','YB_x','roadtype_x','length_x','geometry_x','ExtraInfo_x','V_FREEFLOW_x']+pollutants_total]
FASTRACE.columns=['XA','XB','YA','YB','roadtype','length','geometry','ExtraInfo','V_FREEFLOW']+pollutants_temp+['LV','ZV']
print 'Converted to one directional roads'



# #### Remove roads with too small emissions
print '\n Remove roads with small emissions'
old_length=len(FASTRACE)
FASTRACE = FASTRACE[FASTRACE.NOx > cut_off_nox]
print str(old_length-len(FASTRACE))+' segments with too small emissions observed'
print 'Check NOx (ton/year):  '+str((FASTRACE['NOx']*FASTRACE['length']).sum()*365*24/1e6)
print 'Check PM10 (ton/year): '+str((FASTRACE['PM10']*FASTRACE['length']).sum()*365*24/1e6)  
print 'Final number of line sources '+str(len(FASTRACE))

# #### Add values for other parameters 

factors=pd.read_excel(file_factors)
IFDM=FASTRACE.merge(factors,on='roadtype')
IFDM['NO/NOX ratio']=(IFDM['NOx']-IFDM['NO2'])/IFDM['NOx']
IFDM=IFDM[['XA','YA','XB','YB','width','height','spread']+pollutants+\
  ['temp dependence','time factor','NO/NOX ratio','type','sector','LV','ZV','V_FREEFLOW','geometry','ExtraInfo']]
print '\n Added additional parameters (height, width...)'



# #### Bridges
if 'Bridges' in sheets:
    print '\n Processing Bridges'
    bridges=pd.read_excel(tunnels_bridges,sheetname='Bridges')
    bridges_list=bridges['Naam'].unique()
    columns=IFDM.columns
    for bridge in bridges_list:
        print 'Processing bridge: '+bridge
        bridge_data=bridges[bridges['Naam']==bridge]
        
        # check
        to_change=IFDM.merge(bridge_data,on=['XA','XB','YA','YB'],how='right')
        if to_change['width'].isnull().any():
            print 'WARNING'
            print 'Some road segments are missing'
            print 'Unfound segments:'
            print to_change[to_change['width'].isnull()][['XA','YA','XB','YB']]
            warning_count+=1

        else:
            print 'Bridge OK'
            
        temp=IFDM.merge(bridge_data,on=['XA','XB','YA','YB'],how='left')
        temp.loc[temp['NewHeight'].notnull(),'height']=temp['NewHeight']
        temp.loc[temp['NewHeight'].notnull(),'ExtraInfo']='Bridge'
        IFDM=temp[columns]
        
    print 'Processed bridges'


# #### Screen
if 'Screens' in sheets:
    print '\n Processing screens'
    screens=pd.read_excel(tunnels_bridges,sheetname='Screens')
    screens_list=screens['Naam'].unique()
    columns=IFDM.columns
    for screen in screens_list:
        print 'Processing screen: '+screen
        screen_data=screens[screens['Naam']==screen]
        if screen_data['Type'].isnull().any():
            to_change=IFDM.merge(screen_data,on=['XA','XB','YA','YB'],how='right')
            if to_change['spread'].isnull().any():
                print 'WARNING'
                print 'Some road segments are missing'
                print 'Unfound segments:'
                print to_change[to_change['spread'].isnull()][['XA','YA','XB','YB']]
                warning_count+=1
            else:
                print 'screen OK'
                
                temp=IFDM.merge(screen_data,on=['XA','XB','YA','YB'],how='left')
                temp.loc[temp['NewSpread'].notnull(),'spread']=temp['NewSpread']
                temp.loc[temp['NewSpread'].notnull(),'ExtraInfo']='Screen'
                IFDM=temp[columns]
   
        else:
            location=(IFDM['ExtraInfo']==screen_data.loc[screen_data.index[0],'Type'])\
                &(IFDM['XA']==screen_data.loc[screen_data.index[0],'XA'])
            if not any(location):
                print 'WARNING'
                print 'Some road segments are missing'
                print 'Unfound segments:'
                print screen_data[['XA','YA','XB','YB']]
                warning_count+=1
            else: 
                IFDM.loc[location,'spread']=screen_data.loc[screen_data.index[0],'NewSpread']
    
    print 'Processed screens'



# #### Embedment
if 'Sleuven' in sheets:
    print '\n Processing sleuven'
    sleuven=pd.read_excel(tunnels_bridges,sheetname='Sleuven')
    sleuven_list=sleuven['Naam'].unique()
    columns=IFDM.columns
    for sleuf in sleuven_list:
        print 'Processing sleuf: '+sleuf
        sleuf_data=sleuven[sleuven['Naam']==sleuf]
        
        if sleuf_data['Type'].isnull().any():

            to_change=IFDM.merge(sleuf_data,on=['XA','XB','YA','YB'],how='right')
            if to_change['width'].isnull().any():
                print 'WARNING'
                print 'Some road segments are missing'
                print 'Unfound segments:'
                print to_change[to_change['width'].isnull()][['XA','YA','XB','YB']]
                warning_count+=1

            else:
                print 'sleuf OK'
                
                temp=IFDM.merge(sleuf_data,on=['XA','XB','YA','YB'],how='left')
                temp.loc[temp['NewWidth'].notnull(),'height']=temp['NewWidth']
                temp.loc[temp['NewWidth'].notnull(),'ExtraInfo']='Sleuf'
                IFDM=temp[columns]

        else:
            location=(IFDM['ExtraInfo']==sleuf_data.loc[sleuf_data.index[0],'Type'])\
                &(IFDM['XA']==sleuf_data.loc[sleuf_data.index[0],'XA'])
            if not any(location):
                print 'WARNING'
                print 'Some road segments are missing'
                print 'Unfound segments:'
                print sleuf_data[['XA','YA','XB','YB']]
                warning_count+=1
            else:
                IFDM.loc[location,'spread']=sleuf_data.loc[sleuf_data.index[0],'NewWidth']

            
    print 'Processed sleuven'



# #### Start saving 
IFDM.V_FREEFLOW.fillna(50,inplace=True)
crs='epsg:31370'
geo_df = gpd.GeoDataFrame(IFDM, crs=crs, geometry=IFDM['geometry'])
geo_df.to_file(output_shape)
print 'Saved as Shapefile'


IFDM_save=IFDM.loc[:,['XA','YA','XB','YB','width','height','spread']+pollutants+\
 ['temp dependence','time factor','NO/NOX ratio','type','LV','ZV','V_FREEFLOW','V_FREEFLOW','sector','ExtraInfo']]
IFDM_save.loc[:,['XA','YA','XB','YB']]*=1.0/1000.0
IFDM_save.to_excel(output_excel,index=False,sheet_name='LineSources')
print 'Saved as Excel file'


# #### Some Statistics 
print ''
print 'Some Statistics'
IFDM_save['distance']=np.sqrt((IFDM_save['XA']-IFDM_save['XB'])**2+(IFDM_save['YA']-IFDM_save['YB'])**2)
print 'Total NOx  emissions (ton/year): '+str((IFDM_save['NOx']*IFDM_save['distance']).sum()*365*24/1e3)
print 'Total PM10 emissions (ton/year): '+str((IFDM_save['PM10']*IFDM_save['distance']).sum()*365*24/1e3)
print 'Total PM25 emissions (ton/year): '+str((IFDM_save['PM25']*IFDM_save['distance']).sum()*365*24/1e3)

print '\n \n'
print 'Done.'
if warning_count>0:
    print 'Observed '+str(warning_count)+' warnings during processing!' 

