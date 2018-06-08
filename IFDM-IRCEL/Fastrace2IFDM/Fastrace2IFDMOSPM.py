"""
Conversion of FastRace output to IFDM input
===============================================================

Operations
usage: python Fastrace2IFDMOSPM.py [-h] [-c cut-off] [-o tag] [-m tag] [-t filename]
                            [-f filename]
                            inputfile outputbase

Fastrace2IFDMOSPM.py

positional arguments:
  inputfile             Input file name
  outputbase            Output file base

optional arguments:
  -h, --help            show this help message and exit
  -c cut-off, --cut-off-nox cut-off
                        Cut-off for NOx. Default:0.0
  -o tag, --ospm tag    OSPM: 0 or 1? Default: 0
  -m tag, --multiplication tag
                        Multiplication with factors: 0 or 1? Factors will be
                        requested later on. Default: 0
  -s tag, --shipping tag
                        Add wind speed colums for shipping: 0 or 1. Default: 0
  -t filename, --tunnels_file filename
                        File with tunnels, bridges... Default:
                        TunnelsBridges.xlsx
  -f filename, --conversion filename
                        File with conversion factors. Default:
                        ConversionFactors.xlsx
"""

import sys
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString
pd.options.mode.chained_assignment = None

def main(args):

    """
    =======
    Input
    =======
    """
    # pollutants for output (nota Hans: pas helemaal op het einde in rekening brengen)
    pollutants=['NOx','PM10','PM25','BC']

    # column headers: link names to use with input  [ name as will be used : name as in input]
    column_headers={'XA':'X1',
        'XB':'X2',
        'YA':'Y1',
        'YB':'Y2',
        'NOx':'NOX kg/km',
        'NO2':'NO2 kg/km',
        'PM10':'PM10 kg/km',
        'PM25':'PM25 kg/km',
        'BC':'EC kg/km',                # note: EC-BC conversion occurs 
        'roadtype':'road_type',
        'length':'length (m)',
        'geometry':'geometry',
        'wegID':'Fastrace_I',
        'LV':'LV',
        'ZV':'ZV',
        'V_FREEFLOW':'V_FREEFLOW'}

    # default conversion factor
    conv_emissions=1.0/(24.0*365.0)/1000.0     # required: kg/km/h {this value assumes input in ton/km/year (default fastrace output)]
    conv_counts=1.0/24.0                # required: veh/hour [this values ssumes input in vehicles/day (default fastrace output)] 
    conv_ECBC=1.5
    conv_length=1000                    # required: m [this values assumes input in km]

    """
    =======
    Routine
    ======= 
    """
    # #### Convert input arguments
    ###############################

    input_shape=args.inputfile
    outputbase=args.outputbase
    cut_off_nox=args.cut_off
    mult_factor=bool(args.mult)
    ospm=bool(args.ospm)
    shipping=bool(args.shipping)
    tunnels_bridges=args.tunnels_bridges
    file_factors=args.file_factors

    print 'set-up'
    print 'Inputfile: '+input_shape
    print 'Outputbase: '+outputbase
    print 'OSPM?: '+str(ospm)
    print 'Cut-off NOx: '+str(cut_off_nox)
    print 'Shipping: '+str(shipping)
    print 'Multiplication?: '+str(mult_factor)
    print 'Tunnels file: '+tunnels_bridges
    print 'Conversion factors file: '+file_factors
    print '\n\n'
    
    output_excel=outputbase+'.xlsx'
    output_shape=outputbase+'.shp'



    warning_count=0
    pollutants_withno2=pollutants+['NO2']

    # #### Read Input
    ################
    # check for shape input
    type=input_shape[-3:]
    if type == 'shp':
        print 'Shape input detected'
        FASTRACE_input=gpd.read_file(input_shape)
        
        # correction for OSPM
        if not ospm:
            print 'No ospm, add temporary zero columns'
            FASTRACE_input[column_headers['LV']]=0
            FASTRACE_input[column_headers['ZV']]=0
            FASTRACE_input[column_headers['V_FREEFLOW']]=0
        
        # read input
        FASTRACE_wegendubbel=FASTRACE_input[column_headers.values()]
        FASTRACE_wegendubbel.columns=column_headers.keys()
        
        # conversions
        FASTRACE_wegendubbel[pollutants_withno2]*=conv_emissions
        if 'BC' in pollutants:
            FASTRACE_wegendubbel.BC=FASTRACE_wegendubbel.BC*conv_ECBC
        FASTRACE_wegendubbel[['LV','ZV']]=FASTRACE_wegendubbel[['LV','ZV']].astype(float)
        FASTRACE_wegendubbel[['LV','ZV']]*=conv_counts
        FASTRACE_wegendubbel['length']*=conv_length
        
    else:
        
        print 'invalid input format'
        sys.exit()
        
    # check input    
    check_emissions(FASTRACE_wegendubbel,pollutants_withno2)

    ### Multiplication with factor
    ###############################
    if mult_factor:
        print 'Provide multiplication factors.'
        for pol in pollutants_withno2:     
            var = raw_input('Provide factor for '+pol+': ')
            FASTRACE_wegendubbel[pol]=FASTRACE_wegendubbel[pol]*float(var)
        print 'Multiplied with factor'
        check_emissions(FASTRACE_wegendubbel,pollutants_withno2)


    # #### Find duplicate entries  (roads that appear twice)
    #######################################################
    print '\n Start processing'
    old_length=len(FASTRACE_wegendubbel)
    grouped=FASTRACE_wegendubbel.groupby(['XA','XB','YA','YB'])
    FASTRACE_wegendubbel[pollutants_withno2+['LV','ZV']]=grouped[pollutants_withno2+['LV','ZV']].transform('sum')
    print str(len(FASTRACE_wegendubbel)-len(grouped))+' duplicate line segments removed'
    FASTRACE_wegendubbel.drop_duplicates(subset=['XA','XB','YA','YB'],keep='first',inplace=True)
    check_emissions(FASTRACE_wegendubbel,pollutants_withno2)

    # #### Remove roads with zero length
    old_length=len(FASTRACE_wegendubbel)
    FASTRACE_wegendubbel = FASTRACE_wegendubbel[FASTRACE_wegendubbel.length != 0]
    print str(old_length-len(FASTRACE_wegendubbel))+' zero length segments removed'    
    FASTRACE_wegendubbel['ExtraInfo']=''
    FASTRACE_wegendubbel['ExtraInfo']='none'
    
    # #### Corrections (for IRCEL)
    ##############################
    # read corrections
    xls = pd.ExcelFile(tunnels_bridges)
    sheets = xls.sheet_names    

    if 'Corrections' in sheets:
        print '\n Processing Corrections'
        corrections=pd.read_excel(tunnels_bridges,sheetname='Corrections')
        corr_list=corrections['Number'].unique()

        for corr in corr_list:
            print 'Processing Correction: '+str(corr)
            corr_data=corrections[corrections['Number']==corr]           
            
            ## remove corrected segments
            remove=corr_data[corr_data['Type']=='Remove']
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
               print 'Correction OK'
                
            # actual removal    
            FASTRACE_wegendubbel=FASTRACE_wegendubbel[~FASTRACE_wegendubbel.wegID.isin(remove.wegID)]
            length2=len(FASTRACE_wegendubbel)
            print 'Removed '+str(length1-length2)+' tunnel segments'                
            
            ## add rest
            rest=corr_data[corr_data['Type']=='Rest']
            rest2=rest.merge(remove[pollutants_withno2+['LV','ZV','V_FREEFLOW','Number','roadtype']],on='Number')
            if not rest2.empty:
                rest2['geometry']=rest2.apply(lambda p: LineString([(p['XA'],p['YA']),(p['XB'],p['YB'])]),axis=1)
                rest2['length']=rest2.apply(lambda p:p.geometry.length,axis=1)
                rest2['ExtraInfo']='Correction'+corr
                rest2['wegID']=-10
                merged=FASTRACE_wegendubbel.append(rest2[column_headers.keys()+['ExtraInfo']])
                FASTRACE_wegendubbel=merged
                length3=len(FASTRACE_wegendubbel)
                print 'Added '+str(length3-length2)+' rest segments'       
        
        
    # #### Tunnels
    ###############
    # read tunnels
    if 'Tunnels' in sheets:
        print '\n Processing Tunnels'
        tunnels=pd.read_excel(tunnels_bridges,sheetname='Tunnels')
        tunnel_list=tunnels['Naam'].unique()

        for tunnel in tunnel_list:
            print 'Processing tunnel: '+tunnel
            tunnel_data=tunnels[tunnels['Naam']==tunnel]            
            
            # removal (including check)
            remove=tunnel_data[tunnel_data['Type']=='Remove']
            remove.drop('Length',axis=1,inplace=True)
            remove=FASTRACE_wegendubbel.merge(remove,on=['XA','XB','YA','YB'],how='right')
            length1=len(FASTRACE_wegendubbel)
            if remove['length'].isnull().any():
                print 'WARNING'
                print 'Some road segments are missing'
                print 'Unfound segments:'
                print remove[remove['length'].isnull()][['XA','YA','XB','YB']]
                warning_count+=1
            else:
               print 'Tunnel OK'                
            # actual removal    
            FASTRACE_wegendubbel=FASTRACE_wegendubbel[~FASTRACE_wegendubbel.wegID.isin(remove.wegID)]
            length2=len(FASTRACE_wegendubbel)
            print 'Removed '+str(length1-length2)+' tunnel segments'                
            
            ## add rest
            rest=tunnel_data[tunnel_data['Type']=='Rest']
            rest2=rest.merge(remove[pollutants_withno2+['LV','ZV','V_FREEFLOW','Number','roadtype']],on='Number')
            if not rest2.empty:
                rest2['geometry']=rest2.apply(lambda p: LineString([(p['XA'],p['YA']),(p['XB'],p['YB'])]),axis=1)
                rest2['length']=rest2.apply(lambda p:p.geometry.length,axis=1)
                rest2['ExtraInfo']='TunnelRest'+tunnel
                rest2['wegID']=-10
                merged=FASTRACE_wegendubbel.append(rest2[column_headers.keys()+['ExtraInfo']])
                FASTRACE_wegendubbel=merged
            length3=len(FASTRACE_wegendubbel)
            print 'Added '+str(length3-length2)+' rest segments'            
               
            # add tunnel exit
            uitgang=tunnel_data[tunnel_data['Type']=='Exit']
            if not(uitgang.empty):
                # add rest length
                if not rest2.empty:
                    temp=rest2.groupby('Number',as_index=False).sum()
                    remove_withrest=remove.merge(temp[['Number','length']],on='Number',how='left')            
                    remove_withrest.fillna(0,inplace=True)
                    remove_withrest['LengthRemoved']=remove_withrest['length_x']-remove_withrest['length_y']
                else:
                    remove_withrest=remove
                    remove_withrest['LengthRemoved']=remove['length']
                total_em=pd.Series(index=pollutants_withno2)
                for pol in pollutants_withno2:
                    total_em[pol]=sum(remove_withrest.loc[:,pol]*remove_withrest.loc[:,'LengthRemoved'])                  
                    uitgang.loc[:,pol]=total_em[pol]*uitgang.loc[:,'Percentage']/(uitgang.loc[:,'Length'])
                uitgang=uitgang.merge(remove_withrest[['roadtype','Number','LV','ZV','V_FREEFLOW']],on='Number')
                
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
                uitgang['ExtraInfo']='TunnelExit'+tunnel
                uitgang['wegID']=-10
                
                merged=FASTRACE_wegendubbel.append(uitgang[column_headers.keys()+['ExtraInfo']])
                FASTRACE_wegendubbel=merged                
                length4=len(FASTRACE_wegendubbel)
                print 'Added '+str(length4-length3)+' exit segments'
            
    # #### Go over on single-directional roads
    ###########################################
    print '\n Convert from double to single direction'
    pollutants_x=[w+'_x' for w in pollutants_withno2+['LV','ZV']]
    pollutants_y=[w+'_y' for w in pollutants_withno2+['LV','ZV']]
    pollutants_total=[w+'_total' for w in pollutants_withno2+['LV','ZV']]
    FASTRACE_temp=pd.merge(FASTRACE_wegendubbel,FASTRACE_wegendubbel,left_on=['XA','XB','YA','YB'],
        right_on=['XB','XA','YB','YA'],how='left')
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
    FASTRACE.columns=['XA','XB','YA','YB','roadtype','length','geometry','ExtraInfo','V_FREEFLOW']+pollutants_withno2+['LV','ZV']
    print 'Converted to one directional roads'



    # #### Remove roads with too small emissions
    ############################################
    print '\n Remove roads with small emissions'
    old_length=len(FASTRACE)
    FASTRACE = FASTRACE[FASTRACE.NOx > cut_off_nox]    
    print str(old_length-len(FASTRACE))+' segments with too small emissions observed'
    check_emissions(FASTRACE,pollutants_withno2)

    # #### Add values for other parameters 
    factors=pd.read_excel(file_factors)
    IFDM=FASTRACE.merge(factors,on='roadtype')
    IFDM['NO/NOX ratio']=(IFDM['NOx']-IFDM['NO2'])/IFDM['NOx']
    IFDM=IFDM[['XA','YA','XB','YB','width','height','spread']+pollutants+\
      ['temp dependence','time factor','NO/NOX ratio','type','sector','LV','ZV','V_FREEFLOW','geometry','ExtraInfo']]
    print '\n Added additional parameters (height, width...)'



    # #### Bridges
    ###############
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
                #replacement
                print 'Bridge OK'
                temp=IFDM.merge(bridge_data,on=['XA','XB','YA','YB'],how='left')
                temp.loc[temp['NewHeight'].notnull(),'height']=temp['NewHeight']
                temp.loc[temp['NewHeight'].notnull(),'ExtraInfo']='Bridge'
                IFDM=temp[columns]
            
        print 'Processed bridges'


    # #### Screen
    #############
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
                    temp.loc[temp['NewSpread'].notnull(),'spread']=np.sqrt(temp.loc[temp['NewSpread'].notnull(),'spread']**2.0\
                        +temp.loc[temp['NewSpread'].notnull(),'NewSpread']**2.0)
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
                
                    IFDM.loc[location,'spread']=np.sqrt(screen_data.loc[screen_data.index[0],'NewSpread']**2.0+IFDM.loc[location,'spread']**2.0)
        
        print 'Processed screens'



    # #### Embedment
    #################
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
                    temp.loc[temp['NewWidth'].notnull(),'width']=temp['NewWidth']
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
                    IFDM.loc[location,'width']=sleuf_data.loc[sleuf_data.index[0],'NewWidth']

                
        print 'Processed sleuven'



    # #### Start saving 
    ###################
    IFDM.V_FREEFLOW.fillna(50,inplace=True)
    crs='epsg:31370'
    geo_df = gpd.GeoDataFrame(IFDM, crs=crs, geometry=IFDM['geometry'])
    geo_df.to_file(output_shape)
    print 'Saved as Shapefile'
    
    if shipping:
        IFDM['Wa1']=0
        IFDM['Wa2']=0
        IFDM['Wa3']=0
    
    if ospm:
        if shipping:
            columns_output=['XA','YA','XB','YB','width','height','spread']+pollutants+\
                ['temp dependence','time factor','NO/NOX ratio','type','Wa1','Wa2','Wa3','LV','ZV','V_FREEFLOW','V_FREEFLOW','sector','ExtraInfo']
        else: 
            columns_output=['XA','YA','XB','YB','width','height','spread']+pollutants+\
                ['temp dependence','time factor','NO/NOX ratio','type','LV','ZV','V_FREEFLOW','V_FREEFLOW','sector','ExtraInfo']
    else:
        if shipping:
            columns_output=['XA','YA','XB','YB','width','height','spread']+pollutants+\
                ['temp dependence','time factor','NO/NOX ratio','type','Wa1','Wa2','Wa3','sector','ExtraInfo']
        else:
            columns_output=['XA','YA','XB','YB','width','height','spread']+pollutants+\
                ['temp dependence','time factor','NO/NOX ratio','type','sector','ExtraInfo']
        

    IFDM_save=IFDM.loc[:,columns_output]
    IFDM_save.loc[:,['XA','YA','XB','YB']]*=1.0/1000.0
    IFDM_save.to_excel(output_excel,index=False,sheet_name='LineSources')
    print 'Saved as Excel file'


    # #### Some Statistics 
    print ''
    print 'Some Statistics'
    IFDM_save['length']=np.sqrt(((IFDM_save['XA']-IFDM_save['XB'])**2+(IFDM_save['YA']-IFDM_save['YB'])**2).astype(float))*1e3
    check_emissions(IFDM_save,pollutants)
    print '\n \n'
    print 'Done.'
    if warning_count>0:
        print 'Observed '+str(warning_count)+' warnings during processing!' 

def check_emissions(input,pollutants):
    for pol in pollutants:
        tot_em=(input[pol]*input.length).sum()*365*24/1e6
        print 'Check emissions '+pol+' (ton/year):  '+str(tot_em)
    print 'Check number of roads: ' + str(len(input))
        
        
if __name__ == '__main__':
    
    """"
    ==============
    Input parsing  (nog wat aan te passen)
    ==============
    """    
    
    parser = argparse.ArgumentParser(description='Fastrace2IFDMOSPM.py')
    parser.add_argument('inputfile', metavar='inputfile', type=str, help='Input file name')
    parser.add_argument('outputbase', metavar='outputbase', type=str, help='Output file base')
    parser.add_argument('-c','--cut-off-nox',dest='cut_off',metavar='cut-off', type=float, 
        help='Cut-off for NOx. Default:0.0',default=0.0)
    parser.add_argument('-o','--ospm',dest='ospm',metavar='tag', type=int, help='OSPM: 0 or 1? Default: 0',default=0)
    parser.add_argument('-m','--multiplication',dest='mult',metavar='tag', type=int, 
        help='Multiplication with factors: 0 or 1? Factors will be requested later on. Default: 0',default=0)
    parser.add_argument('-s','--shipping',dest='shipping',metavar='tag', 
        type=str, help='Add wind speed colums for shipping: 0 or 1. Default: 0',default=0)
    parser.add_argument('-t','--tunnels_file',dest='tunnels_bridges',metavar='filename', 
        type=str, help='File with tunnels, bridges... Default: TunnelsBridges.xlsx',default='TunnelsBridges.xlsx')    
    parser.add_argument('-f','--conversion',dest='file_factors',metavar='filename', 
        type=str, help='File with conversion factors. Default: ConversionFactors.xlsx',default='ConversionFactors.xlsx')
    args = parser.parse_args()

    
    main(args)
