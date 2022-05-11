#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                              mergeRoadEmiss_v1.py

This function merges files created by roadEmiss_v1.py which contains the 
spatial disaggregated emission in a specific domain from BRAVES by states. 
Emission input unitis are in ton/year.

Inputs: 
    
    outPath: path to folder with roadDensity_v1.py and roadEmiss_v1.py
        outputs. Folder to outputs from this function.
    
    years: respective years of emission inventories.
    
    IBGE_CODES:  State ID
    
    
Outputs:
    
    BrRoadEmiss_'+name+str(year)+'.csv': file with merged emissions 
    

Last update = 29/09/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br

-------------------------------------------------------------------------------
"""

# Importing packages
import geopandas as gpd
import pandas as pd
from shapely import wkt
import shapely.speedups
shapely.speedups.enable()

# ================================PROCESSING===================================

def mergeRoadEmiss(outPath,years,IBGE_CODES,roadDensPrefix,typeEmiss):
    print('===================STARTING mergeRoadEmiss_v1.py=======================')
    
    basefile = outPath + '/baseGrid_'+roadDensPrefix+'.csv'
    base = pd.read_csv(basefile)
    base = base.set_index('geometry')
    base = base.drop(base.columns[0], axis=1)
    
    names = ['BySource_Light_'+typeEmiss+'_','BySource_ComLight_'+typeEmiss+'_',
             'BySource_Heavy_'+typeEmiss+'_','BySource_Motorcycle_'+typeEmiss+'_']
    
    for name in names:
        
        for year in years:
            kk=0
            # From states to Brazil - merging emissions
            for IBGE_CODE in IBGE_CODES:  
                IBGE_CODE = int(IBGE_CODE)
                # Opening road density file
                #df = pd.read_csv('/home/leohoinaski/Desktop/smartPlatform/roadLength_UF_'+str(IBGE_CODE)+'.csv')        
                file = outPath + '/roadEmiss_'+name +\
                    roadDensPrefix +'_' + str(year)+ '_UF_'+str(IBGE_CODE)+'.csv'
                df = pd.read_csv(file)
                df = df.set_index('geometry')
                df = df.drop(df.columns[0], axis=1)
                df = df.drop(df.columns[0], axis=1)
                mergedDf = base.join(df)
                mergedDf = mergedDf.fillna(0)
                if kk == 0:
                    brEmiss = mergedDf
                else:
                    for jj in range(0,brEmiss.shape[1]):
                        brEmiss[brEmiss.iloc[:,jj].name] = \
                        brEmiss[brEmiss.iloc[:,jj].name] + \
                        mergedDf[brEmiss.iloc[:,jj].name]
                kk = kk + 1
            brEmiss['geometry'] =  brEmiss.index 
            brEmiss['geometry'] = brEmiss['geometry'].apply(wkt.loads) 
            brEmiss = gpd.GeoDataFrame(brEmiss, geometry='geometry')        
            brEmiss.to_csv(outPath+'/mergeRoadEmiss_'+name+\
                roadDensPrefix +'_' + str(year)+'.csv')
        


