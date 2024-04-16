#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 09:07:03 2024

@author: leohoinaski
"""
import numpy as np
import geopandas as gpd

vehicularCategories =['LIGHT','COMMERCIAL-LIGHT','MOTORCYCLES','HEAVY']

def scrappage(vehicularCategory,t):
    
    if vehicularCategory =='LIGHT':
        a = 1.798
        b = -0.137
        s = 1 - np.exp(-np.exp(a+b*t))
    elif vehicularCategory =='COMMERCIAL-LIGHT':
        a = 1.618
        b = -0.141
        s = 1 - np.exp(-np.exp(a+b*t))
    elif vehicularCategory =='MOTORCYCLES':
        if t<=5:
            a = 1.317
            b = -0.175
        else:
            a = 0.923
            b = -0.93
        s = 1 - np.exp(-np.exp(a+b*t))
    elif vehicularCategory =='HEAVY':
        a = 0.100
        t0 = 17.0
        s = (1/(1 + np.exp(a*(t - t0)))) + (1/(1 + np.exp(a*(t + t0))))
    else: 
        print('*******Wrong category********')
        
    return s

def yearModelByVcat(interFolder,filePath,dfVCat,dfYM,
                    shapeFolder,vehicularCategories):
    munShp = gpd.read_file(shapeFolder)
    dfYM[vehicularCategories]=np.nan
    for IBGE_CODE in munShp['CD_MUN']:
        cityYM = dfYM[np.array(dfYM['IBGE_CODE']).astype(int)== int(IBGE_CODE)]
        cityCat = dfVCat[np.array(dfVCat['IBGE_CODE']).astype(int) == int(IBGE_CODE)]
        if np.sum([np.array(dfYM['IBGE_CODE']).astype(int)== int(IBGE_CODE)]) == 0:
            print('City not found')
        else:
            cityPropFleet = cityCat[vehicularCategories]/np.array(cityCat[vehicularCategories].sum(axis=1)).astype(float)[0]
            for vcat in vehicularCategories:
                dfYM[vcat][np.array(dfYM['IBGE_CODE']).astype(int)== int(IBGE_CODE)] = cityYM['N']*np.nanmean(np.array(cityPropFleet[vcat]))
            
        #munShp[munShp['CD_GEOCMU']== IBGE_CODE]
        #dfFuel[dfFuel['IBGE_CODE']== IBGE_CODE]
    dfYM.to_csv(interFolder +'/BRAVES_yearModelByVcat_' + filePath.split('/')[-1]) 
    return dfYM,munShp

# def only_numerics(seq):
#     seq_type= type(seq)
#     return seq_type().join(filter(seq_type.isdigit, seq))

def scrapppageAll(interFolder,filePath,vehicularCategory,dfYM,munShp):
    for vcat in vehicularCategories:
        dfYM[vcat+'_scrapPerc'] = np.nan
    for IBGE_CODE in munShp['CD_MUN']:
        cityData = dfYM[np.array(dfYM['IBGE_CODE']).astype(int)== int(IBGE_CODE)]
        if cityData.shape[0]>0:
            ts = np.nanmax(cityData['yearModel'])-cityData['yearModel']
            for vcat in vehicularCategories:
                print(vcat)
                s = []
                for t in ts:
                    s.append(scrappage(np.array(vcat),t))
                dfYM[vcat+'_scrapPerc'][np.array(dfYM['IBGE_CODE']).astype(int)== int(IBGE_CODE)] = np.array(s)
    dfYM.to_csv(interFolder +'/BRAVES_scrapppageAll_' + filePath.split('/')[-1]) 
    return dfYM     
        
        
    
    
    
    
    

