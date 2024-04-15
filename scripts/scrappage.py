#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 09:07:03 2024

@author: leohoinaski
"""
import numpy as np
import geopandas as gpd

vehicularCategories =['LIGHT','COMMERCIAL-LIGHT','MOTORCYLES','HEAVY']

def scrappage(vehicularCategory,t):
    
    if vehicularCategory =='LIGHT':
        a = 1.798
        b = -0.137
        s = 1 - np.exp(-np.exp(a+b*t))
    elif vehicularCategory =='COMMERCIAL-LIGHT':
        a = 1.618
        b = -0.141
        s = 1 - np.exp(-np.exp(a+b*t))
    elif vehicularCategory =='MOTORCYLES':
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
    for IBGE_CODE in munShp['CD_GEOCMU']:
        cityYM = dfYM[dfYM['IBGE_CODE']== IBGE_CODE]
        cityCat = dfVCat[np.array(dfVCat['IBGE_CODE']).astype(int) == int(IBGE_CODE)]
        cityPropFleet = cityCat[vehicularCategories]/cityCat[vehicularCategories].sum()
        for vcat in vehicularCategories:
            dfYM[vcat][dfYM['IBGE_CODE']== IBGE_CODE] = cityYM['N']*cityPropFleet[vcat]
        
        #munShp[munShp['CD_GEOCMU']== IBGE_CODE]
        #dfFuel[dfFuel['IBGE_CODE']== IBGE_CODE]
    dfYM.to_csv(interFolder +'/BRAVES_' + filePath.split('/')[-1]) 
    return dfYM