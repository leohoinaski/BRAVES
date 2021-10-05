#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                        BRAVES_ChemicalSpec_v1.py

This function speciate the emissions from BRAVES in chemical species. 
These estimates are based on Speciate from US-EPA 
(https://www.epa.gov/air-emissions-modeling/speciate). This function contains 2
subfunction, ChemicalSpeciationLight for light vehicles (light-duty, motorcycle
and commercial-light) and ChemicalSpeciationHeavy for heavy-duty. Each function 
uses a different speciation profile. 


Inputs: 
    
    roadX: 
    
    dfSpc: data from 'BRAVES_speciation.csv'
    
    smm: data from 'CMAQ_speciesMW.csv'
    
    conver: unit conversion factor / conver = 1000000 from ton/year to g/year
    
    
Outputs:
    
    dataEmissX: speciated and spatial disaggregated data from BRAVES.
    

Last update = 29/09/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br

-------------------------------------------------------------------------------
"""
# Importando bibliotecas
import pandas as pd

def ChemicalSpeciationLight(roadX,dfSpc,smm,conver):
    # Creating dataEmiss matrix
    dataEmissX=pd.DataFrame()
    print('===================STARTING BRAVES_ChemicalSpec_v1.py=======================')
    # Filling VOC emissions
    for ii in range(0,dfSpc.shape[0]):
        dataEmissX[dfSpc.ID[ii]] = dfSpc['LightVOC'][ii]*roadX.iloc[:,6]*conver/smm.iloc[ii,1] +\
           dfSpc['LightPM'][ii]*roadX.iloc[:,10]*conver/smm.iloc[ii,1]  +\
               dfSpc['road'][ii]*roadX.iloc[:,18]*conver/smm.iloc[ii,1] +\
                   dfSpc['road'][ii]*roadX.iloc[:,17]*conver/smm.iloc[ii,1] +\
                       dfSpc['brakes'][ii]*roadX.iloc[:,16]*conver/smm.iloc[ii,1] +\
                           dfSpc['tires'][ii]*roadX.iloc[:,16]*conver/smm.iloc[ii,1] +\
                               dfSpc['LightEvap'][ii]*roadX.iloc[:,19]*conver/smm.iloc[ii,1] +\
                                   dfSpc['LightEvap'][ii]*roadX.iloc[:,20]*conver/smm.iloc[ii,1] +\
                                       dfSpc['LightEvap'][ii]*roadX.iloc[:,21]*conver/smm.iloc[ii,1] +\
                                           dfSpc['LightEvap'][ii]*roadX.iloc[:,14]*conver/smm.iloc[ii,1]
    
    # Filling direct emissions
    dataEmissX['ALDX'] = roadX.iloc[:,9]*conver/smm[smm.iloc[:,0]=='ALDX'].MM   
    dataEmissX['CH4_INV'] = roadX.iloc[:,7]*conver
    dataEmissX['CH4'] = roadX.iloc[:,7]*conver/smm[smm.iloc[:,0]=='CH4'].MM
    dataEmissX['CO'] = roadX.iloc[:,4]*conver/smm[smm.iloc[:,0]=='CO'].MM
    dataEmissX['CO2_INV'] = roadX.iloc[:,11]*conver
    dataEmissX['N2O_INV'] = roadX.iloc[:,12]*conver
    dataEmissX['NO'] = roadX.iloc[:,8]*conver*0.495/smm[smm.iloc[:,0]=='NO'].MM # Para motor a diesel... considerei
    dataEmissX['NO2'] = roadX.iloc[:,8]*conver*0.505/smm[smm.iloc[:,0]=='NO2'].MM # Para motor a diesel... considerei
    dataEmissX['PMC'] = roadX.iloc[:,10]*conver
    dataEmissX['PNCOM'] = dataEmissX['POC']*0.25
    dataEmissX['PMOTHR'] = 1-(dataEmissX['POC'] + dataEmissX['PEC'] +
                              dataEmissX['PSO4'] + dataEmissX['PNO3'] +
                              dataEmissX['PNH4'] + dataEmissX['PNCOM'] +
                              dataEmissX['PFE'] + dataEmissX['PAL'] +
                              dataEmissX['PSI'] + dataEmissX['PTI'] +
                              dataEmissX['PCA'] + dataEmissX['PMG'] +
                              dataEmissX['PK'] + dataEmissX['PMN'] +
                              dataEmissX['PNA'] + dataEmissX['PCL'] +
                              dataEmissX['PH2O'])
    dataEmissX['SO2'] = roadX.iloc[:,13]*conver/smm[smm.iloc[:,0]=='SO2'].MM
    dataEmissX['VOC_INV'] = roadX.iloc[:,6]*conver + roadX.iloc[:,14]*conver +\
        roadX.iloc[:,19]*conver + roadX.iloc[:,20]*conver + roadX.iloc[:,21]*conver
    return dataEmissX

def ChemicalSpeciationHeavy(roadX,dfSpc,smm,conver):  
    # Creating dataEmiss matrix
    dataEmissX=pd.DataFrame()
    # Filling VOC emissions
    for ii in range(0,dfSpc.shape[0]):
        dataEmissX[dfSpc.ID[ii]] = dfSpc['HeavyVOC'][ii]*roadX.iloc[:,6]*conver/smm.iloc[ii,1] +\
           dfSpc['HeavyPM'][ii]*roadX.iloc[:,10]*conver/smm.iloc[ii,1]  +\
               dfSpc['road'][ii]*roadX.iloc[:,18]*conver/smm.iloc[ii,1] +\
                   dfSpc['road'][ii]*roadX.iloc[:,17]*conver/smm.iloc[ii,1] +\
                       dfSpc['brakes'][ii]*roadX.iloc[:,16]*conver/smm.iloc[ii,1] +\
                           dfSpc['tires'][ii]*roadX.iloc[:,16]*conver/smm.iloc[ii,1] +\
                               dfSpc['HeavyEvap'][ii]*roadX.iloc[:,19]*conver/smm.iloc[ii,1] +\
                                   dfSpc['HeavyEvap'][ii]*roadX.iloc[:,20]*conver/smm.iloc[ii,1] +\
                                       dfSpc['HeavyEvap'][ii]*roadX.iloc[:,21]*conver/smm.iloc[ii,1] +\
                                           dfSpc['HeavyEvap'][ii]*roadX.iloc[:,14]*conver/smm.iloc[ii,1]
    
    # Filling direct emissions
    dataEmissX['ALDX'] = roadX.iloc[:,9]*conver/smm[smm.iloc[:,0]=='ALDX'].MM  
    dataEmissX['CH4_INV'] = roadX.iloc[:,7]*conver
    dataEmissX['CH4'] = roadX.iloc[:,7]*conver/smm[smm.iloc[:,0]=='CH4'].MM
    dataEmissX['CO'] = roadX.iloc[:,4]*conver/smm[smm.iloc[:,0]=='CO'].MM
    dataEmissX['CO2_INV'] = roadX.iloc[:,11]*conver
    dataEmissX['N2O_INV'] = roadX.iloc[:,12]*conver
    dataEmissX['NO'] = roadX.iloc[:,8]*conver*0.495/smm[smm.iloc[:,0]=='NO'].MM # Para motor a diesel... considerei
    dataEmissX['NO2'] = roadX.iloc[:,8]*conver*0.505/smm[smm.iloc[:,0]=='NO2'].MM # Para motor a diesel... considerei
    dataEmissX['PMC'] = roadX.iloc[:,10]*conver
    dataEmissX['PNCOM'] = dataEmissX['POC']*0.25
    dataEmissX['PMOTHR'] = 1-(dataEmissX['POC'] + dataEmissX['PEC'] +
                              dataEmissX['PSO4'] + dataEmissX['PNO3'] +
                              dataEmissX['PNH4'] + dataEmissX['PNCOM'] +
                              dataEmissX['PFE'] + dataEmissX['PAL'] +
                              dataEmissX['PSI'] + dataEmissX['PTI'] +
                              dataEmissX['PCA'] + dataEmissX['PMG'] +
                              dataEmissX['PK'] + dataEmissX['PMN'] +
                              dataEmissX['PNA'] + dataEmissX['PCL'] +
                              dataEmissX['PH2O'])
    dataEmissX['SO2'] = roadX.iloc[:,13]*conver/smm[smm.iloc[:,0]=='SO2'].MM
    dataEmissX['VOC_INV'] = roadX.iloc[:,6]*conver + roadX.iloc[:,14]*conver +\
        roadX.iloc[:,19]*conver + roadX.iloc[:,20]*conver + roadX.iloc[:,21]*conver
    return dataEmissX