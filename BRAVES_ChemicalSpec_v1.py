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
import numpy as np

def ChemicalSpeciationLight(roadX,dfSpc,smm,conver):
    roadX.columns = roadX.columns.str.replace(' ', '')
    # Creating dataEmiss matrix
    dataEmissX=pd.DataFrame()
    print('===================STARTING BRAVES_ChemicalSpec_v1.py=======================')
    # Filling VOC emissions
    for ii in range(0,dfSpc.shape[0]):
        dataEmissX[dfSpc.ID[ii]] = dfSpc['LightVOC'][ii]*roadX['Exh_EmissNMHC']*conver/smm.iloc[ii,1] +\
           dfSpc['LightPM'][ii]*roadX['Exh_EmissMP2.5']*conver/smm.iloc[ii,1]  +\
               dfSpc['road'][ii]*roadX['wear_EmissMP10']*conver/smm.iloc[ii,1] +\
                   dfSpc['road'][ii]*roadX['Resus_EmissMP10']*conver/smm.iloc[ii,1] +\
                       dfSpc['brakes'][ii]*roadX['break_tier_EmissMP10']*conver/smm.iloc[ii,1] +\
                           dfSpc['tires'][ii]*roadX['break_tier_EmissMP10']*conver/smm.iloc[ii,1] +\
                               dfSpc['LightEvap'][ii]*roadX['Refuel_EmissNMHC']*conver/smm.iloc[ii,1] +\
                                   dfSpc['LightEvap'][ii]*roadX['Evap_D_EmissNMHC']*conver/smm.iloc[ii,1] +\
                                       dfSpc['LightEvap'][ii]*roadX['Evap_H_EmissNMHC']*conver/smm.iloc[ii,1] +\
                                           dfSpc['LightEvap'][ii]*roadX['Evap_R_EmissNMHC']*conver/smm.iloc[ii,1]
    
    # Filling direct emissions
    dataEmissX['ALDX'] = roadX['Exh_EmissRCHO']*conver/np.array(smm[smm.iloc[:,0]=='ALDX'].MM)
    dataEmissX['ALD2_PRIMARY'] = roadX['Exh_EmissRCHO']*conver/np.array(smm[smm.iloc[:,0]=='ALD2_PRIMARY'].MM)
    dataEmissX['CH4_INV'] = roadX['Exh_EmissCH4']*conver
    dataEmissX['CH4'] = roadX['Exh_EmissCH4']*conver/np.array(smm[smm.iloc[:,0]=='CH4'].MM)
    dataEmissX['CO'] = roadX['Exh_EmissCO']*conver/np.array(smm[smm.iloc[:,0]=='CO'].MM)
    dataEmissX['CO2_INV'] = roadX['Exh_EmissCO2']*conver
    dataEmissX['N2O_INV'] = roadX['Exh_EmissN2O']*conver
    dataEmissX['NO'] = roadX['Exh_EmissNOx']*conver*0.495/np.array(smm[smm.iloc[:,0]=='NO'].MM) # Para motor a diesel... considerei
    dataEmissX['NO2'] = roadX['Exh_EmissNOx']*conver*0.505/np.array(smm[smm.iloc[:,0]=='NO2'].MM) # Para motor a diesel... considerei
    dataEmissX['PMC'] = roadX['break_tier_EmissMP10']*conver+\
        roadX['Exh_EmissMP2.5']*conver+\
            roadX['Resus_EmissMP10']*conver+\
                roadX['wear_EmissMP10']*conver
            
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
    dataEmissX['SO2'] = roadX['Exh_EmissSO2']*conver/np.array(smm[smm.iloc[:,0]=='SO2'].MM)
    dataEmissX['VOC_INV'] = roadX['Exh_EmissNMHC']*conver + roadX['Evap_D_EmissNMHC']*conver +\
        roadX['Refuel_EmissNMHC']*conver + roadX['Evap_H_EmissNMHC']*conver + roadX['Evap_R_EmissNMHC']*conver
    dataEmissX['PMFINE'] = roadX['Exh_EmissMP2.5']*conver+ \
        0.4*roadX['break_tier_EmissMP10']*conver + \
            0.17*roadX['Resus_EmissMP10']*conver + \
                0.53*roadX['wear_EmissMP10']*conver
        
        
    return dataEmissX

def ChemicalSpeciationHeavy(roadX,dfSpc,smm,conver):  
    roadX.columns = roadX.columns.str.replace(' ', '')
    # Creating dataEmiss matrix
    dataEmissX=pd.DataFrame()
    # Filling VOC emissions
    for ii in range(0,dfSpc.shape[0]):
        dataEmissX[dfSpc.ID[ii]] = dfSpc['HeavyVOC'][ii]*roadX['Exh_EmissNMHC']*conver/smm.iloc[ii,1] +\
           dfSpc['HeavyPM'][ii]*roadX['Exh_EmissMP2.5']*conver/smm.iloc[ii,1]  +\
               dfSpc['road'][ii]*roadX['wear_EmissMP10']*conver/smm.iloc[ii,1] +\
                   dfSpc['road'][ii]*roadX['Resus_EmissMP10']*conver/smm.iloc[ii,1] +\
                       dfSpc['brakes'][ii]*roadX['break_tier_EmissMP10']*conver/smm.iloc[ii,1] +\
                           dfSpc['tires'][ii]*roadX['break_tier_EmissMP10']*conver/smm.iloc[ii,1] +\
                               dfSpc['HeavyEvap'][ii]*roadX['Refuel_EmissNMHC']*conver/smm.iloc[ii,1] +\
                                   dfSpc['HeavyEvap'][ii]*roadX['Evap_D_EmissNMHC']*conver/smm.iloc[ii,1] +\
                                       dfSpc['HeavyEvap'][ii]*roadX['Evap_H_EmissNMHC']*conver/smm.iloc[ii,1] +\
                                           dfSpc['HeavyEvap'][ii]*roadX['Evap_R_EmissNMHC']*conver/smm.iloc[ii,1]
    
    # Filling direct emissions
    dataEmissX['ALDX'] = roadX['Exh_EmissRCHO']*conver/np.array(smm[smm.iloc[:,0]=='ALDX'].MM)   
    dataEmissX['ALD2_PRIMARY'] = roadX['Exh_EmissRCHO']*conver/np.array(smm[smm.iloc[:,0]=='ALD2_PRIMARY'].MM)
    dataEmissX['CH4_INV'] = roadX['Exh_EmissCH4']*conver
    dataEmissX['CH4'] = roadX['Exh_EmissCH4']*conver/np.array(smm[smm.iloc[:,0]=='CH4'].MM)
    dataEmissX['CO'] = roadX['Exh_EmissCO']*conver/np.array(smm[smm.iloc[:,0]=='CO'].MM)
    dataEmissX['CO2_INV'] = roadX['Exh_EmissCO2']*conver
    dataEmissX['N2O_INV'] = roadX['Exh_EmissN2O']*conver
    dataEmissX['NO'] = roadX['Exh_EmissNOx']*conver*0.495/np.array(smm[smm.iloc[:,0]=='NO'].MM) # Para motor a diesel... considerei
    dataEmissX['NO2'] = roadX['Exh_EmissNOx']*conver*0.505/np.array(smm[smm.iloc[:,0]=='NO2'].MM) # Para motor a diesel... considerei
    dataEmissX['PMC'] = roadX['break_tier_EmissMP10']*conver+\
        roadX['Exh_EmissMP2.5']*conver+\
            roadX['Resus_EmissMP10']*conver+\
                roadX['wear_EmissMP10']*conver
            
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
    dataEmissX['SO2'] = roadX['Exh_EmissSO2']*conver/np.array(smm[smm.iloc[:,0]=='SO2'].MM)
    dataEmissX['VOC_INV'] = roadX['Exh_EmissNMHC']*conver + roadX['Evap_D_EmissNMHC']*conver +\
        roadX['Refuel_EmissNMHC']*conver + roadX['Evap_H_EmissNMHC']*conver + roadX['Evap_R_EmissNMHC']*conver
    dataEmissX['PMFINE'] = roadX['Exh_EmissMP2.5']*conver+ \
        0.4*roadX['break_tier_EmissMP10']*conver + \
            0.17*roadX['Resus_EmissMP10']*conver + \
                0.53*roadX['wear_EmissMP10']*conver
                
    return dataEmissX