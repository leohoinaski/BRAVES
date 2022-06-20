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
References:
    https://fapesp.br/eventos/2016/02/mc/Maria_de_Fatima.pdf
    http://nuance-lapat.iag.usp.br/
    https://www.sciencedirect.com/science/article/pii/S0016236114005535
    https://cfpub.epa.gov/si/si_public_file_download.cfm?p_download_id=527145&Lab=OTAQ
    https://www.epa.gov/sites/default/files/2020-10/documents/emission_factor_documentation_for_ap-42_section_13.2.1_paved_roads_.pdf
    https://cfpub.epa.gov/si/si_public_file_download.cfm?p_download_id=525701


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
        dataEmissX[dfSpc.ID[ii]] = dfSpc['LightVOC'][ii]*roadX['EXH_NMHC']*conver/smm.iloc[ii,1] +\
           dfSpc['LightPM'][ii]*roadX['EXH_PMFINE']*conver/smm.iloc[ii,1]  +\
               dfSpc['road'][ii]*roadX['RW_PMC']*conver/smm.iloc[ii,1] +\
                   dfSpc['road'][ii]*roadX['RDR_PMC']*conver/smm.iloc[ii,1] +\
                       dfSpc['brakes'][ii]*roadX['BT_PMC']*conver/smm.iloc[ii,1] +\
                           dfSpc['tires'][ii]*roadX['BT_PMC']*conver/smm.iloc[ii,1] +\
                               dfSpc['LightEvap'][ii]*roadX['RFUEL_NMHC']*conver/smm.iloc[ii,1] +\
                                   dfSpc['LightEvap'][ii]*roadX['EVD_NMHC']*conver/smm.iloc[ii,1] +\
                                       dfSpc['LightEvap'][ii]*roadX['EVH_NMHC']*conver/smm.iloc[ii,1] +\
                                           dfSpc['LightEvap'][ii]*roadX['EVR_NMHC']*conver/smm.iloc[ii,1]
    
    # Filling direct emissions
    ald = roadX['EXH_RCHO']
    dataEmissX['ALDX'] = ald*0.11*conver/np.array(smm[smm.iloc[:,0]=='ALDX'].MM)
    dataEmissX['ALD2'] =ald*0.50*conver/np.array(smm[smm.iloc[:,0]=='ALD2'].MM)
    dataEmissX['ALD2_PRIMARY'] = ald*0.50*conver/np.array(smm[smm.iloc[:,0]=='ALD2'].MM)
    dataEmissX['ACET'] = ald*0.08*conver/np.array(smm[smm.iloc[:,0]=='ACET'].MM)
    dataEmissX['FORM'] = ald*0.39*conver/np.array(smm[smm.iloc[:,0]=='FORM'].MM)
    dataEmissX['FORM_PRIMARY'] = ald*0.39*conver/np.array(smm[smm.iloc[:,0]=='FORM'].MM)
    
    dataEmissX['CH4_INV'] = roadX['EXH_CH4']*conver
    dataEmissX['CH4'] = roadX['EXH_CH4']*conver/np.array(smm[smm.iloc[:,0]=='CH4'].MM)
    dataEmissX['CO'] = roadX['EXH_CO']*conver/np.array(smm[smm.iloc[:,0]=='CO'].MM)
    dataEmissX['CO2_INV'] = roadX['EXH_CO2']*conver
    dataEmissX['N2O_INV'] = roadX['EXH_N2O']*conver
    dataEmissX['NO'] = roadX['EXH_NOX']*conver*0.495/np.array(smm[smm.iloc[:,0]=='NO'].MM) # Para motor a diesel... considerei
    dataEmissX['NO2'] = roadX['EXH_NOX']*conver*0.505/np.array(smm[smm.iloc[:,0]=='NO2'].MM) # Para motor a diesel... considerei
    dataEmissX['PMC'] = roadX['BT_PMC']*conver+\
        roadX['EXH_PMFINE']*conver+\
            roadX['RDR_PMC']*conver+\
                roadX['RW_PMC']*conver
            
    dataEmissX['PNCOM'] = dataEmissX['POC']*0.25
    dataEmissX['PMOTHR'] = dataEmissX['PMC']-(dataEmissX['POC'] + dataEmissX['PEC'] +
                              dataEmissX['PSO4'] + dataEmissX['PNO3'] +
                              dataEmissX['PNH4'] + dataEmissX['PNCOM'] +
                              dataEmissX['PFE'] + dataEmissX['PAL'] +
                              dataEmissX['PSI'] + dataEmissX['PTI'] +
                              dataEmissX['PCA'] + dataEmissX['PMG'] +
                              dataEmissX['PK'] + dataEmissX['PMN'] +
                              dataEmissX['PNA'] + dataEmissX['PCL'] +
                              dataEmissX['PH2O'])
    dataEmissX['SO2'] = roadX['EXH_SO2']*conver/np.array(smm[smm.iloc[:,0]=='SO2'].MM)
    dataEmissX['VOC_INV'] = roadX['EXH_NMHC']*conver + roadX['EVD_NMHC']*conver +\
        roadX['RFUEL_NMHC']*conver + roadX['EVH_NMHC']*conver + roadX['EVR_NMHC']*conver
    dataEmissX['PMFINE'] = roadX['EXH_PMFINE']*conver+ \
        0.4*roadX['BT_PMC']*conver + \
            0.17*roadX['RDR_PMC']*conver + \
                0.53*roadX['RW_PMC']*conver
     
    etoh = dataEmissX['ETOH']  - \
        dfSpc[dfSpc['ID']=='ETOH']['LightVOC']*\
            roadX['EXHFLEX_NMHC']*conver/np.array(smm[smm.iloc[:,0]=='SO2'].MM) 
    dataEmissX['ETOH'] = np.nansum([etoh,
                                        roadX['EXHFLEX_ETOH']*conver/np.array(smm[smm.iloc[:,0]=='ETOH'].MM)],
                                       axis=0)
        
    return dataEmissX

def ChemicalSpeciationHeavy(roadX,dfSpc,smm,conver):  
    roadX.columns = roadX.columns.str.replace(' ', '')
    # Creating dataEmiss matrix
    dataEmissX=pd.DataFrame()
    # Filling VOC emissions
    for ii in range(0,dfSpc.shape[0]):
        dataEmissX[dfSpc.ID[ii]] = dfSpc['HeavyVOC'][ii]*roadX['EXH_NMHC']*conver/smm.iloc[ii,1] +\
           dfSpc['HeavyPM'][ii]*roadX['EXH_PMFINE']*conver/smm.iloc[ii,1]  +\
               dfSpc['road'][ii]*roadX['RW_PMC']*conver/smm.iloc[ii,1] +\
                   dfSpc['road'][ii]*roadX['RDR_PMC']*conver/smm.iloc[ii,1] +\
                       dfSpc['brakes'][ii]*roadX['BT_PMC']*conver/smm.iloc[ii,1] +\
                           dfSpc['tires'][ii]*roadX['BT_PMC']*conver/smm.iloc[ii,1] +\
                               dfSpc['HeavyEvap'][ii]*roadX['RFUEL_NMHC']*conver/smm.iloc[ii,1] +\
                                   dfSpc['HeavyEvap'][ii]*roadX['EVD_NMHC']*conver/smm.iloc[ii,1] +\
                                       dfSpc['HeavyEvap'][ii]*roadX['EVH_NMHC']*conver/smm.iloc[ii,1] +\
                                           dfSpc['HeavyEvap'][ii]*roadX['EVR_NMHC']*conver/smm.iloc[ii,1]
    
    # Filling direct emissions
    # dataEmissX['ALDX'] = 
    # dataEmissX['ALD2'] = dataEmissX['ALDX']*0.22
    # dataEmissX['ALD2_PRIMARY'] = dataEmissX['ALDX']*0.22
    # dataEmissX['ACET'] = dataEmissX['ALDX']*0.09
    # dataEmissX['FORM'] = dataEmissX['ALDX']*0.69
    # dataEmissX['FORM_PRIMARY'] = dataEmissX['ALDX']*0.69
    
    dataEmissX['CH4_INV'] = roadX['EXH_CH4']*conver
    dataEmissX['CH4'] = roadX['EXH_CH4']*conver/np.array(smm[smm.iloc[:,0]=='CH4'].MM)
    dataEmissX['CO'] = roadX['EXH_CO']*conver/np.array(smm[smm.iloc[:,0]=='CO'].MM)
    dataEmissX['CO2_INV'] = roadX['EXH_CO2']*conver
    dataEmissX['N2O_INV'] = roadX['EXH_N2O']*conver
    dataEmissX['NO'] = roadX['EXH_NOX']*conver*0.495/np.array(smm[smm.iloc[:,0]=='NO'].MM) # Para motor a diesel... considerei
    dataEmissX['NO2'] = roadX['EXH_NOX']*conver*0.505/np.array(smm[smm.iloc[:,0]=='NO2'].MM) # Para motor a diesel... considerei
    dataEmissX['PMC'] = roadX['BT_PMC']*conver+\
        roadX['EXH_PMFINE']*conver+\
            roadX['RDR_PMC']*conver+\
                roadX['RW_PMC']*conver
            
    dataEmissX['PNCOM'] = dataEmissX['POC']*0.25
    dataEmissX['PMOTHR'] = dataEmissX['PMC']-(dataEmissX['POC'] + dataEmissX['PEC'] +
                              dataEmissX['PSO4'] + dataEmissX['PNO3'] +
                              dataEmissX['PNH4'] + dataEmissX['PNCOM'] +
                              dataEmissX['PFE'] + dataEmissX['PAL'] +
                              dataEmissX['PSI'] + dataEmissX['PTI'] +
                              dataEmissX['PCA'] + dataEmissX['PMG'] +
                              dataEmissX['PK'] + dataEmissX['PMN'] +
                              dataEmissX['PNA'] + dataEmissX['PCL'] +
                              dataEmissX['PH2O'])
    dataEmissX['SO2'] = roadX['EXH_SO2']*conver/np.array(smm[smm.iloc[:,0]=='SO2'].MM)
    dataEmissX['VOC_INV'] = roadX['EXH_NMHC']*conver + roadX['EVD_NMHC']*conver +\
        roadX['RFUEL_NMHC']*conver + roadX['EVH_NMHC']*conver + roadX['EVR_NMHC']*conver
    dataEmissX['PMFINE'] = roadX['EXH_PMFINE']*conver+ \
        0.4*roadX['BT_PMC']*conver + \
            0.17*roadX['RDR_PMC']*conver + \
                0.53*roadX['RW_PMC']*conver
                
    return dataEmissX

def ChemicalSpeciationMotorcycle(roadX,dfSpc,smm,conver):
    # Creating dataEmiss matrix
    dataEmissX=pd.DataFrame()
    print('===================STARTING BRAVES_ChemicalSpec_v1.py=======================')
    # Filling VOC emissions
    for ii in range(0,dfSpc.shape[0]):
        dataEmissX[dfSpc.ID[ii]] = dfSpc['LightVOC'][ii]*roadX['EXH_NMHC']*conver/smm.iloc[ii,1] +\
           dfSpc['LightPM'][ii]*roadX['EXH_PMFINE']*conver/smm.iloc[ii,1]  +\
               dfSpc['road'][ii]*roadX['RW_PMC']*conver/smm.iloc[ii,1] +\
                   dfSpc['road'][ii]*roadX['RDR_PMC']*conver/smm.iloc[ii,1] +\
                       dfSpc['brakes'][ii]*roadX['BT_PMC']*conver/smm.iloc[ii,1] +\
                           dfSpc['tires'][ii]*roadX['BT_PMC']*conver/smm.iloc[ii,1] +\
                               dfSpc['LightEvap'][ii]*roadX['RFUEL_NMHC']*conver/smm.iloc[ii,1] +\
                                   dfSpc['LightEvap'][ii]*roadX['EVD_NMHC']*conver/smm.iloc[ii,1] +\
                                       dfSpc['LightEvap'][ii]*roadX['EVH_NMHC']*conver/smm.iloc[ii,1] +\
                                           dfSpc['LightEvap'][ii]*roadX['EVR_NMHC']*conver/smm.iloc[ii,1]

    
    # Filling direct emissions
    # dataEmissX['ALDX'] = roadX.iloc[:,9]*conver/smm[smm.iloc[:,0]=='ALDX'].MM   
    dataEmissX['CH4_INV'] = roadX['EXH_CH4']*conver
    dataEmissX['CH4'] = roadX['EXH_CH4']*conver/np.array(smm[smm.iloc[:,0]=='CH4'].MM)
    dataEmissX['CO'] = roadX['EXH_CO']*conver/np.array(smm[smm.iloc[:,0]=='CO'].MM)
    dataEmissX['CO2_INV'] = roadX['EXH_CO2']*conver
    dataEmissX['N2O_INV'] = roadX['EXH_N2O']*conver
    dataEmissX['NO'] = roadX['EXH_NOX']*conver*0.495/np.array(smm[smm.iloc[:,0]=='NO'].MM) # Para motor a diesel... considerei
    dataEmissX['NO2'] = roadX['EXH_NOX']*conver*0.505/np.array(smm[smm.iloc[:,0]=='NO2'].MM) # Para motor a diesel... considerei
    dataEmissX['PMC'] = roadX['BT_PMC']*conver+\
        roadX['EXH_PMFINE']*conver+\
            roadX['RDR_PMC']*conver+\
                roadX['RW_PMC']*conver
            
    dataEmissX['PNCOM'] = dataEmissX['POC']*0.25
    dataEmissX['PMOTHR'] = dataEmissX['PMC']-(dataEmissX['POC'] + dataEmissX['PEC'] +
                              dataEmissX['PSO4'] + dataEmissX['PNO3'] +
                              dataEmissX['PNH4'] + dataEmissX['PNCOM'] +
                              dataEmissX['PFE'] + dataEmissX['PAL'] +
                              dataEmissX['PSI'] + dataEmissX['PTI'] +
                              dataEmissX['PCA'] + dataEmissX['PMG'] +
                              dataEmissX['PK'] + dataEmissX['PMN'] +
                              dataEmissX['PNA'] + dataEmissX['PCL'] +
                              dataEmissX['PH2O'])
    dataEmissX['SO2'] = roadX['EXH_SO2']*conver/np.array(smm[smm.iloc[:,0]=='SO2'].MM)
    dataEmissX['VOC_INV'] = roadX['EXH_NMHC']*conver + roadX['EVD_NMHC']*conver +\
        roadX['RFUEL_NMHC']*conver + roadX['EVH_NMHC']*conver + roadX['EVR_NMHC']*conver
    dataEmissX['PMFINE'] = roadX['EXH_PMFINE']*conver+ \
        0.4*roadX['BT_PMC']*conver + \
            0.17*roadX['RDR_PMC']*conver + \
                0.53*roadX['RW_PMC']*conver
    return dataEmissX