# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------

                                roadEmiss_v1.py


This function disaggregate the emissions from BRAVES originally in county level 
to grid level. It uses the roadLength_UF....csv file created by roadDensity.py
script from each Brazilian state. BRAVES outputs are in grams per year.

INPUTS:

    outPath = path to folder with roadDensity..csv and outputs from this function 
    
    bravesPath = path to folder with BRAVES outputs
    
    years = years of BRAVES's estimates

    IBGE_CODES = State ids 

OUTPUTS:

    roadEmiss_BySource_....csv = outputs in tons per year.


Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br

Atualização em: 05/05/2021 - implementada a emissão por categoria e tipologia
 
-------------------------------------------------------------------------------    
"""
# Importing packages
import geopandas as gpd
import pandas as pd
from shapely import wkt
import shapely.speedups
shapely.speedups.enable()

#%%
def zeroEmiss (emiss,typeEmiss):
    if typeEmiss =='Exhaust':
        emiss['Refuel_EmissNMHC']=0
        emiss['break_tier_EmissMP10']=0
        emiss['wear_EmissMP10']=0
        emiss['Resus_EmissMP10']=0
        emiss['Evap_D_EmissNMHC']=0
        emiss['Evap_H_EmissNMHC']=0
        emiss['Evap_R_EmissNMHC']=0      
    if typeEmiss =='non-exaust': 
        emiss['Exh_EmissCO']=0
        emiss['Exh_EmissHC']=0
        emiss['Exh_EmissNMHC']=0
        emiss['Exh_EmissCH4']=0
        emiss['Exh_EmissNOx']=0
        emiss['Exh_EmissRCHO']=0
        emiss['Exh_EmissMP2.5']=0
        emiss['Exh_EmissCO2']=0
        emiss['Exh_EmissN2O']=0
        emiss['Exh_EmissSO2']=0
        emiss['Exh_EmissCO2eq']=0
    if typeEmiss =='non-exaustMP': 
        emiss['Exh_EmissCO']=0
        emiss['Exh_EmissHC']=0
        emiss['Exh_EmissNMHC']=0
        emiss['Exh_EmissCH4']=0
        emiss['Exh_EmissNOx']=0
        emiss['Exh_EmissRCHO']=0
        emiss['Exh_EmissMP2.5']=0
        emiss['Exh_EmissCO2']=0
        emiss['Exh_EmissN2O']=0
        emiss['Exh_EmissSO2']=0
        emiss['Refuel_EmissNMHC']=0
        emiss['Exh_EmissCO2eq']=0
        emiss['Evap_D_EmissNMHC']=0
        emiss['Evap_H_EmissNMHC']=0
        emiss['Evap_R_EmissNMHC']=0
    if typeEmiss =='non-exaustMP_no_resusp': 
        #print ('typeEmiss = non-exaustMP_no_resusp')
        emiss['Exh_EmissCO']=0
        emiss['Exh_EmissHC']=0
        emiss['Exh_EmissNMHC']=0
        emiss['Exh_EmissCH4']=0
        emiss['Exh_EmissNOx']=0
        emiss['Exh_EmissRCHO']=0
        emiss['Exh_EmissMP2.5']=0
        emiss['Exh_EmissCO2']=0
        emiss['Exh_EmissN2O']=0
        emiss['Exh_EmissSO2']=0
        emiss['Refuel_EmissNMHC']=0
        emiss['Exh_EmissCO2eq']=0
        emiss['Evap_D_EmissNMHC']=0
        emiss['Evap_H_EmissNMHC']=0
        emiss['Evap_R_EmissNMHC']=0
        emiss['Resus_EmissMP10']=0 
    if typeEmiss =='TOTAL': 
        emiss =  emiss.copy()
    return emiss

def roadEmissBySource(df_emissYearState,df,name,factor,outPath,year,IBGE_CODE):
    roadE = gpd.GeoDataFrame({'geometry':df['geometry']})
    roadE['index'] = df.index.values 
    sumCity = df.sum()
    for jj in range(0, df_emissYearState.shape[1]): 
        emiss = pd.DataFrame()
        res = [i for i, val in enumerate(df.columns=='geometry') if val]
        for ii in range(res[0]+1,df.shape[1]):
            dfec= df_emissYearState[(int(df.columns[ii]) == df_emissYearState.index)]       
            emiss[str(df.columns[ii])] = (df.iloc[:,ii]/sumCity[(str(df.columns[ii]) == sumCity.index)].to_numpy())*dfec.iloc[0,jj]/factor
        sumPOL = emiss.sum(axis=1)
        roadE[df_emissYearState.columns[jj]] = sumPOL                 
    roadE.to_csv(outPath+'/'+ name +'_' + str(year)+ '_UF_'+str(IBGE_CODE)+'.csv')
    return roadE

#%%
def roadEmiss(outPath,bravesPath,years,IBGE_CODES,roadDensPrefix,typeEmiss):
    print('===================STARTING roadEmiss_v1.py=======================')
    factor = 1 # 
    for year in years:
        for IBGE_CODE in IBGE_CODES:
            # Opening road density file
            file = outPath+'/roadDensity_'+roadDensPrefix+'UF_'+str(IBGE_CODE)+'.csv'
            df = pd.read_csv(file)
            #df.set_index(df['index'])
            df = df.drop(df.columns[0], axis=1)
            df.columns = df.columns.str.replace(' ', '')
            df['geometry'] = df['geometry'].apply(wkt.loads)
            roadE = gpd.GeoDataFrame({'geometry':df['geometry']})
            roadE['index'] = df.index.values
        
            # OPEN EMISSION DATA
            df_emiss2 = pd.read_csv(bravesPath+"/EmissCityBRAVES_LightDuty_"+str(int(year))+".csv",
                                    skipinitialspace=True, delimiter=';')          
            df_emiss2.columns = df_emiss2.columns.str.replace(' ', '')
            
            # Selecting emission type
            df_emiss2 = zeroEmiss (df_emiss2,typeEmiss)   
           # Selecting year 
            df_emissYear2 = df_emiss2[df_emiss2.iloc[:,0]==year]
            # Selecting state - IBGE code
            df_emissYearState2 = df_emissYear2[df_emissYear2.iloc[:,1]==IBGE_CODE]
            # Renaming the IBGE code column name 
            df_emissYearState2= df_emissYearState2.rename(columns = {df_emissYearState2.iloc[:,2].name:'IBGEcode'})
            # Set IBGEcode as index 
            df_emissYearState2 = df_emissYearState2.set_index(df_emissYearState2.iloc[:,2].name)
            # Converting index from object to integer
            df_emissYearState2.index.astype(int)
            df_emissYearState2=df_emissYearState2.sort_index(axis = 0) 
            df_emissYearState2 = df_emissYearState2.replace("NaN ", 0)
            
            # ================Commercial-ligh vehicles
            # Reading emissions data 
            df_emiss = pd.read_csv(bravesPath+"/EmissCityBRAVES_CommercialLight_"+str(int(year))+".csv",
                                   skipinitialspace=True, delimiter=';')
            df_emiss.columns = df_emiss.columns.str.replace(' ', '')
            # Selecting emission type
            df_emiss = zeroEmiss (df_emiss,typeEmiss)
            # Selecting year 
            df_emissYear = df_emiss[df_emiss.iloc[:,0]==year]
            # Selecting state - IBGE code
            df_emissYearState1 = df_emissYear[df_emissYear.iloc[:,1]==IBGE_CODE]
            # Renaming the IBGE code column name 
            df_emissYearState1= df_emissYearState1.rename(columns = {df_emissYearState1.iloc[:,2].name:'IBGEcode'})
            # Set IBGEcode as index 
            df_emissYearState1 = df_emissYearState1.set_index(df_emissYearState1.iloc[:,2].name)
            # Converting index from object to integer
            df_emissYearState1.index.astype(int)
            df_emissYearState1=df_emissYearState1.sort_index(axis = 0) 
            df_emissYearState1 = df_emissYearState1.replace("NaN ", 0)
            
            # ================motorcicles vehicles
            # Reading emissions data 
            df_emiss3 = pd.read_csv(bravesPath+"/EmissCityBRAVES_Motorcycle_"+str(int(year))+".csv",
                                    skipinitialspace=True, delimiter=';')
            df_emiss3.columns = df_emiss3.columns.str.replace(' ', '')
            # Selecting emission type
            df_emiss3 = zeroEmiss (df_emiss3,typeEmiss)
            # Selecting year 
            df_emissYear3 = df_emiss3[df_emiss3.iloc[:,0]==year]
            # Selecting state - IBGE code
            df_emissYearState3 = df_emissYear3[df_emissYear3.iloc[:,1]==IBGE_CODE]
            # Renaming the IBGE code column name 
            df_emissYearState3= df_emissYearState3.rename(columns = {df_emissYearState3.iloc[:,2].name:'IBGEcode'})
            # Set IBGEcode as index 
            df_emissYearState3 = df_emissYearState3.set_index(df_emissYearState3.iloc[:,2].name)
            # Converting index from object to integer
            df_emissYearState3.index.astype(int)
            df_emissYearState3.iloc[:,6]= df_emissYearState3.iloc[:,6].replace('        NaN ',0)
            df_emissYearState3=df_emissYearState3.sort_index(axis = 0) 
            df_emissYearState3 = df_emissYearState3.replace("NaN ", 0)
            
            # ================heavy vehicles
            # Reading emissions data 
            df_emiss4 = pd.read_csv(bravesPath+"/EmissCityBRAVES_HeavyDuty_"+str(int(year))+".csv",
                                    skipinitialspace=True, delimiter=';')
            df_emiss4.columns = df_emiss4.columns.str.replace(' ', '')
            # Selecting emission type
            df_emiss4 = zeroEmiss (df_emiss4,typeEmiss)
            # Selecting year 
            df_emissYear4 = df_emiss4[df_emiss4.iloc[:,0]==year]
            # Selecting state - IBGE code
            df_emissYearState4 = df_emissYear4[df_emissYear4.iloc[:,1]==IBGE_CODE]
            # Renaming the IBGE code column name 
            df_emissYearState4= df_emissYearState4.rename(columns = {df_emissYearState4.iloc[:,2].name:'IBGEcode'})
            # Set IBGEcode as index 
            df_emissYearState4 = df_emissYearState4.set_index(df_emissYearState4.iloc[:,2].name)
            df_emissYearState4.iloc[:,6]= df_emissYearState4.iloc[:,6].replace('        NaN ',0)
            df_emissYearState4.iloc[:,8]= df_emissYearState4.iloc[:,8].replace('        NaN',0)
            df_emissYearState4=df_emissYearState4.sort_index(axis = 0)
            df_emissYearState4 = df_emissYearState4.replace(" " , "")
            df_emissYearState4 = df_emissYearState4.replace("NaN ", 0)
            

        
            roadEmissBySource(df_emissYearState1,df,
                              'roadEmiss_BySource_ComLight_'+typeEmiss+'_'+roadDensPrefix,factor,outPath,year,IBGE_CODE)
            roadEmissBySource(df_emissYearState2,df,
                              'roadEmiss_BySource_Light_'+typeEmiss+'_'+roadDensPrefix,factor,outPath,year,IBGE_CODE)
            roadEmissBySource(df_emissYearState4,df,
                              'roadEmiss_BySource_Heavy_'+typeEmiss+'_'+roadDensPrefix,factor,outPath,year,IBGE_CODE)
            roadEmissBySource(df_emissYearState3,df,
                              'roadEmiss_BySource_Motorcycle_'+typeEmiss+'_'+roadDensPrefix,factor,outPath,year,IBGE_CODE)
    

    return roadE       


