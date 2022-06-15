#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 18:57:19 2021

Este script cria a tabela de especiação química do BRAVES com base nas espécies
do CMAQ (CMAQ_species.csv) e espécies coletadas no SPECIATE tool 5.1 da USEPA
 

@author: leohoinaski
"""

# Importando bibliotecas

import pandas as pd
import numpy as np

folderSpec = '/media/leohoinaski/HDD/website_leo/static/inventory/speciation/'

file_path = "CMAQ_species.csv"
dfCMAQspc = pd.read_csv(folderSpec+file_path)

file_path = "SPEC_names.csv"
dfspec = pd.read_csv(folderSpec+file_path,index_col=0)

cmaqSpecies = pd.DataFrame()

#%% ------------------------- ROAD EMISSIONS ----------------------------------
# Paved roads ressuspension 
file_path = "95780.csv"
dfPvdRoad = pd.read_csv(folderSpec+file_path,index_col=0)
dfPvdRoadJoin = dfPvdRoad.join(dfspec)
dfCMAQspc2 = dfCMAQspc.set_index('Formula').join(dfPvdRoadJoin.set_index('Molecular Formula'))
dfCMAQspc2 = dfCMAQspc2.groupby(dfCMAQspc2.index).mean()
dfCMAQspc2 = pd.DataFrame(dfCMAQspc2.iloc[:,0])
dfCMAQspc2 = dfCMAQspc2.rename(columns={dfCMAQspc2.iloc[:,0].name: '95780'})
cmaqSpecies['road'] = dfCMAQspc2['95780']

# Brakes wear
file_path = "95457.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0)
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['95457']=dfCMAQspcX
cmaqSpecies['brakes'] = dfCMAQspc2['95457']


# Tire wear
file_path = "3156.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0)
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['3156']=dfCMAQspcX
cmaqSpecies['tires'] = dfCMAQspc2['3156']

#%% -------------------------LIGHT EMISSIONS ----------------------------------

# --------------------------Light PM exhaust-----------------------------------
file_path = "8993VBS.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0)
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['8993VBS']=dfCMAQspcX
cmaqSpecies['LightPM'] = dfCMAQspc2['8993VBS']

# ------------------------ Light GAS exhaust ----------------------------------
# Média entre os 3 perfis
# Gasoline Exhaust - E10 gasoline - US06 Composite - 75 oF
file_path = "95792.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0)
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['95792']=dfCMAQspcX*1.472151 # TOG to VOC

# Onroad gasoline vehicle cold-start with VBS
file_path = "100VBS.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0)
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['100VBS']=dfCMAQspcX*1.245646

# Gasoline Exhaust - E10 gasoline, winter grade, LA92 cycle composite
file_path = "8908.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0)
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['8908']=dfCMAQspcX*1.109437
cmaqSpecies['LightVOC'] = np.nanmean([dfCMAQspc2['8908'],dfCMAQspc2['100VBS'],
                                      dfCMAQspc2['95792']], axis=0)



#%% ----------------------------HEAVY EMISSIONS -------------------------------

# PM Exhaust -----------------------------------------------------------------
#Diesel Exhaust - Heavy-heavy duty truck - 2007 model year with NCOM, VBS
file_path = "8996VBS.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['8996VBS']=dfCMAQspcX

#Conventional Diesel Exhaust - Idle Cycle, VBS
file_path = "8994VBS.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['8994VBS']=dfCMAQspcX

#Conventional Diesel Exhaust - Transient Cycle, VBS
file_path = "8995VBS.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['8995VBS']=dfCMAQspcX
cmaqSpecies['HeavyPM'] = np.nanmean([dfCMAQspc2['8996VBS'],dfCMAQspc2['8994VBS'],
                                      dfCMAQspc2['8995VBS']], axis=0)

# GAS Exhaust -----------------------------------------------------------------
#Diesel Exhaust Emissions from 2007 Model Year Heavy-Duty Diesel Engines with Controls
file_path = "8775.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['8775']=dfCMAQspcX*2.161389 # TOG to VOC

#Diesel Exhaust - Heavy-heavy duty truck - 2011 model year corrected
file_path = "95335a.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['95335a']=dfCMAQspcX

#Diesel Exhaust - Medium Duty Trucks
file_path = "4674.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['4674']=dfCMAQspcX

#Heavy Duty diesel with DPF, combination of previous measurements with VBS
file_path = "103VBS.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['103VBS']=dfCMAQspcX*2.245104 # TOG to VOC

cmaqSpecies['HeavyVOC'] = np.nanmean([dfCMAQspc2['8775'],dfCMAQspc2['95335a'],
                                      dfCMAQspc2['4674'],dfCMAQspc2['103VBS']], 
                                     axis=0)

#%%-------------------------- EVAPORATIVES ------------------------------------
#Light-Duty Gasoline Vehicles - Evaporative
file_path = "1204.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['1204']=dfCMAQspcX

#Vehicle - Current Fleet (1989) Hot Soak Evaporative
file_path = "2495.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['2495']=dfCMAQspcX

#Vehicle Hot Soak - Atlanta, 1990
file_path = "2567.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['2567']=dfCMAQspcX
	
#Vehicle - Current Fleet (1989) Diurnal Evaporative
file_path = "2493.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['2493']=dfCMAQspcX	
	
cmaqSpecies['LightEvap'] = np.nanmean([dfCMAQspc2['1204'],dfCMAQspc2['2495'],
                                      dfCMAQspc2['2567'],dfCMAQspc2['2493']], 
                                     axis=0)	


# HEAVY VEHICLES 
#Diesel Headspace Vapor Composite
file_path = "DIESEVP.csv"
dfP = pd.read_csv(folderSpec+file_path,index_col=0) 
dfPJoin = dfP.join(dfspec)
dfCMAQspcX = dfCMAQspc.set_index('Formula').join(dfPJoin.set_index('Molecular Formula'))
dfCMAQspcX = dfCMAQspcX.groupby(dfCMAQspcX.index).mean()
dfCMAQspcX = pd.DataFrame(dfCMAQspcX.iloc[:,0])
dfCMAQspc2['DIESEVP']=dfCMAQspcX*1.003551
cmaqSpecies['HeavyEvap'] = dfCMAQspc2['DIESEVP']


#%% Adding description

cmaqSpecies = dfCMAQspc.set_index('Formula').join(cmaqSpecies,on = 'Formula')

cmaqSpecies=cmaqSpecies.fillna(0)

cmaqSpecies.iloc[:,2:]=cmaqSpecies.iloc[:,2:]/100 # converting from percentage to factor

cmaqSpecies.to_csv(folderSpec+'/BRAVES_speciation.csv')
