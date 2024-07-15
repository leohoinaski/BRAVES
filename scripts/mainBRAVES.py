#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 16:02:54 2024

@author: leohoinaski
"""

import name2IdCode as m2id
import scrappage as scrp
import os
import pandas as pd
import geopandas as gpd
import numpy as np 
import BRAVESgridSetup as bgs
import BRAVES_temporalDisag as btd


# ===============================Inputs========================================

# Estabelece a pasta root do BRAVES
rootFolder = os.path.dirname(os.getcwd()) 

 # Estabelce a pasta base 
baseFolder = os.path.dirname(rootFolder)

# Definindo pasta de inputs
inputFolder = rootFolder +'/inputs/'

# Definindo categorias veiculares para uso
vehicularCategories =['LIGHT','COMMERCIAL-LIGHT','MOTORCYCLES','HEAVY']

# Ano Base simulação
year = 2021

# Mes da simulação
month = 1

# Definindo pasta com shapefile de referência
shapeFolder = inputFolder+'shapefiles/BR_Municipios_2022.shp'

# Atributo com o nome dos municipios
atribute = 'CD_MUN'


#Definindo pasta com tabelas
tablePath = inputFolder+'/tables'

# Definindo nome da grade
GDNAM = 'BR_2019'

# Definindo caminho para os arquivos do mcip
mcipPath = baseFolder + '/BR_2019'
metcrod2dPath = mcipPath + '/METCRO2D_BR_2019.nc'
mcipGRIDDOT2DPath = mcipPath+'/GRIDDOT2D_'+GDNAM+'.nc'

#===============================================================================

# Abrindo shapefile de referencia
munShp = gpd.read_file(shapeFolder)

# ----Grid Setup
# Definindo pasta de arquivos intermediarios
interFolder = rootFolder + '/outputs/intermediate/gridFiles'
# Cria pasta dos arquivos intermediarios
os.makedirs(interFolder, exist_ok=True)
# Funcao gridEssential para obter 
baseGrid, datesTime, xlon,ylat, s, cityMat = bgs.gridEssentials(mcipGRIDDOT2DPath,interFolder,GDNAM,munShp,atribute)

# fuelType
interFolder = rootFolder + '/outputs/intermediate/fuelType'
filePath = inputFolder+'fleet/fuelType/fuelType_2021_01.csv' 
if  os.path.isfile(interFolder+'/BRAVES_' + filePath.split('/')[-1]) == False:
    dfFuel = m2id.name2Code(munShp,filePath,interFolder)
else:
    print('YOU ALREADY HAVE THE fuelType file')
    dfFuel = pd.read_csv(interFolder+'/BRAVES_' + filePath.split('/')[-1])

# vehicularCategory
filePath = inputFolder+'fleet/vehicularCategory/vehicularCategory_2021_01.csv'
interFolder = rootFolder + '/outputs/intermediate/vehicularCategory'
if  os.path.isfile(interFolder+'/BRAVES_' + filePath.split('/')[-1]) == False:
    dfVCat = m2id.name2Code(munShp,filePath,interFolder)
else:
    print('YOU ALREADY HAVE THE vehicularCategory file')
    dfVCat = pd.read_csv(interFolder+'/BRAVES_' + filePath.split('/')[-1])
    
# yearModel
filePath = inputFolder+'fleet/yearModel/yearModel_2021_01.csv'
interFolder = rootFolder + '/outputs/intermediate/yearModel'
if  os.path.isfile(interFolder+'/BRAVES_' + filePath.split('/')[-1]) == False:
    dfYM = m2id.name2Code(munShp,filePath,interFolder)
else:
    print('YOU ALREADY HAVE THE yearModel file')
    dfYM = pd.read_csv(interFolder+'/BRAVES_' + filePath.split('/')[-1])

# Scrappage
filePath = inputFolder+'fleet/yearModel/yearModel_2021_01.csv'
interFolder = rootFolder + '/outputs/intermediate/yearModel'
if  os.path.isfile(interFolder+'/BRAVES_scrappage_' + filePath.split('/')[-1]) == False:
    dfYM = scrp.main(munShp,interFolder,filePath,vehicularCategories,
                        dfVCat,dfYM)
else:
    print('YOU ALREADY HAVE THE BRAVES_scrappage file')
    dfYM = pd.read_csv(interFolder+'/BRAVES_scrappage_' + filePath.split('/')[-1])

# fuelType
interFolder = rootFolder + '/outputs/intermediate/fuelType'
filePath = inputFolder+'fuel_consumption/2021.csv' 
dfFuelCons = pd.read_csv(filePath).replace(' -   ','').replace('"','').replace(',','')   
cc = list(np.array(range(1,13),dtype=str))
dfFuelCons[list(np.array(range(1,13),dtype=str))] = dfFuelCons[list(np.array(range(1,13),dtype=str))].apply(pd.to_numeric) 
disagData = bgs.cityData2pixel(dfFuelCons,cityMat,np.array(munShp[atribute],dtype=int),'FUEL',list(np.array(range(1,13),dtype=str)))


# Desagregando consumo de combustível
dataTempo,xX,yY,disvec = btd.BRAVES_temporalDisagFuel(inputFolder,metcrod2dPath,baseGrid)