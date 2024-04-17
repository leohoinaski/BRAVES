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


rootFolder = os.path.dirname(os.getcwd())
inputFolder = rootFolder +'/inputs/'
vehicularCategories =['LIGHT','COMMERCIAL-LIGHT','MOTORCYCLES','HEAVY']
year = 2021
month = 1
shapeFolder = inputFolder+'shapefiles/BR_Municipios_2022.shp'

# fuelType
interFolder = rootFolder + '/outputs/intermediate/fuelType'
filePath = inputFolder+'fleet/fuelType/fuelType_2021_01.csv' 
if  os.path.isfile(interFolder+'/BRAVES_' + filePath.split('/')[-1]) == False:
    dfFuel = m2id.name2Code(shapeFolder,filePath,interFolder)
else:
    print('YOU ALREADY HAVE THE fuelType file')
    dfFuel = pd.read_csv(interFolder+'/BRAVES_' + filePath.split('/')[-1])

# vehicularCategory
filePath = inputFolder+'fleet/vehicularCategory/vehicularCategory_2021_01.csv'
interFolder = rootFolder + '/outputs/intermediate/vehicularCategory'
if  os.path.isfile(interFolder+'/BRAVES_' + filePath.split('/')[-1]) == False:
    dfVCat = m2id.name2Code(shapeFolder,filePath,interFolder)
else:
    print('YOU ALREADY HAVE THE vehicularCategory file')
    dfVCat = pd.read_csv(interFolder+'/BRAVES_' + filePath.split('/')[-1])
    
# yearModel
filePath = inputFolder+'fleet/yearModel/yearModel_2021_01.csv'
interFolder = rootFolder + '/outputs/intermediate/yearModel'
if  os.path.isfile(interFolder+'/BRAVES_' + filePath.split('/')[-1]) == False:
    dfYM = m2id.name2Code(shapeFolder,filePath,interFolder)
else:
    print('YOU ALREADY HAVE THE vehicularCategory file')
    dfYM = pd.read_csv(interFolder+'/BRAVES_' + filePath.split('/')[-1])
    
# Scrappage
dfYM,munShp = scrp.main(shapeFolder,interFolder,filePath,vehicularCategories,
                        dfVCat,dfYM)