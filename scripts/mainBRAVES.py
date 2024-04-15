#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 16:02:54 2024

@author: leohoinaski
"""

import name2IdCode as m2id
import os

rootFolder = os.path.dirname(os.getcwd())
inputFolder = rootFolder +'/inputs/'


# fuelType
shapeFolder = inputFolder+'shapefiles/BR_Municipios_2022.shp'
filePath = inputFolder+'fleet/fuelType/fuelType_2021_01.csv' 
interFolder = rootFolder + '/outputs/intermediate/fuelType'
dfFuel = m2id.name2Code(shapeFolder,filePath,interFolder)

# vehicularCategory
filePath = inputFolder+'fleet/vehicularCategory/vehicularCategory_2021_01.csv'
interFolder = rootFolder + '/outputs/intermediate/vehicularCategory'
dfVCat = m2id.name2Code(shapeFolder,filePath,interFolder)

# yearModel
filePath = inputFolder+'fleet/yearModel/yearModel_2021_01.csv'
interFolder = rootFolder + '/outputs/intermediate/yearModel'
dfYM = m2id.name2Code(shapeFolder,filePath,interFolder)
