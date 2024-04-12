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
shapeFolder = inputFolder+'shapefiles/BRMUE250GC_SIR.shp'
filePath = inputFolder+'fleet/fuelType/fuelType_2021_01.csv' 
dfFuel = m2id.name2Code(shapeFolder, filePath)

# vehicularCategory
filePath = inputFolder+'fleet/vehicularCategory/vehicularCategory_2021_01.csv'
dfVCat = m2id.name2Code(shapeFolder, filePath)

# yearModel
filePath = inputFolder+'fleet/yearModel/yearModel_2021_01.csv'
dfYM = m2id.name2Code(shapeFolder, filePath)
