#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 19:44:29 2024

@author: leohoinaski
"""
import pandas as pd
from unidecode import unidecode
import numpy as np

def name2Code(inputFolder,file,year,month):
    
    inputFolder = '/media/leohoinaski/HDD/BRAVESv2/inputs/' 
    path = inputFolder+'fleet/fuelType/fuelType_2021_01.csv' 
    
    df = pd.read_csv(path)
    df['MUN2'] = np.nan
    
    for ii, mun in enumerate(df['MUN2']):
        try:
            df['MUN2'][ii] = unidecode(df['MUN'][ii])
        except:
            df['MUN2'][ii] = np.nan
            
    
    