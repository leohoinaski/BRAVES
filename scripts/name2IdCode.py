#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 19:44:29 2024

@author: leohoinaski
"""
import pandas as pd
from unidecode import unidecode
import numpy as np
import geopandas as gpd
import difflib

ufs = {
  "UF": ['ACRE', 'ALAGOAS', 'AMAPA', 'AMAZONAS', 'BAHIA', 'CEARA',
       'DISTRITO FEDERAL', 'ESPIRITO SANTO', 'GOIAS', 'MARANHAO',
       'MATO GROSSO', 'MATO GROSSO DO SUL', 'MINAS GERAIS',
       'PARA', 'PARAIBA', 'PARANA', 'PERNAMBUCO',
       'PIAUI', 'RIO DE JANEIRO', 'RIO GRANDE DO NORTE',
       'RIO GRANDE DO SUL', 'RONDONIA', 'RORAIMA', 'SANTA CATARINA',
       'SAO PAULO', 'SERGIPE', 'TOCANTINS'],
  "SIGLA":['AC','AL','AP','AM','BA','CE','DF','ES','GO','MA','MT','MS','MG',
           'PA','PB','PR' ,'PE','PI','RJ','RN','RS','RO','RR','SC','SP','SE','TO'],
  "CODE": [12,27,16,13,29,23,53,32,52,21,51,50,31,15,25,41,26,22,33,24,43,11,14,
           42,35,28,17]
}


def name2Code(inputFolder,file,year,month):
    
    #inputFolder = '/media/leohoinaski/HDD/BRAVESv2/inputs/' 
    inputFolder ='/mnt/sdb1/BRAVESv2/inputs/'
    path = inputFolder+'fleet/fuelType/fuelType_2021_01.csv' 
    shapeFolder = inputFolder+'shapefiles/BRMUE250GC_SIR.shp'
    df = pd.read_csv(path)
    df['MUN2'] = np.nan
    munShp = gpd.read_file(shapeFolder)
    munShp['MUN2'] = np.nan
    munShp['UF'] = ''
    
    for ii, mun in enumerate(df['MUN2']):
        try:
            df['MUN2'][ii] = unidecode(df['MUN'][ii].upper())
        except:
            df['MUN2'][ii] = np.nan
    
    for ii, mun in enumerate(munShp['NM_MUNICIP']):
        try:
            munShp['MUN2'][ii] = unidecode(munShp['NM_MUNICIP'][ii].upper())
            munShp['UF'][ii] = np.array(ufs['UF'])[np.array(ufs['CODE']) == int(str(munShp['CD_GEOCMU'][ii])[0:2])]
                
        except:
            munShp['MUN2'][ii] = np.nan
        
    df['IBGE_CODE'] = np.nan
    uniqueUF = np.unique(munShp['UF'])
    for ii, uuf in enumerate(uniqueUF):
        munUF = df['MUN2'][df['UF']==uuf]
        munShapUF = munShp['MUN2'][munShp['UF']==uuf]
        for jj,mun in enumerate(munUF):
            valMun = difflib.get_close_matches(mun, munShapUF)[0]
            print(valMun)
            munOK = munShp['MUN2'][munShp['UF']==uuf].reset_index()
            dataF = munOK.loc[[np.where((munShp['MUN2'][munShp['UF']==uuf]==valMun))][0][0][0]]
            idx = munOK.loc[[np.where((munShp['MUN2'][munShp['UF']==uuf]==valMun))][0][0][0]][0]
            muni = munOK.loc[[np.where((munShp['MUN2'][munShp['UF']==uuf]==valMun))][0][0][0]][1]
            df['IBGE_CODE'][idx] = muni

