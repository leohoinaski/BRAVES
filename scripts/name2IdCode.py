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
import os
import ismember


ufs = {
  "UF": ['ACRE', 'ALAGOAS', 'AMAPA', 'AMAZONAS', 'BAHIA', 'CEARA',
       'DISTRITO FEDERAL', 'ESPIRITO SANTO', 'GOIAS', 'MARANHAO',
       'MATO GROSSO', 'MATO GROSSO DO SUL', 'MINAS GERAIS',
       'PARA', 'PARAIBA', 'PARANA', 'PERNAMBUCO',
       'PIAUI', 'RIO DE JANEIRO', 'RIO GRANDE DO NORTE',
       'RIO GRANDE DO SUL', 'RONDONIA', 'RORAIMA', 'SANTA CATARINA',
       'SAO PAULO', 'SERGIPE', 'TOCANTINS'],
  "SIGLA":     ['AC','AL','AP','AM','BA','CE',
                'DF','ES','GO','MA',
                'MT','MS','MG',
                'PA','PB','PR' ,'PE',
                'PI','RJ','RN',
                'RS','RO','RR','SC',
                'SP','SE','TO'],
  "IBGE_CODE": [12,27,16,13,29,23,
                53,32,52,21,
                51,50,31,
                15,25,41,26,
                22,33,24,
                43,11,14,42,
                35,28,17]
}



def name2Code(shapeFolder,filePath,interFolder):
    """
    Esta função encontra os códigos IBGE nas planilhas que não possuem este
    dado. A função procura as cidades dentro de cada estado e insere o código 
    IBGE com base em um shapefile de referência (que possui o código e o nome
    das cidades).

    Parameters
    ----------
    shapeFolder : path
        Caminho para o arquivo do shapefile de referencia.
    filePath : path
        Caminho para o arquivo que será inserido o código IBGE

    Returns
    -------
    df : Pandas DataFrame 
        DataFrame com caracteres corrigido e com as colunas de UF e Código IBGE.

    """
    
    df = pd.read_csv(filePath)
    df.columns = df.columns.str.replace(' ', '')
    df['MUN2'] = np.nan
    df['UF_SIGLA'] = ''
    df = df.reset_index(drop=True)
    munShp = gpd.read_file(shapeFolder)
    munShp['MUN2'] = ''
    munShp['UF'] = ''
    munShp = munShp.reset_index(drop=True)
    
    # Loop para arrumar os caracteres do df
    for ii, mun in enumerate(df['MUN2']):
        try:
            df['MUN2'][ii] = unidecode(df['MUN'][ii].upper())
            if len(str(df['UF'][ii]))==2:
                df['UF_SIGLA'] = df['UF'].copy()
                df['UF'][ii] = np.array(ufs['UF'])[np.array(ufs['SIGLA']) == (str(df['UF'][ii]))]
            else:
                df['UF_SIGLA'][ii] = np.array(ufs['SIGLA'])[np.array(ufs['UF']) == (str(df['UF'][ii]))]

        except:
            df['MUN2'][ii] = np.nan
            
    # Loop para arrumar os caracteres do shapefile
    for ii, mun in munShp.iterrows():
        try:
            munShp['MUN2'][ii] = unidecode(munShp['NM_MUN'][ii].upper())
            munShp['UF'][ii] = np.array(ufs['UF'])[np.array(ufs['IBGE_CODE']) == int(str(munShp['CD_MUN'][ii])[0:2])]
        except:
            print('skiping')
            munShp['MUN2'][ii] = np.nan
            
    # Encontrando os códigos com base no arquivo de shapefile    
    df['IBGE_CODE'] = np.nan
    uniqueUF = np.unique(munShp['SIGLA_UF'])
    uniqueUF=['RN']
    # Loop em cada estado
    for ii, uuf in enumerate(uniqueUF):
        munUF = df[df['UF_SIGLA']==str(uuf)]
        munShapUF = munShp['MUN2'][munShp['SIGLA_UF']==str(uuf)]
        print('============'+str(uuf)+'================')
        cityWcode=[]
        cityNot=[]
        munUF['IBGE_CODE']=np.nan
        print(munUF)
        for jj,mun in enumerate(munUF['MUN2']):
            try:
                valMun = difflib.get_close_matches(mun, munShapUF)[0]
                #print(valMun)
                cityWcode.append(munShapUF[munShapUF ==valMun].index[0])
                #print(munShapUF[munShapUF ==valMun].index[0])
            except:
                try:
                    valMun = difflib.get_close_matches(mun.replace('SS','C'), munShapUF)[0]
                    cityWcode.append(munShapUF[munShapUF ==valMun].index[0])
                except:
                    valMun=[]  
                    cityNot.append(mun)
                    cityWcode.append(-999)
        cityWcodeTrue = np.array(cityWcode[np.array(cityWcode)!=-999])
        munUF['IBGE_CODE'][np.array(cityWcode)!=-999]=\
            np.array(munShp['CD_MUN'][cityWcodeTrue]).astype(float)
        df['IBGE_CODE'][df['UF_SIGLA']==str(uuf)] = munUF['IBGE_CODE']
        
    
    # Exportando o dataframe para a pasta do input = subrescreve
    os.makedirs(interFolder, exist_ok=True)
    df.to_csv(interFolder +'/BRAVES_' + filePath.split('/')[-1])
    
    return df

            

