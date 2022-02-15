#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 15:34:23 2022

@author: leohoinaski
"""
#from roadDensity_v1 import roadDensity

import os
import numpy as np
import geopandas as gpd


# Path to functions and BRAVESdatabase_main.py
rootPath = os.path.abspath(os.getcwd())

# Defining folder path with openstreetmaps folder
dirPath = rootPath + "/Shapefiles"
roadFileName = 'gis_osm_roads_free_1.shp'

IBGE_CODES = [11,12,13,14,15,16,17,
              21,22,23,24,25,26,27,28,29,
              31,32,33,35,
              41,42,43,
              50,51,52,53] # include the IBGE code from the states to be considered

# Opening cities shapefile
citySHP = gpd.read_file(dirPath+"/citySHP.shp")
citySHP.columns = ['NM_MUNICIP','CD_GEOCMU','YEAR',
                     'UF','geometry']

# Opening state shapefile
brSHP = gpd.read_file(dirPath+"/Brasil.shp")

for pp in range (0,np.size(IBGE_CODES)):
    print('State number = '+str(IBGE_CODES[pp]))
    ufSHP = brSHP[brSHP['COD_UF']==IBGE_CODES[pp]]
    #ufSHP = brSHP[brSHP['COD_UF']==21]
    ufSHP = ufSHP.dissolve(by='UF')
    ufSHP.crs = "EPSG:4326"
    ufSHP = ufSHP.buffer(0.2, resolution=10)
    print('Cliping roads into State number ' + str(pp) + ' - ' + brSHP[brSHP['COD_UF']==IBGE_CODES[pp]]['UF'].to_numpy()[0])
    #Cliping roads inside city
    print('Selecting states')
    print('Opening roads shapefile')
    if IBGE_CODES[pp]==11 or IBGE_CODES[pp]==12 or IBGE_CODES[pp]==13 or \
        IBGE_CODES[pp]==14 or IBGE_CODES[pp]==15 or IBGE_CODES[pp]==16 or \
            IBGE_CODES[pp]==17:
        roads = gpd.read_file(dirPath+"/norte-latest-free.shp/gis_osm_roads_free_1.shp")
        folder = dirPath+'/norte-latest-free.shp'
    elif IBGE_CODES[pp]==21 or IBGE_CODES[pp]==22 or IBGE_CODES[pp]==23 or \
        IBGE_CODES[pp]==24 or IBGE_CODES[pp]==25 or IBGE_CODES[pp]==26 or \
            IBGE_CODES[pp]==27 or IBGE_CODES[pp]==28 or IBGE_CODES[pp]==29:
        roads = gpd.read_file(dirPath+"/nordeste-latest-free.shp/gis_osm_roads_free_1.shp")
        folder = dirPath+'/nordeste-latest-free.shp'
    elif IBGE_CODES[pp]==31 or IBGE_CODES[pp]==32 or IBGE_CODES[pp]==33 or \
        IBGE_CODES[pp]==34 or IBGE_CODES[pp]==35:
        roads = gpd.read_file(dirPath+"/sudeste-latest-free.shp/gis_osm_roads_free_1.shp")
        folder = dirPath+'/sudeste-latest-free.shp'
    elif IBGE_CODES[pp]==42:
        roads = gpd.read_file(dirPath+"/sul-latest-free.shp/"+roadFileName)
        folder = dirPath+'/sul-latest-free.shp'
    elif IBGE_CODES[pp]==41 or IBGE_CODES[pp]==43:
        roads = gpd.read_file(dirPath+"/sul-latest-free.shp/gis_osm_roads_free_1.shp")
        folder = dirPath+'/sul-latest-free.shp'
    elif IBGE_CODES[pp]==50 or IBGE_CODES[pp]==51 or IBGE_CODES[pp]==52 or \
        IBGE_CODES[pp]==53:
        roads = gpd.read_file(dirPath+"/centro-oeste-latest-free.shp/gis_osm_roads_free_1.shp")
        folder = dirPath+'/centro-oeste-latest-free.shp'                 
    
    highways = roads[((roads['fclass']=='primary') | (roads['fclass']=='secundary') |
          (roads['fclass']=='trunk') | 
          ((roads['maxspeed']>80) & (roads['maxspeed']<160)) | 
          (roads['ref'].isnull()==False))]
    
    
    try:
        roadsUF = gpd.clip(roads['geometry'], ufSHP)
        highwaysUF = gpd.clip(highways['geometry'], ufSHP)
    except:
        roadsUF = gpd.clip(roads, ufSHP)
        highwaysUF = gpd.clip(highways, ufSHP)
        
        # SELECIONANDO VIAS PRIMÁRIAS, SECUNDÁRIAS, ACIMA DE 80 km/h

        
        # Opening cities shapefile
    citySHP = gpd.read_file(dirPath+"/citySHP.shp")
    citySHP.columns = ['NM_MUNICIP','CD_GEOCMU','YEAR',
                         'UF','geometry']
    
    shpSC = citySHP[citySHP.iloc[:,3]==IBGE_CODES[pp]]
    scALL = shpSC.dissolve(by='UF')
    scALL.crs = "EPSG:4326"
    scBuffer = scALL.buffer(0.2, resolution=10)
    scBuffer.crs = "EPSG:4326"
    scBuffer.reset_index(drop=True, inplace=True)
    boundSC = scBuffer.bounds
    valL=[]
    
    
    print('Selecting cells into state')
    for kk in range(0,shpSC.shape[0]):       
        roadCity =[]
        cc=1
        latitudes = []
        longitudes = []
        mlsTest =[]
        print('City number = ' + str(kk) + ' of ' + str(shpSC.shape[0]) +
              '  -'+ shpSC.iloc[kk,0])
        #Cliping roads inside city
        try:
            roadCity = gpd.clip(roadsUF['geometry'], shpSC.iloc[kk,4])
            highwaysCity = gpd.clip(highwaysUF['geometry'], shpSC.iloc[kk,4])
        except:
            roadCity = gpd.clip(roadsUF, shpSC.iloc[kk,4])
            highwaysCity = gpd.clip(highwaysUF, shpSC.iloc[kk,4])
            
        roadCity.to_csv(folder+'/roadCity_'+shpSC['CD_GEOCMU'].values[0]+'.csv')
        highwaysCity.to_csv(folder+'/highwaysCity_'+shpSC['CD_GEOCMU'].values[0]+'.csv')
        
    roadsUF.to_csv(folder+'/roadsUF_'+str(int(IBGE_CODES[pp]))+'.csv')
    highwaysUF.to_csv(folder+'/highwaysUF_'+str(int(IBGE_CODES[pp]))+'.csv')

        
        