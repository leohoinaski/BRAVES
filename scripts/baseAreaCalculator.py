#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  6 20:45:13 2022

@author: leohoinaski
"""
# Importando bibliotecas
import numpy as np
import geopandas as gpd
from shapely import wkt
import pandas as pd


#from BRAVES_temporalDisag import temporalDisagVehicular


def baseAreaCalulator(outPath,roadDensPrefix):
    basefile = outPath + '/baseGrid_'+roadDensPrefix+'.csv'
    base =pd.read_csv(basefile)
    gpdData = gpd.GeoDataFrame(base) 
    gpdData.crs="EPSG:4326"
    gpdData['geometry'] = gpdData['geometry'].apply(wkt.loads) 
    geod = gpdData.crs.get_geod()
    area=[]
    for re in gpdData.geometry:
        a = abs(geod.geometry_area_perimeter(re)[0])
        area.append(a)
    area = np.array(area)/(1000*1000) # area in km2
    return area

#area = baseAreaCalulator(outPath,roadDensPrefix)
