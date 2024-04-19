#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 15:00:54 2024

@author: leohoinaski
"""
from timezonefinder import TimezoneFinder
import netCDF4 as nc4
from shapely.geometry import Polygon
import geopandas as gpd
import numpy as np
import pandas as pd

def baseGrid(mcipGRIDDOT2DPath,outPath,GDNAM):
    print('Extracting MCIP coordinates')
    ds = nc4.Dataset(mcipGRIDDOT2DPath)
    #dataVar = list(ds.variables.keys())
    x = np.unique(ds['LOND'][:])
    y = np.unique(ds['LATD'][:])
    test_naive = pd.date_range('2019-01-01', '2019-04-07', freq='4H')
    tf = TimezoneFinder(in_memory=True)
    #Loop over each cel in x direction
    polygons=[]
    ltc =[]
    print('Creating grid and extraction TimeZones')
    for ii in range(1,x.shape[0]):
        #Loop over each cel in y direction
        for jj in range(1,y.shape[0]):
            #roadClip=[]
            lat_point_list = [y[jj-1], y[jj], y[jj], y[jj-1]]
            lon_point_list = [x[ii-1], x[ii-1], x[ii], x[ii]]
            cel = Polygon(zip(lon_point_list, lat_point_list))
            polygons.append(cel)
            local_time_zone = tf.timezone_at(lng=x[ii], lat=y[jj])
            try:
                ltc.append(float(test_naive.tz_localize(local_time_zone).strftime('%Z')[-1]))
            except:
                ltc.append(np.nan)
    baseGrid = gpd.GeoDataFrame({'geometry':polygons})
    baseGrid['LTZ'] = ltc
    baseGrid.to_csv(outPath+'/baseGrid_'+GDNAM+'.csv')
    return baseGrid