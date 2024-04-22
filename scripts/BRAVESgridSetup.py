#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 15:00:54 2024

@author: leohoinaski
"""
from timezonefinder import TimezoneFinder
import netCDF4 as nc
from shapely.geometry import Polygon
import geopandas as gpd
import numpy as np
import pandas as pd
from datetime import datetime
import pyproj
from shapely.geometry import Point
from ismember import ismember

def baseGridDef(mcipGRIDDOT2DPath,outPath,GDNAM):
    print('Extracting MCIP coordinates')
    ds = nc.Dataset(mcipGRIDDOT2DPath)
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
    return baseGrid,ds

def datePrepCMAQ(ds):
    tf = np.array(ds['TFLAG'][:][:,1,:])
    date=[]
    for ii in range(0,tf.shape[0]):
        date.append(datetime.strptime(tf[:,0].astype(str)[ii] + (tf[:,1]/10000).astype(int).astype(str)[ii], '%Y%j%H').strftime('%Y-%m-%d %H:00:00'))
    date = np.array(date,dtype='datetime64[s]')
    dates = pd.DatetimeIndex(date)
    datesTime=pd.DataFrame()
    datesTime['year'] = dates.year
    datesTime['month'] = dates.month
    datesTime['day'] = dates.day
    datesTime['hour'] = dates.hour
    datesTime['datetime']=dates
    return datesTime

def ioapiCoords(ds):
    # Latlon
    lonI = ds.XORIG
    latI = ds.YORIG
    # Cell spacing 
    xcell = ds.XCELL
    ycell = ds.YCELL
    ncols = ds.NCOLS
    nrows = ds.NROWS
    lon = np.arange(lonI,(lonI+ncols*xcell),xcell)
    lat = np.arange(latI,(latI+nrows*ycell),ycell)
    xv, yv = np.meshgrid(lon,lat)
    return xv,yv,lon,lat

def eqmerc2latlon(ds,xv,yv):
    mapstr = '+proj=merc +a=%s +b=%s +lat_ts=0 +lon_0=%s' % (
              6370000, 6370000, ds.XCENT)
    #p = pyproj.Proj("+proj=merc +lon_0="+str(ds.P_GAM)+" +k=1 +x_0=0 +y_0=0 +a=6370000 +b=6370000 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    p = pyproj.Proj(mapstr)
    xlon, ylat = p(xv, yv, inverse=True)
    return xlon,ylat

def citiesINdomain(xlon,ylat,cities,atribute):
    s = gpd.GeoSeries(map(Point, zip(xlon.flatten(), ylat.flatten())))
    s = gpd.GeoDataFrame(geometry=s)
    s.crs = "EPSG:4326"
    s.to_crs("EPSG:4326")
    pointIn = cities.geometry.clip(s).explode()
    pointIn = gpd.GeoDataFrame({'geometry':pointIn}).reset_index()
    lia, loc = ismember(np.array((s.geometry.x,s.geometry.y)).transpose(),
                        np.array((pointIn.geometry.x,pointIn.geometry.y)).transpose(),'rows')
    s['city']=np.nan
    s.iloc[lia,1]=cities[atribute][pointIn['level_0'][loc]].values
    cityMat = np.reshape(np.array(s.city),(xlon.shape[0],xlon.shape[1])).astype(float)
    return s,cityMat

def dataINcity(aveData,datesTime,cityMat,s,IBGE_CODE):
    #IBGE_CODE=4202404
    cityData = aveData[:,:,cityMat==IBGE_CODE]
    cityDataPoints = s[s.city.astype(float)==IBGE_CODE]
    cityData = cityData[:,0,:]
    matData = aveData.copy()
    matData[:,:,cityMat!=IBGE_CODE]=np.nan
    cityDataFrame=pd.DataFrame(cityData)
    cityDataFrame.columns = cityDataPoints.geometry.astype(str)
    cityDataFrame['Datetime']=datesTime.datetime
    cityDataFrame = cityDataFrame.set_index(['Datetime'])
    return cityData,cityDataPoints,cityDataFrame,matData

def gridEssentials(mcipGRIDDOT2DPath,interFolder,GDNAM,munShp,atribute):
    baseGrid, ds = baseGridDef(mcipGRIDDOT2DPath,interFolder,GDNAM)
    datesTime = datePrepCMAQ(ds)
    xv,yv,lon,lat = ioapiCoords(ds)
    xlon,ylat = eqmerc2latlon(ds,xv,yv)
    s,cityMat = citiesINdomain(xlon,ylat,munShp,atribute)
    return baseGrid, ds, datesTime, xlon,ylat,  s,cityMat


def cityData2pixel(df,cityMat,IBGE_CODES,atribute,cc):
    atr = np.unique(df[atribute])
    disagData = np.zeros([len(atr),len(cc),cityMat.shape[0],cityMat.shape[1]])
    IBGE_CODES = IBGE_CODES[~np.isnan(IBGE_CODES)]
    for jj,at in enumerate(atr):
        fuelData = df[df[atribute]==at]
        for ii, IBGE_CODE in enumerate(IBGE_CODES):
            print(IBGE_CODE)
            cityData = cityMat.copy()
            cityData[cityMat==IBGE_CODE] = 1
            cityData[cityMat!=IBGE_CODE] = 0
            if fuelData[cc][(fuelData[atribute]==at) & (fuelData['IBGE_CODE']==int(IBGE_CODE))].shape[0]>0:
                b = np.repeat(np.array(fuelData[cc][(fuelData[atribute]==at) & (fuelData['IBGE_CODE']==int(IBGE_CODE))])[0],cityMat.shape[0]*cityMat.shape[1]).reshape([len(cc),cityMat.shape[0],cityMat.shape[1]])
                b[np.isnan(b)] = 0
                disagData[jj,:,:,:] = disagData[jj,:,:,:] + cityData*b
            else:
                print('City without ' + at)
    return disagData


