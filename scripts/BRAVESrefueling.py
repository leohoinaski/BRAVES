#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:29:31 2024

https://repositorio.ufpe.br/bitstream/123456789/6883/1/arquivo8070_1.pdf
# A mistura do etanol passa de 22% para 27%, podendo até chegar até 35%
# the RVP of gasolines ranges from 7 to 13
# observou-se que uma mistura com 10% de etanol aumenta a RVP
aproximadamente em 0,8 psi (5,61 kPa), e a fase vapor contém aproximadamente
15% de etanol. 

RVP = https://d35t1syewk4d42.cloudfront.net/file/1410/RVP-Effects-Memo_03_26_12_Final.pdf
# Extrai os dados de RVP do gráfico deste artigo

@author: leohoinaski
"""
import numpy as np
import matplotlib.pyplot as plt 
import netCDF4 as nc
from datetime import datetime
import pandas as pd
import pyproj
from shapely.geometry import Point
import geopandas as gpd
from ismember import ismember
import os

rootFolder = os.path.dirname(os.getcwd())
inputFolder = rootFolder +'/inputs/'
metcrod2dPath = '/media/leohoinaski/HDD/BR_2019/METCRO2D_BR_2019.nc'
shapeFolder = inputFolder+'shapefiles/BR_Municipios_2022.shp'



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

ds = nc.Dataset(metcrod2dPath)
tamb = ds['TEMP2'][:]
atribute = 'CD_MUN'
cities = gpd.read_file(shapeFolder)
cities.crs = "EPSG:4326"
datesTime = datePrepCMAQ(ds)
xv,yv,lon,lat = ioapiCoords(ds)
xlon,ylat = eqmerc2latlon(ds,xv,yv)
s,cityMat = citiesINdomain(xlon,ylat,cities)



#tamb = (45 - 5)*np.random.rand(24*30) + 5
# tConv = (tamb*(9/5)) + 32 # Celsius Farenheit
tConv = (tamb - 273.15)*((9/5)) + 32

rvp = 9 # RVP = Reid Vapor Pressure (psi)

#Td = Dispensed gasoline temperature (degF) = 20.30+0.81*Tamb  
#from MOVES https://www.epa.gov/sites/default/files/2020-11/documents/420r20012.pdf
td = 20.30+0.81*tConv 

#dT = 0.418*TDF -16.6 # From MOVES -
#https://www.epa.gov/sites/default/files/2020-11/documents/420r20012.pdf
deltaT = 0.418*td -16.6   
er = 264.2*(-5.909 - 0.0949*deltaT + 0.084*td + 0.485*rvp)

for ii,IBGE_CODE in enumerate(cities[atribute]):
    print(IBGE_CODE)
    if ii == 0:
        cityData,cityDataPoints,cityDataFrame,matData = dataINcity(
            er,datesTime,cityMat,s,int(IBGE_CODE))
        matData = np.nan_to_num(matData)
    else:
        cityData,cityDataPoints,cityDataFrame,matDataNew = dataINcity(
            er,datesTime,cityMat,s,int(IBGE_CODE))
        matDataNew = np.nan_to_num(matDataNew)
        matData = matData + matDataNew
        


plt.plot(er)
