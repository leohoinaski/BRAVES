#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 17:49:55 2022

@author: leohoinaski
"""
import os
import numpy as np
from datetime import datetime
import pandas as pd
import netCDF4 as nc
from numpy.lib.stride_tricks import sliding_window_view
import pyproj


def dailyAverage (datesTime,data):
    daily = datesTime.groupby(['year','month','day']).count()
    dailyData = np.empty((daily.shape[0],data.shape[1],data.shape[2],data.shape[3]))
    for day in range(0,daily.shape[0]):
        findArr = (datesTime['year'] == daily.index[day][0]) & \
            (datesTime['month'] == daily.index[day][1]) & \
                (datesTime['day'] == daily.index[day][2]) 
        dailyData[day,:,:,:] = data[findArr,:,:,:].mean(axis=0)   
    return dailyData

def monthlyAverage (datesTime,data):
    monthly = datesTime.groupby(['year','month']).count()
    monthlyData = np.empty((monthly.shape[0],data.shape[1],data.shape[2],data.shape[3]))
    for month in range(0,monthly.shape[0]):
        findArr = (datesTime['year'] == monthly.index[month][0]) & \
            (datesTime['month'] == monthly.index[month][1]) 
        monthlyData[month,:,:,:] = data[findArr,:,:,:].mean(axis=0)   

    return monthlyData

def yearlyAverage (datesTime,data):
    yearly = datesTime.groupby(['year']).count()
    yearlyData = np.empty((yearly.shape[0],data.shape[1],data.shape[2],data.shape[3]))
    for year in range(0,yearly.shape[0]):
        if yearly.shape[0]>1:
            findArr = (datesTime['year'] == yearly.index[year][0])
        else:
            findArr = (datesTime['year'] == yearly.index[year])
        yearlyData[year,:,:,:] = data[findArr,:,:,:].mean(axis=0)   
    return yearlyData

def movingAverage (datesTime,data,w):
    daily = datesTime.groupby(['year','month','day']).count()
    mvAveData = np.empty((daily.shape[0],data.shape[1],data.shape[2],data.shape[3]))
    for day in range(0,daily.shape[0]):
        findArr = (datesTime['year'] == daily.index[day][0]) & \
            (datesTime['month'] == daily.index[day][1]) & \
                (datesTime['day'] == daily.index[day][2]) 
        for ii in range(0,findArr.sum()):
            ddData = data[findArr,:,:,:]
            if w+ii<=findArr.sum():
                dataN = ddData[ii:w+ii,:,:,:].mean(axis=0)
                if ii==0:
                    movData=dataN
                else:
                    movData = np.max([movData,dataN],axis=0)                
        mvAveData[day,:,:,:] = movData
    return mvAveData
    
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

def exceedance(data,criteria):
    freqExcd = np.sum(data>criteria,axis=0)
    return freqExcd

def eqmerc2latlon(ds,xv,yv):
    p = pyproj.Proj("+proj=merc +lon_0="+str(ds.P_GAM)+" +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    xlon, ylat = p(xv, yv, inverse=True)
    return xlon,ylat

def trimBorders (data,xv,yv,left,right,top,bottom):
    if np.size(data.shape)==4:
        dataT = data[:,:,bottom:(data.shape[2]-top),left:(data.shape[3]-right)]
    if np.size(data.shape)==3:
        dataT = data[:,bottom:(data.shape[2]-top),left:(data.shape[3]-right)]   
    if np.size(data.shape)==2:
        dataT = data[bottom:(data.shape[2]-top),left:(data.shape[3]-right)]
    
    xvT = xv[bottom:(data.shape[2]-top),left:(data.shape[3]-right)]
    yvT = yv[bottom:(data.shape[2]-top),left:(data.shape[3]-right)]
    
    return dataT,xvT,yvT
              
def getTime(ds,data):
    dd = datePrepCMAQ(ds)
    idx2Remove = np.array(dd.drop_duplicates().index)
    data = data[idx2Remove]
    datesTime = dd.drop_duplicates().reset_index(drop=True)
    return datesTime,data
