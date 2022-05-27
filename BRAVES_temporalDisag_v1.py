#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                           BRAVES_temporalDisag_v1.py

This function provides temporal disaggregated files from the annual speciated and
spatial disaggregated files from BRAVES database.

Inputs:
    
    folder: folter to output files BRAVESdatabaseAnnual....nc
    
    name: output names
         
    month: month for disaggregation
    
    day: days for disaggregation

Outputs:
        
    netdCDF files
      

Last update = 30/09/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br
-------------------------------------------------------------------------------
"""
# Importando bibliotecas
import pandas as pd
import numpy as np
import datetime
import numpy.matlib
import netCDF4 as nc
#from local2UTC import local2UTC
#from netCDFcreator_v1 import createNETCDFtemporal



#%% temporal disagregation


def temporalDisagVehicular(dataNC,year,month,day,hourdis,weekdis,monthdis,ltz):
    # Create the MultiIndex from pollutant and time.
    #year=2013
    print('Temporal disagregation')
    startDate = datetime.datetime(year, month, day, 0, 0)
    endDate = datetime.datetime(year, month, day+1, 0, 0)
    datePfct = np.arange(np.datetime64(startDate),np.datetime64(endDate),3600000000)
    numWeeks = datePfct.shape[0]/(7*24) # Number of weeks
    disvec = pd.DataFrame()
    disvec = disvec.reindex(datePfct, fill_value=np.nan)
    disvec['year'] = disvec.index.year
    disvec['month'] = disvec.index.month
    disvec['day'] = disvec.index.day
    disvec['hour'] = disvec.index.hour
    disvec['weekday'] = disvec.index.weekday # Monday is 0 and Sunday is 6
    disvec['hourdis'] = numpy.matlib.repmat(
        hourdis, 1, int(disvec.shape[0]/24)).transpose() 
   
    disvec['weekdis']=np.zeros([disvec['weekday'].shape[0],1])
    for ii in range(0,disvec['weekday'].shape[0]):
        if disvec['weekday'][ii] == 6:
            disvec['weekdis'][ii] = weekdis[0] 
        if disvec['weekday'][ii] == 0:
            disvec['weekdis'][ii] = weekdis[1]
        if disvec['weekday'][ii] == 1:
            disvec['weekdis'][ii] = weekdis[2] 
        if disvec['weekday'][ii] == 2:
            disvec['weekdis'][ii] = weekdis[3] 
        if disvec['weekday'][ii] == 3:
            disvec['weekdis'][ii] = weekdis[4] 
        if disvec['weekday'][ii] == 4:
            disvec['weekdis'][ii] = weekdis[5] 
        if disvec['weekday'][ii] == 5:
            disvec['weekdis'][ii] = weekdis[6] 
    
    disvec['monthdis']=np.zeros([disvec['month'].shape[0],1])
    for ii in range(0,disvec['month'].shape[0]):
        if disvec['month'][ii] == 1:
            disvec['monthdis'][ii] = monthdis[0] 
        if disvec['month'][ii] == 2:
            disvec['monthdis'][ii] = monthdis[1] 
        if disvec['month'][ii] == 3:
            disvec['monthdis'][ii] = monthdis[2] 
        if disvec['month'][ii] == 4:
            disvec['monthdis'][ii] = monthdis[3] 
        if disvec['month'][ii] == 5:
            disvec['monthdis'][ii] = monthdis[4]     
        if disvec['month'][ii] == 6:
            disvec['monthdis'][ii] = monthdis[5] 
        if disvec['month'][ii] == 7:
            disvec['monthdis'][ii] = monthdis[6] 
        if disvec['month'][ii] == 8:
            disvec['monthdis'][ii] = monthdis[7] 
        if disvec['month'][ii] == 9:
            disvec['monthdis'][ii] = monthdis[8] 
        if disvec['month'][ii] == 10:
            disvec['monthdis'][ii] = monthdis[9]  
        if disvec['month'][ii] == 11:
            disvec['monthdis'][ii] = monthdis[10] 
        if disvec['month'][ii] == 12:
            disvec['monthdis'][ii] = monthdis[11]  
   
    disvec['prod']=disvec['hourdis']*disvec['weekdis']*disvec['monthdis']/numWeeks        
    # converting from hourly to second basis
    disvec['prod'] = disvec['prod']/3600 
    dataTempo = np.zeros([datePfct.shape[0],dataNC.shape[1],
                          dataNC.shape[2],dataNC.shape[3]])
    print(str(dataTempo.shape))
    
    for jj in range(0,dataTempo.shape[1]):
        for ii in range(0,dataTempo.shape[0]):
            utcoffs = numpy.unique(ltz)
            for utcoff in utcoffs:
                idx = ltz==utcoff
                dataTempo[ii,jj,idx]= dataNC[0,jj,idx]* np.roll(disvec['prod'],int(utcoff))[ii]
    
    # #Condition for equal timezone in whole domain - roll vector to UTC time
    # if tag==1:        
    #     for jj in range(0,dataTempo.shape[1]):
    #         for ii in range(0,dataTempo.shape[0]):
    #             utcoffs = numpy.unique(lc2utc)
    #             for utcoff in utcoffs:
    #                 idx = lc2utc==utcoff
    #                 dataTempo[ii,jj,idx.transpose()]= dataNC[0,jj,idx.transpose()]* np.roll(disvec['prod'],int(utcoff))[ii]
    #             #dataTempo[ii,jj,:,:]= dataNC[0,jj,:,:]* np.roll(disvec['prod'],lc2utc)[ii]
    # else:
    #     for jj in range(0,dataTempo.shape[1]):
    #         for ii in range(0,dataTempo.shape[0]):
    #             utcoffs = numpy.unique(lc2utc)
    #             for utcoff in utcoffs:
    #                 idx = lc2utc==utcoff
    #                 dataTempo[ii,jj,idx.transpose()]= dataNC[0,jj,idx.transpose()]* np.roll(disvec['prod'],int(utcoff))[ii]

    return dataTempo,datePfct,disvec


#%%
def BRAVES_temporalDisag(rootPath,outPath,file,month,day):
    print('===================STARTING BRAVES_temporalDisag_v1.py=======================')
    hourdis = list(pd.read_csv(rootPath+'/TemporalAloc/hourdis.csv').iloc[:,1])
    weekdis = list(pd.read_csv(rootPath+'/TemporalAloc/weekdis.csv').iloc[:,1])
    monthdis = list(pd.read_csv(rootPath+'/TemporalAloc/monthdis.csv').iloc[:,1])
    

    prefix = '_'.join(file.split('_')[1:6])  
    
    year = int(file.split('_')[-1][0:4])
    
    ds = nc.Dataset(outPath+'/'+file)
    dataVar = list(ds.variables.keys())
    dataNC = ds[dataVar[5]][:]
    for ii in range(6,np.size(dataVar)):
        
        dataloop = np.array(ds[dataVar[ii]][:])
        print(dataVar[ii])
        dataNC =np.concatenate((dataNC,dataloop), axis=1)
    
    # # Informações das grades 
    # xi = ds.getncattr('XORIG')
    # yi = ds.getncattr('YORIG')
    # dx = ds.getncattr('XCELL')
    # dy = ds.getncattr('YCELL')
    
    # center = {'x': [float(ds.XCENT)],'y': [float(ds.YCENT)] }
    # center = pd.DataFrame(center)

    # lat = np.linspace(yi, yi+(dy*(dataNC.shape[3])), dataNC.shape[3])
    # lon = np.linspace(xi, xi+(dx*(dataNC.shape[2])), dataNC.shape[2])
    # xv, yv = np.meshgrid(lon, lat)
    
    xX = ds['Longitute'][:]
    yY = ds['Latitude'][:]
    area = ds['AREA'][:]
    ltz = ds['LTZ'][:]
    
    # Calling local2UTC
    #lc2utc, tag = local2UTC(xX,yY)
    
    # Calling temporalDisagVehicular
    dataTempo,datePfct,disvec = temporalDisagVehicular(dataNC,year,month,day,
                                                       hourdis,weekdis,
                                                       monthdis,ltz)
    

    
    return dataTempo,xX,yY,disvec,prefix,area