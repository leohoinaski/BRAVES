#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 09:30:01 2022

@author: leohoinaski
"""

import os
#import numpy as np
#import matplotlib.pyplot as plt
#import geopandas as gpd
#from datetime import datetime
import netCDF4 as nc
#import pandas as pd
import temporalStatistics as tst
import GAr_figs as garfig
import matplotlib

#%%
NO2 = {
  "Pollutant": "Nitrogen dioxide ($NO_{2}$)",
  "Criteria": 0.097414,
  "Unit": 'ppm',
  "Criteria_annual": 0.021266,
  "Criteria_average": '1-h average',
  "tag":'NO2',
  "Criteria_ave": 1,
}

CO = {
  "Pollutant": "Carbon monoxide (CO)",
  "Criteria": 9,
  "Unit": 'ppm',
  "Criteria_average": '8-h moving average',
  "tag":'CO',
  "Criteria_ave": 8,
}

O3 = {
  "Pollutant": "Ozone ($O_{3}$)",
  "Criteria":0.050961,
  "Unit": 'ppm',
  "Criteria_average": '8-h moving average',
  "Criteria_ave": 8,
  "tag":'O3'
}

SO2 = {
  "Pollutant": "Sulphur dioxide ($SO_{2}$)",
  "Criteria": 0.007636,
  "Unit": 'ppm',
  "Criteria_annual": 0.007636,
  "Criteria_average": '24-h average',
  "Criteria_ave": 24,
  "tag":'SO2'
}

PM10 = {
  "Pollutant": "Particulate Matter ($PM_{10}$)",
  "Criteria": 50,
  "Unit": '$\u03BCg m_{3}$',
  "Criteria_annual": 20,
  "Criteria_average": '24-h average',
  "tag":'PM10',
  "Criteria_ave": 24,
}

PM25 = {
  "Pollutant": "Particulate Matter ($PM_{2.5}$)",
  "Criteria": 25,
  "Unit": '$\u03BCg m_{3}$',
  "Criteria_annual": 10,
  "Criteria_average": '24-h average',
  "tag":'PM25',
  "Criteria_ave": 24,
}

#%% INPUTS
path = '/media/leohoinaski/Backup'

borderShape = '/media/leohoinaski/Backup/HospDisaggregation/Inputs/shapefiles/Brasil.shp'

fileType='CCTM_CONC'

pollutants = [O3,NO2]

# Trim domain
left = 40
right = 20
top=95
bottom=20




#%% Opening data 

# Moving to dir
os.chdir(path)

# Selecting files and variables
prefixed = sorted([filename for filename in os.listdir(path) if filename.startswith(fileType)])

# Opening netCDF files
ds = nc.MFDataset(prefixed)

for pol in pollutants:
    # Selecting variable
    data = ds[pol['tag']][:]
    
    # Get datesTime and removing duplicates
    datesTime, data = tst.getTime(ds,data)
    
    # Get coordinates from ioapi 
    xv,yv,lon,lat = tst.ioapiCoords(ds)
    
    # Trim borders left/right/bottom/top
    dataT,xvT,yvT= tst.trimBorders(data,xv,yv,left,right,top,bottom)
    
   
    if pol['Criteria_ave']==1:
        aveData = dataT.copy()
    elif pol['Criteria_ave']==8:
        # Daily-maximum 8h-moving average
        aveData = tst.movingAverage(datesTime,dataT,8)
    elif pol['Criteria_ave']==24:
        # Daily averages
        aveData = tst.dailyAverage(datesTime,dataT)
           
    # Monthly averages
    monthlyData = tst.monthlyAverage(datesTime,dataT)
    
    # Frequency of violations
    freqExcd= tst.exceedance(aveData,pol['Criteria'])
    
    # Transforming mercator to latlon/degrees
    xlon, ylat = tst.eqmerc2latlon(ds,xvT,yvT)
    
    # Figures
    # Average
    legend = pol['Pollutant'] +' '+ pol['Criteria_average'] +'\n'+ pol['Unit']
    #cmap = 'YlOrRd'
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightgray","yellow","orange","red"])
    garfig.timeAverageFig(aveData.max(axis=0)[0,:,:],xlon,ylat,legend,cmap,borderShape)
    
    # Exceedence
    legend = pol['Pollutant'] +' '+ pol['Criteria_average'] +'\n'+ 'Number of violations'
    #cmap = 'RdPu'
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["yellow","red","purple"])
    garfig.exceedanceFig(freqExcd[0,:,:],xlon,ylat,legend,cmap,borderShape)
    
    # Criteria
    legend = pol['Pollutant'] +' '+ pol['Criteria_average'] +'\n'+ pol['Unit']
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightgray","yellow","orange","brown"])
    #cmap = 'YlOrBr'     
    garfig.criteriaFig(aveData.max(axis=0)[0,:,:],xlon,ylat,legend,cmap,borderShape,pol['Criteria'])
    
   # Yearly averages
    if "Criteria_annual" in pol:
        yearlyData = tst.yearlyAverage(datesTime,dataT)
        # Frequency of violations
        freqExcdY= tst.exceedance(yearlyData,pol['Criteria_annual'])
        legend = pol['Pollutant'] +' '+ 'Annual average' +'\n'+ pol['Unit']
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightgray","yellow","orange","red"])        
        garfig.timeAverageFig(yearlyData.max(axis=0)[0,:,:],xlon,ylat,legend,cmap,borderShape)
        # Exceedence
        legend = pol['Pollutant'] +' '+ 'Annual average' +'\n'+ 'Number of violations'
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["yellow","red","purple"])
        garfig.exceedanceFig(freqExcdY[0,:,:],xlon,ylat,legend,cmap,borderShape)
        legend = pol['Pollutant'] +' '+ 'Annual average' +'\n'+ pol['Unit']
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightgray","yellow","orange","brown"])    
        garfig.criteriaFig(yearlyData.max(axis=0)[0,:,:],xlon,ylat,legend,cmap,borderShape,pol['Criteria'])
