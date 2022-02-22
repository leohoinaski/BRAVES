#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 18:16:43 2022

@author: leohoinaski
"""

"""
---------------------------------netCDF_ploting--------------------------------

This function plots netCDF file with atmospheric emissions from BRAVES.


INPUTS:
    
    netcdfPath = path to netCDF file
    outputPath = path to outputs
        
    
OUTPUTS:
    
    figures 
    
    

Author: Leonardo Hoinaski,....

Last update: 22/10/2021

"""

# Importing packages

import netCDF4 as nc
from matplotlib import pyplot as plt 
import numpy as np
from mpl_toolkits.basemap import Basemap
import geopandas as gpd
import matplotlib.colors as colors


# ---------------------------------PROCESSING----------------------------------
def BRAVES_ploting(BRAVESPath,outputPath,borderShape,pol):
    # Open netCDF files
    ncData = nc.Dataset(netcdfPath)
    ncData.variables
    # Show netCDF dimensions
    ncData.dimensions
    
    # Latlon
    lonI = ncData.XORIG
    latI = ncData.YORIG
    
    # Cell spacing 
    xcell = ncData.XCELL
    ycell = ncData.YCELL
    ncols = ncData.NCOLS
    nrows = ncData.NROWS
    
    xlon = np.linspace(lonI,(lonI+ncols*xcell),ncols)
    ylat = np.linspace(latI,(latI+nrows*ycell),nrows)
    
    xv, yv = np.meshgrid(xlon, ylat)
    
    
    # Slicing 
    pec = ncData[pol][:]
    
    pec2 = pec[0,0,:,:] 
    fig, ax = plt.subplots()
    m = Basemap(projection='tmerc', 
               lat_0=ylat.mean(),lon_0=xlon.mean(),
               llcrnrlon=xlon[0], 
               llcrnrlat=ylat[0], 
               urcrnrlon=xlon[-1], 
               urcrnrlat=ylat[-1],
               epsg=4326)
    #m.drawcountries()
    #m.drawcoastlines()
    xi, yi = m(xv, yv)
    ax.axis('off')
        
    cs = m.pcolor(xi, yi, np.squeeze(pec), cmap='hot_r', ax=ax,
                  norm=colors.LogNorm(vmin=0.1, vmax=pec.max()))
    
    br = gpd.read_file(borderShape)
    br.boundary.to_crs(epsg=4326).plot(ax=ax,edgecolor='black',linewidth=0.1)
    # Add Colorbar
    cbar = m.colorbar(cs, location='bottom', pad="-5%",shrink=0.1)
    cbar.set_label('Emission')
    #fig.savefig(outputPath+'ncOut_'+pol+'.png')
    plt.show()
    return pec2
    
#%% Running function 
pol ='PAL'

# Inputs    
netcdfPath = '/media/leohoinaski/HDD/BRAVES_database/Outputs/BRAVESdatabaseAnnual_BR_Total_0.05x0.05_2019.nc'
outputPath = '/media/leohoinaski/HDD/BRAVES_database/Outputs'
borderShape = '/media/leohoinaski/HDD/BRAVES_database/Shapefiles/Brasil.shp'
pec2=BRAVES_ploting(netcdfPath,outputPath,borderShape,pol)