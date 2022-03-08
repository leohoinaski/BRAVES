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
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import string

# ---------------------------------PROCESSING----------------------------------
def BRAVES_ploting(netcdfPath,netcdfPath2,netcdfPath3,outputPath,borderShape,pol,colormap):
    

    
    ncData0 = nc.Dataset(netcdfPath+str(2013)+'.nc')
    pec0 = ncData0[pol][:]/(10**6) 
    
    for count, year in enumerate([2013]):
        
        fig, ax = plt.subplots(1,3,gridspec_kw = {'wspace':0, 'hspace':0})
        cm = 1/2.54  # centimeters in inches
        fig.set_size_inches(19*cm, 6*cm)
        
        # Open netCDF files
        ncData = nc.Dataset(netcdfPath+str(year)+'.nc')
        ncData2 = nc.Dataset(netcdfPath2+str(year)+'.nc')
        ncData3 = nc.Dataset(netcdfPath3+str(year)+'.nc')
        
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
        pec = ncData[pol][:]/(10**6)    
    
        pec2 = ncData2[pol][:]/(10**6)   
    
        pec3 = ncData3[pol][:]/(10**6)  

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
        
        ax[0].axis('off')
        ax[1].axis('off')
        ax[2].axis('off')

            
        cs = m.pcolor(xi, yi, np.squeeze(pec)-np.squeeze(pec2), cmap=colormap, ax=ax[0],
                      norm=colors.LogNorm(vmin=pec0[pec0>0].min(), vmax=pec0.max()))
        
        cs2 = m.pcolor(xi, yi, np.squeeze(pec2)-np.squeeze(pec3), cmap=colormap, ax=ax[1],
                      norm=colors.LogNorm(vmin=pec0[pec0>0].min(), vmax=pec0.max()))
        
        cs3 = m.pcolor(xi, yi, np.squeeze(pec3), cmap=colormap, ax=ax[2],
                      norm=colors.LogNorm(vmin=pec0[pec0>0].min(), vmax=pec0.max()))
        
        br = gpd.read_file(borderShape)
        br.boundary.to_crs(epsg=4326).plot(ax=ax[0],edgecolor='black',linewidth=0.1)
        br.boundary.to_crs(epsg=4326).plot(ax=ax[1],edgecolor='black',linewidth=0.1)
        br.boundary.to_crs(epsg=4326).plot(ax=ax[2],edgecolor='black',linewidth=0.1)
        
        # Add Colorbar
        c = plt.colorbar(cs, cax = fig.add_axes([0.27, 0.18, 0.010, 0.32]))
        c.ax.tick_params(labelsize=6)
        c.set_label('Exhaust \n (ton $year^{-1}$)',fontsize=7)
        
        
        # Add Colorbar
        
        c2 = plt.colorbar(cs2, cax = fig.add_axes([0.59, 0.18, 0.010, 0.32]))
        c2.ax.tick_params(labelsize=6)
        c2.set_label('Soil ressusp. \n (ton $year^{-1}$)',fontsize=7 )
        
        
        c3 = plt.colorbar(cs3, cax = fig.add_axes([0.91, 0.18, 0.010, 0.32]))
        c3.ax.tick_params(labelsize=6)
        c3.set_label('Brake and tires \n (ton $year^{-1}$)',fontsize=7)
        

        for index, axs in enumerate(ax):
            axs.text(0.5, 1, '('+string.ascii_lowercase[index]+')',
                    transform=axs.transAxes,
                    size=8, weight='bold')

        fig.subplots_adjust(wspace=0, hspace=0)
        fig.tight_layout()
        
       
        fig.savefig(outputPath+'/ncOut_'+str(year)+'_'+pol+'.png')
        #plt.show()
    return pec2
    
#%% Running function 
pols ='PAL'

pols = ['PAL','PCA', 'PCL', 'PEC', 'PFE', 'PK','PMC','PMG','PMN','PMOTHR','PNA',
       'PNCOM','PNH4','PNO3','POC','PSI','PSO4','PTI']

pols =['PAL']

# Inputs    
netcdfPath = '/media/leohoinaski/HDD/BRAVES_database/Outputs/BRAVESdatabaseAnnual_BR_Total_0.05x0.05_'
netcdfPath3 = '/media/leohoinaski/HDD/BRAVES_database/Outputs/BRAVESdatabaseAnnual_BR_non-exaustMP_no_resusp_Total_0.05x0.05_'
netcdfPath2 = '/media/leohoinaski/HDD/BRAVES_database/Outputs/BRAVESdatabaseAnnual_BR_non-exaustMP_Total_0.05x0.05_'
outputPath = '/media/leohoinaski/HDD/BRAVES_database/Outputs'
borderShape = '/media/leohoinaski/HDD/BRAVES_database/Shapefiles/Brasil.shp'

for pol in pols:
    pec2=BRAVES_ploting(netcdfPath,netcdfPath2,netcdfPath3,outputPath,borderShape,pol,'Spectral_r')
