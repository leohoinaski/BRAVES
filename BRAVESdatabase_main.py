#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                             BRAVESdatabase_main.py

This is the main script form BRAVES database that controls the other functions.

Inputs: 
    
    rootPath: Path to functions and BRAVESdatabase_main.py
    
    dirPath: Defining folder path with openstreetmaps folder
    
    outPath: Path to folder with roadLenght....csv
    
    bravesPath: Path to folder with BRAVES outputs'
    
    folderSpec: Path to folder with chemical speciation profiles
    
    lati: Initial latitude (lower-left)
    
    latf: Final latitude (upper-right)
    
    loni: Initial longitude (lower-left)
    
    lonf: Final longitude (upper-right)
    
    deltaX: Grid resolution/spacing in x direction
    
    deltaY: Grig resolution/spacing in y direction
    
    IBGE_CODES: State ID according to 
        https://www.ibge.gov.br/explica/codigos-dos-municipios.php
    
    years: Base-year for your simulation
    
    roadFileName = 'SC_ROADS.shp' # Only for Santa Catarina State = ID=42
        
    fileId = identification of your output files
    
    file = file to apply the temporal disaggregation
    
    
Outputs:
    
    roadDensity_UF_'+str(IBGE_CODES[pp])+'.csv = csv file with road density by state.

    baseGrid.csv = grid used to calculate the road density.
    
    roadEmiss_BySource_....csv = outputs in tons per year.
    
    BrRoadEmiss_'+name+str(year)+'.csv': file with merged emissions 
    
    BRAVESannualEmiss_ netCDF files
    
    BRAVESdatabase2CMAQ_ ...nc
    
    BRAVESdatabaseTempEmiss_ ...nc
    
    
    
External functions:
    roadDensity_v1, roadEmiss_v1, mergeRoadEmiss_v1, BRAVES2netCDF_v1,
    BRAVES_temporalDisag, createNETCDFtemporalfromNC, createNETCDFtemporalBySpecies
    

Last update = 30/09/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br

-------------------------------------------------------------------------------
"""
#from roadDensity_v1 import roadDensity
from parallelRoadDensity_v1 import roadDensity
from roadEmiss_v1 import roadEmiss
from mergeRoadEmiss_v1 import mergeRoadEmiss
from BRAVES2netCDF_v1 import BRAVES2netCDF
from BRAVES_temporalDisag_v1 import BRAVES_temporalDisag
from netCDFcreator_v1 import createNETCDFtemporalfromNC,createNETCDFtemporalBySpecies
import os
import numpy as np
import pandas as pd

#%%================================INPUTS======================================

#------------------------------Setting folders---------------------------------

# Path to functions and BRAVESdatabase_main.py
rootPath = "/media/leohoinaski/HDD/BRAVES_database"

# Defining folder path with openstreetmaps folder
dirPath = rootPath + "/Shapefiles"

# Path to folder with roadLenght....csv
outPath = rootPath +"/Outputs"

# Path to folder with BRAVES outputs'
bravesPath = rootPath +"/BRAVESoutputs"

# Path to folder with chemical speciation profiles
folderSpec = rootPath +'/ChemicalSpec'

# Road file name - use gis_osm for states other than SC
roadFileName = 'SC_ROADS.shp' # for SC 
#roadFileName='roadFl.shp'
#roadFileName = 'gis_osm_roads_free_1.shp'

#-------------------------Setting grid resolution------------------------------

# Users can change the domain and resolution here.
lati =-30 #lati = int(round(bound.miny)) # Initial latitude (lower-left)

latf = -24 #latf = int(round(bound.maxy)) # Final latitude (upper-right)

loni = -54 #loni = int(round(bound.minx)) # Initial longitude (lower-left)

lonf = -47 #lonf = int(round(bound.maxx)) # Final longitude (upper-right)

deltaX = 0.01 # Grid resolution/spacing in x direction

deltaY = 0.01 # Grig resolution/spacing in y direction

#---------------------------Vehicular emissions--------------------------------

IBGE_CODES = [42] # include the IBGE code from the states to be considered


#---------------------------- Time window--------------------------------------

years=[2013]

months = [1] # Set the month of your simulation

days = [1,2] # Set the day of your simulation


#-------------------Controls and Outputs definition----------------------------

# Run or not road density calculation.
runOrnotRoadDens = 0 #0 for no and 1 for yes

fileId = 'SC' # Code to identify your output files

roadDensPrefix = str(deltaX)+'x'+str(deltaY) # grid definition identification

# Create disaggregated files - temporal, spatial, and one specie
runOrnotTempFiles = 1 # 0 for no and 1 for yes
specs = ['SO2']  # Identification of the specie 

# Create CMAQ emission inputs - temporal, spatial, and all species
runOrnotCMAQemiss = 1 # 0 for no and 1 for yes
files = ['BRAVESdatabaseAnnual_SC_ComLight_2013.nc'] # Define the files to disaggregate


#%%============================PROCESSING========================================

# cd to the main folder
os.chdir(rootPath)

# Creating output directory
if os.path.isdir(outPath)==0:
    os.mkdir(outPath)

# Calling roadDensity function
if runOrnotRoadDens==1:
    roadDensity(dirPath,outPath,IBGE_CODES,lati,latf,loni,lonf,
                deltaX,deltaY,roadFileName,roadDensPrefix)

# Calling roadEmiss function
roadEmiss(outPath,bravesPath,years,IBGE_CODES,roadDensPrefix)

# Calling mergeRoadEmiss function 
mergeRoadEmiss(outPath,years,IBGE_CODES,roadDensPrefix)

# Calling BRAVES2netCDF
BRAVES2netCDF(outPath,folderSpec,outPath,years,fileId)

# Creating input files for CMAQ
if runOrnotCMAQemiss==1:
    for file in files:
        for year in years:
            for month in months:
                for day in days:
                    dataTempo=None
                    dataTempo,xv,yv,lat,lon,center,disvec,prefix=BRAVES_temporalDisag(rootPath,outPath,file,fileId,month,day,deltaX,deltaY)
                    for jj in np.unique(disvec.day):       
                        name = 'BRAVESdatabase2CMAQ_'+fileId+'_'+prefix+'_'+str(year)+'_'+str(month)+'_'+str(jj)+'.nc'
                        dayT = np.where(disvec.day==jj)
                        createNETCDFtemporalfromNC(outPath,name,dataTempo,xv,yv,lat,lon,center,disvec,month)

# Creating temporal files for one specie
if runOrnotTempFiles==1:
    for spec in specs:
        for file in files:
            for year in years:
                for month in months:
                    for day in days:
                        file_path = 'CMAQ_speciesMW.csv'
                        smm = pd.read_csv(folderSpec+'/'+file_path)
                        specIdx=smm[smm.ID==spec].index.to_numpy()[0]
                        dataTempo=None
                        dataTempo,xv,yv,lat,lon,center,disvec,prefix=BRAVES_temporalDisag(rootPath,outPath,file,fileId,month,day)
                        for jj in np.unique(disvec.day):       
                            name = 'BRAVESdatabaseTempEmiss_'+fileId+'_'+prefix+'_'+str(year)+'_'+str(month)+'_'+str(jj)+'.nc'
                            dayT = np.where(disvec.day==jj)
                            createNETCDFtemporalBySpecies(outPath,name,
                                                          dataTempo[:,specIdx,:,:].reshape((dataTempo.shape[0],1,dataTempo.shape[2],dataTempo.shape[3])),
                                                          xv,yv,lat,lon,center,disvec,spec)
