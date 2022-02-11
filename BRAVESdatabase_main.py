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
rootPath = os.path.abspath(os.getcwd())

# Defining folder path with openstreetmaps folder
dirPath = rootPath + "/Shapefiles"

# Path to folder with roadLenght....csv
outPath = rootPath +"/Outputs"

# Path to folder with BRAVES outputs'
bravesPath = rootPath +"/BRAVESoutputs"

# Path to folder with chemical speciation profiles
folderSpec = rootPath +'/ChemicalSpec'

# Road file name - use gis_osm for states other than SC
#roadFileName = 'SC_ROADS.shp' # for SC 
#roadFileName='roadFl.shp'
roadFileName = 'gis_osm_roads_free_1.shp'

#-------------------------Setting grid resolution------------------------------

# Users can change the domain and resolution here.
lati =-40 #lati = int(round(bound.miny)) # Initial latitude (lower-left)

latf = 10 #latf = int(round(bound.maxy)) # Final latitude (upper-right)

loni = -80 #loni = int(round(bound.minx)) # Initial longitude (lower-left)

lonf = -30 #lonf = int(round(bound.maxx)) # Final longitude (upper-right)

deltaX = 0.05 # Grid resolution/spacing in x direction

deltaY = 0.05 # Grig resolution/spacing in y direction

#---------------------------Vehicular emissions--------------------------------

# IBGE_CODES = [11,12,13,14,15,16,17,
#               21,22,23,24,25,26,27,28,29,
#               31,32,33,34,35,
#               41,42,43,
#               50,51,52,53] # include the IBGE code from the states to be considered

IBGE_CODES = [21,22,23,24,25,26,27,28,29,
              31,32,33,34,35,
              41,42,43,
              50,51,52,53] 

#---------------------------- Time window--------------------------------------

years=[2013,2014,2015,2016,2017,2018,2019]

months = [1] # Set the month of your simulation

days = [1,2] # Set the day of your simulation


#-------------------Controls and Outputs definition----------------------------

# Run or not road density calculation.
runOrnotRoadDens = 1 #0 for no and 1 for yes

runOrnotRoadEmiss = 1 # 0 for no and 1 for yes

runOrnotMergeRoadEmiss = 1 # 0 for no and 1 for yes

runOrnotBRAVES2netCDF = 1 # 0 for no and 1 for yes

# Create disaggregated files - temporal, spatial, and one specie
runOrnotTempFiles = 0 # 0 for no and 1 for yes
specs = ['SO2']  # Identification of the specie 

# Create CMAQ emission inputs - temporal, spatial, and all species
runOrnotCMAQemiss = 0 # 0 for no and 1 for yes
files = ['BRAVESdatabaseAnnual_SC_ComLight_2013.nc'] # Define the files to disaggregate


fileId = 'BR' # Code to identify your output files

roadDensPrefix = str(deltaX)+'x'+str(deltaY) # grid definition identification


#%%============================PROCESSING========================================

print('=================== BRAVES database v1 =======================')
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
if runOrnotRoadEmiss==1:
    roadEmiss(outPath,bravesPath,years,IBGE_CODES,roadDensPrefix)

# Calling mergeRoadEmiss function 
if runOrnotMergeRoadEmiss==1:
    mergeRoadEmiss(outPath,years,IBGE_CODES,roadDensPrefix)

# Calling BRAVES2netCDF
if runOrnotBRAVES2netCDF==1:
    BRAVES2netCDF(outPath,folderSpec,outPath,years,fileId,roadDensPrefix)

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
