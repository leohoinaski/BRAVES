#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                             BRAVESdatabase_main.py

This is the main script form BRAVES database that controls the other functions.  
       
    
Outputs:
    
    roadDensity_UF_'+str(IBGE_CODES[pp])+'.csv = csv file with road density by state.

    baseGrid.csv = grid used to calculate the road density.
    
    roadEmiss_BySource_....csv = outputs in tons per year.
    
    BrRoadEmiss_'+name+str(year)+'.csv': file with merged emissions 
    
    BRAVESdatabaseAnnual netCDF files with annual emissions
    
    BRAVESdatabase2CMAQ_ ...nc
    
    BRAVESdatabaseTempEmiss_ ...nc
    
    BRAVESdatabase2WRFCHEM_ ...nc
    
    
Domain suggestions:
    BR - (lati = -36 / latf = 8 / loni = -76/ lonf = -38)
    RO - (lati = -14 / latf = -6 / loni = -68 / lonf = -58)
    AC - (lati = -12 / latf = -6 / loni = -76 / lonf = -64)
    North - (lati = -16 / latf = 8 / loni = -76 / lonf = -44)
    Northeast - (lati = -20 / latf = -2 / loni = -52 / lonf = -32)
    MidWest - (lati = -26 / latf = -6 / loni = -64 / lonf = -42)
    SouthEast - (lati = -28 / latf = -12 / loni = -56 / lonf = -38)
    South - (lati = -36 / latf = -20 / loni = -60 / lonf = -46)          


IBGE code description:                
    IBGE_CODES =  11 RO   12 ACRE  13 AM   14 RR   15 PA   16 AP   17 TO   21 MA
                  22 PI   23 CE    24 RN   25 PB   26 PE   27 AL   28 SE   29 BA
                  31 MG   32 ES    33 RJ   35 SP   41 PR   42 SC   43 RS   50 MS
                  51 MT   52 GO    53 DF    


Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br

-------------------------------------------------------------------------------
"""
#%%============================================================================
#from roadDensity_v1 import roadDensity
from parallelRoadDensity_v2 import roadDensity,roadDensityMCIP,roadDensityWRF
from roadEmiss_v1 import roadEmiss
from mergeRoadEmiss_v1 import mergeRoadEmiss
from BRAVES2netCDF_v1 import BRAVES2netCDF
from BRAVES_temporalDisag_v1 import BRAVES_temporalDisag
from netCDFcreator_v1 import createNETCDFtemporalfromNC
from netCDFcreator_v1 import createNETCDFtemporalBySpecies
from netCDFcreator_v1 import createNETCDFtemporalfromNCforWRFCHEM
import os
import numpy as np
import pandas as pd
import netCDF4 as nc


#================================INPUTS======================================
# userGrid -> 0 for user-defined / 1 for mcip based grid/ 2 for WRF based grid
userGrid = 2

# If userGrid  = 0 - Users can change the domain and resolution here.
lati = -36 #(Brazil) #lati = int(round(bound.miny)) # Initial latitud>

latf = 8 #(Brazil) #latf = int(round(bound.maxy)) # Final latitude (>

loni = -76 #(Brazil) #loni = int(round(bound.minx)) # Initial longit>

lonf = -32 #(Brazil) #lonf = int(round(bound.maxx)) # Final longitu>

deltaX = 0.05 # Grid resolution/spacing in x direction

deltaY = 0.05 # Grig resolution/spacing in y direction

# If userGrid = 1, define the path to mcip grid - GRIDDOT2D file
mcipPath = '/media/leohoinaski/HDD/GRIDDOT2D_SC.nc'

# If userGrid = 2, define the path to WRF geogrid file
wrfPath = '/media/leohoinaski/HDD/wrfout_d02_2019-06-01_00:00:00'
geowrfPath = '/media/leohoinaski/HDD/geo_em.d02.nc'

# -----------------------------Emission type-----------------------------------
# Type of emission to run
typeEmiss = 'TOTAL' 

#'TOTAL' = Total emissions/sum of emissions types
#'Exhaust' = Only exhaust emissions
#'non-exaust' = Only non-exhaust emissions
#'non-exaustMP' = Only Particulate Matter non-exhaust emissions
#'non-exaustMP_no_resusp'= Only Particulate Matter non-exhaust emissions 
#                          excluding road resuspension            
             
#--------------------Brazilian states to incluse in your database--------------
IBGE_CODES = [11,12,13,14,15,16,17,
              21,22,23,24,25,26,27,28,29,
              31,32,33,35,
              41,42,43,
              50,51,52,53] # include the IBGE code from the states to be consid>

#IBGE_CODES = [41,42,43] 

#------------------ Years for running annual BRAVESdatabase--------------------
years=[2013,2014,2015,2016,2017,2018,2019]
#years = [2019]

#----------------------- Output identification---------------------------------
fileId = 'BR_' # Code to identify your output files 

#-------------------Controls and Outputs definition----------------------------
# Run or not road density calculation. If you choose this option, the 
# roadDensity calculation will start. This might take long time if you set 
# a large domain or small detalX/Y
runOrnotRoadDens = 0 #0 for no and 1 for yes

# This option will produce the csv files with emissions from each state
runOrnotRoadEmiss = 0 # 0 for no and 1 for yes

# This option will merge the emissions if you have more than one state
runOrnotMergeRoadEmiss = 0 # 0 for no and 1 for yes

# This option will create annual netCDF files = BRAVESdatabaseAnnual
runOrnotBRAVES2netCDF = 0 # 0 for no and 1 for yes

# -----------------------Temporal disagregation--------------------------------
# If you want temporal disagregated files, you need the BRAVESdatabaseAnnual files
files = ['BRAVESdatabaseAnnual_BR_TOTAL_Total_BR_0.05x0.05_2019.nc'] # Define the files to disaggregate
yearsTempFiles=[2019] # Years to run the temporal files
months = [1] # Set the month of your simulation
days = [1] # Set the day of your simulation

# This option Create disaggregated files - temporal, spatial, and one specie
# You should define the year and specie to create your files
runOrnotTempFiles = 0 # 0 for no and 1 for yes
specs = ['SO2']  # Identification of the specie 

# This option create CMAQ emission inputs - temporal, spatial, and all species
# You should define the annual file to creat the CMAQ inputs
runOrnotCMAQemiss = 1 # 0 for no and 1 for yes

# This option create WRFCHEM emission inputs - temporal, spatial, and all species
# You should define the annual file to creat the WRFCHEM inputs
runOrnotWRFCHEMemiss=0

#%%============================PROCESSING========================================

print('====================== BRAVES database v1 ==========================')
# cd to the main folder
#------------------------------Setting folders---------------------------------

# Path to functions and BRAVESdatabase_main.py
rootPath = os.path.abspath(os.getcwd())
os.chdir(rootPath)

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

# Creating output directory
if os.path.isdir(outPath)==0:
    os.mkdir(outPath)


# THis is your grid identification   
if userGrid == 0:
    roadDensPrefix = fileId+str(deltaX)+'x'+str(deltaY) # grid definition identification
elif userGrid == 1:
    roadDensPrefix = fileId+'MCIPgrid' # grid definition identification
elif userGrid == 2:
    roadDensPrefix = fileId+'WRFgrid' # grid definition identification
    
    
# Calling roadDensity function
if runOrnotRoadDens==1:
    if userGrid == 0:
        roadDensity(dirPath,outPath,IBGE_CODES,lati,latf,loni,lonf,
                    deltaX,deltaY,roadFileName,roadDensPrefix)
    elif userGrid == 1:
        roadDensityMCIP(dirPath,outPath,IBGE_CODES,lati,latf,loni,lonf,
                    deltaX,deltaY,roadFileName,roadDensPrefix,mcipPath)
    elif userGrid == 2:
        roadDensityWRF(dirPath,outPath,IBGE_CODES,lati,latf,loni,lonf,
                    deltaX,deltaY,roadFileName,roadDensPrefix,geowrfPath)


# Calling roadEmiss function
if runOrnotRoadEmiss==1:
    roadEmiss(outPath,bravesPath,years,IBGE_CODES,roadDensPrefix,typeEmiss)


# Calling mergeRoadEmiss function 
if runOrnotMergeRoadEmiss==1:
    mergeRoadEmiss(outPath,years,IBGE_CODES,roadDensPrefix,typeEmiss)


# Calling BRAVES2netCDF
if runOrnotBRAVES2netCDF==1:
    BRAVES2netCDF(outPath,folderSpec,outPath,years,fileId+typeEmiss,roadDensPrefix,typeEmiss)


# Creating input files for CMAQ
if runOrnotCMAQemiss==1:
    for file in files:
        for year in yearsTempFiles:
            for month in months:
                for day in days:
                    dataTempo=None
                    dataTempo,xX,yY,disvec,prefix,area = BRAVES_temporalDisag(rootPath,outPath,file,month,day)
                    for jj in np.unique(disvec.day):       
                        name = 'BRAVESdatabase2CMAQ'+'_'+str(year)+'_'+str(month)+'_'+str(jj)+'.nc'
                        dayT = np.where(disvec.day==jj)
                        createNETCDFtemporalfromNC(outPath,name,dataTempo,xX,yY,disvec,area)

                    
# Creating temporal files for one specie
if runOrnotTempFiles==1:
    for spec in specs:
        for file in files:
            for year in yearsTempFiles:
                for month in months:
                    for day in days:
                        file_path = 'CMAQ_speciesMW.csv'
                        smm = pd.read_csv(folderSpec+'/'+file_path)
                        specIdx=smm[smm.ID==spec].index.to_numpy()[0]
                        dataTempo=None
                        dataTempo,xX,yY,disvec,prefix,area=BRAVES_temporalDisag(rootPath,outPath,file,month,day)
                        for jj in np.unique(disvec.day):       
                            name = 'BRAVESdatabaseTempEmiss_'+spec+'_'+str(year)+'_'+str(month)+'_'+str(jj)+'.nc'
                            dayT = np.where(disvec.day==jj)
                            createNETCDFtemporalBySpecies(outPath,name,
                                                          dataTempo[:,specIdx,:,:].reshape((dataTempo.shape[0],1,dataTempo.shape[2],dataTempo.shape[3])),
                                                          xX,yY,disvec,spec,area)

if runOrnotWRFCHEMemiss==1:
    print('Running BRAVESdatabase2WRFCHEM')
    for file in files:
        ds3 = nc.Dataset(wrfPath)
        time=ds3['Times'][:]       
        years=[]
        months=[]
        days=[]
        for tt in time:
            years.append(np.array("". join(np.array(tt[0:4],dtype=np.str_)), 
                                      dtype=np.int16))
            months.append(np.array("". join(np.array(tt[5:7],dtype=np.str_)), 
                                      dtype=np.int16))       
            days.append(np.array("". join(np.array(tt[8:10],dtype=np.str_)), 
                                      dtype=np.int16))
        
        dates = np.unique(np.array([years,months,days]), axis=1).transpose()
        # years = dates[0,:]
        # months = dates[1,:]
        # days = dates[2,:]

        for date in dates:
            dataTempo,xX,yY,disvec,prefix,area = BRAVES_temporalDisag(rootPath,outPath,file,date[1],date[2])    
            name = 'BRAVESdatabase2WRFCHEM'+'_'+str(date[0])+'_'+str(date[1])+'_'+str(date[2])+'.nc'                      
            createNETCDFtemporalfromNCforWRFCHEM(rootPath,outPath,
                                                 name,dataTempo,
                                                 xX,yY,disvec,area,
                                                 wrfPath)
            del dataTempo
