#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                            BRAVES2netCDF_v1.py

This function creates speciated and spatial disagregated data
ready to convert in netCDF files by the netCDFcreator_v1.py.


Subfunctions:
    
    gridding(lon, lat)
        It creates the meshgrid matrix using lower-left pixel coordinates.
        
    populatingGridMat(dataMat, center, xX, yY)
        It populates the matrix ready to convert the netCDF file
        
    splitnetCDFfiles(dataEmiss,centerX,xX,yY,xv,yv,lat,lon,year,prefix,outPath,fileId)
        This function prepares the data to creat the netCDF file
        
    BRAVES2netCDF(folder,folderSpec,outPath,years,fileId)
        Controling function
    
Inputs:
    
    folder: folder with BrRoadEmiss_BySource_Light_
    
    folderSpec: folder with 'BRAVES_speciation.csv' and 'CMAQ_speciesMW.csv'
            files
            
    outPath: folder to output files
    
    years: respective years of emission inventories
    
    fileId: Code to identify your output files


Outputs:
        
    BRAVESannualEmiss_ 
    
    
External functions:
    
    netCDFcreator_v1 and BRAVES_ChemicalSpec_v1
    

Last update = 29/09/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br
-------------------------------------------------------------------------------
"""

# Importando bibliotecas
import geopandas as gpd
import pandas as pd
import os
from shapely import wkt
import numpy as np
from shapely.geometry import MultiLineString
from shapely.ops import polygonize
from netCDFcreator_v1 import createNETCDFtemporal
#from BRAVES_temporalDisag import temporalDisagVehicular
from BRAVES_ChemicalSpec_v1 import ChemicalSpeciationLight, ChemicalSpeciationHeavy
import datetime
import numpy.matlib

#%% -----------------------------INPUTS---------------------------------------
conver = 1 # from grams to grams  :)

#%% Gridding and populatingGrid functions

def gridding(lon,lat):
    xv, yv = np.meshgrid(lon, lat)
    hlines = [((x1, yi), (x2, yi)) for x1, x2 in zip(lon[:-1], lon[1:]) for yi in lat]
    vlines = [((xi, y1), (xi, y2)) for y1, y2 in zip(lat[:-1], lat[1:]) for xi in lon]
    grids = list(polygonize(MultiLineString(hlines + vlines)))
    grids = gpd.GeoDataFrame(grids) 
    grids.columns =['geometry'] 
    grids['geometry'] = grids['geometry']
    grids.crs = "EPSG:4326"  
    grids['X'] = grids.geometry.centroid.x
    grids['Y'] = grids.geometry.centroid.y
    xX = np.array(grids['X']).reshape((lon.shape[0]-1,lat.shape[0]-1)).transpose()
    yY = np.array(grids['Y']).reshape((lon.shape[0]-1,lat.shape[0]-1)).transpose()
    return grids,xv,yv,xX,yY


def populatingGridMat(dataMat,center,xX,yY):
    print('populatingGridMat')
    dataTempo = np.zeros([dataMat.shape[2],dataMat.shape[1],xX.shape[0], yY.shape[1]])
    for ii in range(0,dataMat.shape[2]):
        for jj in range(0,dataMat.shape[1]):
            dataTempo[ii,jj,:,:] = dataMat[:,jj,ii].reshape(xX.shape[0],yY.shape[1],order='F')
    return dataTempo
            


#%% Creating netCDF file
def splitnetCDFfiles(dataEmiss,centerX,xX,yY,xv,yv,lat,lon,year,prefix,outPath,fileId,roadDensPrefix):
    month = 1
    print('Spliting netCDF files')
       
    dataMat = np.zeros([dataEmiss.shape[0], dataEmiss.shape[1],1])
    
    for ii in range(0,dataMat.shape[1]):
        for jj in range(0,dataMat.shape[2]):
            dataMat[:,ii,jj]= dataEmiss.iloc[:,ii]
    
    startDate = datetime.datetime(year, month, 1, 0, 0)
    endDate = datetime.datetime(year, month, 1, 1, 0)
    datePfct = np.arange(np.datetime64(startDate),np.datetime64(endDate),3600000000)
    disvec = pd.DataFrame()
    disvec = disvec.reindex(datePfct, fill_value=np.nan)    
    disvec['year'] = year
    disvec['month'] = 1
    disvec['day'] = 1
    disvec['hour'] = 00
    
    name = 'BRAVESdatabaseAnnual_'+fileId+'_'+prefix+'_'+roadDensPrefix+'_'+str(year)+'.nc'
    dataTempo = populatingGridMat(dataMat[:,:,:],centerX,xX,yY)
    createNETCDFtemporal(outPath,name,dataTempo,xv,yv,lat,lon,centerX,disvec,month)
    print(name +' is ready')
    
    return dataTempo



#%%--------------Reading Road emission from BRAIN/BRAVES----------------------------------

def BRAVES2netCDF (folder,folderSpec,outPath,years,fileId,roadDensPrefix,typeEmiss):
    print('===================STARTING BRAVES2netCDF_v1.py=======================')
    conver = 1    
    for year in years:
        file_path = [filename for filename in os.listdir(folder) if 
                     filename.startswith("mergeRoadEmiss_BySource_Light_"+typeEmiss+'_'+\
                         roadDensPrefix +'_' + str(year)+'.csv')]
        print(file_path)
        df = pd.read_csv(folder+'/'+file_path[0])
        roadE = gpd.GeoDataFrame(df) 
        roadE['geometry'] = roadE['geometry'].apply(wkt.loads) 
        roadE.crs = "EPSG:4326"     
        roadE=roadE.reset_index(drop=True)
        
        # geod = roadE.crs.get_geod()
        # area=[]
        # for re in roadE.geometry:
        #     a = abs(geod.geometry_area_perimeter(re)[0])
        #     area.append(a)
        # area = np.array(area)/(1000*1000)
        

        
        # Extracting coords
        arrX = []
        arrY = []
        for ii in range(0, roadE.shape[0]):
            arrX.append(np.array(roadE['geometry'][ii].exterior.coords.xy[0]))
            arrY.append(np.array(roadE['geometry'][ii].exterior.coords.xy[1]))
          
        arrX = np.array(arrX)
        arrY = np.array(arrY)
        lon = np.unique(arrX[:])
        lat = np.unique(arrY[:])
        
        #%% Calling gridding functions 
        
        grids,xv,yv,xX,yY = gridding(lon,lat)
        
        # Road emission
        centerRoad = roadE.geometry.centroid
        centerRoad.to_crs("EPSG:4326")
        
        
        #%% Calling speciation function
        
        file_path = 'BRAVES_speciation.csv'
        dfSpc = pd.read_csv(folderSpec+'/'+file_path)
        dfSpc.iloc[:,3:] = dfSpc.iloc[:,3:] # em porcentagem
        
        file_path = 'CMAQ_speciesMW.csv'
        smm = pd.read_csv(folderSpec+'/'+file_path)
        
        #Commercial light
        file_path = [filename for filename in os.listdir(folder) if
                     filename.startswith("mergeRoadEmiss_BySource_ComLight_"+typeEmiss+'_'+\
                         roadDensPrefix +'_' + str(year)+'.csv')]
        df = pd.read_csv(folder+'/'+file_path[0])
        roadX = gpd.GeoDataFrame(df) 
        roadX['geometry'] = roadX['geometry'].apply(wkt.loads) 
        roadX.crs = "EPSG:4326"     
        roadX=roadX.reset_index(drop=True)
        roadX.insert(0, 'A', 0)
        centerX = roadX.geometry.centroid
        centerX.to_crs("EPSG:4326")
        dataEmiss1 = ChemicalSpeciationLight(roadX,dfSpc,smm,conver)
        prefix = 'ComLight'
        splitnetCDFfiles(dataEmiss1,centerX,xX,yY,xv,yv,lat,lon,year,prefix,outPath,fileId,roadDensPrefix)
        
        # Light
        file_path = [filename for filename in os.listdir(folder) if 
                     filename.startswith("mergeRoadEmiss_BySource_Light_"+typeEmiss+'_'+\
                         roadDensPrefix +'_' + str(year)+'.csv')]
        df = pd.read_csv(folder+'/'+file_path[0])
        roadX = gpd.GeoDataFrame(df) 
        roadX['geometry'] = roadX['geometry'].apply(wkt.loads) 
        roadX.crs = "EPSG:4326"     
        roadX=roadX.reset_index(drop=True)
        roadX.insert(0, 'A', 0)
        centerX = roadX.geometry.centroid
        centerX.to_crs("EPSG:4326")
        dataEmiss2 = ChemicalSpeciationLight(roadX,dfSpc,smm,conver)
        prefix = 'Light'
        splitnetCDFfiles(dataEmiss2,centerX,xX,yY,xv,yv,lat,lon,year,prefix,outPath,fileId,roadDensPrefix)
        
        # Motorcycles
        file_path = [filename for filename in os.listdir(folder) if 
                     filename.startswith("mergeRoadEmiss_BySource_Motorcycle_"+typeEmiss+'_'+\
                         roadDensPrefix +'_' + str(year)+'.csv')]
        df = pd.read_csv(folder+'/'+file_path[0])
        roadX = gpd.GeoDataFrame(df) 
        roadX['geometry'] = roadX['geometry'].apply(wkt.loads) 
        roadX.crs = "EPSG:4326"     
        roadX=roadX.reset_index(drop=True)
        roadX.insert(0, 'A', 0)
        centerX = roadX.geometry.centroid
        centerX.to_crs("EPSG:4326")
        dataEmiss3 = ChemicalSpeciationLight(roadX,dfSpc,smm,conver)
        prefix = 'Motorcycles'
        splitnetCDFfiles(dataEmiss3,centerX,xX,yY,xv,yv,lat,lon,year,prefix,outPath,fileId,roadDensPrefix)
        
        # Heavy
        file_path = [filename for filename in os.listdir(folder) if 
                     filename.startswith("mergeRoadEmiss_BySource_Heavy_"+typeEmiss+'_'+\
                         roadDensPrefix +'_' + str(year)+'.csv')]
        df = pd.read_csv(folder+'/'+file_path[0])
        roadX = gpd.GeoDataFrame(df) 
        roadX['geometry'] = roadX['geometry'].apply(wkt.loads) 
        roadX.crs = "EPSG:4326"     
        roadX=roadX.reset_index(drop=True)
        roadX.insert(0, 'A', 0)
        centerX = roadX.geometry.centroid
        centerX.to_crs("EPSG:4326")
        dataEmiss4 = ChemicalSpeciationHeavy(roadX,dfSpc,smm,conver)
        prefix = 'Heavy'
        splitnetCDFfiles(dataEmiss4,centerX,xX,yY,xv,yv,lat,lon,year,prefix,outPath,fileId,roadDensPrefix)
        
        
        #TOTAL EMISSIONS
        dataEmiss = dataEmiss1+dataEmiss2+dataEmiss3+dataEmiss4
        prefix='Total'
        splitnetCDFfiles(dataEmiss,centerX,xX,yY,xv,yv,lat,lon,year,prefix,outPath,fileId,roadDensPrefix)
        print('Your files are ready in: ' + outPath)
    return roadX, dataEmiss
