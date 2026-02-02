#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 18:05:07 2022

@author: leohoinaski
"""

#%%============================================================================
#from roadDensity_v1 import roadDensity
from parallelRoadDensity_v3 import roadDensityMCIP
from roadEmiss_v2 import roadEmiss
from mergeRoadEmiss_v2 import mergeRoadEmiss
from BRAVES2netCDF_v1 import BRAVES2netCDF
from BRAVES_temporalDisag_v1 import BRAVES_temporalDisagMCIP
from netCDFcreator_v1 import createNETCDFtemporalfromNC
import os
import numpy as np
import netCDF4 as nc
import datetime
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', default=0, action='count')
    parser.add_argument('BRAVEShome')
    parser.add_argument('mcipPath')
    parser.add_argument('GDNAM')
    parser.add_argument('YEAR', type=int)
    parser.add_argument('runOrnotRoadDens', type=int)
    parser.add_argument('runOrnotRoadEmiss', type=int)
    parser.add_argument('runOrnotMergeRoadEmiss', type=int)
    parser.add_argument('runOrnotBRAVES2netCDF', type=int)
    parser.add_argument('runOrnotCMAQemiss', type=int)
    args = parser.parse_args()
    
    # Shell inputs
    fileId = args.GDNAM
    mcipPath = args.mcipPath
    year = np.array(args.YEAR)
    rootPath = args.BRAVEShome
    runOrnotRoadDens = args.runOrnotRoadDens #1 #0 for no and 1 for yes
    runOrnotRoadEmiss = args.runOrnotRoadEmiss #1 # 0 for no and 1 for yes
    runOrnotMergeRoadEmiss = args.runOrnotMergeRoadEmiss #1 # 0 for no and 1 for yes
    runOrnotBRAVES2netCDF = args.runOrnotBRAVES2netCDF #1 # 0 for no and 1 for yes
    runOrnotCMAQemiss = args.runOrnotCMAQemiss #1 # 0 for no and 1 for yes

    # BRAVES inputs and paths
    dirPath = rootPath + '/Shapefiles'  
    outPath = rootPath +'/Outputs'
    bravesPath = rootPath +'/BRAVESoutputs'
    folderSpec = rootPath +'/ChemicalSpec'
    roadFileName = 'gis_osm_roads_free_1.shp'
    roadDensPrefix = fileId+'_'+'MCIPgrid' # grid definition identification
    mcipGRIDDOT2DPath = mcipPath+'/GRIDDOT2D_'+fileId+'.nc'
    mcipMETCRO2Dpath = mcipPath+'/METCRO2D_'+fileId+'.nc'
    IBGE_CODES = [41,42,43] 
    typeEmiss = 'TOTAL' 
    fleetEmiss = 'Total'

    
    # Running
    os.chdir(rootPath)
    file = 'BRAVESdatabaseAnnual_'+fileId+'_'+typeEmiss+'_'+fleetEmiss+'_'+fileId+'_MCIPgrid_'+str(year)+'.nc' # Define the files to disaggregate

    if os.path.isdir(outPath)==0:
        os.mkdir(outPath)
    outPath = rootPath +'/Outputs/'+fileId
    if os.path.isdir(outPath)==0:
        os.mkdir(outPath)

    # Calling roadDensity function
    if runOrnotRoadDens==1:
        roadDensityMCIP(dirPath,outPath,IBGE_CODES,roadFileName,roadDensPrefix,mcipGRIDDOT2DPath,'highways')
        roadDensityMCIP(dirPath,outPath,IBGE_CODES,roadFileName,roadDensPrefix,mcipGRIDDOT2DPath,'road')
    
    # Calling roadEmiss function
    if runOrnotRoadEmiss==1:
        roadEmiss(outPath,bravesPath,year,IBGE_CODES,roadDensPrefix,typeEmiss)
    
    # Calling mergeRoadEmiss function 
    if runOrnotMergeRoadEmiss==1:
        mergeRoadEmiss(outPath,year,IBGE_CODES,roadDensPrefix,typeEmiss)
    
    # Calling BRAVES2netCDF
    if runOrnotBRAVES2netCDF==1:
        BRAVES2netCDF(outPath,folderSpec,outPath,year,fileId+'_'+typeEmiss,roadDensPrefix,typeEmiss)

    # Creating input files for CMAQ
    if runOrnotCMAQemiss==1:
        ds3 = nc.Dataset(mcipMETCRO2Dpath)
        time=ds3['TFLAG'][:]       
        dt0 = datetime.datetime.strptime(str(time[:,0,:][:,0][0]),'%Y%j').date()
        dt1 = datetime.datetime.strptime(str(time[:,0,:][:,0][-1]),'%Y%j').date()
        hours = [np.array(time[:,0,:][:,1]/10000)[0],
                 np.array(time[:,0,:][:,1]/10000)[-1]]
        dataTempo=None
        dataTempo,xX,yY,disvec,prefix,area = BRAVES_temporalDisagMCIP(
            rootPath,outPath,file,time)
        name = 'BRAVESdatabase2CMAQ'+\
            '_'+str(dt0.year)+'_'+str(dt0.month).zfill(2)+'_'+str(dt0.day).zfill(2)+'_'+str(int(hours[0])).zfill(2)+'00'+\
                '_to_'+str(dt1.year)+'_'+str(dt1.month).zfill(2)+'_'+str(dt1.day).zfill(2)+'_'+str(int(hours[1])).zfill(2)+'00'+'.nc'
        createNETCDFtemporalfromNC(outPath,name,dataTempo,xX,yY,mcipMETCRO2Dpath)
