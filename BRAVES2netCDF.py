#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------BRAVES2netCDF---------------------------
Created on Wed Jul 15 17:26:47 2020

Esta função converte os dados das emissões veiculares do BRAVES
provenientes do BRAIN em uma matriz para posterior criação do
arquivo em netCDF pela função netCDFcreator.py.
Grade definida de acordo com o arquivo das emissões veiculares. 


Funções associadas:
    gridding(lon, lat)
        Cria a grade com meshgrid dos pontos inferiores e centroides
        dos pixels
    populatingGrid(dataEmiss, center, xX, yY)
        preenche a matriz de dados da grade com os dados de emissão
    populatingGridMat(dataMat, center, xX, yY)
        preenche a matriz com dados de emissão com variabilidade tempora
    temporalDisagVehicular(dataEmissRoad, year, month)
        realiza a disagregação temporal das emissões 
    
Inputs:
    folder: pasta com os arquivos do BRAIN (emissões veiculares em toneladas/ano)
    year: ano para a criação do arquivo, conforme recorte do BRAIN
    month: mês para a criação do arquivo


Outputs:
    em g/s ou mol/s
    BRAVESannualEmiss_ e BRAVEStemporalEmiss_
    
    
Referências:
    https://www.cmascenter.org/speciation_tool/documentation/5.0/Ramboll_sptool_users_guide_V5.pdf


@author: leohoinaski
Atualização em: 04/05/2021

Atualização em: 11/05/2021 - adição da especiação química
---------------------------------------------------------------
"""
# Importando bibliotecas
import geopandas as gpd
import pandas as pd
import os 
from shapely import wkt
import numpy as np
from shapely.geometry import MultiLineString
from shapely.ops import polygonize
from netCDFcreator import createNETCDF, createNETCDFtemporal
import datetime
import numpy.matlib

#============================= INPUTS =========================================
# projection = 'epsg:4326'
# # XORIG, YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK - According to MCIP
# yorig = -4 # ponto inicial da grade em y
# xorig = -48 # ponto inicial da grade em x
# xcell = 500 # espaçamento em x
# ycell = 500 # espaçamento em y
# ncols = 100 # número de colunas da grade em x
# nrols = 100 # número de pontos da grade em x
# startdate = 1 # data de início
# enddate = 10 # data de fim
# tstep = 1 # time step

folder = '/home/leohoinaski/Desktop/BRAIN_output_2021-05-07/'
#year=2013
month = 1

folderSpec = '/home/leohoinaski/website_leo/static/inventory/speciation/'

#%%--------------Reading Road emission from BRAIN/BRAVES----------------------------------

file_path = [filename for filename in os.listdir(folder) if filename.startswith("roadEmissCut_")]
df = pd.read_csv(folder+file_path[0])
roadE = gpd.GeoDataFrame(df) 
roadE['geometry'] = roadE['geometry'].apply(wkt.loads) 
roadE.crs = "EPSG:4326"     
roadE=roadE.reset_index(drop=True)

# Extracting year
year = int(file_path[0].split("_")[1].split(".")[0])

# Extracting coords
arrX = []
arrY = []
for ii in range(0, roadE.shape[0]):
    arrX.append(np.array(roadE.iloc[ii,1].exterior.coords.xy[0]))
    arrY.append(np.array(roadE.iloc[ii,1].exterior.coords.xy[1]))
  
arrX = np.array(arrX)
arrY = np.array(arrY)
lon = np.unique(arrX[:])
lat = np.unique(arrY[:])


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

def populatingGrid(dataEmiss,center,xX,yY):   
    data = np.zeros([np.size(yv,0)-1, np.size(xv,1)-1,dataEmiss.shape[1]])
    xcenter = center.geometry.centroid.x
    ycenter = center.geometry.centroid.y
   
    for ii in range(0,dataEmiss.shape[0]):
        dist = ((xcenter[ii]-xX)**2 + (ycenter[ii]-yY)**2)**(1/2)
        mindist = np.where(dist == np.amin(dist))
        for kk in range (0,dataEmiss.shape[1]):
            data[mindist[0][0],mindist[1][0],kk]= np.nansum([data[mindist[0][0],mindist[1][0],kk],dataEmiss.iloc[ii,kk]])        
    return data

def populatingGridMat(dataMat,center,xX,yY):   
    dataTempo = dataMat[:,:,:].reshape(dataMat.shape[2],dataMat.shape[1],
                                       xX.shape[0],yY.shape[0],order='A')
    dataTempo = np.zeros([dataMat.shape[2],dataMat.shape[1],np.size(yv,0)-1, np.size(xv,1)-1])
    for ii in range(0,dataMat.shape[2]):
        dataTempo[ii,0,:,:] = np.rot90(dataMat[:,:,ii].reshape(xX.shape[0],yY.shape[0]))
    
    # dataTempo = np.zeros([dataMat.shape[2],dataMat.shape[1],np.size(yv,0)-1, np.size(xv,1)-1])
    # xcenter = center.geometry.centroid.x
    # ycenter = center.geometry.centroid.y
   
    # for ii in range(0,dataMat.shape[0]):
    #     dist = ((xcenter[ii]-xX)**2 + (ycenter[ii]-yY)**2)**(1/2)
    #     mindist = np.where(dist == np.amin(dist))
    #     print('cell number = ' + str(ii))
    #     for kk in range (0,dataMat.shape[1]):
    #         dataTempo[:,kk,mindist[0][0],mindist[1][0]]= np.nansum([dataTempo[:,kk,mindist[0][0],mindist[1][0]],dataMat[ii,kk,:]],0)        
    return dataTempo

#data, xv, yv = populatingGrid(dataEmiss,center,lon,lat)
# outProj = Proj(init=projection) 
# inProj = Proj(init='epsg:4326')
# x2,y2 = transform(inProj,outProj,xv[0,:],yv[:,0]) # converting to lat lon  
   
#%% temporal disagregation

def temporalDisagVehicular(dataEmissRoad, year, month):
    # Create the MultiIndex from pollutant and time.
    #year=2013
    startDate = datetime.datetime(year, month, 1, 0, 0)
    endDate = datetime.datetime(year, month+1, 1, 0, 0)
    datePfct = np.arange(np.datetime64(startDate),np.datetime64(endDate),3600000000)
    numWeeks = datePfct.shape[0]/(7*24) # Number of weeks
    disvec = pd.DataFrame()
    disvec = disvec.reindex(datePfct, fill_value=np.nan)
    disvec['year'] = disvec.index.year
    disvec['month'] = disvec.index.month
    disvec['day'] = disvec.index.day
    disvec['hour'] = disvec.index.hour
    disvec['weekday'] = disvec.index.weekday # Monday is 0 and Sunday is 6
    hourdis = [0.014771048744461,0.00738552437223,0.003692762186115,
               0.003692762186115,0.003692762186115,0.00738552437223,
               0.029542097488922,0.062776957163959,0.062776957163955,
               0.059084194977844,0.059084194977844,0.062776957163959,
               0.06056129985229,0.059084194977844,0.059084194977844,
               0.059822747415067,0.059822747415067,0.073855243722304,
               0.073855243722304,0.059084194977844,0.044313146233383,
               0.029542097488922,0.029542097488922,0.014771048744461]   
    disvec['hourdis'] = numpy.matlib.repmat(hourdis,1,int(disvec.shape[0]/24)).transpose()
    
    weekdis=[0.101694915254237,0.169491525423729,0.152542372881356,
             0.169491525423729,0.169491525423729,0.135593220338983,
             0.101694915254237] # Dividing by the number of weeks in a month
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
    
    monthdis = [0.08333333333333333,0.08333333333333333,0.08333333333333333,
                0.08333333333333333,0.08333333333333333,0.08333333333333333,
                0.08333333333333333,0.08333333333333333,0.08333333333333333,
                0.08333333333333333,0.08333333333333333,0.08333333333333333]
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
    dataMat = np.zeros([dataEmissRoad.shape[0], dataEmissRoad.shape[1],
                        datePfct.shape[0]])
    
    for ii in range(0,dataMat.shape[1]):
        for jj in range(0,dataMat.shape[2]):
            dataMat[:,ii,jj]= dataEmissRoad.iloc[:,ii]* disvec['prod'][jj]
    
    return dataMat,datePfct,disvec

#%% Calling gridding functions 

grids,xv,yv,xX,yY = gridding(lon,lat)

# Road emission
conver = 1000000 # from tons to grams 
centerRoad = roadE.geometry.centroid
centerRoad.to_crs("EPSG:4326")
dataEmissRoad = roadE.iloc[:,2:-1]*conver  #em g/ano agora

#%% Chemical speciation 

def ChemicalSpeciationLight(roadX,dfSpc,smm):
    # Creating dataEmiss matrix
    dataEmissX=pd.DataFrame()
    # Filling VOC emissions
    for ii in range(0,dfSpc.shape[0]):
        dataEmissX[dfSpc.ID[ii]] = dfSpc['LightVOC'][ii]*roadX.iloc[:,6]*conver/smm.iloc[ii,1] +\
           dfSpc['LightPM'][ii]*roadX.iloc[:,10]*conver/smm.iloc[ii,1]  +\
               dfSpc['road'][ii]*roadX.iloc[:,18]*conver/smm.iloc[ii,1] +\
                   dfSpc['road'][ii]*roadX.iloc[:,17]*conver/smm.iloc[ii,1] +\
                       dfSpc['brakes'][ii]*roadX.iloc[:,16]*conver/smm.iloc[ii,1] +\
                           dfSpc['tires'][ii]*roadX.iloc[:,16]*conver/smm.iloc[ii,1] +\
                               dfSpc['LightEvap'][ii]*roadX.iloc[:,19]*conver/smm.iloc[ii,1] +\
                                   dfSpc['LightEvap'][ii]*roadX.iloc[:,20]*conver/smm.iloc[ii,1] +\
                                       dfSpc['LightEvap'][ii]*roadX.iloc[:,21]*conver/smm.iloc[ii,1] +\
                                           dfSpc['LightEvap'][ii]*roadX.iloc[:,14]*conver/smm.iloc[ii,1]
    
    # Filling direct emissions
    dataEmissX['ALDX'] = roadX.iloc[:,9]*conver/smm.iloc[4,1]   
    dataEmissX['CH4_INV'] = roadX.iloc[:,7]*conver
    dataEmissX['CH4'] = roadX.iloc[:,7]*conver/smm.iloc[7,1]
    dataEmissX['CO'] = roadX.iloc[:,4]*conver/smm.iloc[10,1]
    dataEmissX['CO2_INV'] = roadX.iloc[:,11]*conver
    dataEmissX['N2O_INV'] = roadX.iloc[:,12]*conver
    dataEmissX['NO'] = roadX.iloc[:,8]*conver*0.495/smm.iloc[28,1] # Para motor a diesel... considerei
    dataEmissX['NO2'] = roadX.iloc[:,8]*conver*0.505/smm.iloc[29,1] # Para motor a diesel... considerei
    dataEmissX['PMC'] = roadX.iloc[:,10]*conver
    dataEmissX['PNCOM'] = dataEmissX['POC']*0.25
    dataEmissX['PMOTHR'] = 1-(dataEmissX['POC'] + dataEmissX['PEC'] +
                              dataEmissX['PSO4'] + dataEmissX['PNO3'] +
                              dataEmissX['PNH4'] + dataEmissX['PNCOM'] +
                              dataEmissX['PFE'] + dataEmissX['PAL'] +
                              dataEmissX['PSI'] + dataEmissX['PTI'] +
                              dataEmissX['PCA'] + dataEmissX['PMG'] +
                              dataEmissX['PK'] + dataEmissX['PMN'] +
                              dataEmissX['PNA'] + dataEmissX['PCL'] +
                              dataEmissX['PH2O'])
    dataEmissX['SO2'] = roadX.iloc[:,13]*conver/smm.iloc[53,1]
    dataEmissX['VOC_INV'] = roadX.iloc[:,6]*conver + roadX.iloc[:,14]*conver +\
        roadX.iloc[:,19]*conver + roadX.iloc[:,20]*conver + roadX.iloc[:,21]*conver
    return dataEmissX

def ChemicalSpeciationHeavy(roadX,dfSpc,smm):  
    # Creating dataEmiss matrix
    dataEmissX=pd.DataFrame()
    # Filling VOC emissions
    for ii in range(0,dfSpc.shape[0]):
        dataEmissX[dfSpc.ID[ii]] = dfSpc['HeavyVOC'][ii]*roadX.iloc[:,6]*conver/smm.iloc[ii,1] +\
           dfSpc['HeavyPM'][ii]*roadX.iloc[:,10]*conver/smm.iloc[ii,1]  +\
               dfSpc['road'][ii]*roadX.iloc[:,18]*conver/smm.iloc[ii,1] +\
                   dfSpc['road'][ii]*roadX.iloc[:,17]*conver/smm.iloc[ii,1] +\
                       dfSpc['brakes'][ii]*roadX.iloc[:,16]*conver/smm.iloc[ii,1] +\
                           dfSpc['tires'][ii]*roadX.iloc[:,16]*conver/smm.iloc[ii,1] +\
                               dfSpc['HeavyEvap'][ii]*roadX.iloc[:,19]*conver/smm.iloc[ii,1] +\
                                   dfSpc['HeavyEvap'][ii]*roadX.iloc[:,20]*conver/smm.iloc[ii,1] +\
                                       dfSpc['HeavyEvap'][ii]*roadX.iloc[:,21]*conver/smm.iloc[ii,1] +\
                                           dfSpc['HeavyEvap'][ii]*roadX.iloc[:,14]*conver/smm.iloc[ii,1]
    
    # Filling direct emissions
    dataEmissX['ALDX'] = roadX.iloc[:,9]*conver/smm.iloc[4,1]   
    dataEmissX['CH4_INV'] = roadX.iloc[:,7]*conver
    dataEmissX['CH4'] = roadX.iloc[:,7]*conver/smm.iloc[7,1]
    dataEmissX['CO'] = roadX.iloc[:,4]*conver/smm.iloc[10,1]
    dataEmissX['CO2_INV'] = roadX.iloc[:,11]*conver
    dataEmissX['N2O_INV'] = roadX.iloc[:,12]*conver
    dataEmissX['NO'] = roadX.iloc[:,8]*conver*0.495/smm.iloc[28,1] # Para motor a diesel... considerei
    dataEmissX['NO2'] = roadX.iloc[:,8]*conver*0.505/smm.iloc[29,1] # Para motor a diesel... considerei
    dataEmissX['PMC'] = roadX.iloc[:,10]*conver
    dataEmissX['PNCOM'] = dataEmissX['POC']*0.25
    dataEmissX['PMOTHR'] = 1-(dataEmissX['POC'] + dataEmissX['PEC'] +
                              dataEmissX['PSO4'] + dataEmissX['PNO3'] +
                              dataEmissX['PNH4'] + dataEmissX['PNCOM'] +
                              dataEmissX['PFE'] + dataEmissX['PAL'] +
                              dataEmissX['PSI'] + dataEmissX['PTI'] +
                              dataEmissX['PCA'] + dataEmissX['PMG'] +
                              dataEmissX['PK'] + dataEmissX['PMN'] +
                              dataEmissX['PNA'] + dataEmissX['PCL'] +
                              dataEmissX['PH2O'])
    dataEmissX['SO2'] = roadX.iloc[:,13]*conver/smm.iloc[53,1]
    dataEmissX['VOC_INV'] = roadX.iloc[:,6]*conver + roadX.iloc[:,14]*conver +\
        roadX.iloc[:,19]*conver + roadX.iloc[:,20]*conver + roadX.iloc[:,21]*conver
    return dataEmissX

#%% Creating netCDF file
def splitnetCDFfiles(dataEmiss,centerX,xX,yY,year,month,prefix,folder):
    dataMat,datePfct,disvec = temporalDisagVehicular(dataEmiss, year, month)
    
    for jj in np.unique(disvec.day):
        name = 'BRAVEStemporalEmiss_'+prefix+'_'+str(year)+'_'+str(month)+'_'+str(jj)+'.nc'
        dayT = np.where(disvec.day==jj)
        dataTempo = populatingGridMat(dataMat[:,:,dayT[0][0]:dayT[0][-1]],centerX,xX,yY)
        createNETCDFtemporal(folder,name,dataTempo,xv,yv,lat,lon,centerX,disvec.iloc[dayT[0][0]:dayT[0][-1],:],month)
    # name = 'BRAVEStemporalEmiss_'+prefix+'_'+str(year)+'_'+str(month)+'_1to5.nc'
    # dataTempo = populatingGridMat(dataMat[:,:,0:24*5],centerX,xX,yY)
    # createNETCDFtemporal(folder,name,dataTempo,xv,yv,lat,lon,centerX,disvec.iloc[0:24*5,:],month)
    # print(name +' is ready')
    # name = 'BRAVEStemporalEmiss_'+prefix+'_'+str(year)+'_'+str(month)+'_5to10.nc'
    # dataTempo = populatingGridMat(dataMat[:,:,24*5:24*10],centerX,xX,yY)
    # createNETCDFtemporal(folder,name,dataTempo,xv,yv,lat,lon,centerX,disvec.iloc[24*5:24*10,:],month)
    # print(name +' is ready')
    # name = 'BRAVEStemporalEmiss_'+prefix+'_'+str(year)+'_'+str(month)+'_10to15.nc'
    # dataTempo = populatingGridMat(dataMat[:,:,24*10:24*15],centerX,xX,yY)
    # createNETCDFtemporal(folder,name,dataTempo,xv,yv,lat,lon,centerX,disvec.iloc[24*10:24*15,:],month)
    # print(name +' is ready')
    # name = 'BRAVEStemporalEmiss_'+prefix+'_'+str(year)+'_'+str(month)+'_15to20.nc'
    # dataTempo = populatingGridMat(dataMat[:,:,24*15:24*20],centerX,xX,yY)
    # createNETCDFtemporal(folder,name,dataTempo,xv,yv,lat,lon,centerX,disvec.iloc[24*15:24*20,:],month)
    # print(name +' is ready')
    # name = 'BRAVEStemporalEmiss_'+prefix+'_'+str(year)+'_'+str(month)+'_20to25.nc'
    # dataTempo = populatingGridMat(dataMat[:,:,24*20:24*25],centerX,xX,yY)
    # createNETCDFtemporal(folder,name,dataTempo,xv,yv,lat,lon,centerX,disvec.iloc[24*20:24*25,:],month) 
    # print(name +' is ready')
    # name = 'BRAVEStemporalEmiss_'+prefix+'_'+str(year)+'_'+str(month)+'_25toEnd.nc'
    # dataTempo = populatingGridMat(dataMat[:,:,24*25:],centerX,xX,yY)
    # createNETCDFtemporal(folder,name,dataTempo,xv,yv,lat,lon,centerX,disvec.iloc[24*25:,:],month)
    return dataTempo
#%% Calling speciation function

file_path = 'BRAVES_speciation.csv'
dfSpc = pd.read_csv(folderSpec+file_path)
dfSpc.iloc[:,3:] = dfSpc.iloc[:,3:]

file_path = 'CMAQ_speciesMW.csv'
smm = pd.read_csv(folderSpec+file_path)

#Commercial light
file_path = [filename for filename in os.listdir(folder) if filename.startswith("roadEmissCutBySource_ComLight_")]
df = pd.read_csv(folder+file_path[0])
roadX = gpd.GeoDataFrame(df) 
roadX['geometry'] = roadX['geometry'].apply(wkt.loads) 
roadX.crs = "EPSG:4326"     
roadX=roadX.reset_index(drop=True)
centerX = roadX.geometry.centroid
centerX.to_crs("EPSG:4326")
dataEmiss1 = ChemicalSpeciationLight(roadX,dfSpc,smm)
prefix = 'ComLight'
#splitnetCDFfiles(dataEmiss1,centerX,xX,yY,year,month,prefix,folder)


# Light
file_path = [filename for filename in os.listdir(folder) if filename.startswith("roadEmissCutBySource_Light_")]
df = pd.read_csv(folder+file_path[0])
roadX = gpd.GeoDataFrame(df) 
roadX['geometry'] = roadX['geometry'].apply(wkt.loads) 
roadX.crs = "EPSG:4326"     
roadX=roadX.reset_index(drop=True)
centerX = roadX.geometry.centroid
centerX.to_crs("EPSG:4326")
dataEmiss2 = ChemicalSpeciationLight(roadX,dfSpc,smm)
prefix = 'Light'
#splitnetCDFfiles(dataEmiss2,centerX,xX,yY,year,month,prefix,folder)


# Motorcycles
file_path = [filename for filename in os.listdir(folder) if filename.startswith("roadEmissCutBySource_Moto")]
df = pd.read_csv(folder+file_path[0])
roadX = gpd.GeoDataFrame(df) 
roadX['geometry'] = roadX['geometry'].apply(wkt.loads) 
roadX.crs = "EPSG:4326"     
roadX=roadX.reset_index(drop=True)
centerX = roadX.geometry.centroid
centerX.to_crs("EPSG:4326")
dataEmiss3 = ChemicalSpeciationLight(roadX,dfSpc,smm)
prefix = 'Motorcycles'
#splitnetCDFfiles(dataEmiss3,centerX,xX,yY,year,month,prefix,folder)


# Heavy
file_path = [filename for filename in os.listdir(folder) if filename.startswith("roadEmissCutBySource_Heavy")]
df = pd.read_csv(folder+file_path[0])
roadX = gpd.GeoDataFrame(df) 
roadX['geometry'] = roadX['geometry'].apply(wkt.loads) 
roadX.crs = "EPSG:4326"     
roadX=roadX.reset_index(drop=True)
centerX = roadX.geometry.centroid
centerX.to_crs("EPSG:4326")
dataEmiss4 = ChemicalSpeciationHeavy(roadX,dfSpc,smm)
prefix = 'Heavy'
#splitnetCDFfiles(dataEmiss4,centerX,xX,yY,year,month,prefix,folder)


#TOTAL EMISSIONS
dataEmiss = dataEmiss1+dataEmiss2+dataEmiss3+dataEmiss4
prefix='Total'
splitnetCDFfiles(dataEmiss,centerX,xX,yY,year,month,prefix,folder)
#%% Calling functions

# # Creating annual basis inventory - results in grams per year
# name = 'BRAVESannualEmiss_'+str(year)+'_'+str(month)+'.nc'
# dataRoad = populatingGrid(dataEmissRoad,centerRoad,xX,yY)
# createNETCDF(folder,name,dataRoad,xv,yv,lat,lon,centerRoad,year,month)

# # Creating hourly basis inventory - results in grams per seconds
# name = 'BRAVEStemporalEmiss_'+str(year)+'_'+str(month)+'.nc'
# dataMat,datePfct,disvec = temporalDisagVehicular(dataEmissRoad, year, month)
# dataTempo = populatingGridMat(dataMat,centerRoad,xX,yY)
# datePfct = disvec
# createNETCDFtemporal(folder,name,dataTempo,xv,yv,lat,lon,centerRoad,disvec,month)
    

     
     
