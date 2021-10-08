# -*- coding: utf-8 -*-
"""
---------------------------------------------------------------------------------
                               roadDensity_v1.py

This function calculates the road density using openstretmap shapefiles, Brazilian 
municipalities shapefile, and a user defined grid domain and resolution. The road
density is estimated in each cell, considering the total road length of the city
which contains the respective cell. Openstreetmap files could be download at 
https://download.geofabrik.de/south-america/brazil.html. 

INPUTS:

    dirPath = path with openstreetmap shapefiles, Brazil territory shapefile,
                and Brazilian municipalities shapefiles.
    
    outPath = path to outputs of this function

    lati = initial (lower-left) latitude in degrees

    latf = final (upper-right) longitude in degrees

    loni = initial (lower-left) longitude in degrees

    lonf = initial (upper-right) longitude in degrees

    deltaX = grid spacing in X direction in degrees

    deltaY = grid spacing in Y direction in degrees

    IBGE_CODES = Id codes from Brazilian states
    
    roadFileName = specific name of roads shapefile from Santa Catarina state

OUTPUTS:

    roadDensity_UF_'+str(IBGE_CODES[pp])+'.csv = csv file with road density by state.

    baseGrid.csv = grid used to calculate the road density.


Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br
---------------------------------------------------------------------------------
"""
# Importing packages
import geopandas as gpd
import pandas as pd
#import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString
import numpy as np
import warnings; warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)


def roadDensity (dirPath,outPath,IBGE_CODES,lati,latf,loni,lonf,
                 deltaX,deltaY,roadFileName,roadDensPrefix):

    # Opening cities shapefile
    citySHP = gpd.read_file(dirPath+"/citySHP.shp")
    citySHP.columns = ['NM_MUNICIP','CD_GEOCMU','YEAR',
                         'UF','geometry']

    #==========================START PROCESSING==================================
    #Open shapefile with Brazilian states
    print('===================STARTING roadDensity_v1.py=======================')
    print('')
    print('Author: Leonardo Hoinaski')
    print('Last update: 30/09/2021')
    print(' ')
    print(' ')
    print('Opening shapefiles ')

    # -------------------------- Domain borders ----------------------------------
    # Grind covering entire country 
    print('Setting domain borders')
    x = np.linspace(loni, lonf, int((lonf-loni)/deltaX))
    y = np.linspace(lati, latf, int((latf-lati)/deltaY))

    #Loop over each cel in x direction
    polygons=[]
    for ii in range(1,x.shape[0]):
        #Loop over each cel in y direction
        for jj in range(1,y.shape[0]):
            roadClip=[]
            lat_point_list = [y[jj-1], y[jj], y[jj], y[jj-1]]
            lon_point_list = [x[ii-1], x[ii-1], x[ii], x[ii]]
            cel = Polygon(zip(lon_point_list, lat_point_list))
            polygons.append(cel)

    baseGrid = gpd.GeoDataFrame({'geometry':polygons})
    baseGrid.to_csv(outPath+'/baseGrid_'+roadDensPrefix+'.csv')


    #========================== Selecting cities ==========================
    #for IBGE_CODE in IBGE_CODES:
    for pp in range (0,np.size(IBGE_CODES)):
        print('Selecting states')
        print('Opening roads shapefile')
        if IBGE_CODES[pp]==11 or IBGE_CODES[pp]==12 or IBGE_CODES[pp]==13 or \
            IBGE_CODES[pp]==14 or IBGE_CODES[pp]==15 or IBGE_CODES[pp]==16 or \
                IBGE_CODES[pp]==17:
            roads = gpd.read_file(dirPath+"/norte-latest-free.shp/gis_osm_roads_free_1.shp")
        elif IBGE_CODES[pp]==21 or IBGE_CODES[pp]==22 or IBGE_CODES[pp]==23 or \
            IBGE_CODES[pp]==24 or IBGE_CODES[pp]==25 or IBGE_CODES[pp]==26 or \
                IBGE_CODES[pp]==27 or IBGE_CODES[pp]==28 or IBGE_CODES[pp]==29:
            roads = gpd.read_file(dirPath+"/nordeste-latest-free.shp/gis_osm_roads_free_1.shp")
        elif IBGE_CODES[pp]==31 or IBGE_CODES[pp]==32 or IBGE_CODES[pp]==33 or \
            IBGE_CODES[pp]==34 or IBGE_CODES[pp]==35:
            roads = gpd.read_file(dirPath+"/sudeste-latest-free.shp/gis_osm_roads_free_1.shp")
        elif IBGE_CODES[pp]==42:
            roads = gpd.read_file(dirPath+"/sul-latest-free.shp/"+roadFileName)
        elif IBGE_CODES[pp]==41 or IBGE_CODES[pp]==43:
                roads = gpd.read_file(dirPath+"/sul-latest-free.shp/gis_osm_roads_free_1.shp")
        elif IBGE_CODES[pp]==50 or IBGE_CODES[pp]==51 or IBGE_CODES[pp]==52 or \
            IBGE_CODES[pp]==53:
            roads = gpd.read_file(dirPath+"/centro-oeste-latest-free.shp/gis_osm_roads_free_1.shp")                 
        
        
        print('State number = '+str(IBGE_CODES[pp]))
        shpSC = citySHP[citySHP.iloc[:,3]==IBGE_CODES[pp]]
        scALL = shpSC.dissolve(by='UF')
        scBuffer = scALL.buffer(0.2, resolution=10)
        scBuffer.crs = "EPSG:4326"
        scBuffer.reset_index(drop=True, inplace=True)
        boundSC = scBuffer.bounds
        valL=[]
        print('Selecting cells into state')
        for poly in polygons:
            maxX = max(poly.exterior.coords.xy[0])
            minX = min(poly.exterior.coords.xy[0])
            maxY = max(poly.exterior.coords.xy[1])
            minY = min(poly.exterior.coords.xy[1])
            if maxY<boundSC.maxy[0]:
                if minY>boundSC.miny[0]: 
                    if maxX<boundSC.maxx[0]:
                        if minX>boundSC.minx[0]: 
                            val = True
                        else:
                            val = False
                    else:
                        val = False
                else:
                    val = False
            else:
                val = False
            valL.append(val)
        
        indX = np.arange(np.size(polygons))
        
        arr = np.array([x==True for x in valL])
        polAr = np.array(polygons) 
        polUsed = polAr[arr]
        indXUsed = indX[arr]
        
        
        # ====================== Calculating road density =========================== 
        
        #grid = gpd.GeoDataFrame()
        gridLength = gpd.GeoDataFrame()
        
        #---------------------- Loop over each city
        print('Start processing')
        for kk in range(0,shpSC.shape[0]):
        #for kk in range(0,2):
            # Seting initial values
            roadCity =[]
            cc=1
            latitudes = []
            longitudes = []
            mlsTest =[]
            print('City number = ' + str(kk) + ' of ' + str(shpSC.shape[0]) +
                  '  -'+ shpSC.iloc[kk,0])
            # Cliping roads inside city
            try:
                roadCity = gpd.clip(roads['geometry'], shpSC.iloc[kk,4])
            except:
                roadCity = gpd.clip(roads, shpSC.iloc[kk,4])
                
            # Converting to geodataframe
            if roadCity.shape[0]>0 and roadCity[roadCity!=None].shape[0] >0 and roadCity.empty==False:
                #print('City with roads')
                roadCity = gpd.GeoDataFrame({'geometry':roadCity})
                # Reset index of geodataframe
                roadCity=roadCity.reset_index()    
                #Loop over each road
                roadCityNEW =gpd.GeoDataFrame()
                roadLine=[]    
                for rr in range(0,roadCity.shape[0]):
                    roadLenght = []
                    #polygons = []
                    mlsTest = roadCity.iloc[rr,1].wkt
                    mlsTest = mlsTest.split('(')
                    #print('Extracting coordinates')
                    # Check if geometry is multilinestring               
                    if mlsTest[0] == 'MULTILINESTRING ' or mlsTest[0] == 'MULTILINESTRING' :
                        #roadCity = roadCity.drop(rr)
                        #print('Multilinestring detected')
                        # Replacing multiline for linestring
                        for lin in range(2,len(mlsTest)):
                            mlsTest[lin] = mlsTest[lin].replace(',', '')
                            mlsTest[lin] = mlsTest[lin].replace(')', '')
                            lineStr= mlsTest[2].split(' ')
                            lineStr= np.array(lineStr[0:len(lineStr)-1])
                            lineStr = lineStr.astype(np.float).reshape((-1, 2))
                            lineStr = LineString(lineStr)
                            lineStrGPD = gpd.GeoSeries(lineStr)
                            lineStrGPD = gpd.GeoDataFrame({'geometry':lineStrGPD})
                            roadCityNEW= pd.concat([roadCityNEW, lineStrGPD])
                            # Get coordinates from roadcity to determine road limits
                            coordC = lineStr.coords[:]
                            # Loop over each coordinate to extract latitudes        
                            for coordi in coordC:
                                latitudes.append(coordi[1])
                                longitudes.append(coordi[0])
                    else:             
                        # Get coordinates from roadcity to determine road limits
                        coordC = roadCity.iloc[rr,1].coords[:]
                        roadLine.append(roadCity.iloc[rr,1])
                        # Loop over each coordinate to extract latitudes        
                        for coordi in coordC:
                            latitudes.append(coordi[1])
                            longitudes.append(coordi[0])
                roadLine2 = gpd.GeoDataFrame({'geometry':roadLine}) 
                roadCityNEW= roadCityNEW.reset_index()
                # Check if there is no multilines in city
                if roadCityNEW.empty:
                    roadCity = roadLine2
                else:        
                    roadCityNEW= gpd.GeoDataFrame({'geometry':roadCityNEW.iloc[:,1]}) 
                    roadCity = pd.concat([roadLine2,roadCityNEW])            
                # Converting roadCity to geodataframe 
                roadCity = gpd.GeoDataFrame({'geometry':roadCity['geometry']})       
                roadCity= roadCity.reset_index()
                # Converting latitudes and longitudes to numpy array
                latMaxF = np.asarray(latitudes).max()
                latMinF = np.asarray(latitudes).min() 
                lonMaxF = np.asarray(longitudes).max()
                lonMinF = np.asarray(longitudes).min() 
                # Setting initial values
            
                #Loop over each cel in x direction
                for ii in range(0,polUsed.shape[0]):
                    celgpd = gpd.GeoDataFrame(geometry=[polUsed[ii]])
                    roadClip=[]
                    #print('City number = ' + str(kk) + ' of ' + str(shpSC.shape[0]) + '  -'+ shpSC.iloc[kk,0]+ '  cel number = '+ str(cc) + ' of ' +str(polUsed.shape[0]))
                    #print(str(max(lat_point_list)) +' and ' + str(min(lat_point_list)))
                    # Check roads close to domain
                    maxX = max(polUsed[ii].exterior.coords.xy[0])
                    minX = min(polUsed[ii].exterior.coords.xy[0])
                    maxY = max(polUsed[ii].exterior.coords.xy[1])
                    minY = min(polUsed[ii].exterior.coords.xy[1])
                    if (latMaxF+0.5) > maxY and (latMinF-0.5) < minY and (lonMaxF+0.5) > maxX and (lonMinF-0.5) < minX :
                        try:
                            roadClip = gpd.clip(roadCity['geometry'], celgpd)  
                        except:
                            roadClip = gpd.clip(roadCity, celgpd)  
                                              
                        #roadCel.append(roadClip)            
                        sumCel = roadClip.length
                        # Check if value is empty
                        if sumCel.empty==False:
                            roadLenght.append(sumCel.to_numpy())
                            #print('ROAD INSIDE City number = ' + str(kk) + ' of ' + str(shpSC.shape[0]) )
                        else:
                            roadLenght.append(float(0))
                    # Values is empty        
                    else:
                        roadLenght.append(float(0))
                    
                    cc=cc+1
            else:
                roadLenght = np.zeros(polUsed.shape[0])
                #print('City with roads')
            # Appending road density
            print('APPENDING ROAD DENSITY')
            rdensity=[]       
            for rl in range(0,np.size(roadLenght)):
                rdensity.append(np.sum(roadLenght[rl]))   
            # Appending values of roaddensity to geodataframe
            #grid['geometry'] = polUsed
            #grid[shpSC.iloc[kk,1]] = rdensity/sum(rdensity)
            gridLength['index'] = indXUsed
            gridLength['geometry'] = polUsed
            gridLength[shpSC.iloc[kk,1]] = rdensity
            
        
        gridLength.to_csv(outPath+'/roadDensity_'+roadDensPrefix+'UF_'+str(IBGE_CODES[pp])+'.csv')
        return gridLength







