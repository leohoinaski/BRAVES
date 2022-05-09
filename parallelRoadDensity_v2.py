
# -*- coding: utf-8 -*-
"""
---------------------------------------------------------------------------------
                        parallelRoadDensity_v2.py

This function calculates the road density using openstretmap shapefiles, Brazilian 
municipalities shapefile, and a user defined grid domain and resolution. The road
density is estimated in each cell, considering the total road length of the city
which contains the respective cell. Openstreetmap files could be download at 
https://download.geofabrik.de/south-america/brazil.html. This is a parallel version
of road_density_v1.py. This function uses

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
import multiprocessing as mp
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
import shapely.wkt
#import os

#%%
def roadDensCity (shpSC,folder,polUsed,indXUsed): 
    # kk in range(0,shpSC.shape[0]):
    # Seting initial values
    #print('Start serial processing - number of cities small than 2')
    
    gridLength = gpd.GeoDataFrame()
    for kk in range(0,shpSC.shape[0]): 
        roads = pd.read_csv(folder+'/roadCity_'+shpSC['CD_GEOCMU'].values[kk]+'.csv')
        roads = roads.drop(roads.columns[0], axis=1)
        roads['geometry'] = roads['geometry'].map(shapely.wkt.loads)
        roads = gpd.GeoDataFrame(roads, crs="EPSG:4326", geometry=roads['geometry'])
        roadCity =[]
        cc=1
        latitudes = []
        longitudes = []
        mlsTest =[]
        
        if shpSC.shape[0]==1:
            roadCity = roads.copy()
            print('You have already cliped')
            
        else:
            print('City number = ' + str(kk) + ' of ' + str(shpSC.shape[0]) +
                  '  -'+ shpSC.iloc[kk,0])
            #Cliping roads inside city
            print('Cliping roads inside city')
            try:
                roadCity = gpd.clip(roads['geometry'], shpSC.iloc[kk,4])
            except:
                roadCity = gpd.clip(roads, shpSC.iloc[kk,4])
            print('Roads inside city already cliped')
            
        # Converting to geodataframe
        if roadCity.shape[0]>0 and roadCity[roadCity!=None].shape[0] >0 and roadCity.empty==False:
            #print('City with roads')
            if shpSC.shape[0]>1:
                roadCity = gpd.GeoDataFrame({'geometry':roadCity})
                roadCity=roadCity.reset_index()  
                # Reset index of geodataframe
            else:
                roadCity=roadCity.reset_index()    
            #Loop over each road
            roadCityNEW =gpd.GeoDataFrame()
            roadLine=[]    
            for rr in range(0,roadCity.shape[0]):
                roadLenght = []
                #polygons = []
                
                mlsTest = roadCity['geometry'][rr].wkt
                mlsTest = mlsTest.split('(')
                #print('Extracting coordinates')
                # Check if geometry is multilinestring               
                if mlsTest[0] == 'MULTILINESTRING ' or mlsTest[0] == 'MULTILINESTRING' :
                    #roadCity = roadCity.drop(r
                    #print('Multilinestring detected')
                    # Replacing multiline for linestring
                    for lin in range(2,len(mlsTest)):
                        mlsTest[lin] = mlsTest[lin].replace(',', '')
                        mlsTest[lin] = mlsTest[lin].replace(')', '')
                        lineStr= mlsTest[2].split(' ')
                        lineStr= np.array(lineStr[0:len(lineStr)-1])
                        lineStr = lineStr.astype(np.float64).reshape((-1, 2))
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
                    coordC = roadCity['geometry'][rr].coords[:]
                    roadLine.append(roadCity['geometry'][rr])
                    # Loop over each coordinate to extract latitudes        
                    for coordi in coordC:
                        latitudes.append(coordi[1])
                        longitudes.append(coordi[0])
            roadLine2 = gpd.GeoDataFrame({'geometry':roadLine}) 
            roadCityNEW= roadCityNEW.reset_index()
            # Check if there is no multilines in city
            if roadCityNEW.empty:
                roadCity = roadLine2.copy()
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
        try: 
            del rdensity,roadLenght,sumCel,roadClip,roadCity,roadCityNEW,\
            roadLine,roadLine2
        except:
            rdensity=[]
            roadLenght=[]
            sumCel=[]
            roadClip=[]
            roadCity=[]
            roadCityNEW=[]
            roadLine=[]
            roadLine2=[]
    return gridLength 



#%% Main function

def roadDensity (dirPath,outPath,IBGE_CODES,lati,latf,loni,lonf,
                 deltaX,deltaY,roadFileName,roadDensPrefix):

    # Opening cities shapefile
    citySHP = gpd.read_file(dirPath+"/citySHP.shp")
    citySHP.columns = ['NM_MUNICIP','CD_GEOCMU','YEAR',
                         'UF','geometry']
    
    
    #==========================START PROCESSING==================================
    #Open shapefile with Brazilian states
    print('===================STARTING roadDensity_v2.py=======================')
    print('')
    print('Author: Leonardo Hoinaski')
    print('Last update: 25/02/2022')
    print(' ')
    print(' ')
    print('Opening shapefiles ')

    # -------------------------- Domain borders ----------------------------------
    # Grind covering entire country 
    print('Setting domain borders')
    x = np.arange(loni, lonf+2*deltaX, deltaX)
    y = np.arange(lati, latf+2*deltaY, deltaY)

    #Loop over each cel in x direction
    polygons=[]
    for ii in range(1,x.shape[0]):
        #Loop over each cel in y direction
        for jj in range(1,y.shape[0]):
            #roadClip=[]
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
            #roads = gpd.read_file(dirPath+"/norte-latest-free.shp/gis_osm_roads_free_1.shp")
            folder = dirPath+'/norte-latest-free.shp'
        elif IBGE_CODES[pp]==21 or IBGE_CODES[pp]==22 or IBGE_CODES[pp]==23 or \
            IBGE_CODES[pp]==24 or IBGE_CODES[pp]==25 or IBGE_CODES[pp]==26 or \
                IBGE_CODES[pp]==27 or IBGE_CODES[pp]==28 or IBGE_CODES[pp]==29:
            #roads = gpd.read_file(dirPath+"/nordeste-latest-free.shp/gis_osm_roads_free_1.shp")
            folder = dirPath+'/nordeste-latest-free.shp'
        elif IBGE_CODES[pp]==31 or IBGE_CODES[pp]==32 or IBGE_CODES[pp]==33 or \
            IBGE_CODES[pp]==34 or IBGE_CODES[pp]==35:
            #roads = gpd.read_file(dirPath+"/sudeste-latest-free.shp/gis_osm_roads_free_1.shp")
            folder = dirPath+'/sudeste-latest-free.shp'
        elif IBGE_CODES[pp]==42:
            #roads = gpd.read_file(dirPath+"/sul-latest-free.shp/"+roadFileName)
            folder = dirPath+'/sul-latest-free.shp'
        elif IBGE_CODES[pp]==41 or IBGE_CODES[pp]==43:
            #roads = gpd.read_file(dirPath+"/sul-latest-free.shp/gis_osm_roads_free_1.shp")
            folder = dirPath+'/sul-latest-free.shp'
        elif IBGE_CODES[pp]==50 or IBGE_CODES[pp]==51 or IBGE_CODES[pp]==52 or IBGE_CODES[pp]==53:
            #roads = gpd.read_file(dirPath+"/centro-oeste-latest-free.shp/gis_osm_roads_free_1.shp")
            folder = dirPath+'/centro-oeste-latest-free.shp'                 
    
        
        print('State number = '+str(IBGE_CODES[pp]))
        
        # if os.path.isfile(folder+'/roadsUF_'+str(int(IBGE_CODES[pp]))+'.csv'):
        #     print('Using cliped roads in states')
        #     roads = pd.read_csv(folder+'/roadsUF_'+str(int(IBGE_CODES[pp]))+'.csv')
        #     roads = roads.drop(roads.columns[0], axis=1)
        #     roads['geometry'] = roads['geometry'].map(shapely.wkt.loads)
        #     roads = gpd.GeoDataFrame(roads, crs="EPSG:4326", geometry=roads['geometry'])
        
        # else:
        #     # Opening state shapefile
        #     print('Using roads in region - it might take a bit longer')
        #     brSHP = gpd.read_file(dirPath+"/Brasil.shp")
        #     roads = gpd.read_file(folder+"/gis_osm_roads_free_1.shp")           
        #     ufSHP = brSHP[brSHP['COD_UF']==IBGE_CODES[pp]]
        #     ufSHP = brSHP[brSHP['COD_UF']==21]
        #     ufSHP = ufSHP.dissolve(by='UF')
        #     ufSHP.crs = "EPSG:4326"
        #     ufSHP = ufSHP.buffer(0.2, resolution=10)
        #     print('Cliping roads into State number ' + str(pp) + ' - ' + brSHP[brSHP['COD_UF']==IBGE_CODES[pp]]['UF'].to_numpy()[0])
        #     #Cliping roads inside city
        #     try:
        #         roads = gpd.clip(roads['geometry'], ufSHP)
        #     except:
        #         roads = gpd.clip(roads, ufSHP)

      
        shpSC = citySHP[citySHP.iloc[:,3]==IBGE_CODES[pp]]
        scALL = shpSC.dissolve(by='UF')
        scALL.crs = "EPSG:4326"
        scBuffer = scALL.buffer(0.2, resolution=10)
        scBuffer.crs = "EPSG:4326"
        scBuffer.reset_index(drop=True, inplace=True)
        boundSC = scBuffer.bounds
        
        valL=[]
        #tuple_list=[]     
        print('Selecting cells into state')
        for poly in polygons:
            maxX = max(poly.exterior.coords.xy[0])
            minX = min(poly.exterior.coords.xy[0])
            maxY = max(poly.exterior.coords.xy[1])
            minY = min(poly.exterior.coords.xy[1])
        #     tuple_list.append([minX,maxX,minY,maxY])      
        # val =  [_ for _ in tuple_list if (_[0] > boundSC.minx[0]) and
        #        (_[1] < boundSC.maxx[0]) and 
        #        (_[2] > boundSC.miny[0]) and 
        #        (_[3] < boundSC.maxy[0])]
           
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
        del polAr,arr,valL,scALL,scBuffer,boundSC
        
        
        # ====================== Calculating road density =========================== 
        
        #grid = gpd.GeoDataFrame()
        gridLength = gpd.GeoDataFrame()
        
        #---------------------- Loop over each city
        if shpSC.shape[0] < 2:

            gridLength = roadDensCity (shpSC,folder,polUsed,indXUsed)
            
        else:
        
            print('Start parallell processing')
            cpus = mp.cpu_count()-2
            #cpus = 4
            cityChunks = np.array_split(shpSC, cpus)
            pool = mp.Pool(processes=cpus)  
            chunk_processes = [pool.apply_async(roadDensCity, 
                                                args=(chunk,folder,polUsed,indXUsed)) for chunk in cityChunks]                    
            #new section
            pool.close()
            pool.join()  
            pool.terminate()
            #end new section
            del polUsed,indXUsed            
            roadDensChunks = [chunk.get() for chunk in chunk_processes]
            #gridLength=roadDensCity (shpSC,roads,polUsed,indXUsed)
            gridLength= baseGrid.copy()
            for gg in range(0,len(roadDensChunks)):
                if roadDensChunks[gg].shape[0]>0:
                   try:    
                      ggdf= roadDensChunks[gg].set_index('index')
                      gridLength = pd.merge(gridLength, ggdf, on="geometry")
                   except:
                      print('No roadDensChunks inside')            
            #gpd.GeoDataFrame(pd.concat(roadDensChunks,axis=1), crs=shpSC.crs)
            del roadDensChunks, chunk_processes   
        try:   
            del polUsed,indXUsed 
        except:
            roads=[]
            polUsed=[]
            indXUsed=[]
                
        gridLength.to_csv(outPath+'/roadDensity_'+roadDensPrefix+'UF_'+str(IBGE_CODES[pp])+'.csv')
    return gridLength







