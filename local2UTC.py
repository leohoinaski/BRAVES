#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 20:38:04 2021

@author: leohoinaski
"""

from timezonefinder import TimezoneFinder
import pandas as pd
import numpy as np
#import multiprocessing as mp

def local2UTC(xX,yY):
    lc2utc = np.zeros([xX.shape[1],xX.shape[0]])
    test_naive = pd.date_range('2019-01-01', '2019-04-07', freq='4H')
    tf = TimezoneFinder(in_memory=True)
    
    ltz0 = tf.timezone_at(lng=xX[0,0], lat=yY[0,0])
    ltz0 = float(test_naive.tz_localize(ltz0).strftime('%Z')[-1])
    ltzn = tf.timezone_at(lng=xX[0,-1], lat=yY[0,0])
    ltzn = float(test_naive.tz_localize(ltzn).strftime('%Z')[-1])
    
    if ltz0==ltzn: 
        #lc2utc=ltz0
        lc2utc =np.ones([xX.shape[1],xX.shape[0]])*ltz0
        tag = 1
        
    else: 
        for ii in range(0,xX.shape[1]):       
            print("Longitude " +str(ii))
            local_time_zone0 = tf.timezone_at(lng=xX[0,ii], lat=yY[0,ii])
            t0 = float(test_naive.tz_localize(local_time_zone0).strftime('%Z')[-1])
            local_time_zone5 = tf.timezone_at(lng=xX[round(xX.shape[0]/2),ii], 
                                              lat=yY[round(xX.shape[0]/2),ii])
            t1= float(test_naive.tz_localize(local_time_zone5).strftime('%Z')[-1])
            local_time_zone10 = tf.timezone_at(lng=xX[-1,ii], lat=yY[-1,ii])
            t2= float(test_naive.tz_localize(local_time_zone10).strftime('%Z')[-1])
            
            if (t0==t1) and (t1==t2):                      
                lc2utc[:,ii] = float(test_naive.tz_localize(local_time_zone0).strftime('%Z')[-1])
                
            else:
                for jj in range(0,xX.shape[1]):
                    local_time_zone = tf.timezone_at(lng=xX[jj,ii], lat=yY[jj,ii])
                    lc2utc[jj,ii] = float(test_naive.tz_localize(local_time_zone).strftime('%Z')[-1])
            #lc2utc[ii,jj] = float(re.split('GMT',local_time_zone)[1])
        tag=0
    return lc2utc, tag


#lc2utc = local2UTC(xv,yv)
#np.roll(hourdis,-3)

# cpus = mp.cpu_count()
# xChunks = np.array_split(yv, cpus)
# pool = mp.Pool(processes=cpus)  
# chunk_processes = [pool.apply_async(local2UTC, 
#                                     args=(xv,chunk)) for chunk in xChunks]        
# utcChunks = [chunk.get() for chunk in chunk_processes]
# #gridLength=roadDensCity (shpSC,roads,polUsed,indXUsed)


# gridLength= baseGrid
# for gg in range(0,len(roadDensChunks)):
#     ggdf= roadDensChunks[gg].set_index('index')
#     gridLength = pd.merge(gridLength, ggdf, on="geometry")