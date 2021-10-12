#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 20:38:04 2021

@author: leohoinaski
"""

from timezonefinder import TimezoneFinder
import pandas as pd
import numpy as np


def local2UTC(xv,yv):
    lc2utc = np.zeros([xv.shape[0],xv.shape[1]])
    test_naive = pd.date_range('2019-01-01', '2019-04-07', freq='4H')
    tf = TimezoneFinder(in_memory=True)
    for ii in range(0,xv.shape[0]):
        for jj in range(0,xv.shape[1]):
            local_time_zone = tf.timezone_at(lng=xv[ii,jj], lat=yv[ii,jj])
            lc2utc[ii,jj] = float(test_naive.tz_localize(local_time_zone).strftime('%Z')[-1])
            #lc2utc[ii,jj] = float(re.split('GMT',local_time_zone)[1])
