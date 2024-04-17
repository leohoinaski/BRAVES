#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:29:31 2024

https://repositorio.ufpe.br/bitstream/123456789/6883/1/arquivo8070_1.pdf
# A mistura do etanol passa de 22% para 27%, podendo até chegar até 35%
# the RVP of gasolines ranges from 7 to 13
# observou-se que uma mistura com 10% de etanol aumenta a RVP
aproximadamente em 0,8 psi (5,61 kPa), e a fase vapor contém aproximadamente
15% de etanol. 

RVP = https://d35t1syewk4d42.cloudfront.net/file/1410/RVP-Effects-Memo_03_26_12_Final.pdf
# Extrai os dados de RVP do gráfico deste artigo

@author: leohoinaski
"""
import numpy as np
import matplotlib.pyplot as plt 

tamb = (45 - 5)*np.random.rand(24*30) + 5

tConv = (tamb*(9/5)) + 32


rvp = 9 # RVP = Reid Vapor Pressure (psi)

#Td = Dispensed gasoline temperature (degF) = 20.30+0.81*Tamb  
#from MOVES https://www.epa.gov/sites/default/files/2020-11/documents/420r20012.pdf
td = 20.30+0.81*tConv 

#dT = 0.418*TDF -16.6 # From MOVES -
#https://www.epa.gov/sites/default/files/2020-11/documents/420r20012.pdf
deltaT = 0.418*td -16.6   


er = 264.2*(-5.909 - 0.0949*deltaT + 0.084*td + 0.485*rvp)

plt.plot(er)
