#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                            netCDFcreator_v1.py

This function creates the netCDF files ready to use in CMAQ.

Inputs:
    
    folder: folter to output files
    
    name: output names
    
    data: matrix with data ready to convert in netCDF
    
    xv, yv: meshgrid outputs - grid definition
    
    center: pixel centroids 
    
    year: respective years of emission inventories
    
    month: respective month of emission inventories

Outputs:
        
    netdCDF files
      

Last update = 29/09/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br
-------------------------------------------------------------------------------
"""
import netCDF4 as nc4
#from pyproj import Proj, transform
import numpy as np
import datetime
import pandas as pd
#import tarfile
#import os




#%% Creating a dataset
def createNETCDFannual(folder,name,data,xX,yY,dates,area,ltz):
    print('===================STARTING netCDFcreator_v1.py=======================')
    cdate = datetime.datetime.now()
    cdateStr = int(str(cdate.year)+str(cdate.timetuple().tm_yday))
    ctime = int(str(cdate.hour)+str(cdate.minute)+str(cdate.second))
    # if month<10:
    #     sdate = int(str(dates['year'][0])+'00'+str(month)) 
    # else:
    #     sdate = int(str(dates['year'][0])+'0'+str(month)) 
    #dates['TSTEP'] = np.array(range(0,dates.shape[0]))
    #datesString = dates.index.strftime('%Y-%m-%d %H:%M:%S')
    #datesUsed= np.array([np.array(dates['TSTEP'].astype(str)),datesString]).transpose()
    #datesString= np.repeat(datesUsed[:, :, np.newaxis], data.shape[1], axis=2)
    tflag = np.empty([dates.shape[0],data.shape[1],2],dtype='i4')
    for ii in range(0,dates.shape[0]):
        tflag[ii,:,0]=int(dates['year'][0]*1000 + dates.index[ii].timetuple().tm_yday)
        tflag[ii,:,1]=int(str(dates['hour'][ii])+'0000')
    
    sdate =  dates['year'][0]*1000 + dates.index[0].timetuple().tm_yday            

    # Reshape area array
    area1 = area.reshape((np.size(xX,1),np.size(yY,0))).transpose()
    ltz = ltz.reshape((np.size(xX,1),np.size(yY,0))).transpose()
    
    f2 = nc4.Dataset(folder+'/'+name,'w', format='NETCDF4_CLASSIC') #'w' stands for write    
    #Add global attributes
    #f2.IOAPI_VERSION ='$Id: @(#) ioapi library version 3.1 $'
    #f2.EXEC_ID = '???????????????'
    # f2.FTYPE =  1
    f2.CDATE= cdateStr
    f2.CTIME= ctime
    f2.WDATE= cdateStr
    f2.WTIME= ctime
    f2.SDATE= sdate
    f2.STIME= 0
    f2.TSTEP= 10000
    f2.NTHIK= 1
    f2.NCOLS= data.shape[3]
    print('NCOLS=' +str(data.shape[3]))
    f2.NROWS= data.shape[2]
    print('NROWS=' +str(data.shape[2]))
    f2.NLAYS= 1
    f2.NVARS= data.shape[1]
    # f2.GDTYP= 1
    f2.P_ALP= -10
    f2.P_BET= 0
    f2.P_GAM= round(xX.mean(),6)
    f2.XCENT= round(xX.mean(),6)
    print('XCENT='+ str(round(xX.mean(),6)))
    f2.YCENT= round(yY.mean(),6)
    print('YCENT='+ str(round(yY.mean(),6)))
    f2.XORIG= np.round((xX.min() - (xX[0,1] - xX[0,0])/2),6)
    print('XORIG = ' + str(np.round((xX.min() - (xX[0,1] - xX[0,0])/2),6))) 
    f2.YORIG= np.round((yY.min() - (yY[1,0] - yY[0,0])/2),4)
    print('YORIG = ' + str(np.round((yY.min() - (yY[1,0] - yY[0,0])/2),6))) 
    f2.XCELL= round((xX[0,1] - xX[0,0]), 6)
    print('XCELL = '+str(round(xX[0,1] - xX[0,0],6)))
    f2.YCELL= round((yY[1,0] - yY[0,0]), 6)
    print('YCELL = '+str(round(yY[1,0] - yY[0,0],6)))
    # f2.VGTYP= -1
    # f2.VGTOP= 0.0
    # f2.VGLVLS= [0,0]
    
    #f2.GDNAM= 'SE53BENCH'       
    #f2.UPNAM= 'M3WNDW'   
    strVAR = ' ACET            ACROLEIN        ALD2           \n\
    ALD2_PRIMARY    ALDX            BENZ            BUTADIENE13\n\
    CH4             CH4_INV         CL2             CO\n\
    CO2_INV         ETH             ETHA            ETHY\n\
    ETOH            FORM            FORM_PRIMARY    HCL\n\
    HONO            IOLE            ISOP            KET\n\
    MEOH            N2O_INV         NAPH            NH3\n\
    NH3_FERT        NO              NO2             NVOL\n\
    OLE             PAL             PAR             PCA\n\
    PCL             PEC             PFE             PH2O\n\
    PK              PMC             PMG             PMN\n\
    PMOTHR          PNA             PNCOM           PNH4\n\
    PNO3            POC             PRPA            PSI\n\
    PSO4            PTI             SO2             SOAALK\n\
    SULF            TERP            TOL             UNK\n\
    UNR             VOC_INV         XYLMN           PMFINE'        
    f2.VAR_LIST=strVAR
    f2.FILEDESC= 'BRAVES database ANNUAL vehicular emissions'
    f2.HISTORY ='' 
       
    # # Specifying dimensions
    #tempgrp = f.createGroup('vehicularEmissions_data')
    f2.createDimension('TSTEP', dates.shape[0])
    f2.createDimension('DATE-TIME', 2)
    f2.createDimension('LAY', 1)
    f2.createDimension('VAR', data.shape[1])
    f2.createDimension('ROW', np.size(yY,0))
    f2.createDimension('COL', np.size(xX,1))
    
    # Building variables
    TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    LON = f2.createVariable('Longitude', 'f4', ( 'ROW','COL'))
    LAT = f2.createVariable('Latitude', 'f4', ( 'ROW','COL'))
    AREA = f2.createVariable('AREA', 'f4', ( 'ROW','COL'))
    LTZ = f2.createVariable('LTZ', 'f4', ( 'ROW','COL'))
    #TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP'))
    ACET = f2.createVariable('ACET', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ACROLEIN = f2.createVariable('ACROLEIN', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ALD2 = f2.createVariable('ALD2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ALD2_PRIMARY = f2.createVariable('ALD2_PRIMARY', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ALDX = f2.createVariable('ALDX', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    BENZ = f2.createVariable('BENZ', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    BUTADIENE13 = f2.createVariable('BUTADIENE13', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    CH4 = f2.createVariable('CH4', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    CH4_INV = f2.createVariable('CH4_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    CL2 = f2.createVariable('CL2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    CO = f2.createVariable('CO', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    CO2_INV = f2.createVariable('CO2_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ETH = f2.createVariable('ETH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ETHA = f2.createVariable('ETHA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ETHY = f2.createVariable('ETHY', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ETOH = f2.createVariable('ETOH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    FORM = f2.createVariable('FORM', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    FORM_PRIMARY = f2.createVariable('FORM_PRIMARY', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    HCL = f2.createVariable('HCL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    HONO = f2.createVariable('HONO', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    IOLE = f2.createVariable('IOLE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ISOP = f2.createVariable('ISOP', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    KET = f2.createVariable('KET', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    MEOH = f2.createVariable('MEOH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    N2O_INV = f2.createVariable('N2O_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NAPH = f2.createVariable('NAPH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NH3 = f2.createVariable('NH3', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NH3_FERT = f2.createVariable('NH3_FERT', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NO = f2.createVariable('NO', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NO2 = f2.createVariable('NO2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NVOL = f2.createVariable('NVOL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    OLE = f2.createVariable('OLE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PAL = f2.createVariable('PAL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PAR = f2.createVariable('PAR', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PCA = f2.createVariable('PCA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PCL = f2.createVariable('PCL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PEC = f2.createVariable('PEC', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PFE = f2.createVariable('PFE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PH2O = f2.createVariable('PH2O', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PK = f2.createVariable('PK', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PMC = f2.createVariable('PMC', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PMG = f2.createVariable('PMG', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PMN = f2.createVariable('PMN', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PMOTHR = f2.createVariable('PMOTHR', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PNA = f2.createVariable('PNA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PNCOM = f2.createVariable('PNCOM', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PNH4 = f2.createVariable('PNH4', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PNO3 = f2.createVariable('PNO3', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    POC = f2.createVariable('POC', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PRPA = f2.createVariable('PRPA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PSI = f2.createVariable('PSI', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PSO4 = f2.createVariable('PSO4', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PTI = f2.createVariable('PTI', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    SO2 = f2.createVariable('SO2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    SOAALK = f2.createVariable('SOAALK', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    SULF = f2.createVariable('SULF', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    TERP = f2.createVariable('TERP', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    TOL = f2.createVariable('TOL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    UNK = f2.createVariable('UNK', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    UNR = f2.createVariable('UNR', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    VOC_INV = f2.createVariable('VOC_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    XYLMN = f2.createVariable('XYLMN', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PMFINE = f2.createVariable('PMFINE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))

    
    # Passing data into variables
    TFLAG[:,:,:] = tflag
    #TFLAG[:] = dates.index
    print(yY.shape)
    print(xX.shape)
    LAT[:,:] =  yY
    LON[:,:] = xX
    print(area1.shape)
    AREA[:,:] = area1
    LTZ[:,:] = ltz
    print(str(ACET.shape) + str(data.shape))
    ACET[:,:,:,:] =  data[:,0,:,:]
    ACROLEIN[:,:,:,:] = data[:,1,:,:]
    ALD2[:,:,:,:] = data[:,2,:,:]
    ALD2_PRIMARY[:,:,:,:] = data[:,3,:,:]
    ALDX[:,:,:,:] = data[:,4,:,:]
    BENZ[:,:,:,:] = data[:,5,:,:]
    BUTADIENE13[:,:,:,:] = data[:,6,:,:]
    CH4[:,:,:,:] = data[:,7,:,:]
    CH4_INV[:,:,:,:] = data[:,8,:,:]
    CL2[:,:,:,:] =  data[:,9,:,:]
    CO[:,:,:,:] = data[:,10,:,:]
    CO2_INV[:,:,:,:] =  data[:,11,:,:]
    ETH[:,:,:,:] =  data[:,12,:,:]
    ETHA[:,:,:,:] =  data[:,13,:,:]
    ETHY[:,:,:,:] =  data[:,14,:,:]
    ETOH[:,:,:,:] = data[:,15,:,:]
    FORM[:,:,:,:] = data[:,16,:,:]
    FORM_PRIMARY[:,:,:,:] = data[:,17,:,:]
    HCL[:,:,:,:] = data[:,18,:,:]
    HONO[:,:,:,:] = data[:,19,:,:]
    IOLE[:,:,:,:] = data[:,20,:,:]
    ISOP[:,:,:,:] = data[:,21,:,:]
    KET[:,:,:,:] = data[:,22,:,:]
    MEOH[:,:,:,:] = data[:,23,:,:]
    N2O_INV[:,:,:,:] = data[:,24,:,:]
    NAPH[:,:,:,:] = data[:,25,:,:]
    NH3[:,:,:,:] = data[:,26,:,:]
    NH3_FERT[:,:,:,:] = data[:,27,:,:]
    NO[:,:,:,:] = data[:,28,:,:]
    NO2[:,:,:,:] = data[:,29,:,:]
    NVOL[:,:,:,:] = data[:,30,:,:]
    OLE[:,:,:,:] = data[:,31,:,:]
    PAL[:,:,:,:] = data[:,32,:,:]
    PAR[:,:,:,:] = data[:,33,:,:]
    PCA[:,:,:,:] = data[:,34,:,:]
    PCL[:,:,:,:] = data[:,35,:,:]
    PEC[:,:,:,:] = data[:,36,:,:]
    PFE[:,:,:,:] = data[:,37,:,:]
    PH2O[:,:,:,:] = data[:,38,:,:]
    PK[:,:,:,:] = data[:,39,:,:]
    PMC[:,:,:,:] = data[:,40,:,:]
    PMG[:,:,:,:] = data[:,41,:,:]
    PMN[:,:,:,:] = data[:,42,:,:]
    PMOTHR[:,:,:,:] = data[:,43,:,:]
    PNA[:,:,:,:] = data[:,44,:,:]
    PNCOM[:,:,:,:] = data[:,45,:,:]
    PNH4[:,:,:,:] = data[:,46,:,:]
    PNO3[:,:,:,:] = data[:,47,:,:]
    POC[:,:,:,:] = data[:,48,:,:]
    PRPA[:,:,:,:] = data[:,49,:,:]
    PSI[:,:,:,:] = data[:,50,:,:]
    PSO4[:,:,:,:] = data[:,51,:,:]
    PTI[:,:,:,:] = data[:,52,:,:]
    SO2[:,:,:,:] = data[:,53,:,:]
    SOAALK[:,:,:,:] = data[:,54,:,:]
    SULF[:,:,:,:] = data[:,55,:,:]
    TERP[:,:,:,:] = data[:,56,:,:]
    TOL[:,:,:,:] = data[:,57,:,:]
    UNK[:,:,:,:]= data[:,58,:,:]
    UNR[:,:,:,:] = data[:,59,:,:]
    VOC_INV[:,:,:,:] = data[:,60,:,:]
    XYLMN[:,:,:,:] = data[:,61,:,:]
    PMFINE[:,:,:,:] = data[:,62,:,:]
    
    
    #Add local attributes to variable instances
    TFLAG.units = '<YYYYDDD,HHMMSS>'
    ACET.units = 'moles/year '
    ACROLEIN.units = 'moles/year '
    ALD2.units = 'moles/year '
    ALD2_PRIMARY.units = 'moles/year '
    ALDX.units = 'moles/year '
    BENZ.units = 'moles/year '
    BUTADIENE13.units = 'moles/year '
    CH4.units = 'moles/year '
    CH4_INV.units = 'g/year '
    CL2.units = 'moles/year ' 
    CO.units = 'moles/year '
    CO2_INV.units = 'g/year '
    ETH.units = 'moles/year '
    ETHA.units = 'moles/year '
    ETHY.units = 'moles/year '
    ETOH.units = 'moles/year '
    FORM.units = 'moles/year '
    FORM_PRIMARY.units = 'moles/year '
    HCL.units = 'moles/year '
    HONO.units = 'moles/year '
    IOLE.units = 'moles/year '
    ISOP.units = 'moles/year '
    KET.units = 'moles/year '
    MEOH.units = 'moles/year '
    N2O_INV.units = 'g/year '
    NAPH.units = 'moles/year '
    NH3.units = 'moles/year '
    NH3_FERT.units = 'moles/year '
    NO.units = 'moles/year '
    NO2.units = 'moles/year '
    NVOL.units = 'moles/year '
    OLE.units = 'moles/year '
    PAL.units = 'moles/year '
    PAR.units = 'moles/year '
    PCA.units = 'g/year '
    PCL.units = 'g/year '
    PEC.units = 'g/year '
    PFE.units = 'g/year '
    PH2O.units = 'g/year '
    PK.units = 'g/year '
    PMC.units = 'g/year ' 
    PMG.units = 'g/year '
    PMN.units = 'g/year ' 
    PMOTHR.units = 'g/year ' 
    PNA.units = 'g/year ' 
    PNCOM.units = 'g/year ' 
    PNH4.units = 'g/year ' 
    PNO3.units = 'g/year ' 
    POC.units = 'g/year ' 
    PRPA.units = 'moles/year '
    PSI.units = 'g/year ' 
    PSO4.units = 'g/year ' 
    PTI.units = 'g/year ' 
    SO2.units = 'moles/year ' 
    SOAALK.units = 'moles/year ' 
    SULF.units = 'moles/year ' 
    TERP.units = 'moles/year ' 
    TOL.units = 'moles/year ' 
    UNK.units = 'moles/year ' 
    UNR.units = 'moles/year ' 
    VOC_INV.units = 'g/year ' 
    XYLMN.units = 'moles/year '
    PMFINE.units = 'g/year '
    LON.units = 'degrees '
    LAT.units = 'degrees '
    AREA.units = 'km2 '
    LTZ.units = 'hours '
    
    f2.close()
    print('Your annual BRAVESdatabase netCDF file is ready!')
    
    #%% Creating a dataset
def createNETCDFtemporalfromNC(folder,name,data,xX,yY,mcipPath):
    print('===================STARTING netCDFcreator_v1.py=======================')
    # cdate = datetime.datetime.now()
    # cdateStr = int(str(cdate.year)+str(cdate.timetuple().tm_yday))
    # ctime = int(str(cdate.hour)+str(cdate.minute)+str(cdate.second))

    # tflag = np.empty([dates.shape[0],data.shape[1],2],dtype='i4')
    # for ii in range(0,dates.shape[0]):
    #     tflag[ii,:,0]=int(dates['year'][0]*1000 + dates.index[ii].timetuple().tm_yday)
    #     tflag[ii,:,1]=int(str(dates['hour'][ii])+'0000')
    
    # sdate =  dates['year'][0]*1000 + dates.index[0].timetuple().tm_yday            
    f2 = nc4.Dataset(folder+'/'+name,'w', format='NETCDF3_CLASSIC') #'w' stands for write   
    #Add global attributes
    f2.IOAPI_VERSION ='$Id: @(#) ioapi library version 3.1 $'
    f2.EXEC_ID = '???????????????'
    f2.FTYPE =  1
    ds3 = nc4.Dataset(mcipPath)
    for gatr in ds3.ncattrs() :
        setattr(f2, gatr, ds3.getncattr(gatr))
    # f2.VGTYP= -9999
    f2.VGTOP= 0.0
    f2.VGLVLS= [0,0]
    f2.NLAYS = 1 
    # f2.NCOLS = f2.NCOLS
    # f2.NROWS = f2.NROWS



    strVAR ='ACET            ACROLEIN        ALD2            ALD2_PRIMARY    ALDX            BENZ            BUTADIENE13     CH4             CH4_INV         CL2             CO              CO2_INV         ETH             ETHA            ETHY            ETOH            FORM            FORM_PRIMARY    HCL             HONO            IOLE            ISOP            KET             MEOH            N2O_INV         NAPH            NH3             NH3_FERT        NO              NO2             NVOL            OLE             PAL             PAR             PCA             PCL             PEC             PFE             PH2O            PK              PMC             PMG             PMN             PMOTHR          PNA             PNCOM           PNH4            PNO3            POC             PRPA            PSI             PSO4            PTI             SO2             SOAALK          SULF            TERP            TOL             UNK             UNR             VOC_INV         XYLMN           PMFINE          '
    setattr(f2, 'VAR-LIST', strVAR)
    f2.FILEDESC= 'BRAVES database vehicular emissions ready for CMAQ'
    f2.HISTORY =''
    # # Specifying dimensions
    #tempgrp = f.createGroup('vehicularEmissions_data')
    f2.createDimension('TSTEP', None )
    f2.createDimension('DATE-TIME', 2)
    f2.createDimension('LAY', 1)
    f2.createDimension('VAR', data.shape[1])
    f2.createDimension('ROW', np.size(yY,0))
    f2.createDimension('COL', np.size(xX,1))
    
    # Building variables
    TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))

    ACET = f2.createVariable('ACET', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    ACROLEIN = f2.createVariable('ACROLEIN', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    ALD2 = f2.createVariable('ALD2', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    ALD2_PRIMARY = f2.createVariable('ALD2_PRIMARY', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    ALDX = f2.createVariable('ALDX', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    BENZ = f2.createVariable('BENZ', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    BUTADIENE13 = f2.createVariable('BUTADIENE13', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    CH4 = f2.createVariable('CH4', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    CH4_INV = f2.createVariable('CH4_INV', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    CL2 = f2.createVariable('CL2', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    CO = f2.createVariable('CO', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    CO2_INV = f2.createVariable('CO2_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ETH = f2.createVariable('ETH', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    ETHA = f2.createVariable('ETHA', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    ETHY = f2.createVariable('ETHY', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    ETOH = f2.createVariable('ETOH', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    FORM = f2.createVariable('FORM', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    FORM_PRIMARY = f2.createVariable('FORM_PRIMARY', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    HCL = f2.createVariable('HCL', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    HONO = f2.createVariable('HONO', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    IOLE = f2.createVariable('IOLE', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    ISOP = f2.createVariable('ISOP', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    KET = f2.createVariable('KET', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    MEOH = f2.createVariable('MEOH', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    N2O_INV = f2.createVariable('N2O_INV', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    NAPH = f2.createVariable('NAPH',np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    NH3 = f2.createVariable('NH3', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    NH3_FERT = f2.createVariable('NH3_FERT', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    NO = f2.createVariable('NO', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    NO2 = f2.createVariable('NO2', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    NVOL = f2.createVariable('NVOL', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    OLE = f2.createVariable('OLE', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PAL = f2.createVariable('PAL', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PAR = f2.createVariable('PAR', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PCA = f2.createVariable('PCA', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PCL = f2.createVariable('PCL', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PEC = f2.createVariable('PEC', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PFE = f2.createVariable('PFE', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PH2O = f2.createVariable('PH2O', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PK = f2.createVariable('PK', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PMC = f2.createVariable('PMC', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PMG = f2.createVariable('PMG', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PMN = f2.createVariable('PMN', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PMOTHR = f2.createVariable('PMOTHR', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PNA = f2.createVariable('PNA', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PNCOM = f2.createVariable('PNCOM', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PNH4 = f2.createVariable('PNH4', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PNO3 = f2.createVariable('PNO3', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    POC = f2.createVariable('POC', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PRPA = f2.createVariable('PRPA', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PSI = f2.createVariable('PSI', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PSO4 = f2.createVariable('PSO4', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PTI = f2.createVariable('PTI', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    SO2 = f2.createVariable('SO2', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    SOAALK = f2.createVariable('SOAALK', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    SULF = f2.createVariable('SULF', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    TERP = f2.createVariable('TERP', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    TOL = f2.createVariable('TOL', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    UNK = f2.createVariable('UNK', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    UNR = f2.createVariable('UNR', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    VOC_INV = f2.createVariable('VOC_INV', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    XYLMN = f2.createVariable('XYLMN', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
    PMFINE = f2.createVariable('PMFINE', np.float32, ('TSTEP', 'LAY', 'ROW','COL'))

    
    # Passing data into variables
    TFLAG[:,:,:] = np.repeat(ds3['TFLAG'][:,0,:][:, np.newaxis,:], 
                             data.shape[1], axis=1)
    print(yY.shape)
    print(xX.shape)
    print(str(ACET.shape) + str(data.shape))

    ACET[:,:,:,:] =  data[:,0,:,:]
    ACROLEIN[:,:,:,:] = data[:,1,:,:]
    ALD2[:,:,:,:] = data[:,2,:,:]
    ALD2_PRIMARY[:,:,:,:] = data[:,3,:,:]
    ALDX[:,:,:,:] = data[:,4,:,:]
    BENZ[:,:,:,:] = data[:,5,:,:]
    BUTADIENE13[:,:,:,:] = data[:,6,:,:]
    CH4[:,:,:,:] = data[:,7,:,:]
    CH4_INV[:,:,:,:] = data[:,8,:,:]
    CL2[:,:,:,:] =  data[:,9,:,:]
    CO[:,:,:,:] = data[:,10,:,:]
    CO2_INV[:,:,:,:] =  data[:,11,:,:]
    ETH[:,:,:,:] =  data[:,12,:,:]
    ETHA[:,:,:,:] =  data[:,13,:,:]
    ETHY[:,:,:,:] =  data[:,14,:,:]
    ETOH[:,:,:,:] = data[:,15,:,:]
    FORM[:,:,:,:] = data[:,16,:,:]
    FORM_PRIMARY[:,:,:,:] = data[:,17,:,:]
    HCL[:,:,:,:] = data[:,18,:,:]
    HONO[:,:,:,:] = data[:,19,:,:]
    IOLE[:,:,:,:] = data[:,20,:,:]
    ISOP[:,:,:,:] = data[:,21,:,:]
    KET[:,:,:,:] = data[:,22,:,:]
    MEOH[:,:,:,:] = data[:,23,:,:]
    N2O_INV[:,:,:,:] = data[:,24,:,:]
    NAPH[:,:,:,:] = data[:,25,:,:]
    NH3[:,:,:,:] = data[:,26,:,:]
    NH3_FERT[:,:,:,:] = data[:,27,:,:]
    NO[:,:,:,:] = data[:,28,:,:]
    NO2[:,:,:,:] = data[:,29,:,:]
    NVOL[:,:,:,:] = data[:,30,:,:]
    OLE[:,:,:,:] = data[:,31,:,:]
    PAL[:,:,:,:] = data[:,32,:,:]
    PAR[:,:,:,:] = data[:,33,:,:]
    PCA[:,:,:,:] = data[:,34,:,:]
    PCL[:,:,:,:] = data[:,35,:,:]
    PEC[:,:,:,:] = data[:,36,:,:]
    PFE[:,:,:,:] = data[:,37,:,:]
    PH2O[:,:,:,:] = data[:,38,:,:]
    PK[:,:,:,:] = data[:,39,:,:]
    PMC[:,:,:,:] = data[:,40,:,:]
    PMG[:,:,:,:] = data[:,41,:,:]
    PMN[:,:,:,:] = data[:,42,:,:]
    PMOTHR[:,:,:,:] = data[:,43,:,:]
    PNA[:,:,:,:] = data[:,44,:,:]
    PNCOM[:,:,:,:] = data[:,45,:,:]
    PNH4[:,:,:,:] = data[:,46,:,:]
    PNO3[:,:,:,:] = data[:,47,:,:]
    POC[:,:,:,:] = data[:,48,:,:]
    PRPA[:,:,:,:] = data[:,49,:,:]
    PSI[:,:,:,:] = data[:,50,:,:]
    PSO4[:,:,:,:] = data[:,51,:,:]
    PTI[:,:,:,:] = data[:,52,:,:]
    SO2[:,:,:,:] = data[:,53,:,:]
    SOAALK[:,:,:,:] = data[:,54,:,:]
    SULF[:,:,:,:] = data[:,55,:,:]
    TERP[:,:,:,:] = data[:,56,:,:]
    TOL[:,:,:,:] = data[:,57,:,:]
    UNK[:,:,:,:]= data[:,58,:,:]
    UNR[:,:,:,:] = data[:,59,:,:]
    VOC_INV[:,:,:,:] = data[:,60,:,:]
    XYLMN[:,:,:,:] = data[:,61,:,:]
    PMFINE[:,:,:,:] = data[:,62,:,:]
       
    #Add local attributes to variable instances
    TFLAG.units = '<YYYYDDD,HHMMSS>'
    TFLAG.var_desc = 'Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS'

    ACET.units = 'moles/s '
    ACET.var_desc ='ACET[1]' 
    ACROLEIN.var_desc ='ACROLEIN[1] '
    ACROLEIN.units = 'moles/s '
    ALD2.units = 'moles/s '
    ALD2.var_desc = 'ALD2[1]'
    ALD2_PRIMARY.units = 'moles/s '
    ALD2_PRIMARY.var_desc = 'ALD2_PRIMARY[1]'
    ALDX.units = 'moles/s '
    ALDX.var_desc = 'ALDX[1]'
    BENZ.units = 'moles/s '
    BENZ.var_desc = 'BENZ[1]'
    BUTADIENE13.units = 'moles/s '
    BUTADIENE13.var_desc = 'BUTADIENE13[1]'
    CH4.units = 'moles/s '
    CH4.var_desc = 'CH4[1]'
    CH4_INV.units = 'g/s '
    CH4_INV.var_desc = 'CH4_INV[1]'
    CL2.units = 'moles/s '
    CL2.var_desc = 'CL2[1]' 
    CO.units = 'moles/s '
    CO.var_desc = 'CO[1]'
    CO2_INV.units = 'g/s '
    CO2_INV.var_desc = 'CO2_INV[1]'
    ETH.units = 'moles/s '
    ETH.var_desc = 'ETH[1]'
    ETHA.units = 'moles/s '
    ETHA.var_desc = 'ETHA[1]'
    ETHY.units = 'moles/s '
    ETHY.var_desc = 'ETHY[1]'
    ETOH.units = 'moles/s '
    ETOH.var_desc = 'ETOH[1]'
    FORM.units = 'moles/s '
    FORM.var_desc = 'FORM[1]'
    FORM_PRIMARY.units = 'moles/s '
    FORM_PRIMARY.var_desc = 'FORM_PRIMARY[1]'
    HCL.units = 'moles/s '
    HCL.var_desc = 'HCL[1]'
    HONO.units = 'moles/s '
    HONO.var_desc = 'HONO[1]'
    IOLE.units = 'moles/s '
    IOLE.var_desc = 'IOLE[1]'
    ISOP.units = 'moles/s '
    ISOP.var_desc = 'ISOP[1]'
    KET.units = 'moles/s '
    KET.var_desc = 'KET[1]'
    MEOH.units = 'moles/s '
    MEOH.var_desc = 'MEOH[1]'
    N2O_INV.units = 'g/s '
    N2O_INV.var_desc = 'N2O_INV[1]'
    NAPH.units = 'moles/s '
    NAPH.var_desc = 'NAPH[1]'
    NH3.units = 'moles/s '
    NH3.var_desc = 'NH3[1]'
    NH3_FERT.units = 'moles/s '
    NH3_FERT.var_desc = 'NH3_FERT[1]'
    NO.units = 'moles/s '
    NO.var_desc = 'NO[1]'
    NO2.units = 'moles/s '
    NO2.var_desc = 'NO2[1]'
    NVOL.units = 'moles/s '
    NVOL.var_desc = 'NVOL[1]'
    OLE.units = 'moles/s '
    OLE.var_desc = 'OLE[1]'
    PAL.units = 'moles/s '
    PAL.var_desc = 'PAL[1]'
    PAR.units = 'moles/s '
    PAR.var_desc = 'PAR[1]'
    PCA.units = 'g/s '
    PCA.var_desc = 'PCA[1]'
    PCL.units = 'g/s '
    PCL.var_desc = 'PCL[1]'
    PEC.units = 'g/s '
    PEC.var_desc = 'PEC[1]'
    PFE.units = 'g/s '
    PFE.var_desc = 'PFE[1]'
    PH2O.units = 'g/s '
    PH2O.var_desc = 'PH2O[1]'
    PK.units = 'g/s '
    PK.var_desc = 'PK[1]'
    PMC.units = 'g/s ' 
    PMC.var_desc = 'PMC[1]' 
    PMG.units = 'g/s '
    PMG.var_desc = 'PMG[1]'
    PMN.units = 'g/s ' 
    PMN.var_desc = 'PMN[1]'
    PMOTHR.units = 'g/s ' 
    PMOTHR.var_desc = 'PMOTHR[1]'
    PNA.units = 'g/s ' 
    PNA.var_desc = 'PNA[1]'
    PNCOM.units = 'g/s ' 
    PNCOM.var_desc = 'PNCOM[1]'
    PNH4.units = 'g/s ' 
    PNH4.var_desc = 'PNH4[1]'
    PNO3.units = 'g/s ' 
    PNO3.var_desc = 'PNO3[1]'
    POC.units = 'g/s ' 
    POC.var_desc = 'POC[1]'
    PRPA.units = 'moles/s '
    PRPA.var_desc = 'PRPA[1]'
    PSI.units = 'g/s ' 
    PSI.var_desc = 'PSI[1]'
    PSO4.units = 'g/s ' 
    PSO4.var_desc = 'PSO4[1]'
    PTI.units = 'g/s '
    PTI.var_desc = 'PTI[1]'
    SO2.units = 'moles/s ' 
    SO2.var_desc = 'SO2[1]'
    SOAALK.units = 'moles/s '
    SOAALK.var_desc = 'SOAALK[1]' 
    SULF.units = 'moles/s ' 
    SULF.var_desc = 'SULF[1]'
    TERP.units = 'moles/s '
    TERP.var_desc = 'TERP[1]' 
    TOL.units = 'moles/s '
    TOL.var_desc = 'TOL[1]'
    UNK.units = 'moles/s ' 
    UNK.var_desc = 'UNK[1]'
    UNR.units = 'moles/s ' 
    UNR.var_desc = 'UNR[1]'
    VOC_INV.units = 'g/s ' 
    VOC_INV.var_desc = 'VOC_INV[1]'
    XYLMN.units = 'moles/s '
    XYLMN.var_desc = 'XYLMN[1]'
    PMFINE.units = 'g/s '
    PMFINE.var_desc = 'PMFINE[1]'
   
    f2.close()
    print('Your BRAVESdatabase netCDF file is ready for CMAQ!')

#%%
def createNETCDFtemporalBySpecies(folder,name,data,xX,yY,dates,specie,area):
    cdate = datetime.datetime.now()
    cdateStr = int(str(cdate.year)+str(cdate.timetuple().tm_yday))
    ctime = int(str(cdate.hour)+str(cdate.minute)+str(cdate.second))

    tflag = np.empty([dates.shape[0],data.shape[1],2],dtype='i4')
    for ii in range(0,dates.shape[0]):
        tflag[ii,:,0]=int(dates['year'][0]*1000 + dates.index[ii].timetuple().tm_yday)
        tflag[ii,:,1]=int(str(dates['hour'][ii])+'0000')
    
    sdate =  dates['year'][0]*1000 + dates.index[0].timetuple().tm_yday            
    
    f2 = nc4.Dataset(folder+'/'+name,'w', format='NETCDF4_CLASSIC') #'w' stands for write    
    #Add global attributes
    #f2.IOAPI_VERSION ='$Id: @(#) ioapi library version 3.1 $'
    #f2.EXEC_ID = '???????????????'
    f2.FTYPE =  1
    f2.CDATE= cdateStr
    f2.CTIME= ctime
    f2.WDATE= cdateStr
    f2.WTIME= ctime
    f2.SDATE= sdate
    f2.STIME= 0
    f2.TSTEP= 10000
    f2.NTHIK= 1
    f2.NCOLS= data.shape[3]
    print('NCOLS=' +str(data.shape[3]))
    f2.NROWS= data.shape[2]
    print('NROWS=' +str(data.shape[2]))
    f2.NLAYS= 1
    f2.NVARS= data.shape[1]
    f2.GDTYP= 1
    f2.P_ALP= -10
    f2.P_BET= 0
    f2.P_GAM= round(xX.mean(),6)
    f2.XCENT= round(xX.mean(),6)
    print('XCENT='+ str(round(xX.mean(),6)))
    f2.YCENT= round(yY.mean(),6)
    print('YCENT='+ str(round(yY.mean(),6)))
    f2.XORIG= np.round((xX.min() - (xX[0,1] - xX[0,0])/2),6)
    print('XORIG = ' + str(np.round((xX.min() - (xX[0,1] - xX[0,0])/2),6))) 
    f2.YORIG= np.round((yY.min() - (yY[1,0] - yY[0,0])/2),6)
    print('YORIG = ' + str(np.round((yY.min() - (yY[1,0] - yY[0,0])/2),6))) 
    f2.XCELL= round((xX[0,1] - xX[0,0]), 6)
    print('XCELL = '+str(round(xX[0,1] - xX[0,0],6)))
    f2.YCELL= round((yY[1,0] - yY[0,0]), 6)
    print('YCELL = '+str(round(yY[1,0] - yY[0,0],6)))
    f2.VGTYP= -1
    f2.VGTOP= 0.0
    f2.VGLVLS= [0,0]
    
    #f2.GDNAM= 'SE53BENCH'       
    #f2.UPNAM= 'M3WNDW'   
    strVAR= specie       
    f2.VAR_LIST=strVAR
    f2.FILEDESC= 'BRAVES temporal emissions from 1 specie'
    f2.HISTORY ='' 
       
    # # Specifying dimensions
    #tempgrp = f.createGroup('vehicularEmissions_data')
    f2.createDimension('TSTEP', dates.shape[0])
    f2.createDimension('DATE-TIME', 2)
    f2.createDimension('LAY', 1)
    f2.createDimension('VAR', data.shape[1])
    f2.createDimension('ROW', np.size(yY,0))
    f2.createDimension('COL', np.size(xX,1))
    
    # Building variables
    TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    #TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP'))
    TFLAG[:,:,:] = tflag
    TFLAG.units = '<YYYYDDD,HHMMSS>'
    
    LON = f2.createVariable('Longitute', 'f4', ( 'ROW','COL'))
    LAT = f2.createVariable('Latitude', 'f4', ( 'ROW','COL'))  
    AREA = f2.createVariable('AREA', 'f4', ( 'ROW','COL'))
    LAT[:,:,] =  yY
    LON[:,:,] = xX
    AREA[:,:] = area
    LAT.units = 'degrees '
    LON.units = 'degrees '
    AREA.units = 'km2 '
    
    if specie=='ACET':
        ACET = f2.createVariable('ACET', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        ACET[:,:,:,:] =  data[:,0,:,:]
        ACET.units = 'moles/s '
    if specie=='ACROLEIN':
        ACROLEIN = f2.createVariable('ACROLEIN', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        ACROLEIN[:,:,:,:] =  data[:,0,:,:]
        ACROLEIN.units = 'moles/s '
    if specie=='ALD2':
        ALD2 = f2.createVariable('ALD2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        ALD2[:,:,:,:] =  data[:,0,:,:]
        ALD2.units = 'moles/s '
    if specie=='ALD2_PRIMARY':
        ALD2_PRIMARY = f2.createVariable('ALD2_PRIMARY', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        ALD2_PRIMARY[:,:,:,:] =  data[:,0,:,:]
        ALD2_PRIMARY.units = 'moles/s '
    if specie=='ALDX':
        ALDX = f2.createVariable('ALDX', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        ALDX[:,:,:,:] =  data[:,0,:,:]
        ALDX.units = 'moles/s '
    if specie=='BENZ':
        BENZ = f2.createVariable('BENZ', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        BENZ[:,:,:,:] =  data[:,0,:,:]
        BENZ.units = 'moles/s '        
    if specie=='BUTADIENE13':
        BUTADIENE13 = f2.createVariable('BUTADIENE13', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        BUTADIENE13[:,:,:,:] =  data[:,0,:,:]
        BUTADIENE13.units = 'moles/s '
    if specie=='CH4':
        CH4 = f2.createVariable('CH4', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        CH4[:,:,:,:] =  data[:,0,:,:]
        CH4.units = 'moles/s '
    if specie=='CH4_INV':
        CH4_INV = f2.createVariable('CH4_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        CH4_INV[:,:,:,:] =  data[:,0,:,:]
        CH4_INV.units = 'g/s '       
    if specie=='CL2':
        CL2 = f2.createVariable('CL2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        CL2[:,:,:,:] =  data[:,0,:,:]
        CL2.units = 'moles/s '   
    if specie=='CO':
        CO = f2.createVariable('CO', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        CO[:,:,:,:] =  data[:,0,:,:]
        CO.units = 'moles/s '   
    if specie=='CO2_INV':
        CO2_INV = f2.createVariable('CO2_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        CO2_INV[:,:,:,:] =  data[:,0,:,:]
        CO2_INV.units = 'g/s '           
    if specie=='ETH':
        ETH = f2.createVariable('ETH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        ETH[:,:,:,:] =  data[:,0,:,:]
        ETH.units = 'moles/s ' 
    if specie=='ETHA':
        ETHA = f2.createVariable('ETHA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        ETHA[:,:,:,:] =  data[:,0,:,:]
        ETHA.units = 'moles/s ' 
    if specie=='ETHY':
        ETHY = f2.createVariable('ETHY', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        ETHY[:,:,:,:] =  data[:,0,:,:]
        ETHY.units = 'moles/s '  
    if specie=='ETOH':
        ETOH = f2.createVariable('ETOH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        ETOH[:,:,:,:] =  data[:,0,:,:]
        ETOH.units = 'moles/s '  
    if specie=='FORM':
        FORM = f2.createVariable('FORM', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        FORM[:,:,:,:] =  data[:,0,:,:]
        FORM.units = 'moles/s '  
    if specie=='FORM_PRIMARY':
        FORM_PRIMARY = f2.createVariable('FORM_PRIMARY', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        FORM_PRIMARY[:,:,:,:] =  data[:,0,:,:]
        FORM_PRIMARY.units = 'moles/s '  
    if specie=='HCL':
        HCL = f2.createVariable('HCL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        HCL[:,:,:,:] =  data[:,0,:,:]
        HCL.units = 'moles/s '  
    if specie=='HONO':
        HONO = f2.createVariable('HONO', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        HONO[:,:,:,:] =  data[:,0,:,:]
        HONO.units = 'moles/s '  
    if specie=='IOLE':
        IOLE = f2.createVariable('IOLE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        IOLE[:,:,:,:] =  data[:,0,:,:]
        IOLE.units = 'moles/s '  
    if specie=='ISOP':
        ISOP = f2.createVariable('ISOP', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        ISOP[:,:,:,:] =  data[:,0,:,:]
        ISOP.units = 'moles/s '  
    if specie=='KET':
        KET = f2.createVariable('KET', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        KET[:,:,:,:] =  data[:,0,:,:]
        KET.units = 'moles/s '  
    if specie=='MEOH':
        MEOH = f2.createVariable('MEOH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        MEOH[:,:,:,:] =  data[:,0,:,:]
        MEOH.units = 'moles/s '  
    if specie=='N2O_INV':
        N2O_INV = f2.createVariable('N2O_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        N2O_INV[:,:,:,:] =  data[:,0,:,:]
        N2O_INV.units = 'g/s ' 
    if specie=='NAPH':
        NAPH = f2.createVariable('NAPH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        NAPH[:,:,:,:] =  data[:,0,:,:]
        NAPH.units = 'moles/s '  
    if specie=='NH3':
        NH3 = f2.createVariable('NH3', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        NH3[:,:,:,:] =  data[:,0,:,:]
        NH3.units = 'moles/s '
    if specie=='NH3_FERT':
        NH3_FERT = f2.createVariable('NH3_FERT', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        NH3_FERT[:,:,:,:] =  data[:,0,:,:]
        NH3_FERT.units = 'moles/s '  
    if specie=='NO':
        NO = f2.createVariable('NO', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        NO[:,:,:,:] =  data[:,0,:,:]
        NO.units = 'moles/s '  
    if specie=='NO2':
        NO2 = f2.createVariable('NO2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        NO2[:,:,:,:] =  data[:,0,:,:]
        NO2.units = 'moles/s '       
    if specie=='NVOL':
        NVOL = f2.createVariable('NVOL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        NVOL[:,:,:,:] =  data[:,0,:,:]
        NVOL.units = 'moles/s '  
    if specie=='OLE':
        OLE = f2.createVariable('OLE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        OLE[:,:,:,:] =  data[:,0,:,:]
        OLE.units = 'moles/s '  
    if specie=='PAL':
        PAL = f2.createVariable('PAL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PAL[:,:,:,:] =  data[:,0,:,:]
        PAL.units = 'moles/s '
    if specie=='PAR':
        PAR = f2.createVariable('PAR', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PAR[:,:,:,:] =  data[:,0,:,:]
        PAR.units = 'moles/s '  
    if specie=='PCA':
        PCA = f2.createVariable('PCA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PCA[:,:,:,:] =  data[:,0,:,:]
        PCA.units = 'g/s ' 
    if specie=='PCL':
        PCL = f2.createVariable('PCL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PCL[:,:,:,:] =  data[:,0,:,:]
        PCL.units = 'g/s '
    if specie=='PEC':
        PEC = f2.createVariable('PEC', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PEC[:,:,:,:] =  data[:,0,:,:]
        PEC.units = 'g/s '  
    if specie=='PFE':
        PFE = f2.createVariable('PFE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PFE[:,:,:,:] =  data[:,0,:,:]
        PFE.units = 'g/s '  
    if specie=='PH2O':
        PH2O = f2.createVariable('PH2O', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PH2O[:,:,:,:] =  data[:,0,:,:]
        PH2O.units = 'g/s '
    if specie=='PK':
        PK = f2.createVariable('PK', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PK[:,:,:,:] =  data[:,0,:,:]
        PK.units = 'g/s '  
    if specie=='PMC':
        PMC = f2.createVariable('PMC', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PMC[:,:,:,:] =  data[:,0,:,:]
        PMC.units = 'g/s '  
    if specie=='PMG':
        PMG = f2.createVariable('PMG', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PMG[:,:,:,:] =  data[:,0,:,:]
        PMG.units = 'g/s '
    if specie=='PMN':
        PMN = f2.createVariable('PMN', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PMN[:,:,:,:] =  data[:,0,:,:]
        PMN.units = 'g/s ' 
    if specie=='PMOTHR':
        PMOTHR = f2.createVariable('PMOTHR', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PMOTHR[:,:,:,:] =  data[:,0,:,:]
        PMOTHR.units = 'g/s ' 
    if specie=='PNA':
        PNA = f2.createVariable('PNA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PNA[:,:,:,:] =  data[:,0,:,:]
        PNA.units = 'g/s '     
    if specie=='PNCOM':
        PNCOM = f2.createVariable('PNCOM', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PNCOM[:,:,:,:] =  data[:,0,:,:]
        PNCOM.units = 'g/s ' 
    if specie=='PNH4':
        PNH4 = f2.createVariable('PNH4', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PNH4[:,:,:,:] =  data[:,0,:,:]
        PNH4.units = 'g/s ' 
    if specie=='PNO3':
        PNO3 = f2.createVariable('PNO3', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PNO3[:,:,:,:] =  data[:,0,:,:]
        PNO3.units = 'g/s '        
    if specie=='POC':
        POC = f2.createVariable('POC', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        POC[:,:,:,:] =  data[:,0,:,:]
        POC.units = 'g/s ' 
    if specie=='PRPA':
        PRPA = f2.createVariable('PRPA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PRPA[:,:,:,:] =  data[:,0,:,:]
        PRPA.units = 'moles/s '
    if specie=='PSI':
        PSI = f2.createVariable('PSI', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PSI[:,:,:,:] =  data[:,0,:,:]
        PSI.units = 'g/s '         
    if specie=='PSO4':
        PSO4 = f2.createVariable('PSO4', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PSO4[:,:,:,:] =  data[:,0,:,:]
        PSO4.units = 'g/s ' 
    if specie=='PTI':
        PTI = f2.createVariable('PTI', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PTI[:,:,:,:] =  data[:,0,:,:]
        PTI.units = 'g/s ' 
    if specie=='SO2':
        SO2 = f2.createVariable('SO2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        print(str(SO2.shape) + str(data.shape))
        SO2[:,:,:,:] =  data[:,0,:,:]
        SO2.units = 'moles/s '
    if specie=='SOAALK':
        SOAALK = f2.createVariable('SOAALK', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        SOAALK[:,:,:,:] =  data[:,0,:,:]
        SOAALK.units = 'moles/s '
    if specie=='SULF':
        SULF = f2.createVariable('SULF', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        SULF[:,:,:,:] =  data[:,0,:,:]
        SULF.units = 'moles/s '
    if specie=='TERP':
        TERP = f2.createVariable('TERP', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        TERP[:,:,:,:] =  data[:,0,:,:]
        TERP.units = 'moles/s '       
    if specie=='TOL':
        TOL = f2.createVariable('TOL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        TOL[:,:,:,:] =  data[:,0,:,:]
        TOL.units = 'moles/s '  
    if specie=='UNK':
        UNK = f2.createVariable('UNK', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        UNK[:,:,:,:] =  data[:,0,:,:]
        UNK.units = 'moles/s '  
    if specie=='UNR':
        UNR = f2.createVariable('UNR', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        UNR[:,:,:,:] =  data[:,0,:,:]
        UNR.units = 'moles/s '       
    if specie=='VOC_INV':
        VOC_INV = f2.createVariable('VOC_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        VOC_INV[:,:,:,:] =  data[:,0,:,:]
        VOC_INV.units = 'g/s '  
    if specie=='XYLMN':
        XYLMN = f2.createVariable('XYLMN', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        XYLMN[:,:,:,:] =  data[:,0,:,:]
        XYLMN.units = 'moles/s '      
    if specie=='PMFINE':
        PMFINE = f2.createVariable('PMFINE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
        PMFINE[:,:,:,:] =  data[:,0,:,:]
        PMFINE.units = 'g/s '
        
    
    f2.close()
    print('Your BRAVESdatabase netCDF file is ready!')


    #%% Creating a dataset
def createNETCDFtemporalfromNCforWRFCHEM(rootPath,folder,name,data,xX,yY,dates,area,wrfPath):
    #https://github.com/wrf-model/WRF/blob/master/Registry/registry.chem
    # import netCDF4 as nc
    # path='/media/leohoinaski/HDD/wrfchemi_00_d02/wrfchemi_00z_d02'
    # ds = nc.Dataset(path)
    # dataVar = list(ds.variables.keys())
    # dataNC = ds[dataVar[5]][:]
    
    print('===================STARTING netCDFcreator_v1.py=======================')
    cdate = datetime.datetime.now()
    cdateStr = int(str(cdate.year)+str(cdate.timetuple().tm_yday))
    ctime = int(str(cdate.hour)+str(cdate.minute)+str(cdate.second))
    
    sdate =  dates['year'][0]*1000 + dates.index[0].timetuple().tm_yday            
    
    
    # Converting emissions in g(mol)/s to g(mol)/h.km2 - WRF-CHEM
    print('Converting to emission density in ug/h.km2')
    convPath = rootPath + '/ChemicalSpec/WRFCHEM_speciesConv.csv'
    conv = pd.read_csv(convPath)
    # data = np.zeros([dataTempo.shape[0],dataTempo.shape[1],
    #                       dataTempo.shape[2],dataTempo.shape[3]])

    for ii in range(0,data.shape[1]):
        print(str(ii))
        data[:,ii,:,:] = (data[:,ii,:,:]/area)*conv['CONV'][ii]
        
    
    f2 = nc4.Dataset(folder+'/'+name,'w', format='NETCDF4_CLASSIC') #'w' stands for write    
    # ff = '/media/leohoinaski/HDD/wrfchemi_00_d02/wrfchemi_00z_d02'
    # ds = nc.Dataset(ff) 
    ds3 = nc4.Dataset(wrfPath)
    for gatr in ds3.ncattrs() :
        setattr(f2, gatr, ds3.getncattr(gatr))

    
    f2.TITLE= ' BRAVES emissions for WRFCHEM'
    f2.CDATE= cdateStr
    f2.CTIME= ctime
    f2.WDATE= cdateStr
    f2.WTIME= ctime
    f2.SDATE= sdate
    f2.STIME= 0
    f2.TSTEP= 10000
    f2.NTHIK= 1
    f2.NCOLS= data.shape[3]
    print('NCOLS=' +str(data.shape[3]))
    f2.NROWS= data.shape[2]
    print('NROWS=' +str(data.shape[2]))
    f2.NLAYS= 1
    f2.NVARS= data.shape[1]
    f2.GDTYP= 1
    f2.P_ALP= -10
    f2.P_BET= 0
    f2.P_GAM= round(xX.mean(),6)
    f2.XCENT= round(xX.mean(),6)
    print('XCENT='+ str(round(xX.mean(),6)))
    f2.YCENT= round(yY.mean(),6)
    print('YCENT='+ str(round(yY.mean(),6)))
    f2.XORIG= np.round((xX.min() - (xX[0,1] - xX[0,0])/2),6)
    print('XORIG = ' + str(np.round((xX.min() - (xX[0,1] - xX[0,0])/2),6))) 
    f2.YORIG= np.round((yY.min() - (yY[1,0] - yY[0,0])/2),6)
    print('YORIG = ' + str(np.round((yY.min() - (yY[1,0] - yY[0,0])/2),6))) 
    f2.XCELL= round((xX[0,1] - xX[0,0]), 6)
    print('XCELL = '+str(round(xX[0,1] - xX[0,0],6)))
    f2.YCELL= round((yY[1,0] - yY[0,0]), 6)
    print('YCELL = '+str(round(yY[1,0] - yY[0,0],6)))
    f2.VGTYP= -1
    f2.VGTOP= 0.0
    f2.VGLVLS= [0,0]


    #f2.FILEDESC= 'BRAVES database vehicular emissions ready for CMAQ'
    f2.HISTORY ='' 
       
    # # Specifying dimensions
    #tempgrp = f.createGroup('vehicularEmissions_data')
    f2.createDimension('DateStrLen', len(ds3['Times'][0,:]))
    f2.createDimension('Time', data.shape[0])
    f2.createDimension('emissions_zdim_stag', 1)
    #f2.createDimension('VAR', data.shape[1])
    f2.createDimension('south_north', np.size(yY,0))
    f2.createDimension('west_east', np.size(xX,1))
    
    # Building variables
    #https://ruc.noaa.gov/wrf/wrf-chem/Emission_guide.pdf
    #TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    # DateStrLen = f2.createVariable('DateStrLen', 'i4', ('DateStrLen'))  
    # Time = f2.createVariable('Time', 'i4', ('Time'))
    # west_east = f2.createVariable('west_east', 'i4', ( 'west_east')) 
    # south_north = f2.createVariable('south_north', 'i4', ( 'south_north'))
    # emissions_zdim_stag = f2.createVariable('emissions_zdim_stag', 'i4', ( 'emissions_zdim_stag'))
    

     
    dtes = np.array(dates.index.strftime("%Y-%m-%d_%H:%M:%S").str.split('',expand=True))
    dts=np.empty((dtes.shape[0],19), dtype=object)
    for ii in range(0,dtes.shape[0]):
        for jj in range(0,19):
            dts[ii,jj] = dtes[ii][jj+1]
    
    
    Times = f2.createVariable('Times', 'S1', ('Time','DateStrLen'))
    Times[:,:] = dts
    
    XLAT = f2.createVariable('XLAT', 'f4', ( 'south_north','west_east'))
    print(ds3['XLAT'][0,:,:].shape)
    XLAT[:,:] = ds3['XLAT'][0,:,:]
    XLAT.units = 'degree north '
    
    XLONG = f2.createVariable('XLONG', 'f4', ( 'south_north','west_east'))
    XLONG[:,:] = ds3['XLONG'][0,:,:]
    XLONG.units = 'degree_east '
    
    e_acet = f2.createVariable('e_acet', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_acet[:,:,:,:] =  data[:,0,:,:]
    e_acet.units = 'mol km^-2 hr^-1'
          
    # Sum of acrolein and Butadien13
    e_macr = f2.createVariable('e_macr', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_macr[:,:,:,:] =  data[:,1,:,:] + data[:,6,:,:]
    e_macr.units = 'mol km^-2 hr^-1'
    
    e_ald2 = f2.createVariable('e_ald2', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_ald2[:,:,:,:] =  data[:,2,:,:]
    e_ald2.units = 'mol km^-2 hr^-1'
    
    e_ald = f2.createVariable('e_ald', 'f4',('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_ald[:,:,:,:] =  data[:,3,:,:]
    e_ald.units = 'mol km^-2 hr^-1'
     
    # Higher aldehydes
    e_aldx = f2.createVariable('e_aldx', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_aldx[:,:,:,:] =  data[:,4,:,:]
    e_aldx.units = 'mol km^-2 hr^-1'
    
    e_benzene = f2.createVariable('e_benzene', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_benzene[:,:,:,:] =  data[:,5,:,:]
    e_benzene.units = 'mol km^-2 hr^-1'
    
    
    e_ch4 = f2.createVariable('e_ch4', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_ch4[:,:,:,:] =  data[:,7,:,:]
    e_ch4.units = 'mol km^-2 hr^-1'
    
    #CH4_INV = f2.createVariable('CH4_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    
    e_clc = f2.createVariable('e_clc', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_clc[:,:,:,:] =  data[:,9,:,:]
    e_clc.units = 'ug/m3 m/s'
    
    e_co = f2.createVariable('e_co', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_co[:,:,:,:] =  data[:,10,:,:]
    e_co.units = 'mol km^-2 hr^-1'
    
    e_co2 = f2.createVariable('e_co2', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_co2[:,:,:,:] =  data[:,11,:,:]
    e_co2.units = 'mol km^-2 hr^-1'
    
    e_eth = f2.createVariable('e_eth', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_eth[:,:,:,:] =  data[:,12,:,:]
    e_eth.units = 'mol km^-2 hr^-1'
    
    e_etha = f2.createVariable('e_etha', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_etha[:,:,:,:] =  data[:,13,:,:]
    e_etha.units = 'mol km^-2 hr^-1'
    
    # ETHY
    e_c2h2 = f2.createVariable('e_c2h2', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_c2h2[:,:,:,:] =  data[:,14,:,:]
    e_c2h2.units = 'mol km^-2 hr^-1'
    
    e_etoh = f2.createVariable('e_etoh', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_etoh[:,:,:,:] =  data[:,15,:,:]
    e_etoh.units = 'mol km^-2 hr^-1'
    
    e_c2h5oh = f2.createVariable('e_c2h5oh', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_c2h5oh[:,:,:,:] =  data[:,15,:,:]
    e_c2h5oh.units = 'mol km^-2 hr^-1'
    
    e_form = f2.createVariable('e_form', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_form[:,:,:,:] =  data[:,16,:,:]
    e_form.units = 'mol km^-2 hr^-1'
    
   #FORM_PRIMARY = f2.createVariable('FORM_PRIMARY', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
   
    e_hcl = f2.createVariable('e_hcl', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_hcl[:,:,:,:] =  data[:,18,:,:]
    e_hcl.units = 'mol km^-2 hr^-1'
    
    e_hono = f2.createVariable('e_hono', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_hono[:,:,:,:] =  data[:,19,:,:]
    e_hono.units = 'mol km^-2 hr^-1'   
    
    e_iole = f2.createVariable('e_iole', 'f4',  ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_iole[:,:,:,:] =  data[:,20,:,:]
    e_iole.units = 'mol km^-2 hr^-1'       
    
    e_isop = f2.createVariable('e_isop', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_isop[:,:,:,:] =  data[:,21,:,:]
    e_isop.units = 'mol km^-2 hr^-1'     
    
    e_ket = f2.createVariable('e_ket', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_ket[:,:,:,:] =  data[:,22,:,:]
    e_ket.units = 'mol km^-2 hr^-1'
    
    e_meoh = f2.createVariable('e_meoh', 'f4',('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_meoh[:,:,:,:] =  data[:,23,:,:]
    e_meoh.units = 'mol km^-2 hr^-1'
    
    e_n2o = f2.createVariable('e_n2o', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_n2o[:,:,:,:] =  data[:,24,:,:]
    e_n2o.units = 'mol km^-2 hr^-1'    
    
    e_nh3 = f2.createVariable('e_nh3', 'f4',('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_nh3[:,:,:,:] =  data[:,26,:,:]
    e_nh3.units = 'mol km^-2 hr^-1'      
        
    e_no = f2.createVariable('e_no', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_no[:,:,:,:] =  data[:,28,:,:]
    e_no.units = 'mol km^-2 hr^-1'    
    
    e_no2 = f2.createVariable('e_no2', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_no2[:,:,:,:] =  data[:,29,:,:]
    e_no2.units = 'mol km^-2 hr^-1' 
        
    e_ole = f2.createVariable('e_ole', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_ole[:,:,:,:] =  data[:,31,:,:]
    e_ole.units = 'mol km^-2 hr^-1' 
    
    e_par = f2.createVariable('e_par', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_par[:,:,:,:] =  data[:,33,:,:]
    e_par.units = 'mol km^-2 hr^-1' 
           
    e_cli = f2.createVariable('e_cli', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_cli[:,:,:,:] =  data[:,35,:,:]
    e_cli.units = 'ug/m3 m/s'
    
    e_clj = f2.createVariable('e_clj', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_clj[:,:,:,:] =  data[:,35,:,:]
    e_clj.units = 'ug/m3 m/s'
    
    # e_clc = f2.createVariable('e_clc', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    # e_clc[:,:,:,:] =  data[:,36,:,:]
    # e_clc.units = 'ug/m3 m/s' 
    
    e_eci = f2.createVariable('e_eci', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_eci[:,:,:,:] =  data[:,36,:,:]
    e_eci.units = 'ug/m3 m/s'   
    
    e_ecj = f2.createVariable('e_ecj', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_ecj[:,:,:,:] =  data[:,36,:,:]
    e_ecj.units = 'ug/m3 m/s'   
    
    e_ecc = f2.createVariable('e_ecc', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_ecc[:,:,:,:] =  data[:,36,:,:]
    e_ecc.units = 'ug/m3 m/s'   
    
    e_pm_10 = f2.createVariable('e_pm_10', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_pm_10[:,:,:,:] =  data[:,40,:,:]
    e_pm_10.units = 'ug/m3 m/s' 
    
    e_nac = f2.createVariable('e_nac', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_nac[:,:,:,:] =  data[:,44,:,:]
    e_nac.units = 'ug/m3 m/s' 
    
    e_nh4c = f2.createVariable('e_nh4c', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_nh4c[:,:,:,:] =  data[:,46,:,:]
    e_nh4c.units = 'ug/m3 m/s' 
   
    e_no3c = f2.createVariable('e_no3c', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_no3c[:,:,:,:] =  data[:,47,:,:]
    e_no3c.units = 'ug/m3 m/s' 

    e_oc = f2.createVariable('e_oc', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_oc[:,:,:,:] =  data[:,48,:,:]
    e_oc.units = 'ug/m3 m/s' 

    e_c3h8 = f2.createVariable('e_c3h8', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_c3h8[:,:,:,:] =  data[:,49,:,:]
    e_c3h8.units = 'mol km^-2 hr^-1' 

    e_so4c = f2.createVariable('e_so4c', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_so4c[:,:,:,:] =  data[:,51,:,:]
    e_so4c.units = 'ug/m3 m/s'
    
    e_so4i = f2.createVariable('e_so4i', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_so4i[:,:,:,:] =  data[:,51,:,:]
    e_so4i.units = 'ug/m3 m/s'
    
    e_so4j = f2.createVariable('e_so4j', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_so4j[:,:,:,:] =  data[:,51,:,:]
    e_so4i.units = 'ug/m3 m/s'

    e_so2 = f2.createVariable('e_so2', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_so2[:,:,:,:] =  data[:,53,:,:]
    e_so2.units = 'mol km^-2 hr^-1' 

    e_terp = f2.createVariable('e_terp', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_terp[:,:,:,:] =  data[:,56,:,:]
    e_terp.units = 'mol km^-2 hr^-1' 
    
    e_tol = f2.createVariable('e_tol', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_tol[:,:,:,:] =  data[:,57,:,:]
    e_tol.units = 'mol km^-2 hr^-1' 

    e_xyl = f2.createVariable('e_xyl', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_xyl[:,:,:,:] =  data[:,61,:,:]
    e_xyl.units = 'mol km^-2 hr^-1'  
    
    e_pm25 = f2.createVariable('e_pm25', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_pm25[:,:,:,:] =  data[:,62,:,:]
    e_pm25.units = 'ug/m3 m/s' 
    
    e_pm25j = f2.createVariable('e_pm25j', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_pm25j[:,:,:,:] =  data[:,62,:,:]
    e_pm25j.units = 'ug/m3 m/s' 
    
    e_pm25i = f2.createVariable('e_pm25i', 'f4', ('Time', 'emissions_zdim_stag','south_north', 'west_east'))
    e_pm25i[:,:,:,:] =  data[:,62,:,:]
    e_pm25i.units = 'ug/m3 m/s' 
       

   
    f2.close()
    print('Your BRAVESdatabase netCDF file is ready for WRFCHEM!')
