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
    f2.P_GAM= round(xX.mean(),5)
    f2.XCENT= round(xX.mean(),5)
    print('XCENT='+ str(round(xX.mean(),5)))
    f2.YCENT= round(yY.mean(),5)
    print('YCENT='+ str(round(yY.mean(),5)))
    f2.XORIG= np.round((xX.min() - (xX[0,1] - xX[0,0])/2),4)
    print('XORIG = ' + str(np.round((xX.min() - (xX[0,1] - xX[0,0])/2),4))) 
    f2.YORIG= np.round((yY.min() - (yY[1,0] - yY[0,0])/2),4)
    print('YORIG = ' + str(np.round((yY.min() - (yY[1,0] - yY[0,0])/2),4))) 
    f2.XCELL= round((xX[0,1] - xX[0,0]), 4)
    print('XCELL = '+str(round(xX[0,1] - xX[0,0],4)))
    f2.YCELL= round((yY[1,0] - yY[0,0]), 4)
    print('YCELL = '+str(round(yY[1,0] - yY[0,0],4)))
    f2.VGTYP= -1
    f2.VGTOP= 0.0
    f2.VGLVLS= [0,0]
    
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
    LON = f2.createVariable('LON', 'f4', ( 'ROW','COL'))
    LAT = f2.createVariable('LAT', 'f4', ( 'ROW','COL'))
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
    print('Your BRAVESdatabase netCDF file is ready!')
    
    #%% Creating a dataset
def createNETCDFtemporalfromNC(folder,name,data,xX,yY,dates,area):
    print('===================STARTING netCDFcreator_v1.py=======================')
    cdate = datetime.datetime.now()
    cdateStr = int(str(cdate.year)+str(cdate.timetuple().tm_yday))
    ctime = int(str(cdate.hour)+str(cdate.minute)+str(cdate.second))

    tflag = np.empty([dates.shape[0],data.shape[1],2],dtype='i4')
    for ii in range(0,dates.shape[0]):
        tflag[ii,:,0]=int(dates['year'][0]*1000 + dates.index[ii].timetuple().tm_yday)
        tflag[ii,:,1]=int(str(dates['hour'][ii])+'0000')
    
    sdate =  dates['year'][0]*1000 + dates.index[0].timetuple().tm_yday            
    
    # Reshape area array
    
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
    f2.P_GAM= round(xX.mean(),5)
    f2.XCENT= round(xX.mean(),5)
    print('XCENT='+ str(round(xX.mean(),5)))
    f2.YCENT= round(yY.mean(),5)
    print('YCENT='+ str(round(yY.mean(),5)))
    f2.XORIG= np.round((xX.min() - (xX[0,1] - xX[0,0])/2),4)
    print('XORIG = ' + str(np.round((xX.min() - (xX[0,1] - xX[0,0])/2),4))) 
    f2.YORIG= np.round((yY.min() - (yY[1,0] - yY[0,0])/2),4)
    print('YORIG = ' + str(np.round((yY.min() - (yY[1,0] - yY[0,0])/2),4))) 
    f2.XCELL= round((xX[0,1] - xX[0,0]), 4)
    print('XCELL = '+str(round(xX[0,1] - xX[0,0],4)))
    f2.YCELL= round((yY[1,0] - yY[0,0]), 4)
    print('YCELL = '+str(round(yY[1,0] - yY[0,0],4)))
    f2.VGTYP= -1
    f2.VGTOP= 0.0
    f2.VGLVLS= [0,0]
    
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
    f2.FILEDESC= 'BRAVES database vehicular emissions ready for CMAQ'
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
    LON = f2.createVariable('LON', 'f4', ( 'ROW','COL'))
    LAT = f2.createVariable('LAT', 'f4', ( 'ROW','COL'))
    AREA = f2.createVariable('AREA', 'f4', ( 'ROW','COL'))
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
    print(area.shape)
    AREA[:,:] = area
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
    ACET.units = 'moles/s '
    ACROLEIN.units = 'moles/s '
    ALD2.units = 'moles/s '
    ALD2_PRIMARY.units = 'moles/s '
    ALDX.units = 'moles/s '
    BENZ.units = 'moles/s '
    BUTADIENE13.units = 'moles/s '
    CH4.units = 'moles/s '
    CH4_INV.units = 'g/s '
    CL2.units = 'moles/s ' 
    CO.units = 'moles/s '
    CO2_INV.units = 'g/s '
    ETH.units = 'moles/s '
    ETHA.units = 'moles/s '
    ETHY.units = 'moles/s '
    ETOH.units = 'moles/s '
    FORM.units = 'moles/s '
    FORM_PRIMARY.units = 'moles/s '
    HCL.units = 'moles/s '
    HONO.units = 'moles/s '
    IOLE.units = 'moles/s '
    ISOP.units = 'moles/s '
    KET.units = 'moles/s '
    MEOH.units = 'moles/s '
    N2O_INV.units = 'g/s '
    NAPH.units = 'moles/s '
    NH3.units = 'moles/s '
    NH3_FERT.units = 'moles/s '
    NO.units = 'moles/s '
    NO2.units = 'moles/s '
    NVOL.units = 'moles/s '
    OLE.units = 'moles/s '
    PAL.units = 'moles/s '
    PAR.units = 'moles/s '
    PCA.units = 'g/s '
    PCL.units = 'g/s '
    PEC.units = 'g/s '
    PFE.units = 'g/s '
    PH2O.units = 'g/s '
    PK.units = 'g/s '
    PMC.units = 'g/s ' 
    PMG.units = 'g/s '
    PMN.units = 'g/s ' 
    PMOTHR.units = 'g/s ' 
    PNA.units = 'g/s ' 
    PNCOM.units = 'g/s ' 
    PNH4.units = 'g/s ' 
    PNO3.units = 'g/s ' 
    POC.units = 'g/s ' 
    PRPA.units = 'moles/s '
    PSI.units = 'g/s ' 
    PSO4.units = 'g/s ' 
    PTI.units = 'g/s ' 
    SO2.units = 'moles/s ' 
    SOAALK.units = 'moles/s ' 
    SULF.units = 'moles/s ' 
    TERP.units = 'moles/s ' 
    TOL.units = 'moles/s ' 
    UNK.units = 'moles/s ' 
    UNR.units = 'moles/s ' 
    VOC_INV.units = 'g/s ' 
    XYLMN.units = 'moles/s '
    PMFINE.units = 'g/s '
    LON.units = 'degrees '
    LAT.units = 'degrees '
    AREA.units = 'km2 '
   
    f2.close()
    print('Your BRAVESdatabase netCDF file is ready!')

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
    f2.P_GAM= round(xX.mean(),5)
    f2.XCENT= round(xX.mean(),5)
    print('XCENT='+ str(round(xX.mean(),5)))
    f2.YCENT= round(yY.mean(),5)
    print('YCENT='+ str(round(yY.mean(),5)))
    f2.XORIG= np.round((xX.min() - (xX[0,1] - xX[0,0])/2),4)
    print('XORIG = ' + str(np.round((xX.min() - (xX[0,1] - xX[0,0])/2),4))) 
    f2.YORIG= np.round((yY.min() - (yY[1,0] - yY[0,0])/2),4)
    print('YORIG = ' + str(np.round((yY.min() - (yY[1,0] - yY[0,0])/2),4))) 
    f2.XCELL= round((xX[0,1] - xX[0,0]), 4)
    print('XCELL = '+str(round(xX[0,1] - xX[0,0],4)))
    f2.YCELL= round((yY[1,0] - yY[0,0]), 4)
    print('YCELL = '+str(round(yY[1,0] - yY[0,0],4)))
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
    
    LON = f2.createVariable('LON', 'f4', ( 'ROW','COL'))
    LAT = f2.createVariable('LAT', 'f4', ( 'ROW','COL'))  
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


    #%% Creating a dataset
def createNETCDFtemporalfromNCforWRFCHEM(folder,name,data,xX,yY,dates,area):
    print('===================STARTING netCDFcreator_v1.py=======================')
    cdate = datetime.datetime.now()
    cdateStr = int(str(cdate.year)+str(cdate.timetuple().tm_yday))
    ctime = int(str(cdate.hour)+str(cdate.minute)+str(cdate.second))

    tflag = np.empty([dates.shape[0],data.shape[1],2],dtype='i4')
    for ii in range(0,dates.shape[0]):
        tflag[ii,:,0]=int(dates['year'][0]*1000 + dates.index[ii].timetuple().tm_yday)
        tflag[ii,:,1]=int(str(dates['hour'][ii])+'0000')
    
    sdate =  dates['year'][0]*1000 + dates.index[0].timetuple().tm_yday            
    
    
    # Converting emissions in g(mol)/s to g(mol)/h.km2 - WRF-CHEM
    data = (data/area)*3600
    
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
    f2.P_GAM= round(xX.mean(),5)
    f2.XCENT= round(xX.mean(),5)
    print('XCENT='+ str(round(xX.mean(),5)))
    f2.YCENT= round(yY.mean(),5)
    print('YCENT='+ str(round(yY.mean(),5)))
    f2.XORIG= np.round((xX.min() - (xX[0,1] - xX[0,0])/2),4)
    print('XORIG = ' + str(np.round((xX.min() - (xX[0,1] - xX[0,0])/2),4))) 
    f2.YORIG= np.round((yY.min() - (yY[1,0] - yY[0,0])/2),4)
    print('YORIG = ' + str(np.round((yY.min() - (yY[1,0] - yY[0,0])/2),4))) 
    f2.XCELL= round((xX[0,1] - xX[0,0]), 4)
    print('XCELL = '+str(round(xX[0,1] - xX[0,0],4)))
    f2.YCELL= round((yY[1,0] - yY[0,0]), 4)
    print('YCELL = '+str(round(yY[1,0] - yY[0,0],4)))
    f2.VGTYP= -1
    f2.VGTOP= 0.0
    f2.VGLVLS= [0,0]
    
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
    f2.FILEDESC= 'BRAVES database vehicular emissions ready for CMAQ'
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
    LON = f2.createVariable('LON', 'f4', ( 'ROW','COL'))
    LAT = f2.createVariable('LAT', 'f4', ( 'ROW','COL'))
    AREA = f2.createVariable('AREA', 'f4', ( 'ROW','COL'))
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
    print(area.shape)
    AREA[:,:] = area
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
    ACET.units = 'moles/(km2.h) '
    ACROLEIN.units = 'moles/(km2.h) '
    ALD2.units = 'moles/(km2.h) '
    ALD2_PRIMARY.units = 'moles/(km2.h) '
    ALDX.units = 'moles/(km2.h) '
    BENZ.units = 'moles/(km2.h) '
    BUTADIENE13.units = 'moles/(km2.h) '
    CH4.units = 'moles/(km2.h) '
    CH4_INV.units = 'moles/(km2.h) '
    CL2.units = 'moles/(km2.h) '
    CO.units = 'moles/(km2.h) '
    CO2_INV.units = 'g/(km2.h) '
    ETH.units = 'moles/(km2.h) '
    ETHA.units = 'moles/(km2.h) '
    ETHY.units = 'moles/(km2.h) '
    ETOH.units = 'moles/(km2.h) '
    FORM.units = 'moles/(km2.h) '
    FORM_PRIMARY.units = 'moles/(km2.h) '
    HCL.units = 'moles/(km2.h) '
    HONO.units = 'moles/(km2.h) '
    IOLE.units = 'moles/(km2.h) '
    ISOP.units = 'moles/(km2.h) '
    KET.units = 'moles/(km2.h) '
    MEOH.units = 'moles/(km2.h) '
    N2O_INV.units = 'moles/(km2.h) '
    NAPH.units = 'moles/(km2.h) '
    NH3.units = 'moles/(km2.h) '
    NH3_FERT.units = 'moles/(km2.h) '
    NO.units = 'moles/(km2.h) '
    NO2.units = 'moles/(km2.h) '
    NVOL.units = 'moles/(km2.h) '
    OLE.units = 'moles/(km2.h) '
    PAL.units = 'moles/(km2.h) '
    PAR.units = 'moles/(km2.h) '
    PCA.units = 'g/(km2.h) '
    PCL.units = 'g/(km2.h) '
    PEC.units = 'g/(km2.h) '
    PFE.units = 'g/(km2.h) '
    PH2O.units = 'g/(km2.h) '
    PK.units = 'g/(km2.h) '
    PMC.units = 'g/(km2.h) '
    PMG.units = 'g/(km2.h) '
    PMN.units = 'g/(km2.h) '
    PMOTHR.units = 'g/(km2.h) '
    PNA.units = 'g/(km2.h) '
    PNCOM.units = 'g/(km2.h) '
    PNH4.units = 'g/(km2.h) '
    PNO3.units = 'g/(km2.h) '
    POC.units = 'g/(km2.h) ' 
    PRPA.units = 'moles/(km2.h) '
    PSI.units = 'g/(km2.h) '
    PSO4.units = 'g/(km2.h) '
    PTI.units = 'g/(km2.h) '
    SO2.units = 'moles/(km2.h) '
    SOAALK.units = 'moles/(km2.h) '
    SULF.units = 'moles/(km2.h) '
    TERP.units = 'moles/(km2.h) '
    TOL.units = 'moles/(km2.h) '
    UNK.units = 'moles/(km2.h) '
    UNR.units = 'moles/(km2.h) '
    VOC_INV.units = 'g/(km2.h) '
    XYLMN.units = 'moles/(km2.h) '
    PMFINE.units = 'g/(km2.h) '
    LON.units = 'degrees '
    LAT.units = 'degrees '
    AREA.units = 'km2 '
   
    f2.close()
    print('Your BRAVESdatabase netCDF file is ready!')
