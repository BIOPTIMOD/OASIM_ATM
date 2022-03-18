import netCDF4 as NC
import numpy as np
from commons.mask import Mask

def readfile(filename, var):
    '''
    Generic file reader
    '''
    D = NC.Dataset(filename, "r")
    VAR = np.array(D.variables[var])
    D.close()
    return VAR

    
hourly_file="/g100_scratch/userexternal/gbolzon0/OASIM_ATM/out/atm.20190101-18:00:00.nc"
clim_file  ="/g100_scratch/userexternal/gbolzon0/OASIM_ATM/out/atm.yyyy0101-00:00:00.nc"
TheMask    ="/g100_scratch/userexternal/gbolzon0/DEBUG/OGSTM/wrkdir/MODEL/meshmask.nc"


sp     = readfile(hourly_file,  'sp')
msl    = readfile(hourly_file, 'mls')
d2m    = readfile(hourly_file, 'd2m')
t2m    = readfile(hourly_file, 't2m')
tcc    = readfile(hourly_file, 'tcc')
ws10   = readfile(hourly_file, 'w10')
tclw   = readfile(clim_file,  'tclw')
tco3   = readfile(clim_file,  'tco3')
cdrem  = readfile(clim_file, 'cdrem')
cldtcm = readfile(clim_file,'cldtcm')
lon = TheMask.xlevels()
lat = TheMask.ylevels()


j=100
result = oasim_lib(year, month, day, sec_start, sec_end, lon[j,:], lat[j,:], sp[j,:],msl[j,:],d2m[j,:],t2m[j,:],tcc[j,:],ws10[j,:],tclw[j,:],tco3[j,:],cdrem[j,:],cldtcm[j,:], aerosol)
#result = oasim_lib(time, lon, lat, sp,msl,d2m,t2m,tcc,ws10,tclw,tco3,cdrem,cldtcm,aerosol)
