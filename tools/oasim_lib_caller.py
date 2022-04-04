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
aero_file  ="/g100_scratch/userexternal/gbolzon0/OASIM_ATM/out/aero.20190101-00:00:00.nc"
TheMask    = Mask("/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/V8C/meshmask.nc")


sp     = readfile(hourly_file,  'sp')
msl    = readfile(hourly_file, 'msl')
d2m    = readfile(hourly_file, 'd2m')
t2m    = readfile(hourly_file, 't2m')
tcc    = readfile(hourly_file, 'tcc')
ws10   = readfile(hourly_file, 'w10')

tclw   = readfile(clim_file,  'tclw')
tco3   = readfile(clim_file,  'tco3')
cdrem  = readfile(clim_file, 'cdrem')
cldtcm = readfile(clim_file,'cldtcm')

taua   = readfile(aero_file, 'taua')[ :33,:]
asymp  = readfile(aero_file, 'asymp')[:33,:]
ssalb  = readfile(aero_file, 'ssalb')[:33,:]

lon = TheMask.xlevels
lat = TheMask.ylevels
year=2019
month=1
day=1
sec_start=66600 # 18:30
sec_end=70200   # 19:30

i=388
j=306
print("lon " , lon[j,i])
print("lat " , lat[j,i])
print("sp " , sp[j,i])
print("msl " , msl[j,i])
print("d2m " , d2m[j,i])
print("t2m " , t2m[j,i])
print("tcc " , tcc[j,i])
print("ws10 " , ws10[j,i])
print("tclw " , tclw[j,i])
print("tco3 " , tco3[j,i])
print("cdrem " , cdrem[j,i])
print("cldtcm " , cldtcm[j,i])
print("taua " , taua[:,j,i])
print("asymp " , asymp[:,j,i])
print("ssalb " , ssalb[:,j,i])

#result = oasim_lib(year, month, day, sec_start, sec_end, lon[j,:], lat[j,:], sp[j,:],msl[j,:],d2m[j,:],t2m[j,:],tcc[j,:],ws10[j,:],tclw[j,:],tco3[j,:],cdrem[j,:],cldtcm[j,:], taua[:,j,:],asymp[:,j,:], ssalb[:,j,:])
#result = oasim_lib(year, month, day, sec_start, sec_end, lon, lat, sp,msl,d2m,t2m,tcc,ws10,tclw,tco3,cdrem,cldtcm, taua,asymp, ssalb)
