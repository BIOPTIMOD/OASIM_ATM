import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates daily files for OASIM by reading ECMWF file provided by CMCC
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputfile', '-i',
                                type = str,
                                required = True,
                                help = ''' Input CMCC file'''
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = ''' mask filename'''
                                )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' path of the output optical dir '''
                                )
    return parser.parse_args()


args = argument()


import os
from datetime import datetime, timedelta
from commons import netcdf4
import numpy as np
from scipy import interpolate
from commons.mask import Mask
import netCDF4
from commons.utils import addsep



TheMask=Mask(args.maskfile)
jpk,jpj,jpi = TheMask.shape
filename=args.inputfile
OUTDIR = addsep(args.outdir)


ncIN=netCDF4.Dataset(filename)
nframes_in_day = ncIN.dimensions['time'].size
ncIN.close()
deltaH = int(24/nframes_in_day)
yyyymmdd=os.path.basename(filename)[:8]



lon= netcdf4.readfile(filename, 'lon')
lat= netcdf4.readfile(filename, 'lat')



xMin = TheMask.xlevels[0,0]
xMax = TheMask.xlevels[0,-1]
yMin = TheMask.ylevels[0,0]
yMax = TheMask.ylevels[-1,0]


I_start = np.argmin(np.abs(lon - xMin))
I_end   = np.argmin(np.abs(lon - xMax))
J_start = np.argmin(np.abs(lat - yMax))
J_end   = np.argmin(np.abs(lat - yMin))

lon = lon[I_start:I_end]
lat = lat[J_start:J_end]
lat = lat[-1::-1] # reverse


def readframe(filename,var, timeframe):
    A=netcdf4.readfile(filename, var)[timeframe,J_start:J_end,I_start:I_end]
    return A[-1::-1,:]

def interp(Min):
    f = interpolate.interp2d(lon, lat, Min, kind='linear')
    Mout  = f(TheMask.xlevels[0,:], TheMask.ylevels[:,0])
    return Mout

def getframe(filename,var, timeframe):
    M2d_orig = readframe(filename, var, timeframe)
    return  interp(M2d_orig)

def dumpfile(filename, maskObj, sp,msl, t2m,d2m, tcc,w10):
    ncOUT   = netCDF4.Dataset(filename,"w");

    ncOUT.createDimension('lon',jpi);
    ncOUT.createDimension('lat',jpj);

    ncvar=ncOUT.createVariable('lon','f',('lon',))
    setattr(ncvar,'units','degrees')
    setattr(ncvar,'long_name'    ,'longitude')
    ncvar[:]=maskObj.xlevels[0,:]
    ncvar=ncOUT.createVariable('lat','f',('lat',))
    setattr(ncvar,'units','degrees')
    setattr(ncvar,'long_name'    ,'latitude')
    ncvar[:]=maskObj.ylevels[:,0]


    ncvar = ncOUT.createVariable('sp','f',('lat','lon'))
    ncvar[:]=sp
    setattr(ncvar, 'long_name',  'Surface pressure' )
    setattr(ncvar, 'units','Pa' )
    setattr(ncvar, 'code', 34 )
    

    ncvar = ncOUT.createVariable('msl','f',('lat','lon'))
    setattr(ncvar, 'units','Pa' )
    setattr(ncvar,'long_name','Mean sea-level pressure')
    setattr(ncvar,'code',151)
    ncvar[:] = msl

    ncvar = ncOUT.createVariable('d2m','f',('lat','lon'))
    setattr(ncvar,'units','K')
    setattr(ncvar,'long_name','2 metre dewpoint temperature')
    setattr(ncvar, 'code', 168)
    ncvar[:] = d2m

    ncvar = ncOUT.createVariable('t2m','f',('lat','lon'))
    setattr(ncvar,'units','K')
    setattr(ncvar,'long_name','2 metre temperature')
    setattr(ncvar, 'code', 167)
    ncvar[:] = t2m

    ncvar = ncOUT.createVariable('tcc','f',('lat','lon'))
    setattr(ncvar,'units','%')
    setattr(ncvar,'long_name','Total cloud cover')
    setattr(ncvar, 'code', 164)
    ncvar[:] = tcc

    ncvar = ncOUT.createVariable('w10','f',('lat','lon'))
    setattr(ncvar,'units','m/s')
    setattr(ncvar,'long_name','10 metre wind speed module')
    setattr(ncvar,'code', '165 and 166')
    ncvar[:] = w10

    setattr(ncOUT, 'input_file', args.inputfile)
    ncOUT.close()


for iframe in range(nframes_in_day):
    d=datetime.strptime(yyyymmdd,'%Y%m%d') + timedelta(hours = deltaH*iframe)
    outfile = OUTDIR + d.strftime("atm.%Y%m%d-%H:%M:%S.nc")
    print(outfile)
    
    
    msl = getframe(filename,'MSL' , iframe)
    sp =  getframe(filename,'SP'  , iframe)
    u10 = getframe(filename,'U10M', iframe)
    v10 = getframe(filename,'V10M', iframe)
    t2m = getframe(filename,'T2M' , iframe)
    d2m = getframe(filename,'D2M' , iframe)
    tcc = getframe(filename,'TCC' , iframe)

    w10 = np.sqrt(u10**2 + v10**2)
    dumpfile(outfile, TheMask, sp,msl, t2m,d2m, tcc,w10)
    

# a1 = 611.21 # Pascal
# a3 = 17.502 # dimensionless
# a4 = 32.19  # Kelvin
# To = 273.16 # Kelvin
# b1 = 0.14 * 0.01  # cm/Pascal
# b2 = 0.21   # cm
# 
# 
# for t in range(nTimes):
#     for j in range(nLat):
#         for i  in range(nLon):
#              T            = t2m[t,nLat-1-j,i]
#              Td           = d2m[t,nLat-1-j,i]
#              es_Td        = a1 * np.exp( a3 * (Td-To)/(Td-a4) )
#              es_T         = a1 * np.exp( a3 *  (T-To)/(T-a4) )
#              rhorg[t,j,i] = 100. * es_Td /es_T
#              wvorg[t,j,i] = b1 * es_Td * sporg[t,j,i]/slporg[t,j,i] + b2
# 

#
# computation of precipitable water [absorption by water vapour] from saturation water vapour pressure es_T as in Garrison and Adler (1990) and then used in Gregg and Casey (1990)
#
# note that saturation water pressure is computed in 3 ways:
#
# 1) in Garrison and Adler (1990) by the Tabata relation (Tabata, 1973)
#    es_T = 10^(8.42926609 - 1827.17843/T - 71208.271/T^2) in mb
#
# 2) in Gregg and Carder (1990) they use the method of Lowe (Lowe, 1977)
#    es_T = 1013.25*exp(13.3185*t - 1.9760 *t^2 - 0.6445*t^3 - 0.1299*t^4) in mb, where t=1-373.16/T
#
# 3) in ECMWF model (https://www.ecmwf.int/sites/default/files/elibrary/2015/9211-part-iv-physical-processes.pdf) by the Teten's formula (in Pa, 1 Pa = 0.01 mb) which is the one then used to compute rhorg above
#
# however, the three relations give the same results, as shown in saturation-water-vapour-pressure.png
#
# !!!! OCCORRE IMPORTARE ANCHE LA SURFACE PRESSURE sfcpr da ECMWF
# !!!! CONTROLLARE SE SIA NECESSARIO DIVIDERE PER 100 IN QUANTO in GarAdl90 si definisce
# e = h * es dove e=vapour pressure (ovvero es_Td) e h=relative humidity
