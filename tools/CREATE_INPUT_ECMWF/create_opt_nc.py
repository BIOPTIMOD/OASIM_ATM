import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates opt file for OASIM by reading ECMWF file
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputfile', '-i',
                                type = str,
                                required = True,
                                help = ''' Input ECMWF file'''
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


import sys
from datetime import datetime
from commons import netcdf4
import numpy as np
from scipy import interpolate
from commons.mask import Mask
import netCDF4
from commons.utils import addsep
    
def set_nan_values(M, fill_val):
    bad=np.isnan(M) | M < 0
    M[bad] = fill_val


nframes_in_day=4

#d=datetime.strptime(args.date,"%Y%m%d")
#frame_start = ( int(d.strftime("%j")) - 1 ) * nframes_in_day
#frame_end   = frame_start + nframes_in_day



TheMask=Mask(args.maskfile)
jpk,jpj,jpi = TheMask.shape
filename=args.inputfile
OUTDIR = addsep(args.outdir)

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

def dumpfile(filename, sp,si10,t2m,d2m,msl):
    ncOUT   = netCDF4.Dataset(filename,"w");

    ncOUT.createDimension('lon',jpi);
    ncOUT.createDimension('lat',jpj);


    ncvar = ncOUT.createVariable('sp','f',('lat','lon'))
    ncvar[:]=sp
    setattr(ncvar, 'long_name',  'surface pressure' )
    setattr(ncvar, 'units',  '[mb]' )
    
    ncvar = ncOUT.createVariable('wsmorg','f',('lat','lon'))
    ncvar[:] = si10
    setattr(ncvar, 'long_name',  '10 m wind speed' )
    setattr(ncvar, 'units',  '[m.s-1]' )
    
    ncvar = ncOUT.createVariable('msl','f',('lat','lon'))
    ncvar[:] = msl
    
    ncvar = ncOUT.createVariable('d2m','f',('lat','lon'))
    ncvar[:] = d2m

    ncvar = ncOUT.createVariable('t2m','f',('lat','lon'))
    ncvar[:] = t2m

    setattr(ncOUT, 'input_file', args.inputfile)
    ncOUT.close()


for iframe in range(nframes_in_day):
    outfile = OUTDIR + "ave.20190101-00:00:00.nc"
    
    
    msl = getframe(filename,'MSL', iframe)
    sp =  getframe(filename,'SP', iframe)
    u10 = getframe(filename,'U10M', iframe)
    v10 = getframe(filename,'V10M', iframe)
    t2m = getframe(filename,'T2M', iframe)
    d2m = getframe(filename,'D2M', iframe)

    si10 = np.sqrt(u10**2 + v10**2)
    dumpfile(outfile, sp, si10, t2m, d2m,msl)
    
    import sys
    sys.exit()

    
    
    
    
    
# ncin=NC4.Dataset(args.inputfile,"r")
# 
# 
# 
# 
# msl =  ncin.variables['msl'][frame_start:frame_end,1:181,:] # Mean sea level pressurea Pascal
# sp  =  ncin.variables['sp' ][frame_start:frame_end,1:181,:]  # surface pressure Pascal
# u10 =  ncin.variables['u10'][frame_start:frame_end,1:181,:] # 10 metres u wind speed
# v10 =  ncin.variables['v10'][frame_start:frame_end,1:181,:] # 10 metres v wind speed
# # relative humidity and
# # precipitable water (absorption by water vapour)
# t2m  = ncin.variables['t2m'][frame_start:frame_end,1:181,:] # 2m temperature K
# d2m =  ncin.variables['d2m'][frame_start:frame_end,1:181,:] # dew temperature K
# ozo =  ncin.variables['tco3'][frame_start:frame_end,1:181,:] # total column ozone
# 
# 
# slporg=np.zeros(msl.shape)
# sporg =np.zeros(sp.shape)
# 
# nTimes, nLat, nLon  = slporg.shape
# 
# for t in range(nTimes):
#     for j in range(nLat):
#         for i  in range(nLon):
#              slporg[t,j,i] = msl[t,nLat-1-j,i] / 100. # convert to mb
#              sporg[t,j,i]  = sp[t,nLat-1-j,i]  / 100. # convert to mb
# 
# 
# set_nan_values(slporg, -1.0)
# set_nan_values(sporg, -1.0)
# 
# 
# 
# 
# si10 = np.sqrt(u10*u10 + v10*v10)
# wsmorg=np.zeros(si10.shape)
# 
# nTimes, nLat, nLon  = wsmorg.shape
# 
# for t in range(nTimes):
#     for j in range(nLat):
#         for i  in range(nLon):
#              wsmorg[t,j,i] = si10[t,nLat-1-j,i]
# 
# set_nan_values(wsmorg, -1.0)
# 
# 
# 
# a1 = 611.21 # Pascal
# a3 = 17.502 # dimensionless
# a4 = 32.19  # Kelvin
# To = 273.16 # Kelvin
# b1 = 0.14 * 0.01  # cm/Pascal
# b2 = 0.21   # cm
# 
# rhorg=np.zeros(si10.shape)
# wvorg=np.zeros(si10.shape)
# 
# nTimes, nLat, nLon  = rhorg.shape
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
# set_nan_values(rhorg, -1.0)
# set_nan_values(wvorg, -1.0)
# #
# # computation of precipitable water [absorption by water vapour] from saturation water vapour pressure es_T as in Garrison and Adler (1990) and then used in Gregg and Casey (1990)
# #
# # note that saturation water pressure is computed in 3 ways:
# #
# # 1) in Garrison and Adler (1990) by the Tabata relation (Tabata, 1973)
# #    es_T = 10^(8.42926609 - 1827.17843/T - 71208.271/T^2) in mb
# #
# # 2) in Gregg and Carder (1990) they use the method of Lowe (Lowe, 1977)
# #    es_T = 1013.25*exp(13.3185*t - 1.9760 *t^2 - 0.6445*t^3 - 0.1299*t^4) in mb, where t=1-373.16/T
# #
# # 3) in ECMWF model (https://www.ecmwf.int/sites/default/files/elibrary/2015/9211-part-iv-physical-processes.pdf) by the Teten's formula (in Pa, 1 Pa = 0.01 mb) which is the one then used to compute rhorg above
# #
# # however, the three relations give the same results, as shown in saturation-water-vapour-pressure.png
# #
# # !!!! OCCORRE IMPORTARE ANCHE LA SURFACE PRESSURE sfcpr da ECMWF
# # !!!! CONTROLLARE SE SIA NECESSARIO DIVIDERE PER 100 IN QUANTO in GarAdl90 si definisce
# # e = h * es dove e=vapour pressure (ovvero es_Td) e h=relative humidity
# #
# 
# # ozone
# 
# 
# 
# ozorg=np.zeros(si10.shape)
# nTimes, nLat, nLon  = ozorg.shape
# 
# for t in range(nTimes):
#     for j in range(nLat):
#         for i  in range(nLon):
#              ozorg[t,j,i] = ozo[t,nLat-1-j,i] / 2.1414 / 0.00001
# 
# set_nan_values(ozorg, -1.0)
# 
# 
# 
# ncOUT   = NC4.Dataset(args.outfile,"w");
# 
# ncOUT.createDimension('lon',nLon);
# ncOUT.createDimension('lat',nLat);
# ncOUT.createDimension('time',nTimes);
# 
# #   ncvar = ncOUT.createVariable('slporg','f',('time','lat','lon')); ncvar[:] = sporg; 
# 
# ncvar = ncOUT.createVariable('slporg','f',('lat','lon'))
# ncvar[:] = sporg.mean(axis=0)
# setattr(ncvar, 'missing_value',-1.0 )
# setattr(ncvar, 'long_name',  'surface pressure' )
# setattr(ncvar, 'units',  '[mb]' )
# 
# ncvar = ncOUT.createVariable('wsmorg','f',('lat','lon'))
# ncvar[:] = wsmorg.mean(axis=0)
# setattr(ncvar, 'missing_value',-1.0 )
# setattr(ncvar, 'long_name',  '10 m wind speed' )
# setattr(ncvar, 'units',  '[m.s-1]' )
# 
# ncvar = ncOUT.createVariable('rhorg','f',('lat','lon'))
# ncvar[:] =  rhorg.mean(axis=0)
# setattr(ncvar, 'missing_value',-1.0 )
# setattr(ncvar, 'long_name',  'relative humidity');
# setattr(ncvar, 'units',  '[%]' );
# 
# ncvar = ncOUT.createVariable('ozorg' ,'f',('lat','lon'))
# ncvar[:] =  ozorg.mean(axis=0)
# setattr(ncvar, 'missing_value',-1.0 )
# setattr(ncvar, 'long_name',  'ozone');
# setattr(ncvar, 'unit',  '[DU]' );
# 
# 
# ncvar = ncOUT.createVariable('wvorg','f',('lat','lon'))
# ncvar[:] = wvorg.mean(axis=0)
# setattr(ncvar, 'missing_value',-1.0 )
# setattr(ncvar, 'long_name',  'precipitable water')
# setattr(ncvar, 'unit',  '[cm]' )
# 
# 
# setattr(ncOUT, 'input_file', args.inputfile)
# ncOUT.close()


