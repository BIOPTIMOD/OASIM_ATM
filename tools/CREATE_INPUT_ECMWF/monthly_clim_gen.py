from commons import netcdf4
import netCDF4
from commons.mask import Mask
import numpy as np
from scipy import interpolate

inputfile="/g100_scratch/userexternal/plazzari/OASIM/MODIS_DATA/modcld0000.nc"
input_clim="/g100_scratch/userexternal/gbolzon0/OASIM_ATM/tools/CREATE_INPUT_ECMWF/ERA5_Med.nc"
TheMask=Mask('/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/V8C/meshmask.nc')
OUTDIR="/g100_scratch/userexternal/gbolzon0/OASIM_ATM/out/"
jpk,jpj,jpi = TheMask.shape

xMin = TheMask.xlevels[0,0]
xMax = TheMask.xlevels[0,-1]
yMin = TheMask.ylevels[0,0]
yMax = TheMask.ylevels[-1,0]

lon = netcdf4.readfile(inputfile, 'lon')
lat = netcdf4.readfile(inputfile, 'lat')

I_start = np.argmin(np.abs(lon - xMin))
I_end   = np.argmin(np.abs(lon - xMax))
J_start = np.argmin(np.abs(lat - yMin))
J_end   = np.argmin(np.abs(lat - yMax))

lon = lon[I_start:I_end]
lat = lat[J_start:J_end]

def interp(Min):
    f = interpolate.interp2d(lon, lat, Min, kind='linear')
    Mout  = f(TheMask.xlevels[0,:], TheMask.ylevels[:,0])
    return Mout


def getframe(filename,var, timeframe):
    M2d_orig=netcdf4.readfile(filename, var)[timeframe,J_start:J_end,I_start:I_end]    
    return  interp(M2d_orig)

def dumpfile(filename, maskObj, cdrem,cldtcm,tclw,tco3):
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


    ncvar = ncOUT.createVariable('cdrem','f',('lat','lon'))
    ncvar[:]=cdrem
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )
    setattr(ncvar, 'inputfile', inputfile )

    ncvar = ncOUT.createVariable('cldtcm','f',('lat','lon'))
    ncvar[:]=cdrem
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )
    setattr(ncvar, 'inputfile', inputfile)
    
    ncvar = ncOUT.createVariable('tclw','f',('lat','lon'))
    ncvar[:]=tclw
    setattr(ncvar, 'long_name',  'Total column cloud liquid water' )
    setattr(ncvar, 'units','kg m**-2' )
    setattr(ncvar, 'inputfile', input_clim)

    ncvar = ncOUT.createVariable('tco3','f',('lat','lon'))
    ncvar[:]=tco3
    setattr(ncvar, 'long_name',  'Total column ozone' )
    setattr(ncvar, 'units','kg m**-2' )
    setattr(ncvar, 'inputfile', input_clim)

    ncOUT.close()


tclw = np.ones((jpj,jpi),dtype=np.float32)*0.05
tco3 = np.ones((jpj,jpi),dtype=np.float32)*0.0065



for iframe in range(12):
    month='%02d' %iframe 
    outfile = "%satm.yyyy%02d01-00:00:00.nc" %(OUTDIR,iframe) 
    cdrem = getframe(inputfile, "cdrem", iframe)
    cldtcm = getframe(inputfile, "cldtcm", iframe)
    dumpfile(outfile,TheMask,cdrem,cldtcm, tclw,tco3)



