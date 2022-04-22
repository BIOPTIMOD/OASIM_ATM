from commons import netcdf4
import netCDF4
from commons.mask import Mask
import numpy as np
from scipy import interpolate
from commons.Timelist import TimeList

INPUTDIR="/g100_scratch/userexternal/gbolzon0/OASIM/modaer/NC/"
OUTDIR="/g100_scratch/userexternal/gbolzon0/OASIM_ATM/out/"


TheMask=Mask('/g100_scratch/userexternal/gbolzon0/V9C/DEV_OASIM_INPUTS/wrkdir/MODEL/meshmask.nc')
TL=TimeList.fromfilenames(None, INPUTDIR, "modaer201901*", prefix="modaer", dateformat="%Y%m")

jpk,jpj,jpi = TheMask.shape


xMin = TheMask.xlevels[0,0]
xMax = TheMask.xlevels[0,-1]
yMin = TheMask.ylevels[0,0]
yMax = TheMask.ylevels[-1,0]

inputfile=TL.filelist[0]
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


def readfile(filename,maskObj,var):
    mask0=TheMask.mask_at_level(0)
    Mout = np.ones((jpk,jpj,jpi),np.float32)*1.e+20
    M3d_orig=netcdf4.readfile(filename, var)[:,J_start:J_end,I_start:I_end]
    nl,_,_= M3d_orig.shape
    for k in range(nl):
        M2d = interp(M3d_orig[k,:,:])
        M2d[~mask0] = 1.e+20
        sea=M2d[mask0]
        bad = (sea<0) | (sea > 2) | np.isnan(sea)
        if bad.any(): raise ValueError('ferma qua')
        Mout[k,:,:] = M2d
    return Mout

def dumpfile(filename, maskObj, taua,asymp,ssalb):
    ncOUT   = netCDF4.Dataset(filename,"w");

    ncOUT.createDimension('lon',jpi);
    ncOUT.createDimension('lat',jpj);
    ncOUT.createDimension('depth',jpk);
    

    ncvar=ncOUT.createVariable('lon','f',('lon',))
    setattr(ncvar,'units','degrees')
    setattr(ncvar,'long_name'    ,'longitude')
    ncvar[:]=maskObj.xlevels[0,:]
    ncvar=ncOUT.createVariable('lat','f',('lat',))
    setattr(ncvar,'units','degrees')
    setattr(ncvar,'long_name'    ,'latitude')
    ncvar[:]=maskObj.ylevels[:,0]


    ncvar = ncOUT.createVariable('taua','f',('depth','lat','lon'), zlib=True, fill_value=1.0e+20)
    ncvar[:]=taua
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )

    ncvar = ncOUT.createVariable('asymp','f',('depth','lat','lon'), zlib=True, fill_value=1.0e+20)
    ncvar[:]=asymp
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )

    ncvar = ncOUT.createVariable('ssalb','f',('depth','lat','lon'), zlib=True, fill_value=1.0e+20)
    ncvar[:]=ssalb
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )


    ncOUT.close()







for i, inputfile in enumerate(TL.filelist):
    
    outfile = "%saero.%s.nc" %(OUTDIR,TL.Timelist[i].strftime('%Y%m%d-%H:%M:%S') )
    print(outfile)
    taua  = readfile(inputfile,TheMask,'taua')
    asymp = readfile(inputfile,TheMask,'asymp')
    ssalb = readfile(inputfile,TheMask,'ssalb')

    dumpfile(outfile,TheMask,taua,asymp,ssalb)
