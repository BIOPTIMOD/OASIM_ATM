import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates monthly clim files for OASIM by cutting and interpolating:
    - ERA5 data files for tclw tco3
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--modis', '-M',
                                type = str,
                                required = True,
                                help = '''MODCLD NetCDF dir '''
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

from commons import netcdf4
import netCDF4
from commons.mask import Mask
import numpy as np
from scipy import interpolate
from commons.utils import addsep
from commons import genUserDateList as DL
from commons.Timelist import TimeList
from commons.utils import Time_Interpolation
from commons.mask import Mask

INPUT_MODCLD = addsep(args.modis)
OUTDIR       = addsep(args.outdir)
TheMask      = Mask(args.maskfile)

jpk,jpj,jpi = TheMask.shape

xMin = TheMask.lon[0]
xMax = TheMask.lon[-1]
yMin = TheMask.lat[0]
yMax = TheMask.lat[-1]



def interp(Min,lon,lat):
#   lon2,lat2=np.meshgrid(lon, lat)
    f = interpolate.interp2d(lon, lat, Min, kind='linear')
    Mout  = f(TheMask.lon, TheMask.lat)
    return Mout


def getMap(filename,var,lon,lat):
    I_start = np.argmin(np.abs(lon - xMin))
    I_end   = np.argmin(np.abs(lon - xMax))
    if lat[0]< lat[-1] : #increasing, MODCLD_files
        J_start = np.argmin(np.abs(lat - yMin))
        J_end   = np.argmin(np.abs(lat - yMax))
        M2d_orig=netcdf4.readfile(filename, var)[0, J_start:J_end,I_start:I_end]
        local_lon = lon[I_start:I_end]
        local_lat = lat[J_start:J_end]
    else: # decreasing, ERA5 files
        J_start = np.argmin(np.abs(lat - yMax))
        J_end   = np.argmin(np.abs(lat - yMin))
        M2d_orig=netcdf4.readfile(filename, var)[0,J_start:J_end,I_start:I_end]
        local_lon = lon[I_start:I_end]
        local_lat = lat[J_start:J_end]
        M2d_orig = M2d_orig[-1::-1,:]
        local_lat = local_lat[-1::-1]

    return  interp(M2d_orig,local_lon, local_lat)

def dumpfile(filename, maskObj, cdrem,cldtcm):
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
    setattr(ncvar, 'orig', 'MODCLD' )

    ncvar = ncOUT.createVariable('cldtcm','f',('lat','lon'))
    ncvar[:]=cldtcm
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )
    setattr(ncvar, 'orig', 'MODCLD')

    ncOUT.close()




TL = TimeList.fromfilenames(None, INPUT_MODCLD, "*nc", prefix="MODIS", dateformat="%Y%m%d-%H:%M:%S")

modisfile0=TL.filelist[0]
lon_modis = netcdf4.readfile(modisfile0, 'lon')
lat_modis = netcdf4.readfile(modisfile0, 'lat')
print(lon_modis)
print(lat_modis)

for it,inputfile in enumerate(TL.filelist[:10]):

    print('Processing... ' + inputfile)
    cdrem  = getMap(inputfile, "cdrem" , lon_modis, lat_modis)
    cldtcm = getMap(inputfile, "cldtcm", lon_modis, lat_modis)

    current_date=TL.Timelist[it].strftime('%Y%m%d-%H:%M:%S')
    outfile = "%sMODIS_MED.%s.nc" %(OUTDIR,current_date)

    dumpfile(outfile,TheMask, cdrem,cldtcm)



