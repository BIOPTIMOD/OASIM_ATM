import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates monthly clim files for OASIM by cutting and interpolating:
    - ERA5 data files for tclw tco3
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--era5', '-e',
                                type = str,
                                required = True,
                                help = '''Clim ERA5 NetCDF dir'''
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

INPUT_ERA5   = addsep(args.era5)
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

def dumpfile(filename, maskObj, sp,msl,u10,v10,tco3,t2m,d2m,tcc,tclw):
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
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )
    setattr(ncvar, 'orig', 'ERA5' )

    ncvar = ncOUT.createVariable('msl','f',('lat','lon'))
    ncvar[:]=msl
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )
    setattr(ncvar, 'orig', 'ERA5')

    ncvar = ncOUT.createVariable('u10','f',('lat','lon'))
    ncvar[:]=u10
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','m s*-1' )
    setattr(ncvar, 'orig', 'ERA5')

    ncvar = ncOUT.createVariable('v10','f',('lat','lon'))
    ncvar[:]=v10
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','m s*-1' )
    setattr(ncvar, 'orig', 'ERA5')
    
    ncvar = ncOUT.createVariable('tclw','f',('lat','lon'))
    ncvar[:]=tclw
    setattr(ncvar, 'long_name',  'Total column cloud liquid water' )
    setattr(ncvar, 'units','kg m**-2' )
    setattr(ncvar, 'orig', 'ERA5')

    ncvar = ncOUT.createVariable('tco3','f',('lat','lon'))
    ncvar[:]=tco3
    setattr(ncvar, 'long_name',  'Total column ozone' )
    setattr(ncvar, 'units','kg m**-2' )
    setattr(ncvar, 'orig', 'ERA5')

    ncvar = ncOUT.createVariable('t2m','f',('lat','lon'))
    ncvar[:]=t2m
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )
    setattr(ncvar, 'orig', 'ERA5')

    ncvar = ncOUT.createVariable('d2m','f',('lat','lon'))
    ncvar[:]=d2m
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )
    setattr(ncvar, 'orig', 'ERA5')

    ncvar = ncOUT.createVariable('tcc','f',('lat','lon'))
    ncvar[:]=tcc
    setattr(ncvar, 'long_name',  'Total cloud cover' )
    setattr(ncvar, 'units','[-]' )
    setattr(ncvar, 'orig', 'ERA5')

    ncOUT.close()




TL = TimeList.fromfilenames(None, INPUT_ERA5, "*nc", prefix="ERA5_MED", dateformat="%Y%m%d-%H:%M:%S")

era5file0=TL.filelist[0]
lon_era5 = netcdf4.readfile(era5file0, 'longitude')
lat_era5 = netcdf4.readfile(era5file0, 'latitude')
print(lon_era5)
print(lat_era5)

for it,inputfile in enumerate(TL.filelist[:10]):

    print('Processing... ' + inputfile)
    sp   = getMap(inputfile, 'sp',lon_era5,lat_era5)
    msl  = getMap(inputfile, 'msl',lon_era5,lat_era5)
    u10  = getMap(inputfile, 'u10',lon_era5,lat_era5)
    v10  = getMap(inputfile, 'v10',lon_era5,lat_era5)
    tco3 = getMap(inputfile, 'tco3',lon_era5,lat_era5)
    t2m  = getMap(inputfile, 't2m',lon_era5,lat_era5)
    d2m  = getMap(inputfile, 'd2m',lon_era5,lat_era5)
    tcc  = getMap(inputfile, 'tcc',lon_era5,lat_era5)
    tclw = getMap(inputfile, 'tclw',lon_era5,lat_era5)

    current_date=TL.Timelist[it].strftime('%Y%m%d-%H:%M:%S')
    outfile = "%sERA5_MED.%s.nc" %(OUTDIR,current_date)

    dumpfile(outfile,TheMask,sp,msl,u10,v10,tco3,t2m,d2m,tcc,tclw)



