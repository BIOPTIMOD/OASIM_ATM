import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates monthly clim files for OASIM by cutting and interpolating:
    - ERA5 climatology files for tclw tco3
    - MODCLD climatology files for cldtcm cdrem
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--clim_era5', '-e',
                                type = str,
                                required = True,
                                help = '''Clim ERA5 NetCDF dir'''
                                )
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = ''' mask filename'''
                                )

    parser.add_argument(   '--clim_mod', '-M',
                                type = str,
                                required = True,
                                help = '''Clim MODCLD NetCDF dir '''
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

INPUT_ERA5   = addsep(args.clim_era5)
INPUT_MODCLD = addsep(args.clim_mod)
OUTDIR       = addsep(args.outdir)
TheMask      = Mask(args.maskfile)

jpk,jpj,jpi = TheMask.shape

xMin = TheMask.lon[0]
xMax = TheMask.lon[-1]
yMin = TheMask.lat[0]
yMax = TheMask.lat[-1]



def interp(Min,lon,lat):
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
        M2d_orig=netcdf4.readfile(filename, var)[J_start:J_end,I_start:I_end]
        local_lon = lon[I_start:I_end]
        local_lat = lat[J_start:J_end]
        M2d_orig = M2d_orig[-1::-1,:]
        local_lat = local_lat[-1::-1]

    return  interp(M2d_orig,local_lon, local_lat)

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
    setattr(ncvar, 'orig', 'MODCLD' )

    ncvar = ncOUT.createVariable('cldtcm','f',('lat','lon'))
    ncvar[:]=cdrem
    setattr(ncvar, 'long_name',  'TODO' )
    setattr(ncvar, 'units','TODO' )
    setattr(ncvar, 'orig', 'MODCLD')
    
    ncvar = ncOUT.createVariable('tclw','f',('lat','lon'))
    ncvar[:]=tclw
    setattr(ncvar, 'long_name',  'Total column cloud liquid water' )
    setattr(ncvar, 'units','kg m**-2' )
    setattr(ncvar, 'orig', 'ERA5_clim')

    ncvar = ncOUT.createVariable('tco3','f',('lat','lon'))
    ncvar[:]=tco3
    setattr(ncvar, 'long_name',  'Total column ozone' )
    setattr(ncvar, 'units','kg m**-2' )
    setattr(ncvar, 'orig', 'ERA5_clim')

    ncOUT.close()




era5file="%s%02d.%s.nc"  %(INPUT_ERA5,1,'tco3')
lon_era5 = netcdf4.readfile(era5file, 'lon')
lat_era5 = netcdf4.readfile(era5file, 'lat')

modcldfile="%smodcld%02d.nc"  %(INPUT_MODCLD,1)
lon_cld = netcdf4.readfile(modcldfile, 'lon')
lat_cld = netcdf4.readfile(modcldfile, 'lat')


for iframe in range(1,13):
    inputfile = "%s%02d.%s.nc"  %(INPUT_ERA5,iframe,'tco3')
    tco3= getMap(inputfile, 'tco3',lon_era5,lat_era5)
    inputfile = "%s%02d.%s.nc"  %(INPUT_ERA5,iframe,'tclw')
    tclw= getMap(inputfile, 'tclw',lon_era5,lat_era5)

    inputfile = "%smodcld%02d.nc"  %(INPUT_MODCLD,iframe)

    outfile = "%sclimatm.yyyy%02d01-00:00:00.nc" %(OUTDIR,iframe)
    cdrem  = getMap(inputfile, "cdrem" , lon_cld, lat_cld)
    cldtcm = getMap(inputfile, "cldtcm", lon_cld, lat_cld)
    dumpfile(outfile,TheMask,cdrem,cldtcm, tclw,tco3)



