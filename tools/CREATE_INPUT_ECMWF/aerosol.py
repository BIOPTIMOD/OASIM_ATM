import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates monthly files for OASIM by cutting and interpolating
    NASA modaer files for: taua asymp ssalb
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '''Clim Nasa NetCDF dir'''
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
    parser.add_argument(   '--clim',
                                action = "store_true",
                                help = ''' To work in monthly climatological mode'''
                                )

    return parser.parse_args()


args = argument()


from commons.Timelist import TimeList
from commons import netcdf4
import netCDF4
from commons.mask import Mask
import numpy as np
from scipy import interpolate
from commons.interpolators import shift
from commons.utils import addsep

INPUTDIR=addsep(args.inputdir)
OUTDIR  =addsep(args.outdir)


TheMask=Mask(args.maskfile)

jpk,jpj,jpi = TheMask.shape

def smoother(M2d,mask2d):
    NsmoothingCycles = 20
    jpj,jpi = mask2d.shape
    M2d[~mask2d] = np.nan
    auxs = np.zeros((5,jpj,jpi),np.float32)
    for _ in range(NsmoothingCycles):
        auxs[0,:,:] = M2d
        auxs[1,:,:] = shift(M2d,2,'r')
        auxs[2,:,:] = shift(M2d,2,'l')
        auxs[3,:,:] = shift(M2d,2,'u')
        auxs[4,:,:] = shift(M2d,2,'l')
        M2d = np.nanmean(auxs,axis=0)
    return M2d

def uniform(M2d,mask2d):
    ii = mask2d & (M2d>0)
    meanvalue = M2d[ii].mean()
    return np.ones_like(M2d) * meanvalue


xMin = TheMask.xlevels[0,0]
xMax = TheMask.xlevels[0,-1]
yMin = TheMask.ylevels[0,0]
yMax = TheMask.ylevels[-1,0]

if args.clim:
    inputfile="%smodaer%02d.nc" %(INPUTDIR, 1)
else:
    TL = TimeList.fromfilenames(None, INPUTDIR, "*nc", prefix="modaer", dateformat="%Y%m")
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
    M3d_orig[M3d_orig < -99] = 0.5
    nl,_,_= M3d_orig.shape
    for k in range(nl):
        M2d = interp(M3d_orig[k,:,:])
        #M2d = smoother(M2d, mask0)
        #M2d = uniform(M2d, mask0)
        sea=M2d[mask0]
        negative = sea<0
        if negative.any():
            sea[negative]=0.5
            print(filename, var, 'jl=', k)
            print(negative.sum(), "negative values corrected with 0.5")
        nans = np.isnan(sea)
        if nans.any(): raise ValueError('Nans: stop')
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
    setattr(ncvar, 'long_name',  'aerosol optical thickness' )
    setattr(ncvar, 'units','[-]' )
    setattr(ncvar, 'orig', 'MODIS_AEROSOL' )

    ncvar = ncOUT.createVariable('asymp','f',('depth','lat','lon'), zlib=True, fill_value=1.0e+20)
    ncvar[:]=asymp
    setattr(ncvar, 'long_name',  'aerosol asymmetry parameter' )
    setattr(ncvar, 'units','[-]' )
    setattr(ncvar, 'orig', 'MODIS_AEROSOL' )

    ncvar = ncOUT.createVariable('ssalb','f',('depth','lat','lon'), zlib=True, fill_value=1.0e+20)
    ncvar[:]=ssalb
    setattr(ncvar, 'long_name',  'aerosol single scattering albedo' )
    setattr(ncvar, 'units','[-]' )
    setattr(ncvar, 'orig', 'MODIS_AEROSOL' )


    ncOUT.close()



if args.clim:
    for iframe in range(1,13):
        inputfile = "%smodaer%02d.nc" %(INPUTDIR, iframe)
        outfile = "%saero.yyyy%02d15-00:00:00.nc" %(OUTDIR,iframe)
        print(outfile)
        taua  = readfile(inputfile,TheMask,'taua')
        asymp = readfile(inputfile,TheMask,'asymp')
        ssalb = readfile(inputfile,TheMask,'ssalb')

        dumpfile(outfile,TheMask,taua,asymp,ssalb)
else:

    dateformat="%Y%m%d-%H:%M:%S"
    for it, inputfile in enumerate(TL.filelist):
        outfile = "%saero.%s.nc" %(OUTDIR, TL.Timelist[it].strftime(dateformat))
        print(outfile)
        taua  = readfile(inputfile,TheMask,'taua')
        asymp = readfile(inputfile,TheMask,'asymp')
        ssalb = readfile(inputfile,TheMask,'ssalb')

        dumpfile(outfile,TheMask,taua,asymp,ssalb)
