from commons import genUserDateList as DL
from commons.Timelist import TimeList, TimeInterval
from commons import timerequestors
from commons import netcdf4
import netCDF4
import numpy as np
import glob


DIR="/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/ERA5/"
OUTDIR="/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/CLIM_ERA5/"
firstfile = glob.glob(DIR + "*nc")[0]
lon = netcdf4.readfile(firstfile, 'longitude')
lat = netcdf4.readfile(firstfile, 'latitude')
jpi=len(lon)
jpj=len(lat)

def dumpfile(filename, M, varname):
    ncOUT   = netCDF4.Dataset(filename,"w");

    ncOUT.createDimension('lon',jpi);
    ncOUT.createDimension('lat',jpj);

    ncvar=ncOUT.createVariable('lon','f',('lon',))
    ncvar[:]=lon
    ncvar=ncOUT.createVariable('lat','f',('lat',))
    ncvar[:]=lat


    ncvar = ncOUT.createVariable(varname,'f',('lat','lon'))
    ncvar[:]=M

    ncOUT.close()



datestart="20100101-00:00:00"
date__end="20211231-23:00:00"


datelist=DL.getTimeList(datestart, date__end, hours=1)

TL = TimeList(datelist)
nFrames = TL.nTimes

M = np.zeros((nFrames,jpj,jpi),np.float32)

YEARLY_REQS=TL.getYearlist()

for var in ["tco3","tclw"] :
    counter=0
    # queueing all matrices
    for req in YEARLY_REQS:
        ii, w = TL.select(req)
        if len(ii)<300: continue
        filename=DIR + "ERA5_Med_" + req.string+ ".nc"
        print("reading " + filename)

        A=netcdf4.readfile(filename, var)
        yearly_frames, _,_ = A.shape
        M[counter:counter+yearly_frames,:] = A
        counter=counter+yearly_frames
    
    
    for month in range(1,13):
        req = timerequestors.Clim_month(month)
        ii,w = TL.select(req)
        Monthy_Mean =  M[ii,:].mean(axis=0)
        outfile = OUTDIR + req.string + "." + var + ".nc"
        print(outfile)
        dumpfile(outfile, Monthy_Mean, var)
        
