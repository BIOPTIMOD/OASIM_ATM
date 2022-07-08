#!/usr/bin/python3
import sys
import ctypes as ct
import numpy as np
import datetime
from commons import genUserDateList as DL
from commons.Timelist import TimeList
from commons.utils import Time_Interpolation
from commons.mask import Mask
import pandas as pd
import netCDF4 as NC

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False

print("Init ...", flush=True)
print(isParallel, flush=True)

class oasim_lib:
    def __init__(self, lib_path, config_path, lat, lon):
        self._lib = ct.cdll.LoadLibrary(lib_path)
        self._init_lib = self._lib.py_init_lib
        self._init_lib.restype = ct.c_void_p
        self._init_lib.argtypes = [ct.c_int,
                                   ct.c_int,
                                   ct.c_char_p,
                                   np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                   np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS")]

        self._finalize_lib = self._lib.py_finalize_lib
        self._finalize_lib.restype = None
        self._finalize_lib.argtypes = [ct.c_void_p]

        error = ct.c_bool()
        if lat.shape != lon.shape or len(lat.shape) != 1:
            raise ValueError("The lat and lon mesh are incompatible")

        self._ptr = self._init_lib(len(config_path),
                                   lat.shape[0],
                                   config_path.encode("utf-8"),
                                   lat,
                                   lon)

        if not self._ptr:
            raise FileNotFoundError("Corrupted file: {}".format(config_path))
    

    def __del__(self):
        if self._ptr:
            self._finalize_lib(self._ptr)
        


class calc_unit:
    def __init__(self, p_size, lib):
        self._lib = lib
        self._init_calc = self._lib._lib.py_init_calc
        self._init_calc.restype = ct.c_void_p
        self._init_calc.argtypes = [ct.c_int,
                                    ct.c_void_p]

        self._monrad = self._lib._lib.py_monrad
        self._monrad.restype = ct.c_bool
        self._monrad.argtypes = [ct.c_void_p,
                                 ct.c_int,
                                 ct.c_int,
                                 np.ctypeslib.ndpointer(ct.c_int, flags="F_CONTIGUOUS"),
                                 ct.c_int,
                                 ct.c_int,
                                 ct.c_double,
                                 ct.c_double,
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS")]
                                 

        self._finalize_calc = self._lib._lib.py_finalize_calc
        self._finalize_calc.restype = None
        self._finalize_calc.argtypes = [ct.c_void_p]
        
        self._ptr = self._init_calc(p_size, self._lib._ptr)


    def monrad(self, points, iyr, iday, sec_b, sec_e,
               sp, msl, ws10, tco3, t2m, d2m, tcc, tclw, cdrem,
               taua, asymp, ssalb):
        p_size, rows = taua.shape
        edout = np.zeros((p_size, rows), dtype=np.float64, order="F")
        esout = np.zeros((p_size, rows), dtype=np.float64, order="F")

        error = self._monrad(self._ptr, p_size, rows, 
                             np.asfortranarray(points, dtype=np.int32),
                             iyr, iday, sec_b, sec_e,
                             np.asfortranarray(sp, dtype=np.float64), 
                             np.asfortranarray(msl, dtype=np.float64), 
                             np.asfortranarray(ws10, dtype=np.float64), 
                             np.asfortranarray(tco3, dtype=np.float64), 
                             np.asfortranarray(t2m, dtype=np.float64), 
                             np.asfortranarray(d2m, dtype=np.float64), 
                             np.asfortranarray(tcc, dtype=np.float64),
                             np.asfortranarray(tclw, dtype=np.float64), 
                             np.asfortranarray(cdrem, dtype=np.float64),
                             np.asfortranarray(taua, dtype=np.float64), 
                             np.asfortranarray(asymp, dtype=np.float64), 
                             np.asfortranarray(ssalb, dtype=np.float64), 
                             edout, esout)
        
        if error:
            raise ValueError("Error in monrad computation.")

        return edout, esout


    def __del__(self):
        if self._ptr:
            self._finalize_calc(self._ptr)

def Load_Data(DATADIR, TimeList,before,after,prefix,var ):
    Before_date17 = TimeList[before].strftime(dateFormat)
    After__date17 = TimeList[after ].strftime(dateFormat)
    Before_File = prefix + Before_date17 + ".nc"
    After__File = prefix + After__date17 + ".nc"
    ncB = NC.Dataset(DATADIR + Before_File,'r')
    ncA = NC.Dataset(DATADIR + After__File,'r')

    if prefix == 'ERA5_MED.':
        Before_Data = ncB.variables[var][:,:].copy().astype(np.float32)
        After__Data = ncA.variables[var][:,:].copy().astype(np.float32)
    elif prefix == 'MODIS_CLD_MED.':
        Before_Data = ncB.variables[var][:,:].copy().astype(np.float32)
        After__Data = ncA.variables[var][:,:].copy().astype(np.float32)
    elif prefix == 'MODIS_AEROSOL_MED.':
        Before_Data = ncB.variables[var][:,:,:].copy().astype(np.float32)
        After__Data = ncA.variables[var][:,:,:].copy().astype(np.float32)

    ncA.close()
    ncB.close()
    return Before_Data, After__Data

dateFormat="%Y%m%d-%H:%M:%S"

#import timelist from data
datadir='/g100_scratch/userexternal/plazzari/TEST_ERA5/READY_FOR_MODEL/'
TL_ERA5 = TimeList.fromfilenames(None, datadir, "ERA5_MED*nc", prefix="ERA5_MED.", dateformat="%Y%m%d-%H:%M:%S")
TL_MODIS_AER = TimeList.fromfilenames(None, datadir, "MODIS_AEROSOL_MED*nc", prefix="MODIS_AEROSOL_MED.", dateformat="%Y%m%d-%H:%M:%S")
TL_MODIS_CLD = TimeList.fromfilenames(None, datadir, "MODIS_CLD_MED*nc", prefix="MODIS_CLD_MED.", dateformat="%Y%m%d-%H:%M:%S")

#output time iframes
dtnew = 60*60. #1 hour in seconds
OUT_Timelist = DL.getTimeList("20190115-00:00:00","20211215-00:00:00", seconds=dtnew) #datelist to be interpolated

#import wavelenghts
wl = pd.read_csv('../data/bin.txt', delim_whitespace=True, header=None).to_numpy()
wl = np.mean(wl,1).astype(int)

location='BOUSSOLE'

if location == 'BOUSSOLE':
    lat = np.array([43.367], dtype=np.float64, order='F') # Boussole
    lon = np.array([7.9], dtype=np.float64, order='F')    # Boussole
elif location == 'LAMPEDUSA':
    lat = np.array([35.49], dtype=np.float64, order='F') # Lampedusa
    lon = np.array([12.47], dtype=np.float64, order='F') # Lampedusa

olib = oasim_lib("../../OASIMlib/liboasim-py.so", "config.yaml", lat, lon)
cunit = calc_unit(1, olib)
dt = (OUT_Timelist[1]-OUT_Timelist[0]).total_seconds() #seconds

#spatial mash
from commons.mask import Mask
TheMask = Mask('/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc')
ji_sel, jj_sel = TheMask.convert_lon_lat_to_indices(lon[0], lat[0])
print(ji_sel,jj_sel)

points = np.array([1])
npoints=points.shape[0]
nframe = len(OUT_Timelist)
print(nframe)
nwavelengths = len(wl)
edout = np.zeros((nframe,nwavelengths))
esout = np.zeros((nframe,nwavelengths))
year = np.zeros(nframe, dtype=int)
month = np.zeros(nframe, dtype=int)
day = np.zeros(nframe, dtype=int)
hour = np.zeros(nframe, dtype=int)
minute = np.zeros(nframe, dtype=int)

INPUT_DIR='/g100_scratch/userexternal/plazzari/TEST_ERA5/READY_FOR_MODEL/'

sp=np.zeros(nframe)
msl=np.zeros(nframe)
u10=np.zeros(nframe)
v10=np.zeros(nframe)
ws10=np.zeros(nframe)
tco3=np.zeros(nframe)
t2m=np.zeros(nframe)
d2m=np.zeros(nframe)
tcc=np.zeros(nframe)
tclw=np.zeros(nframe)
cdrem=np.zeros(nframe)
#ttlist=np.zeros(nframe, dtype=str)
ttlist=["" for x in range(nframe)]
taua=np.zeros((npoints,nframe,33))
asymp=np.zeros((npoints,nframe,33))
ssalb=np.zeros((npoints,nframe,33))


for it,tt in enumerate(OUT_Timelist[:]):
   i = it % nranks
   if i == rank :
       ttlist[it] = str(tt)

########### ERA5    
       beforeERA5,afterERA5,T_interpERA5 = Time_Interpolation(tt,TL_ERA5.Timelist)

#sp
#sp    = f.variables['sp'][:,:,:]#np.array([102377.13])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_ERA5.Timelist, beforeERA5, afterERA5,"ERA5_MED.","sp")
       sp[it] = (1-T_interpERA5)*Before_DATA[jj_sel,ji_sel] + T_interpERA5*After__DATA[jj_sel,ji_sel]
#msl
#msl   = f.variables['msl'][:,:,:]#np.array([102410.06])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_ERA5.Timelist, beforeERA5, afterERA5,"ERA5_MED.","msl")
       msl[it] = (1-T_interpERA5)*Before_DATA[jj_sel,ji_sel] + T_interpERA5*After__DATA[jj_sel,ji_sel]
#u10
#u10   = f.variables['u10'][:,:,:]#np.array([10.0])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_ERA5.Timelist, beforeERA5, afterERA5,"ERA5_MED.","u10")
       u10[it] = (1-T_interpERA5)*Before_DATA[jj_sel,ji_sel] + T_interpERA5*After__DATA[jj_sel,ji_sel]
#v10
#v10   = f.variables['v10'][:,:,:]#np.array([10.0])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_ERA5.Timelist, beforeERA5, afterERA5,"ERA5_MED.","v10")
       v10[it] = (1-T_interpERA5)*Before_DATA[jj_sel,ji_sel] + T_interpERA5*After__DATA[jj_sel,ji_sel]
#ws10
       ws10[it]  = np.sqrt(u10[it]*u10[it]+v10[it]*v10[it])
#tco3
#tco3  = f.variables['tco3'][:,:,:]#np.array([0.0065])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_ERA5.Timelist, beforeERA5, afterERA5,"ERA5_MED.","tco3")
       tco3[it] = (1-T_interpERA5)*Before_DATA[jj_sel,ji_sel] + T_interpERA5*After__DATA[jj_sel,ji_sel]
#t2m
#t2m   = f.variables['t2m'][:,:,:]#np.array([286.96487])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_ERA5.Timelist, beforeERA5, afterERA5,"ERA5_MED.","t2m")
       t2m[it] = (1-T_interpERA5)*Before_DATA[jj_sel,ji_sel] + T_interpERA5*After__DATA[jj_sel,ji_sel]
#d2m
#d2m   = f.variables['d2m'][:,:,:]#np.array([284.10672])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_ERA5.Timelist, beforeERA5, afterERA5,"ERA5_MED.","d2m")
       d2m[it] = (1-T_interpERA5)*Before_DATA[jj_sel,ji_sel] + T_interpERA5*After__DATA[jj_sel,ji_sel]
#tcc
#tcc   = f.variables['tcc'][:,:,:]#np.array([12.03943])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_ERA5.Timelist, beforeERA5, afterERA5,"ERA5_MED.","tcc")
       tcc[it] = 100.*((1-T_interpERA5)*Before_DATA[jj_sel,ji_sel] + T_interpERA5*After__DATA[jj_sel,ji_sel])
#tclw
#tclw  = f.variables['tclw'][:,:,:]#np.array([0.05])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_ERA5.Timelist, beforeERA5, afterERA5,"ERA5_MED.","tclw")
       tclw[it] = (1-T_interpERA5)*Before_DATA[jj_sel,ji_sel] + T_interpERA5*After__DATA[jj_sel,ji_sel]

########### MODIS CLOUD

       beforeMODIS_CLD,afterMODIS_CLD,T_interpMODIS_CLD = Time_Interpolation(tt,TL_MODIS_CLD.Timelist)
#cdrem
#cdrem = np.array([10.042508])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_MODIS_CLD.Timelist, beforeMODIS_CLD, afterMODIS_CLD,"MODIS_CLD_MED.","cdrem")
       cdrem[it] = (1-T_interpMODIS_CLD)*Before_DATA[jj_sel,ji_sel] + T_interpMODIS_CLD*After__DATA[jj_sel,ji_sel]

########### MODIS AEROSOL

       beforeMODIS_AER,afterMODIS_AER,T_interpMODIS_AER = Time_Interpolation(tt,TL_MODIS_AER.Timelist)
#taua
#      taua  = np.array([[0.12861817,0.11693237,0.11303711,0.10914185,0.10524658,0.10135132,
#          0.09745605,0.0934308,0.08906767,0.08519087,0.08188277,0.079188,
#          0.07697173,0.07506739,0.07332103,0.07159076,0.06982587,0.06804233,
#          0.0645105,0.05960191,0.05457829,0.05089346,0.0479044,0.04534442,
#          0.04311129,0.04106918,0.03914813,0.0373441,0.03567282,0.03332167,
#          0.02901247,0.01838189,0.00623265]])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_MODIS_AER.Timelist, beforeMODIS_AER, afterMODIS_AER,"MODIS_AEROSOL_MED.","taua")
       taua[0,it,:] = (1-T_interpMODIS_AER)*Before_DATA[:,jj_sel,ji_sel] + T_interpMODIS_AER*After__DATA[:,jj_sel,ji_sel]
#asymp
#      asymp = np.array([[0.78112185,0.76577055,0.76065344,0.7555364,0.75041926,0.74530214,
#          0.7401851,0.7352586,0.7308154,0.72576404,0.71973574,0.7124992,
#          0.7045832,0.69671303,0.68963724,0.68404067,0.6799671,0.6770385,
#          0.6733678,0.6695107,0.65929943,0.64360243,0.62599105,0.6091321,
#          0.59477234,0.58193725,0.5697906,0.5578721,0.54593724,0.5280586,
#          0.49310532,0.41066745,0.3164527]])
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_MODIS_AER.Timelist, beforeMODIS_AER, afterMODIS_AER,"MODIS_AEROSOL_MED.","asymp")
       asymp[0,it,:] = (1-T_interpMODIS_AER)*Before_DATA[:,jj_sel,ji_sel] + T_interpMODIS_AER*After__DATA[:,jj_sel,ji_sel]
#ssalb
#      ssalb = np.array([[0.9653854,0.9685926,0.96966165,0.9707307,0.9717998,0.97286886,
#          0.9739379,0.9750599,0.97628933,0.97734255,0.97821194,0.9788913,
#          0.9793933,0.97973263,0.9799186,0.9799566,0.9798687,0.9796886,
#          0.9791781,0.9784721,0.9785546,0.9787444,0.97678363,0.96918947,
#          0.95303714,0.9326887,0.91252154,0.897008,0.88944095,0.88677096,
#          0.89569014,0.89257544,0.88901585]])    
       Before_DATA, After__DATA = Load_Data(INPUT_DIR, TL_MODIS_AER.Timelist, beforeMODIS_AER, afterMODIS_AER,"MODIS_AEROSOL_MED.","ssalb")
       ssalb[0,it,:] = (1-T_interpMODIS_AER)*Before_DATA[:,jj_sel,ji_sel] + T_interpMODIS_AER*After__DATA[:,jj_sel,ji_sel]


       iyr = tt.year
       iday = datetime.date(tt.year,tt.month,tt.day).timetuple().tm_yday 
       sec   = tt.second+tt.minute*60+tt.hour*3600 
       sec_b = sec - dt/2
       sec_e = sec + dt/2
       year[it]=tt.year
       month[it]=tt.month
       day[it]=tt.day
       hour[it]=tt.hour
       minute[it]=tt.minute
       
       edout[it,:], esout[it,:] = cunit.monrad(points, iyr, iday, sec_b, sec_e,
                                   sp[it], msl[it], ws10[it], tco3[it], t2m[it], d2m[it], tcc[it], tclw[it], cdrem[it],
                                   taua[:,it,:], asymp[:,it,:], ssalb[:,it,:])

if rank == 0:
    for it,tt in enumerate(OUT_Timelist[:]):
        i = it % nranks
        if i != 0 :
            dic = comm.recv(source=i, tag=1)
            year[int(dic['idx'])] = dic['year']
            month[int(dic['idx'])] = dic['month']
            day[int(dic['idx'])] = dic['day']
            hour[int(dic['idx'])] = dic['hour']
            minute[int(dic['idx'])] = dic['min']
            ttlist[int(dic['idx'])] = dic['tt']
            edout[int(dic['idx']),:] = dic['edout']
            esout[int(dic['idx']),:] = dic['esout']
            sp[int(dic['idx'])] = dic['sp']
            msl[int(dic['idx'])] = dic['msl']
            ws10[int(dic['idx'])] = dic['ws10']
            tco3[int(dic['idx'])] = dic['tco3']
            t2m[int(dic['idx'])] = dic['t2m']
            d2m[int(dic['idx'])] = dic['d2m']
            tcc[int(dic['idx'])] = dic['tcc']
            tclw[int(dic['idx'])] = dic['tclw']
            cdrem[int(dic['idx'])] = dic['cdrem']
            taua[:,int(dic['idx']),:] = dic['taua']
            asymp[:,int(dic['idx']),:] = dic['asymp']
            ssalb[:,int(dic['idx']),:] = dic['ssalb']
else :
    for it,tt in enumerate(OUT_Timelist[:]):
        i = it % nranks
        dic = {'idx':it, 'year':year[it], 'month':month[it], 'day':day[it], 'hour':hour[it], 'min':minute[it], 'tt':ttlist[it],
                'edout':edout[it,:], 'esout':esout[it,:],
                'sp':sp[it], 'msl':msl[it], 'ws10':ws10[it], 'tco3':tco3[it],
                't2m':t2m[it],'d2m':d2m[it],'tcc':tcc[it], 'tclw':tclw[it], 'cdrem':cdrem[it],
                'taua':taua[:,it,:], 'asymp':asymp[:,it,:], 'ssalb':ssalb[:,it,:]}
        if i == rank :
            comm.send( dic,dest=0,tag=1 )
comm.Barrier()
#write data
if rank == 0 :
    d={'Year':year, 'Month':month, 'Day':day, 'Hour':hour, 'Min':minute, 'tt':ttlist}
    for iw,w in enumerate(wl):
        d["ed_"+str(w)]=edout[:,iw]
    for iw,w in enumerate(wl):
        d["es_"+str(w)]=esout[:,iw]
    d['sp'] =sp
    d['msl']=msl
    d['ws10']=ws10
    d['tco3']=tco3
    d['t2m']=t2m
    d['d2m']=d2m
    d['tcc']=tcc
    d['tclw']=tclw
    d['cdrem']=cdrem
    
    for iw,w in enumerate(wl):
        d["taua"+str(w)]   = taua[0,:,iw]
    for iw,w in enumerate(wl):
        d["asymp"+str(w)]  = asymp[0,:,iw]
    for iw,w in enumerate(wl):
        d["ssalb"+str(w)]  = ssalb[0,:,iw]
    
    df=pd.DataFrame(d)
    outfile=location + '2019-2021.txt'
    df.to_csv(outfile, index=False)
    

del cunit
del olib
