#! /bin/bash

INDIR=/g100_work/OGS21_PRACE_P/NASA
OUTDIR=/g100_scratch/userexternal/plazzari/OASIM/modcld/NC

mkdir -p ${OUTDIR}

for YEAR in `seq 2000 2021`; do
    ./dat2NetCDF.xx ${INDIR}/modcld${YEAR}.dat ${OUTDIR}/modcld${YEAR}.nc ${YEAR}01-00:00:00  ${YEAR}1231-23:59:59
done

#./dat2NetCDF.xx modcld2001.dat modcld2001.nc  20010101-00:00:00 20011231-23:59:59
#./dat2NetCDF.xx modcld2002.dat modcld2002.nc  20020101-00:00:00 20021231-23:59:59
#./dat2NetCDF.xx modcld2003.dat modcld2003.nc  20030101-00:00:00 20031231-23:59:59
#./dat2NetCDF.xx modcld2004.dat modcld2004.nc  20040101-00:00:00 20041231-23:59:59
#./dat2NetCDF.xx modcld2005.dat modcld2005.nc  20050101-00:00:00 20051231-23:59:59
#./dat2NetCDF.xx modcld2006.dat modcld2006.nc  20060101-00:00:00 20061231-23:59:59
