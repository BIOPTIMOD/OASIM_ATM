#! /bin/bash
INDIR=/g100_work/OGS21_PRACE_P/NASA
OUTDIR=/g100_scratch/userexternal/gbolzon0/OASIM/modaer/NC

mkdir -p ${OUTDIR}

for YEAR in `seq 2017 2020`; do 
 for month in `seq 1 12 `; do 
    MONTH=`printf %02d $month`
    ./dat2NetCDF.xx ${INDIR}/modaer${YEAR}${MONTH}.dat ${OUTDIR}/modaer${YEAR}${MONTH}.nc ${YEAR}${MONTH}  ${YEAR}${MONTH}
 done
done
