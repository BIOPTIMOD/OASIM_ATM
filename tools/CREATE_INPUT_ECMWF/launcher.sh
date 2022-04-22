#! /bin/bash



MASKFILE=/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/V8C/meshmask.nc

INPUTDIR=/g100_scratch/userexternal/gbolzon0/V9C/DEV_OASIM_INPUTS/ORIG_ATM_CMCC/

OUTDIR=/g100_scratch/userexternal/gbolzon0/V9C/DEV_OASIM_INPUTS/daily_atm_6h

for INPUTFILE in $INPUTDIR/2019012* ; do

#INPUTFILE=/g100_scratch/userexternal/gbolzon0/V9C/DEV_OASIM_INPUTS/ORIG_ATM_CMCC/20190101-ECMWF---AM0100-MEDATL-b20190102_an-fv12.00.nc





python daily_oasim_input_gen.py -i $INPUTFILE -o $OUTDIR  -m $MASKFILE
done

#python3  create_cloud_nc.py  ${yyyymmdd}


