#! /bin/bash



MASKFILE=/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/V8C/meshmask.nc


INPUTFILE=/g100_scratch/userexternal/gbolzon0/OASIM/20190101-ECMWF---AM0100-MEDATL-b20190102_an-fv12.00.nc
OUTDIR=/g100_scratch/userexternal/gbolzon0/OASIM_ATM/out




python daily_oasim_input_gen.py -i $INPUTFILE -o $OUTDIR  -m $MASKFILE
#python3  create_cloud_nc.py  ${yyyymmdd}


