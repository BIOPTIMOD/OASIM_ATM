#! /bin/bash



MASKFILE=/g100_scratch/userexternal/gbolzon0/DEBUG/OGSTM/wrkdir/MODEL/meshmask.nc


INPUTFILE=/g100_scratch/userexternal/gbolzon0/OASIM/20190101-ECMWF---AM0100-MEDATL-b20190102_an-fv12.00.nc
OUTDIR=/g100_scratch/userexternal/gbolzon0/OASIM_ATM/out




echo python create_opt_nc.py -i $INPUTFILE -o $OUTDIR  -m $MASKFILE
#python3  create_cloud_nc.py  ${yyyymmdd}


