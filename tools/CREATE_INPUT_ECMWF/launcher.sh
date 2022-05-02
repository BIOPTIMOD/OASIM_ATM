#! /bin/bash

MASKFILE=/gss/gss_work/DRES_OGS_BiGe/gbolzon/masks/eas/V8C/meshmask.nc

INPUTDIR=/g100_scratch/userexternal/gbolzon0/V9C/DEV_OASIM_INPUTS/ORIG_ATM_CMCC/

OUTDIR=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/READY_FOR_MODEL/

for INPUTFILE in $INPUTDIR/2019012* ; do
    python daily_oasim_input_gen.py -i $INPUTFILE -o $OUTDIR  -m $MASKFILE
done

ORIG_ERA5=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/ERA5/ # downloaded by plazzari
CLIM_ERA5=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/CLIM_ERA5/

python era5_climatology.py -i $ORIG_ERA5 -o $CLIM_ERA5

## MODCLD SECTION ###
# DIR with modcldyyyymm.dat --> conversion to NETCDF --> ncea to have modcldmm.nc
CLIM_MODCLD=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/MODCLD/clim/monthly

## Putting together  ##
python monthly_clim_gen.py -e $CLIM_ERA5 -M $CLIM_MODCLD -o $OUTDIR -m $MASKFILE

## AEROSOL SECTION ###
# DIR with modaer.dat --> conversion to NetCDF --> ncea to have modaermm.nc
NASA_CLIM=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/modaer_NASA/clim/monthly

python aerosol.py -i $NASA_CLIM -m $MASKFILE -o $OUTDIR




