#! /bin/bash

#SBATCH --job-name=OPTICS-pre
#SBATCH -N1
#SBATCH --ntasks-per-node=16
#SBATCH --time=01:30:00
#SBATCH --mem=300gb
#SBATCH --account=OGS21_PRACE_P
#SBATCH --partition=g100_usr_prod
#SBATCH --qos=g100_qos_dbg

cd $SLURM_SUBMIT_DIR

. ./profile.inc

module load autoload
module load intel/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary
module load mkl/oneapi-2021--binary
module load cmake/3.18.4--gcc--10.2.0
module load netcdf/4.7.4--oneapi--2021.2.0-ifort
module load netcdff/4.5.3--oneapi--2021.2.0-ifort
source /g100_work/OGS20_PRACE_P_2/COPERNICUS/py_env_3.6.8/bin/activate
export PYTHONPATH=$PYTHONPATH:/g100_work/OGS20_PRACE_P_2/COPERNICUS/bit.sea3

unset I_MPI_PMI_LIBRARY
export UCX_TLS=ib
export SLURM_PMIX_DIRECT_CONN_UCX=false

MASKFILE=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc

INPUTDIR=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/ORIG_ATM_CMCC

OUTDIR=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/READY_FOR_MODEL/

my_prex_or_die "mpirun python daily_oasim_input_gen.py -i $INPUTDIR -o $OUTDIR  -m $MASKFILE"

ORIG_ERA5=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/ERA5/ # downloaded by plazzari
CLIM_ERA5=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/CLIM_ERA5/

my_prex_or_die "python era5_climatology.py -i $ORIG_ERA5 -o $CLIM_ERA5"

## MODCLD SECTION ###
# DIR with modcldyyyymm.dat --> conversion to NETCDF --> ncea to have modcldmm.nc
CLIM_MODCLD=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/MODCLD/clim/monthly

## Putting together  ##
my_prex_or_die "python monthly_clim_gen.py -e $CLIM_ERA5 -M $CLIM_MODCLD -o $OUTDIR -m $MASKFILE"

## AEROSOL SECTION ###
# DIR with modaer.dat --> conversion to NetCDF --> ncea to have modaermm.nc
NASA_CLIM=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/OPTICS/modaer_NASA/clim/monthly

my_prex_or_die "python aerosol.py -i $NASA_CLIM -m $MASKFILE -o $OUTDIR"




