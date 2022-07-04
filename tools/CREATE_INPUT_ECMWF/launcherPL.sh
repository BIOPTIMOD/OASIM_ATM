#! /bin/bash

#SBATCH --job-name=OPTICS-pre
#SBATCH -N1
#SBATCH --ntasks-per-node=16
#SBATCH --time=01:30:00
#SBATCH --mem=300gb
#SBATCH --account=OGS21_PRACE_P
#SBATCH --partition=g100_usr_prod
#SBATCH --qos=g100_qos_dbg

#cd $SLURM_SUBMIT_DIR

. ./profile.inc

module load autoload
module load intel/oneapi-2021--binary
module load intelmpi/oneapi-2021--binary
module load mkl/oneapi-2021--binary
#module load netcdf/4.7.4--oneapi--2021.2.0-ifort
#module load netcdff/4.5.3--oneapi--2021.2.0-ifort
source /g100_work/OGS20_PRACE_P_2/COPERNICUS/py_env_3.6.8/bin/activate
export PYTHONPATH=$PYTHONPATH:/g100_work/OGS20_PRACE_P_2/COPERNICUS/bit.sea3

unset I_MPI_PMI_LIBRARY
export UCX_TLS=ib
export SLURM_PMIX_DIRECT_CONN_UCX=false

MASKFILE=/g100_work/OGS_prod100/OPA/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc


ORIG_ERA5=/g100_scratch/userexternal/gocchipi/ERA5_Med # downloaded by plazzari
OUT_ERA5=/g100_scratch/userexternal/plazzari/TEST_ERA5/READY_FOR_MODEL

#my_prex_or_die "python ERA5_indata_gen.py -m $MASKFILE  -e $ORIG_ERA5 -o $OUT_ERA5"

## MODCLD SECTION ###
# DIR with modcldyyyymm.dat --> conversion to NETCDF --> ncea to have modcldmm.nc
ORIG_MODIS=/g100_scratch/userexternal/plazzari/PROVA # downloaded by plazzari
OUT_MODIS=/g100_scratch/userexternal/plazzari/TEST_ERA5/READY_FOR_MODEL

#my_prex_or_die "python MODIS_indata_gen.py -m $MASKFILE  -M $ORIG_MODIS -o $OUT_MODIS"


## AEROSOL SECTION ###
# DIR with modaer.dat --> conversion to NetCDF --> ncea to have modaermm.nc
ORIG_NASA=/g100_scratch/userexternal/plazzari/OASIM/modaer/NC
OUT_NASA=/g100_scratch/userexternal/plazzari/TEST_ERA5/READY_FOR_MODEL

my_prex_or_die "python AEROSOL_indata_gen.py -A $ORIG_NASA -m $MASKFILE -o $OUT_NASA"




