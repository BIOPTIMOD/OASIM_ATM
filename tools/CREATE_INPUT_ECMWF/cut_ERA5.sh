#! /bin/bash

module load autoload nco

filein=/g100_scratch/userexternal/gocchipi/ERA5_Med.nc
outdir=/g100_scratch/userexternal/plazzari/PROVA
year=2019
echo $year
startDate="2019-01-01 00:00:00"

for hour in `seq 1 8760`; do

  yyyymmdd_hhmmss=$(date -d "$startDate ${hour} hours" +%Y%m%d-%H:%M:%S)

  fileout=${outdir}/ERA5_MED${yyyymmdd_hhmmss}.nc

  echo $fileout

  ncks -O -F -d time,$hour,$hour $filein $fileout

done
