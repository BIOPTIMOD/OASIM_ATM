#! /bin/bash

module load autoload nco
indir=/g100_scratch/userexternal/plazzari/OASIM/modcld/NC
outdir=/g100_scratch/userexternal/plazzari/PROVA
mkdir -p $outdir
year=2019
echo $year

 for year in `seq 2019 2021`; do
    startDate="${year}-01-01 00:00:00"
    for month in `seq 0 11`; do

      filein=$indir/modcld${year}.nc
    
      yyyymmdd_hhmmss=$(date -d "$startDate ${month} month" +%Y%m15-12:%M:%S)
    
      fileout=${outdir}/MODIS${yyyymmdd_hhmmss}.nc
    
      echo $fileout
    
      ncks -O  -d time,$month,$month $filein $fileout
    
    done
done
