#! /bin/bash


module purge
conda activate python_conda

while read yyyymmdd ;
do

echo "$yyyymmdd"

INPUTFILE=/g100_scratch/userexternal/plazzari/OASIM/ECMWF/ERAINTERIM_2017.nc
OUTFILE=/g100_scratch/userexternal/OASIM/output_${yyyymmdd}.nc




python3 create_opt_nc.py -i $INPUTFILE -o $OUTFILE -d ${yyyymmdd} -m mask.nc
python3  create_cloud_nc.py  ${yyyymmdd}

done < lista_date.txt
