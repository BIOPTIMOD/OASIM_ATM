#! /bin/bash


module purge
conda activate python_conda

while read yyyymmdd ;
do

echo "$yyyymmdd"

python3  create_opt_nc.py    ${yyyymmdd}
python3  create_cloud_nc.py  ${yyyymmdd}

done < lista_date.txt
