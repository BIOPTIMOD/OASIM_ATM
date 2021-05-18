#!/bin/bash  
IRRDIR=/gpfs/scratch/userexternal/plazzari/OASIM_HF/OUTPUT/
yyyymmdd=`cat bcs/yr.dat`
YR=$( echo $yyyymmdd | cut -c 1-4)
MON=$( echo $yyyymmdd | cut -c 5-6)
DAY=$( echo $yyyymmdd | cut -c 7-8)
JDAY=$(echo $(date --date=$yyyymmdd +%j))
echo ' yyyymmdd ' $yyyymmdd
echo ' YR ' $YR
echo ' MON ' $MON
echo ' DAY ' $DAY
echo 'JDAY ' $JDAY
IRRDIRYR=$IRRDIR/$yyyymmdd
if [ ! -d $IRRDIRYR ] ; then mkdir $IRRDIRYR ; fi
cd src
rm *.o
#make monrad
./compile.sh
cp monrad ../
cd ../
#  
# Set up input data files
for AFLE in  atmo25b abw25b acbc25b clddays slingo
do
 if [ -f bcs/$AFLE.dat.gz ]
 then
  echo 'Uncompressing ' $AFLE.dat.gz
  gunzip bcs/$AFLE.dat.gz
 fi
done

ln -fs /gpfs/scratch/userexternal/plazzari/OASIM_HF/CREATE_INPUT_ECMWF/output/opt${yyyymmdd}_ECMWF.nc opt.nc
ln -fs /gpfs/scratch/userexternal/plazzari/OASIM_HF/CREATE_INPUT_ECMWF/output/clouds${yyyymmdd}_ECMWF.nc clouds.nc

#
#  Run
echo $JDAY > day.dat
echo $MON > month.dat
echo $YR  > year.dat
if [ ! -f $IRRDIRYR/swr$yyyymmdd.dat.gz ]
 then
   cp -p atmdata/modaer$YR$MON.dat.gz modaer.dat.gz
   gunzip modaer.dat.gz
  fi
  ./monrad > out
  mv out $IRRDIRYR/out$yyyymmdd
  if [ -f modaer.dat ] ; then rm modaer.dat ; fi

  mv rad rad$yyyymmdd.dat
  gzip rad$yyyymmdd.dat
  mv rad$yyyymmdd.dat.gz $IRRDIRYR
  mv rad_0p.nc rad_0p$yyyymmdd.nc
  mv rad_0m.nc rad_0m$yyyymmdd.nc
  gzip rad_0p$yyyymmdd.nc
#  gzip rad_0m$yyyymmdd.nc
  mv rad_0p$yyyymmdd.nc.gz $IRRDIRYR
# mv rad_0m$yyyymmdd.nc.gz $IRRDIRYR
  mv rad_0m$yyyymmdd.nc    $IRRDIRYR

  mv avgirr.dat swr$yyyymmdd.dat
  gzip swr$yyyymmdd.dat
  mv swr$yyyymmdd.dat.gz $IRRDIRYR
  mv swr.nc swr$yyyymmdd.nc
  gzip swr$yyyymmdd.nc
  mv swr$yyyymmdd.nc.gz $IRRDIRYR

  mv avgparq.dat parq$yyyymmdd.dat
  gzip parq$yyyymmdd.dat
  mv parq$yyyymmdd.dat.gz $IRRDIRYR
  mv parq.nc parq$yyyymmdd.nc
  gzip parq$yyyymmdd.nc
  mv parq$yyyymmdd.nc.gz $IRRDIRYR

  mv eds.dat eds$yyyymmdd.dat
  gzip eds$yyyymmdd.dat
  mv eds$yyyymmdd.dat.gz $IRRDIRYR
  mv eds.nc eds$yyyymmdd.nc
  gzip eds$yyyymmdd.nc
  mv eds$yyyymmdd.nc.gz $IRRDIRYR

fi
#
