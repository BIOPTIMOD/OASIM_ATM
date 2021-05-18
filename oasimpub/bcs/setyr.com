#!/bin/sh
if [ $# -lt 1 ]
then
 echo 'Enter year to setup'
 exit 1
fi
YR=$1
if [ -f yr.dat ] ; then rm yr.dat ; fi
echo $YR >yr.dat
