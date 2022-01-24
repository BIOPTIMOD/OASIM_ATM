#!/bin/bash

WRKDIR=$PWD

while read yyyymmdd; do
  echo "$yyyymmdd"

 cd bcs

 ./setyr.com $yyyymmdd

  cd $WRKDIR

 ./monradyr.com >monradyr.doc

done <lista_date.txt
