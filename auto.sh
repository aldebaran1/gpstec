#!/usr/bin/env bash

D=$1
ODIR=$2
FLAG=""
if [ $3 ]; then
  FLAG=$FLAG" "$3
fi

year=$(date -d $D +'%Y')
dBy=$(date -d $D +'%m%d')
dBy_l="${dBy,,}"

d1=$(date -d $D +'%m%dT%H%M')
d2=$(date -d $D +'%m%dT%H%M' -d "$D + 1 day")

ODIR2=$ODIR$year'$dBy_l'/

python dltec.py $D $D $ODIR $FLAG
sleep 1
python convert.py $ODIR2
