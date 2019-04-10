#!/usr/bin/env bash

D=$1
FLAG=""
if [ $2 ]; then
  FLAG=$FLAG" --mode "$2
fi
if [ $3 ] && [ $4 ]; then
  FLAG=$FLAG" --clim "$3" "$4
fi


year=$(date -d $D +'%Y')
dBy=$(date -d $D +'%d%B%y')
dBy_l="${dBy,,}"

d1=$(date -d $D +'%m%dT%H%M')
d2=$(date -d $D +'%m%dT%H%M' -d "$D + 1 day")

ODIR=/media/smrak/figures/gpstec/
ODIR2=$ODIR$year'/gps/'$dBy_l'/'
filename='conv_'$d1'-'$d2'.h5'

python dltec.py $D $D $ODIR
python convert.py $ODIR2
python plottec.py $ODIR2$filename $FLAG
