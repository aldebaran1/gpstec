#!/usr/bin/env bash

FLAG=$1
echo 'Converting los 2 grid'
python los2grid.py $FLAG
echo 'Done'
