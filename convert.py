#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 15:24:36 2019

@author: smrak
"""

from glob import glob
import os
from gpstec import gpstec
from datetime import datetime
from dateutil import parser
import subprocess
from argparse import ArgumentParser

def convert(root:str = None, 
            date:str = None,
            tlim:str = None,
            ofn:str = None,
            force:bool = False):
    
    if date is not None:
        datedt = parser.parse(date)
        dateday = datedt.strftime('%j')
    else:
        fdate = root.split(os.sep)
        if fdate[-1] == '':
            fdate.pop(-1)
        dateyear = fdate[-2]
        datedt = datetime.strptime("{}{}".format(dateyear,fdate[-1]), '%Y%m%d')
        dateday = datedt.strftime('%j')
    if tlim is None:
        if (int(dateday) + 1) > 365:
            year1 = int(datedt.year) + 1
            day1 = 1
        else:
            year1 = int(datedt.year)
            day1 = int(dateday) + 1
        tlim = [datetime.strptime("{} {}".format(str(dateyear), fdate[-1]), "%Y %m%d"), 
                datetime.strptime("{} {}".format(str(year1), str(day1)), "%Y %j")]
    else:
        tlim = parser.parse(tlim)
        
        
    if os.path.isfile(root) and root.endswith('.hdf5'):
        fn = root
    else:
        fn = sorted(glob(root + os.sep +  '*.hdf5'))[0]
        
    if ofn is None:
        ofn = root
    if os.path.isdir(ofn):
        ofn = os.path.join(ofn, 'conv_' + tlim[0].strftime('%m%dT%H%M') + '-' + tlim[1].strftime('%m%dT%H%M') + '.h5')

    if os.path.isfile(ofn):
        if not os.path.splitext(ofn)[1] in ['.h5', '.hdf5']:
            ofn = os.path.splitext(ofn)[0] + '.h5'
    
    if not os.path.exists(os.path.split(ofn)[0]):
        
        subprocess.call('mkdir "{}"'.format(os.path.split(ofn)[0]), shell=True, timeout=5)
    
    if os.path.exists(ofn) and not force:
        print ('File already exist')
    else:
        print ('Loading and rearranging data ...')
        D = gpstec.returnGlobalTEC(datafolder=fn, timelim=tlim)
        print ('Saving data to ... {}'.format(ofn))
        gpstec.save2HDF(D['time'], D['xgrid'], D['ygrid'], D['tecim'], ofn)

if __name__ == '__main__':
    
    p = ArgumentParser()
    
    p.add_argument('folder', type = str)
    p.add_argument('-d', '--date', type = str, help='date in YYYY-mm-dd format', default=None)
    p.add_argument('-o', '--ofn', help = 'Destination folder, if None-> the same as input folder', default = None)
    p.add_argument('--tlim', help = 'set time limints for the file to cenvert', nargs = 2)
    p.add_argument('--force', help='Do you want to override existing files?', action='store_true')
    P = p.parse_args()
    
    convert(root = P.folder, date = P.date, tlim = P.tlim, ofn = P.ofn, force=P.force)