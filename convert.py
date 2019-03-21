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

months = {'jan': 1, 'feb': 2, 'mar':3, 'apr':4, 'may':5, 'jun':6, 'jul':7, 'aug':8, 'sep':9, 'oct':10, 'nov':11, 'dec':12}

def convert(root:str = None, 
            date:str = None,
            tlim:str = None,
            ofn:str = None):
    
    if date is not None:
        datedt = parser.parse(date)
    else:
        try:
            path = os.path.expanduser(root)
            parts = path.split(os.sep)
            assert len(parts[-2]) == 7
            part = parts[-2]
            month = str(months[part[2:5]])
            dlist = [part[:2], str(month), part[-2:]]
            datedt = datetime.strptime('-'.join(dlist), '%d-%m-%y')
        except Exception as e:
                raise (e)
    if tlim is None:
        tlim = [datetime(datedt.year, datedt.month, datedt.day, 0, 0, 0), 
                datetime(datedt.year, datedt.month, datedt.day  + 1, 0, 0, 0)]
    else:
        tlim = parser.parse(tlim)
    if os.path.splitext(root[1]) in ['.h5', '.hdf5']:
        fn = root
        folder = os.path.split(fn)[0]
    else:
        flist = sorted(glob(root + '*.hdf5'))
        if len(flist) > 0:
            fn = flist[0]
            folder = os.path.split(root)[0]
        else:
            try:
                d = str(datedt.day) if len(str(datedt.day)) == 2 else '0' + str(datedt.day)
                m = str(datedt.month) if len(str(datedt.month)) == 2 else '0' + str(datedt.month)
                ddir = d + tlim[0].strftime("%B")[:3].lower() + str(datedt.year)[2:]
                folder = os.path.join(root, ddir)
                pattern = 'gps*{}{}{}g.*.hdf5'.format(str(datedt.year)[2:], m, d)
                fn = sorted(glob(os.path.join(folder, pattern)))[0]
            except Exception as e:
                raise (e)
        
        

    if ofn is None:
        ofn = folder
    if os.path.isdir(ofn):
        ofn = os.path.join(ofn, 'conv_' + tlim[0].strftime('%m%dT%H%M') + '-' + tlim[1].strftime('%m%dT%H%M') + '.h5')

    if os.path.isfile(ofn):
        if not os.path.splitext(ofn)[1] in ['.h5', '.hdf5']:
            ofn = os.path.splitext(ofn)[0] + '.h5'
    
    if not os.path.exists(folder):
        subprocess.call('mkdir -p {}'.format(folder+'/'), shell=True, timeout=5)
    
    if os.path.exists(ofn):
        answer = str(input("{} already exists. Do you want to override it? Y/n".format(ofn)))
        if answer in ['y', 'Y']:
            print ('Loading and rearranging data ...')
            D = gpstec.returnGlobalTEC(datafolder=fn, timelim=tlim)
            print ('Saving data to ... {}'.format(ofn))
            gpstec.save2HDF(D['time'], D['xgrid'], D['ygrid'], D['tecim'], ofn)
        else:
            pass
    else:
        print ('Loading and rearranging data ...')
        D = gpstec.returnGlobalTEC(datafolder=fn, timelim=tlim)
        print ('Saving data to ... {}'.format(ofn))
        gpstec.save2HDF(D['time'], D['xgrid'], D['ygrid'], D['tecim'], ofn)

if __name__ == '__main__':
    
    p = ArgumentParser()
    
    p.add_argument('folder', type = str)
    p.add_argument('-d', '--date', type = str, help='date in YYYY-mm-dd format')
    p.add_argument('-o', '--ofn', help = 'Destination folder, if None-> the same as input folder', default = None)
    p.add_argument('--tlim', help = 'set time limints for the file to cenvert', nargs = 2)
    P = p.parse_args()
    
    convert(root = P.folder, date = P.date, tlim = P.tlim, ofn = P.ofn)