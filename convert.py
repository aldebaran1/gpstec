#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 15:24:36 2019

@author: smrak
"""

from glob import glob
import os, platform
import numpy as np
from gpstec import gpstec
from datetime import datetime
from dateutil import parser
import subprocess
from argparse import ArgumentParser

def convert(root:str = None, 
            date:str = None,
            tlim:str = None,
            force:bool = False):
    
    if os.path.isdir(root):
        files = np.array(sorted(glob(root + 'gps*.hdf5')))
        file_dates = np.array([datetime.strptime(os.path.split(f)[1][3:9], "%y%m%d") for f in files])

        if date is not None:
            datedt = parser.parse(date)
            idt = np.isin(file_dates, datedt)
        else:
            idt = np.ones(file_dates.size, dtype=bool)
    
        files_convert = files[idt]
        fd_convert = file_dates[idt]
        for i,f in enumerate(files_convert):
            if tlim is None:
                dateday = fd_convert[i].strftime("%j")
                dateyear = fd_convert[i].strftime("%Y")
                if (int(dateday) + 1) > 365:
                    year1 = int(dateyear) + 1
                    day1 = 1
                else:
                    year1 = int(dateyear)
                    day1 = int(dateday) + 1
                tl = [datetime.strptime("{} {}".format(str(dateyear), str(dateday)), "%Y %j"), 
                        datetime.strptime("{} {}".format(str(year1), str(day1)), "%Y %j")]
            else:
                tl = parser.parse(tlim)
            
            ofn = os.path.join(root, 'conv_' + tl[0].strftime('%Y%m%dT%H%M') + '-' + tl[1].strftime('%Y%m%dT%H%M') + '.h5')
    
            if not os.path.exists(os.path.split(ofn)[0]):    
                if platform.system() in ('Darwin', 'Linux'):
                    subprocess.call('mkdir -p "{}"'.format(os.path.split(ofn)[0]), shell=True, timeout=5)
                else:
                    subprocess.call('mkdir "{}"'.format(os.path.split(ofn)[0]), shell=True, timeout=5)
        
            if os.path.exists(ofn) and not force:
                print ('File already exist')
            else:
                print ('Loading and rearranging data ...')
                D = gpstec.returnGlobalTEC(datafolder=f, timelim=tl)
                print ('Saving data to ... {}'.format(ofn))
                gpstec.save2HDF(D['time'], D['xgrid'], D['ygrid'], D['tecim'], ofn)
    else:
        assert(os.path.isfile(root))
        if os.path.endswith('.hdf5'):
            print ('Loading and rearranging data ...')
            D = gpstec.returnGlobalTEC(datafolder=f, timelim=tl)
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
    
    convert(root = P.folder, date = P.date, tlim = P.tlim, force=P.force)