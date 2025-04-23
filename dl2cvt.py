# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 12:15:36 2022

@author: smrak@bu.edu
"""

from datetime import timedelta, datetime
from dateutil import parser
from argparse import ArgumentParser
import subprocess
import os

def main(start: str = None, stop : str = None, odir: str = None,  los : bool = False, v : bool = False):
    
    timeout = 1900 if los else 180
    
    startdt, stopdt = parser.parse(start), parser.parse(stop)
    t = startdt
    dtlist = []
    
    while t <= stopdt:
        dtlist.append(t)
        t += timedelta(days=1)
        
    for date in dtlist:
        d = date.strftime('%Y-%m-%d')
        odir2 = odir + os.sep + str(date.year) + os.sep + date.strftime('%m%d') + os.sep
        line = f'./auto.sh {d} {odir} --los' if los else f'./auto.sh {d} {odir}'
        t0 = datetime.now()
        subprocess.call(line, shell = True, timeout=timeout)
        if v:
            print (f'It took {datetime.now()-t0} to download.')
        if los:
            line_l2g = f'./auto_l2g.sh {odir2}'
            t0 = datetime.now()
            subprocess.call(line_l2g, shell=True)
            if v:
                print ('It took {datetime.now()-t0} to re-grid')

if __name__ == '__main__':
    
    p = ArgumentParser()
    
    p.add_argument('startdate', type = str)
    p.add_argument('enddate', type = str)
    p.add_argument('odir', type=str)
    p.add_argument('--los', action='store_true')
    p.add_argument('-v', help='Verbose?', action='store_true')
    P = p.parse_args()
    
    main(start = P.startdate, stop = P.enddate, odir = P.odir, los = P.los, v = P.v)
