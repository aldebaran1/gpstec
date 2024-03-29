# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 12:15:36 2022

@author: smrak@bu.edu
"""

from datetime import timedelta, datetime
from dateutil import parser
from argparse import ArgumentParser
import subprocess

def main(start: str = None, stop : str = None, odir: str = None,  los : bool = False):
    
    timeout = 900 if los else 180
    
    startdt, stopdt = parser.parse(start), parser.parse(stop)
    t = startdt
    dtlist = []
    
    while t <= stopdt:
        dtlist.append(t.strftime('%Y-%m-%d'))
        t += timedelta(days=1)
        
    for d in dtlist:
        line = f'./auto.sh {d} {odir} --los' if los else f'./auto.sh {d} {odir}'
        t0 = datetime.now()
        print (line)
        subprocess.call(line, shell = True, timeout=timeout)
        print (f'It took {datetime.now()-t0} to download.')
        
if __name__ == '__main__':
    
    p = ArgumentParser()
    
    p.add_argument('startdate', type = str)
    p.add_argument('enddate', type = str)
    p.add_argument('odir', type=str)
    p.add_argument('--los', type = str)
    
    P = p.parse_args()
    
    main(start = P.startdate, stop = P.enddate, odir=P.odir, los = P.los)
