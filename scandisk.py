#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 19:10:08 2019

@author: smrak
"""
import os
from glob import glob
import subprocess
from datetime import datetime, timedelta

months = {'jan': 1, 'feb': 2, 'mar':3, 'apr':4, 'may':5, 'jun':6, 'jul':7, 'aug':8, 'sep':9, 'oct':10, 'nov':11, 'dec':12}
FN = '/media/smrak/figures/gpstec/'

if __name__ == "__main__":

    years = glob(FN+'*')
    for year  in years:
        diryr = os.path.join(year,'gps/*')
        days = glob(diryr+'*')
        for day in days:
            print (os.path.split(day)[1])
            fnlist = glob(day+'/*')
            names = [os.path.split(f)[1][:3] for f in fnlist]
            odir = day + '/'
            date = os.path.split(day)[1]
            mnth = months[date[2:-2]]
            dt = datetime.strptime("{}-{}-{}".format(date[:2], mnth, date[-2:]), "%d-%m-%y")
            m = dt.month if len(str(dt.month)) == 2 else '0' + str(dt.month)
            d = dt.day if len(str(dt.day)) == 2 else '0' + str(dt.day)
            FMT = "{}-{}-{}".format(dt.year, m, d)
            if 'gps' not in names:
                print ('Downloading raw data...')
                subprocess.call("python dltec.py {} {} {} --fixpath".format(FMT, FMT, odir), 
                            shell=True)
            t0 = dt.strftime("%m%dT0000")
            t1 = (dt + timedelta(days=1)).strftime("%m%dT0000")
            CVTFN = 'conv_' + t0 + '-' + t1 + '.h5'
            CFN = os.path.join(odir, CVTFN)
            if not os.path.exists(CFN):
                print ('Converting raw data ...')
                subprocess.call("python convert.py {}".format(odir), 
                                shell=True)
            if 'lin' not in names:
                print ('Plotting the data (lin)...')
                subprocess.call("python plottec.py {} --mode lin".format(CFN), shell=True)
            if 'log' not in names:
                print ('Plotting the data (log)...')
                subprocess.call("python plottec.py {} --mode log".format(CFN), shell=True)
