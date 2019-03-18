# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:32:25 2019

@author: smrak
"""
from glob import glob
from gpstec import gpstec
from datetime import datetime

import numpy as np

mode = ['lin', 'log']
year = 2017
month = 9
day = 11
tlim = [datetime(year,month,day,0,0,0), datetime(year,month,day+1,0,0,0)]
d = str(day) if len(str(day)) == 2 else '0' + str(day)
m = str(month) if len(str(month)) == 2 else '0' + str(month)
ddir = d + tlim[0].strftime("%B")[:3].lower() + str(year)[2:]
folder = 'C:\\Users\\smrak\\Documents\\LWSI\\' + ddir + '\\'
fn = glob(folder + 'gps*{}{}{}g.*.hdf5'.format(str(year)[2:], m, d))[0]

D = gpstec.returnGlobaTEC(datafolder=fn, timelim=tlim)

c = 0

for mode in mode:
    clim = [0,40] if mode == 'lin' else [0,1.5]
    for i in range(D['time'].shape[0]):
        figname = folder + '\\' + mode +'\\' + D['time'][i].strftime('%m%d_%H%M')
        z = np.log10(D['tecim'][i]) if mode == 'log' else D['tecim'][i]
        fig = gpstec.plotTECmap(D['xgrid'], D['ygrid'], z, 
                                latlim=[10,70],title = D['time'][i],
                                clim = clim, projection='lambert',cmap='jet',
                                savefn=figname)
    #    break
    if c == 0:
        fns = str(tlim[0].day) + str(tlim[0].hour) + '-' + str(tlim[1].day) + str(tlim[1].hour)
        gpstec.save2HDF(D['time'], D['xgrid'], D['ygrid'], D['tecim'], folder + fns + 'ut.h5')
    c += 1

