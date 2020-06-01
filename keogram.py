# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:35:51 2020

@author: smrak@bu.edu
"""

from pyGnss import dm
from gpstec import gpstec
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib import dates
from dateutil import parser
import h5py, os
import numpy as np

d = '2015-4-17'

folder = 'G:\\My Drive\\scintillation_data\\tid\\{}\\'.format(parser.parse(d).strftime("%Y%m%d"))

fn = folder + 'conv_{}T0000-{}T0000.h5'.format(parser.parse(d).strftime("%Y%m%d"), 
           (parser.parse(d)+timedelta(days=1)).strftime("%Y%m%d") )

fntid = folder + '{}_{}T0000-{}T0000_all{}.yaml_10el_60s_altkm_350_res_03.h5'.format(
            parser.parse(d).year, parser.parse(d).strftime("%m%d"), 
            (parser.parse(d)+timedelta(days=1)).strftime("%m%d"), 
            parser.parse(d).strftime("%m%d"))
if not os.path.exists(fntid):
    fntid = folder + '{}_{}T0000-{}T0000_all{}.yaml_20el_300s_altkm_350_res_03.h5'.format(
            parser.parse(d).year, parser.parse(d).strftime("%m%d"), 
            (parser.parse(d)+timedelta(days=1)).strftime("%m%d"), 
            parser.parse(d).strftime("%m%d"))

fnscint = folder + 'ix_{}_{}T0000-{}T0000_reduced{}_d1_r2_yaml_30el_1s_350km.h5'.format(
            parser.parse(d).year, parser.parse(d).strftime("%m%d"), 
            (parser.parse(d)+timedelta(days=1)).strftime("%m%d"), 
            parser.parse(d).strftime("%m%d"))

tlim = [parser.parse(d), parser.parse(d)+timedelta(hours=10)]
Xt = -80
Yt=[10, 60]

dtec_clim = [-0.02, 0.02]

TEC = gpstec.readFromHDF(fn)
xgrid = TEC['xgrid']
ygrid = TEC['ygrid']
im = TEC['tecim']

D = dm.keogram(time=TEC['time'], xgrid=xgrid, ygrid=ygrid, im=im,
                      tlim=tlim, Xt=Xt, Yt=Yt, Xn=2, fillPixel=2)

fig = plt.figure(figsize=[10, 5])
plt.title('Lon0 = {}'.format(Xt))
plt.pcolormesh(D['time'], D['Y'], D['keo'].T, cmap='nipy_spectral', vmin=5, vmax=30)
plt.colorbar()
formatter = '%H:%M'
ax=plt.gca()
ax.xaxis.set_major_formatter(dates.DateFormatter(formatter))

TID = h5py.File(fntid, 'r')
time = np.array([datetime.utcfromtimestamp(t) for t in TID['data/time'][:]])
x = TID['data/xgrid'][:]
y = TID['data/ygrid'][:]

try:
    TD = dm.keogram(time=time, xgrid=x, ygrid=y, im=TID['data/im'][:][:][:],
                    tlim=tlim, Xt=Xt, Yt=Yt, Xn=5, imskip=1, 
                    fillPixel=2, filter='median')
    xg, yg = np.meshgrid(np.arange(TD['time'].size), 
                         np.arange(TD['Y'].size))
    T00 = gpstec.interpolateTEC(im=D['keo'], 
                         xl = TD['time'].size, yl=TD['Y'].size,
#                         xgrid=xg, ygrid=yg,
                         method='linear')
    
    fig = plt.figure(figsize=[10,5])
    plt.title('Lon0 = {}'.format(Xt))
    plt.pcolormesh(TD['time'], TD['Y'], (TD['keo']/T00).T, cmap='bwr', 
                vmin=dtec_clim[0], vmax=dtec_clim[1])
#    plt.contour(TD['time'], TD['Y'], TD['keo'].T, cmap='bwr', 
#                levels=np.arange(-2, 2, 0.5))
    plt.colorbar()
    formatter = '%H:%M'
    ax=plt.gca()
    ax.xaxis.set_major_formatter(dates.DateFormatter(formatter))
except Exception as e:
    TID.close()
    print (e)
TID.close()

# TEC + dTEC
k_m = (TD['keo']/T00)
k_m[abs(k_m) < 0.001] = np.nan

fig = plt.figure(figsize=[10,5])
plt.title('Lon0 = {}'.format(Xt))
tidim = plt.pcolormesh(TD['time'], TD['Y'], k_m.T, cmap='bwr', 
                vmin=dtec_clim[0], vmax=dtec_clim[1])
tecim = plt.contour(D['time'], D['Y'], D['keo'].T, levels=np.arange(0,51,5), 
            colors='k', extend='both')
formatter = '%H:%M'
ax=plt.gca()
ax.xaxis.set_major_formatter(dates.DateFormatter(formatter))

clabels = ax.clabel(tecim, tecim.levels, inline=True, fmt='%i', fontsize=10)
#[txt.set_backgroundcolor('white') for txt in clabels]

posn = ax.get_position()
cax = fig.add_axes([posn.x0+posn.width+0.01, posn.y0, 0.02, posn.height])
fig.colorbar(tidim, cax=cax, label='$\delta TEC / TEC$ [%]')
#cax = fig.add_axes([posn.x0+posn.width+0.12, posn.y0, 0.02, posn.height])
#fig.colorbar(tecim, cax=cax, label='TEC [TECu]')

# SCINT
SCINT = h5py.File(fnscint, 'r')
scint_time = np.array([datetime.utcfromtimestamp(t) for t in SCINT['data/time'][:]])
tgrid = D['time'][1:]
ygrid = D['Y']
st_keo = np.copy((D['keo'].shape[0] - 1, D['keo'].shape[1]))
snr_keo = np.copy((D['keo'].shape[0] -1, D['keo'].shape[1]))
for it, t in enumerate(tgrid):
    idt = abs(scint_time - t).argmin()
    glon = SCINT['data/ipp'][idt-300 : idt+1, :, :, 1]
    idx = (glon >= Xt-2) & (glon <= Xt+2)
    glat = SCINT['data/ipp'][idt-300 : idt+1, :, :, 0][idx]
    
    snr_row = np.nan * np.zeros(ygrid.size)
    st_row = np.nan * np.zeros(ygrid.size)
    for j in range(ygrid.size):
        idy = abs(glat - ygrid[j]) <= 1
        if np.sum(idy) > 0:
#            print ('sss')
            snr_row[j] = np.nanmedian(SCINT['data/snr4'][idt-300 : idt+1, :, :][idx][idy])
            st_row[j] = np.nanmedian(SCINT['data/sigma_tec'][idt-300 : idt+1, :, :][idx][idy])
    if np.sum(np.isfinite(st_row)) > 0:
        st_keo[it] = st_row
    if np.sum(np.isfinite(snr_row)) > 0:
        snr_keo[it] = snr_row
    
SCINT.close()

plt.pcolormesh(tgrid, ygrid, st_keo.T, cmap='jet')

