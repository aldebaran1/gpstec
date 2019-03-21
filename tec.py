# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:32:25 2019

@author: smrak
"""
from glob import glob
from gpstec import gpstec
from datetime import datetime
import numpy as np
import os
import yaml
import subprocess
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartomap import geogmap as gm

mcfg = '/home/smrak/Documents/cartomap/map/conus.yaml'

def _round(x, base=5):
    r = base * round(x/base)
    if r < x:
        return r + 5
    else:
        return r
    return 

#mode = ['lin', 'log']
year = 2017
month = 5
day = 28
tlim = [datetime(year,month,day,0,0,0), datetime(year,month,day+1,0,0,0)]
d = str(day) if len(str(day)) == 2 else '0' + str(day)
m = str(month) if len(str(month)) == 2 else '0' + str(month)
ddir = d + tlim[0].strftime("%B")[:3].lower() + str(year)[2:]
#folder = 'C:\\Users\\smrak\\Documents\\LWSI\\' + ddir + '\\'
folder = '/media/smrak/gnss/gpstec/2017/gps/' + ddir + '/'
fn = glob(folder + 'gps*{}{}{}g.*.hdf5'.format(str(year)[2:], m, d))[0]
# Save to:
save = False
odir = '/media/smrak/figures/gpstec/' + str(year) + '/' + ddir
# HDF FN
fnhdf = folder + '270-280ut.h5'
D = gpstec.readFromHDF(fnhdf)
# Load Data from original database
#D = gpstec.returnGlobalTEC(datafolder=fn, timelim=tlim)
# Save to reshuffled format ??
savehdf = 0
# Map details
#vmax = 3 * _round(np.nanmedian(D['tecim']), base=5)
#vmin = 0
#clim = [vmin, vmax]
#colorbar = 1
#
#projection = 'stereo'
#cmap = 'jet'
#nightshade = False
#latlim = [15, 55]
#lonlim = [-130, -60]
#
#parallels = np.arange(0,70,10)
#meridians = np.arange(-160,-10,20)
mst = yaml.load(open(mcfg, 'r'))
projection = mst.get('projection')
cmap = mst.get('cmap')
nightshade = mst.get('nightshade')
latlim = mst.get('latlim')
lonlim = mst.get('lonlim')
paralim = mst.get('parallels')
meridilim = mst.get('meridians')
parallels = np.arange(paralim[0], paralim[1], paralim[2])
meridians = np.arange(meridilim[0], meridilim[1], meridilim[2])

mode = mst.get('mode') if mst.get('mode') is not None else 'lin'
colorbar = mst.get('colorbar')
cbar_label = mst.get('cbar_label') if mst.get('cbar_label') is not None else 'TEC [TECu]'
clim = mst.get('clim')
if clim is None:
    vmax = 3 * _round(np.nanmedian(D['tecim']), base=5)
    if mode == 'log':
        vmax = np.round(np.log10(vmax), 1)
    vmin = 0
    clim = [vmin, vmax]

#vmax = 3 * _round(np.nanmedian(D['tecim']), base=5)
#vmin = 0
#clim = [vmin, vmax]
#colorbar = 1

c = 1
for mode in mode:
    savedir = os.path.join(odir, mode)
    if mode == 'log': clim[1] = np.round(np.log10(vmax), 1)
    for i in range(D['time'].shape[0]):
        
        
        
        z = np.log10(D['tecim'][i]) if mode == 'log' else D['tecim'][i]
        fig, ax = gm.plotCartoMap(figsize=(12,7), projection=projection,
                              title=D['time'][i],
                              latlim=latlim,lonlim=lonlim,
                              meridians=meridians, parallels=parallels,
                              background_color='grey', grid_color='white',
                              grid_linewidth=1,
                              figure=True, 
                              nightshade=nightshade, ns_dt=D['time'][i])
        im = plt.pcolormesh(D['xgrid'],D['ygrid'], z.T, cmap=cmap, transform=ccrs.PlateCarree())
        plt.clim(clim)
        
        if colorbar:
            cax = fig.add_axes([0, 0, 0.1, 0.1])
            cbar = plt.colorbar(im, cax=cax, label='TEC [TECu]')
            
            axp = fig.gca()
            posn = ax.get_position()
            axp.set_position([posn.x0 + posn.width + 0.04, posn.y0,
                          0.02, posn.height])
            fig.canvas.draw()
        if save:
            figdir =os.path.join(odir,mode)
            if not os.path.exists(figdir):
                print ('Creating a new directory: ', savedir)
                subprocess.call('mkdir -p ', figdir, shell=True, timeout=5)
            fn = D['time'][i].strftime('%m%d_%H%M')
            figname = fn + str('.png')
            plt.savefig(os.path.join(savedir, figname), dpi=100)
            plt.close(fig=fig)
        else:
            plt.show()
        break
if savehdf and c == 0:
    fns = 'cnv_' + tlim[0].strftime('%m%dT%H%M') + '-' + tlim[1].strftime('%m%dT%H%M')
    print ('Saving data to ... {}'.format(fns))
    gpstec.save2HDF(D['time'], D['xgrid'], D['ygrid'], D['tecim'], folder + fns + 'ut.h5')
c += 1
#    break