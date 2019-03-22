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
from argparse import ArgumentParser

mcfg = '/home/smrak/Documents/cartomap/map/conus.yaml'

def _round(x, base=5):
    r = base * round(x/base)
    if r < x:
        return r + 5
    else:
        return r
    return 

def plot(fnhdf:str = None,
         show:bool = False,
         odir:str = None,
         cfg:str = None):

    assert fnhdf is not None and os.path.splitext(fnhdf)[1] in ['.h5', '.hdf5']
    # Converted file
    D = gpstec.readFromHDF(fnhdf)
    if cfg is None:
        folder = os.path.expanduser(os.path.join(os.getcwd(), 'map'))
        mcfg = sorted(glob(os.path.join(folder, '*.yaml')))[0]
    else:
        assert os.path.splitext(cfg)[1] in ['.yaml', 'yml']
        mcfg = cfg
    # Map
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
    if isinstance(mode, str):
        mode = [mode]
    colorbar = mst.get('colorbar')
    cbar_label = mst.get('cbar_label') if mst.get('cbar_label') is not None else 'TEC [TECu]'
    clim = mst.get('clim')
    if clim is None:
        vmax = 3 * _round(np.nanmedian(D['tecim']), base=5)
        vmin = 0
        clim = [vmin, vmax]
    
    figsize = mst.get('figsize') if mst.get('figsize') is not None else [12,8]
    background_color = mst.get('background_color') if mst.get('background_color') is not None else 'grey'
    grid_color = mst.get('grid_color') if mst.get('grid_color') is not None else 'white'
    grid_linewidth = mst.get('grid_linewidth') if mst.get('grid_linewidth') is not None else 1
    grid_linestyle = mst.get('grid_linestyle') if mst.get('grid_linestyle') is not None else '--'
    
    dpi = mst.get('dpi') if mst.get('dpi') is not None else 100

    # Plot
    for mode in mode:
        if mode == 'log': clim[1] = np.log10(vmax)
        else: clim[1] = vmax
        for i in range(D['time'].shape[0]):
            z = np.log10(D['tecim'][i]) if mode == 'log' else D['tecim'][i]
            fig, ax = gm.plotCartoMap(figsize=figsize, projection=projection,
                                  title=D['time'][i],
                                  latlim=latlim,lonlim=lonlim,
                                  meridians=meridians, parallels=parallels,
                                  background_color=background_color, 
                                  grid_color=grid_color,grid_linewidth=grid_linewidth,
                                  grid_linestyle=grid_linestyle,
                                  figure=True, 
                                  nightshade=nightshade, ns_dt=D['time'][i])
            im = plt.pcolormesh(D['xgrid'],D['ygrid'], z.T, cmap=cmap, transform=ccrs.PlateCarree())
            plt.clim(clim)
            
            if colorbar:
                cax = fig.add_axes([0, 0, 0.1, 0.1])
                plt.colorbar(im, cax=cax, label=cbar_label)
                
                axp = fig.gca()
                posn = ax.get_position()
                axp.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                              0.02, posn.height])
                fig.canvas.draw()
            if show:
                plt.show()
            else:
                figdir = os.path.expanduser(os.path.join(odir,mode))
                if not os.path.exists(figdir):
                    print ('Creating a new directory: ', figdir)
                    subprocess.call('mkdir -p {}'.format(figdir), shell=True, timeout=5)
                fn = D['time'][i].strftime('%m%d_%H%M')
                figname = fn + str('.png')
                plt.savefig(os.path.join(figdir, figname), dpi=dpi)
                plt.close(fig=fig)

if __name__ == '__main__':
    
    p = ArgumentParser()
    
    p.add_argument('fn', type = str, help='converted hdf5 file')
    p.add_argument('-c', '--cfg', type = str, help='path to donfig.yaml file for map configuration. Default = ~/map/conus.yaml')
    p.add_argument('-o', '--odir', help = 'Destination folder, if None-> the same as input folder', default = None)
    p.add_argument('--show', help = 'set time limints for the file to cenvert', action='store_true')
    P = p.parse_args()
    
    plot(fnhdf = P.fn, cfg = P.cfg, odir = P.odir, show = P.show)