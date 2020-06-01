# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 15:32:25 2019

@author: smrak
"""
from glob import glob
from gpstec import gpstec
#from datetime import datetime
import numpy as np
import os, platform
import yaml
import subprocess
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartomap import geogmap as gm
from argparse import ArgumentParser

#mcfg = '/home/smrak/Documents/cartomap/map/conus.yaml'

def _round(x, base=5):
    r = base * round(x/base)
    if r < x:
        return r + 5
    else:
        return r
    return 

def plot(fnhdf:str = None,
         odir:str = None,
         cfg:str = None,
         clim:list = None, mode:str = None, average:int=None, cmap=None,
         projection=None, lonlim=None, latlim=None, terminator=None):
    C0 = clim
    fnhdf = os.path.expanduser(fnhdf.strip())
    assert fnhdf is not None 
    assert os.path.splitext(fnhdf)[1] in ['.h5', '.hdf5']
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
    if projection is None:
        projection = mst.get('projection')
    if cmap is None:
        cmap = mst.get('cmap')
    if lonlim is None:
        lonlim = mst.get('lonlim')
    if latlim is None:
        latlim = mst.get('latlim')
    try:
        nightshade = mst.get('nightshade')
    except:
        pass
    try:
        paralim = mst.get('parallels')
        meridilim = mst.get('meridians')
        parallels = np.arange(paralim[0], paralim[1], paralim[2])
        meridians = np.arange(meridilim[0], meridilim[1], meridilim[2])
    except:
        parallels = []
        meridians = []
    
    try:
        apex = True
        mlatlim = mst.get('mlat')
        mltlim = mst.get('mlt')
        mlat_levels = np.arange(mlatlim[0], mlatlim[1]+.1, mlatlim[2])
        mlon_levels = np.arange(mltlim[0], mltlim[1]+.1, mltlim[2])
    except:
        apex = False
        mlat_levels = None
        mlon_levels = None
    
    
    if mode is None:
        mode = 'lin'
    if isinstance(mode, str):
        mode = [mode]
        
    colorbar = mst.get('colorbar')
    cbar_label = mst.get('cbar_label') if mst.get('cbar_label') is not None else 'TEC [TECu]'
    if clim is None:
        clim = mst.get('clim')
    
    if clim is None:
        vmax = 3 * _round(np.nanmedian(D['tecim']), base=5)
        vmin = 0
        clim = [vmin, vmax]
    else:
        vmax = float(clim[1])
        vmin = float(clim[0])
        clim = [vmin, vmax]
        
    figsize = mst.get('figsize') if mst.get('figsize') is not None else [12,8]
    background_color = mst.get('background_color') if mst.get('background_color') is not None else 'grey'
    grid_color = mst.get('grid_color') if mst.get('grid_color') is not None else 'white'
    grid_linewidth = mst.get('grid_linewidth') if mst.get('grid_linewidth') is not None else 1
    grid_linestyle = mst.get('grid_linestyle') if mst.get('grid_linestyle') is not None else '--'
    
    dpi = int(mst.get('DPI')) if mst.get('DPI') is not None else 50

    # Plot
    for mode in mode:
        if mode == 'log': clim[1] = np.log10(vmax)
        else: clim[1] = vmax
        iterate = np.arange(0, D['time'].size, average)
        for i in iterate: #range(D['time'].shape[0]):
            t0 = D['time'][i]
            z = np.log10(np.nanmean(D['tecim'][i:i+average], axis=0)) if mode == 'log' else np.nanmean(D['tecim'][i:i+average], axis=0)
            fig, ax = gm.plotCartoMap(latlim=latlim, lonlim=lonlim, 
                                  projection='merc', #lon0 = -90,#glons[idmidnight],
                                  title = t0, 
                                  #lat0 = 40,
                                  meridians=None, parallels=None, figsize=figsize,
                                  background_color='gray', border_color='k', states=0,
                                  apex=True,mlat_labels=0,mlon_labels=0,
                                  mlat_levels=mlat_levels,
                                  mlon_levels=mlon_levels,
                                  date=t0, 
                                  mlon_colors='w', mlat_colors='w',mlon_cs='mlt',
                                  terminator=True, ter_color='r', terminator_altkm=350
                                  )
            
            print ("Plotting {}".format(D['time'][i]))
            im = plt.pcolormesh(D['xgrid'], D['ygrid'], z.T, cmap=cmap, transform=ccrs.PlateCarree())
            plt.clim(clim)
            
            if colorbar:
                posn = ax.get_position()
                cax = fig.add_axes([posn.x0+posn.width+0.01, posn.y0, 0.02, posn.height])
                fig.colorbar(im, cax=cax, label='TEC [TECu]')
            
            if odir is None:
                odir = os.path.join(os.path.split(fnhdf)[0], D['time'][i].strftime("%Y%m%d"))
            svdir = projection + '_' + mode + '_' + str(int(clim[0])) + '-' + str(int(clim[1]))
            figdir = os.path.expanduser(os.path.join(odir, svdir))
            if not os.path.exists(figdir):
                print ('Creating a new directory: ', figdir)
                if platform.system == 'Linux':
                    subprocess.call('mkdir -p {}'.format(figdir), shell=True, timeout=5)
                else:
                    subprocess.call('mkdir "{}"'.format(figdir), shell=True, timeout=5)
            fn = D['time'][i].strftime('%m%d_%H%M')
            figname = fn + str('.png')
            plt.savefig(os.path.join(figdir, figname), dpi=dpi)
            plt.close(fig=fig)

if __name__ == '__main__':
    
    p = ArgumentParser()
    
    p.add_argument('fn', type = str, help='converted hdf5 file')
    p.add_argument('-c', '--cfg', type = str, help='path to donfig.yaml file for map configuration. Default = ~/map/conus.yaml')
    p.add_argument('-o', '--odir', help = 'Destination folder, if None-> the same as input folder', default = None)
    p.add_argument('--clim', nargs=2, default=None)
    p.add_argument('--projection', default=None, type=str)
    p.add_argument('--lonlim', nargs=2, type=float)
    p.add_argument('--latlim', nargs=2, type=float)
    p.add_argument('--mode', type=str, default=None)
    p.add_argument('--average', type=int, default=1)
    p.add_argument('--cmap', type=str, default='jet')
    p.add_argument('--terminator', action='store_true')
    
    P = p.parse_args()

    plot(fnhdf = P.fn, cfg = P.cfg, odir = P.odir, mode=P.mode,
         clim=P.clim, cmap=P.cmap, projection=P.projection, average=P.average,
         lonlim=P.lonlim, latlim=P.latlim, terminator=P.terminator)