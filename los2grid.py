#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 10:17:12 2024

@author: mraks1
"""

import numpy as np
import h5py
import os, glob
from datetime import datetime
import matplotlib.pyplot as plt
import pyGnss
import pandas as pd
import xarray as xr
from scipy import ndimage
from astropy.convolution import convolve, Gaussian2DKernel
# grid

def getNeighbours(image,i,j,N=3):
    """
    Return an array of <=9 neighbour pixel of an image with a center at (i,j)
    """
    nbg = []
    m = int(np.floor(N/2))
    M = int(np.ceil(N/2))
    for k in np.arange(i-m, i+M):
        for l in np.arange(j-m, j+M):
            try:
                nbg.append(image[k,l])
            except:
                pass
    return np.array(nbg)

def fillPixels(im, N=1):
    """
    Fill in the dead pixels. If a dead pixel has a least 4 finite neighbour
    pixel, than replace the center pixel with a mean valuse of the neighbours
    """
    X = im.shape[0]-1
    Y = im.shape[1]-1
    imcopy = np.copy(im)
    for n in range(N):
        skip = int(np.floor((3+n)/2))
        starti = 0
        startj = 0
        forwardi = int(np.floor(0.7*X))
        backwardi = int(np.floor(0.3*X))
        if n%2 == 0:
            for i in np.arange(starti, forwardi, skip):
                for j in np.arange(startj, Y, skip):
                    # Check if th epixel is dead, i.e. empty
                    if np.isnan(im[i,j]):
                        # Get its neighbours as a np array
                        nbg = getNeighbours(imcopy,i,j,N=(3+n))
                        # If there are at leas 4 neighbours, replace the value with a mean
                        if sum(np.isfinite(nbg)) >= 4:
                            ix = np.where(np.isfinite(nbg))[0]
                            avg = np.mean(nbg[ix])
                            im[i,j] = avg
            for i in np.arange(X, backwardi, -skip):
                for j in np.arange(Y, 0, -skip):
                    # Check if th epixel is dead, i.e. empty
                    if np.isnan(im[i,j]):
                        # Get its neighbours as a np array
                        nbg = getNeighbours(imcopy,i,j,N=(3+n))
                        # If there are at leas 4 neighbours, replace the value with a mean
                        if sum(np.isfinite(nbg)) >= 4:
                            ix = np.where(np.isfinite(nbg))[0]
                            avg = np.mean(nbg[ix])
                            im[i,j] = avg
        else:
            for j in np.arange(startj, Y, skip):
                for i in np.arange(starti, forwardi, skip):
                    # Check if th epixel is dead, i.e. empty
                    if np.isnan(im[i,j]):
                        # Get its neighbours as a np array
                        nbg = getNeighbours(imcopy,i,j,N=(3+n))
                        # If there are at leas 4 neighbours, replace the value with a mean
                        if sum(np.isfinite(nbg)) >= 4:
                            ix = np.where(np.isfinite(nbg))[0]
                            avg = np.mean(nbg[ix])
                            im[i,j] = avg

            for j in np.arange(Y, 0, -skip):
                for i in np.arange(X, backwardi, -skip):
                    # Check if th epixel is dead, i.e. empty
                    if np.isnan(im[i,j]):
                        # Get its neighbours as a np array
                        nbg = getNeighbours(imcopy,i,j,N=(3+n))
                        # If there are at leas 4 neighbours, replace the value with a mean
                        if sum(np.isfinite(nbg)) >= 4:
                            ix = np.where(np.isfinite(nbg))[0]
                            avg = np.mean(nbg[ix])
                            im[i,j] = avg
    return im

def ImageNew(glon, glat, tid, 
             latlim=[-90,90], lonlim=[-180,180], res=None,
             filter_type='gaussian', sigma=2, filter_size=5):
    """
    """
    xgrid, ygrid = np.meshgrid(np.arange(lonlim[0], lonlim[1]+.01, res),
                               np.arange(latlim[0], latlim[1]+.01, res))
    im = np.empty(xgrid.shape, dtype=object)
    # Fill out the image pixels
    for i in range(glon.size):
        idx = abs(xgrid[0, :] - glon[i]).argmin() if abs(xgrid[0, :] - glon[i]).min() < 3*res else np.nan
        idy = abs(ygrid[:, 0] - glat[i]).argmin() if abs(ygrid[:, 0] - glat[i]).min() < 3*res else np.nan
        # If image indexes are valid
        if np.isfinite(idx) and np.isfinite(idy):
            # Assign the value to the pixel
            if im[idy,idx] is None:
                im[idy,idx] = [tid[i]]
            # If this is not the first value to assign, assign a
            # mean of both values
            else:
                im[idy,idx].append(tid[i])
    
    imout = np.nan * np.empty(xgrid.shape)
    for i in range(xgrid.shape[0]):
        for j in range(xgrid.shape[1]):
            if im[i,j] is not None:
                imout[i,j] = np.nanmedian(im[i,j])
    if filter_type == 'median':
        imout = fillPixels(imout)
        imout = ndimage.median_filter(imout, filter_size)
    elif filter_type == 'gaussian':
        kernel = Gaussian2DKernel(x_stddev=sigma, y_stddev=sigma, x_size=filter_size, y_size=filter_size)
        imout = convolve(imout, kernel)
        
        imout[:filter_size, :] = np.nan
        imout[:, :filter_size] = np.nan
        imout[-filter_size:, :] = np.nan
        imout[:, -filter_size:] = np.nan

    del im
    return xgrid, ygrid, imout

def main(f, el_mask, altkm, out):
    if out is None:
        out = os.path.split(f)[0] + os.sep
    D = h5py.File(f, 'r')
    M = np.array([a[0].decode('ascii') for a in D['Metadata/Data Parameters'][:]])
    S = D['Data/Table Layout'].size
    iter = np.arange(0, S, 10000e3, dtype=int)
    iter = np.append(iter, S-1)
    ut1 = np.zeros(D['Data/Table Layout'].size, dtype=int)
    for ii, i in enumerate(iter[1:]):
        ut1[iter[ii]:i] = np.asarray(pd.DataFrame(D['Data/Table Layout'][iter[ii]:i])['ut1_unix'])
    
    utu = np.unique(np.delete(ut1, ut1==0))
    dt = np.array([datetime.utcfromtimestamp(t) for t in utu])
    
    for i, ut in enumerate(utu):
        print (dt[i])
        idt = np.isin(ut1, ut)
        az, el = np.array(pd.DataFrame(D['Data/Table Layout'][idt])['azm']), np.array(pd.DataFrame(D['Data/Table Layout'][idt])['elm']) 
        if el_mask is not None:
            el[el<el_mask] = np.nan
            az[el<el_mask] = np.nan
            idf = np.logical_and(np.isfinite(el), np.isfinite(az))
        if 'DALTR' in M: 
            rxp =  np.array(pd.DataFrame(D['Data/Table Layout'][idt])['gdlatr']), np.array(pd.DataFrame(D['Data/Table Layout'][idt])['gdlonr']), np.array(pd.DataFrame(D['Data/Table Layout'][idt])['galtr'])
        else:
            rxp =  np.array(pd.DataFrame(D['Data/Table Layout'][idt])['gdlatr']), np.array(pd.DataFrame(D['Data/Table Layout'][idt])['gdlonr']), np.zeros(np.sum(idt))
        ipp_lla = pyGnss.aer2ipp(az[idf], el[idf], np.asanyarray(rxp)[:,idf].T, H=altkm)
        xg, yg, im = ImageNew(ipp_lla[1], ipp_lla[0], np.array(pd.DataFrame(D['Data/Table Layout'][idt])['tec'])[idf], res=1, filter_type=None)
        if i == 0:
            IM = np.copy(im)
        else:
            IM = np.dstack((IM, im))
    
    X = xr.Dataset(coords={'time': dt, 'glon': xg[0,:], 'glat': yg[:,0]})
    X['tec'] = (('glat', 'glon', 'time'), IM)
    X['elmask'] = el_mask
    X['altkm'] = altkm
    X.to_netcdf(out+f'{os.sep}conv_from_los_{dt[0].strftime("%Y%m%d")}_{altkm}km_{el_mask}el.nc')
        
    
if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('ifn', help='input file')
    p.add_argument('--odir', help='output dirrectory', default=None)
    p.add_argument('--elmask', default=30, type=float)
    p.add_argument('--altkm',  default=350, type=float)
    P = p.parse_args()
    
    
    main(P.ifn, el_mask=P.elmask, out=P.odir, altkm=P.altkm)
    
    
    
    
    