# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 20:25:26 2020

@author: smrak@bu.edu
"""

from gpstec import gpstec
from datetime import datetime
import numpy as np
from typing import Union
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

def interpolateTEC(im: Union[list,np.ndarray] = None,
                          x0 = None, y0 = None,
                          xgrid = None, ygrid = None,
                          N: int = 512, res=1,
                          method: str = 'linear'):
    assert im is not None, 'Invalid input argument. Has to be a list or np.ndarray with a length of at least 1'
    if x0 is None or y0 is None:
        x0, y0 = np.meshgrid(np.arange(im.shape[0]),
                             np.arange(im.shape[1]))
    x0 = x0.T
    y0 = y0.T
    mask = np.ma.masked_invalid(im)
    x0 = x0[~mask.mask]
    y0 = y0[~mask.mask]
    X = im[~mask.mask]
    if xgrid is None or ygrid is None:
        xgrid, ygrid = np.meshgrid(np.arange(0, im.shape[0], res), 
                                   np.arange(0, im.shape[1], res))
    xgrid = xgrid.T
    ygrid = ygrid.T
    z = griddata((x0,y0), X.ravel(), (xgrid, ygrid), method=method, fill_value=np.nan)
    return z

def interp(values, vtx, wts, fill_value=np.nan):
    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    ret[np.any(wts < 0, axis=1)] = fill_value
    return ret

def _interpWeights(xyz, uvw, d=2):
    tri = Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def _interpolate(values, vtx, wts, fill_value=np.nan):
        ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
        ret[np.any(wts < 0, axis=1)] = fill_value
        return ret

tlim = [datetime(2016,5,8,0), datetime(2016,5,8,5)]
lonlim = [-140, -55]
latlim = [0, 60]
res = 0.3

fnhdf = 'G:\\My Drive\\scintillation_data\\tid\\20160508\\conv_20160508T0000-20160509T0000.h5'

D = gpstec.readFromHDF(fnhdf)
idt = (D['time'] >= tlim[0]) & (D['time'] <= tlim[1])
idx = (D['xgrid'] >= lonlim[0]) & (D['xgrid'] <= lonlim[1])
idy = (D['ygrid'] >= latlim[0]) & (D['ygrid'] <= latlim[1])

xgt, ygt = np.meshgrid(np.arange(D['xgrid'][idx][0], D['xgrid'][idx][-1]+0.1, 1),
                       np.arange(D['ygrid'][idy][0], D['ygrid'][idy][-1]+0.1, 1))

T0t = D['tecim'][idt]
T0x = T0t[:, idx, :]
T0 = T0x[:, :, idy]

xg, yg = np.meshgrid(np.arange(D['xgrid'][idx][0], D['xgrid'][idx][-1]+0.1, res),
                       np.arange(D['ygrid'][idy][0], D['ygrid'][idy][-1]+0.1, res))

z = T0[0].T
T00 = interpolateTEC(im=T0[0], x0=xgt, y0=ygt, xgrid=xg, ygrid=yg, method='linear')

plt.figure()
plt.pcolormesh(xgt.T, ygt.T, T0[0], vmin=0, vmax=30, cmap='jet')
plt.figure()
plt.pcolormesh(xg.T, yg.T, T00, vmin=0, vmax=30, cmap='jet')