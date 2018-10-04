# -*- coding: utf-8 -*-
"""
Created on Sun May 20 20:56:20 2018

@author: smrak
"""
from numpy import ma, isnan
from datetime import datetime
from scipy import ndimage
from gpstec import gpstec as gtec
import os

force = True
saveh5 = True
date = '2008-03-26'
folder = 'C:\\Users\\smrak\\Google Drive\\BU\\Projects\\steve\\data\\gps\\'
hdffn = folder + 'struct_' + date.replace('-','') + '.h5'
timelim = [datetime(2008,3,26,0,0,0), datetime(2008,3,26,20,0,0)]

# 
savefolder = 'C:\\Users\\smrak\\Google Drive\\BU\\Projects\\steve\\images\\gpstec\\fill1\\'

if os.path.exists(hdffn) and force == False:
    t,xgrid,ygrid,im = gtec.readFromHDF(hdffn)
    gtec.save2HDF(t=t,lon=xgrid,lat=ygrid,images=im,h5fn=hdffn)
else:
    t,xgrid,ygrid,im = gtec.returnGlobaTEC(date=date,datafolder=folder,timelim=timelim)
    if saveh5:
        gtec.save2HDF(t=t,lon=xgrid,lat=ygrid,images=im,h5fn=hdffn)

######## MAP #########
projection='lambert'
figsize=(10,6)
latlim = [20,75]
lonlim= [-150,-70]
meridians = [-200,-180,-160,-140,-120,-100,-80,-60,-40,-20,0]
parallels = [-20,20,40,60,80]


for i in range(t.shape[0]):
    z = im[i]
    z = gtec.fillPixels(z,1)
    z = ndimage.median_filter(z, 3)
    image = ma.masked_where(isnan(z),z)
    title = 'GPSTEC: {} UT'.format(t[i])
    savefn = savefolder + datetime.strftime(t[i],'%H%M%S') + '.png'
    fig = gtec.plotTECmap(xgrid,ygrid,z,title=title,clim=[4,8],
                          figsize=figsize,projection=projection,
                          lonlim=lonlim,latlim=latlim,
                          meridians=meridians,parallels=parallels,
                          tight=False,savefn=savefn,DPI=200)
#    break