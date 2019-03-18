# -*- coding: utf-8 -*-
"""
Created on Sun May 20 20:56:20 2018

@author: smrak
"""
from numpy import ma, isnan
from datetime import datetime
from scipy import ndimage
from gpstec import returnGlobalTEC
#force = True
#saveh5 = True
date = '2010-10-18'
folder = 'C:\\Users\\smrak\\Google Drive\\BU\\Projects\\steve\\data\\gps\\'
folder = '/home/smrak/Documents/precipitation/data/'
hdffn = folder + 'struct_' + date.replace('-','')[2:] + '.h5'
timelim = [datetime(2008,3,26,6,0,0), datetime(2008,3,26,12,0,0)]

#
t,xgrid,ygrid,im = returnGlobalTEC(date=date,datafolder=folder)
t,xgrid,ygrid,im = gtec.readFromHDF(hdffn)
gtec.save2HDF(t=t,lon=xgrid,lat=ygrid,images=im,h5fn=hdffn)
#else:
#    
#    if saveh5:
#        gtec.save2HDF(t=t,lon=xgrid,lat=ygrid,images=im,h5fn=hdffn)

t,xgrid,ygrid,im = gtec.returnGlobaTEC(date=date,datafolder=folder,timelim=timelim)
gtec.save2HDF(t=t,lon=xgrid,lat=ygrid,images=im,h5fn=hdffn)

######## MAP #########
projection='lambert'
figsize=(10,6)
latlim = [58,70]
lonlim= [-165,-135]
meridians = [-200,-180,-160,-140,-120,-100,-80,-60,-40,-20,0]
parallels = [40,55,60,65,70,75,80]

# Plot
savefolder = '/home/smrak/Documents/steve/plots/gps/'
for i in range(t.shape[0]):
    z = im[i]
#    z = gtec.fillPixels(z,1)
    z = ndimage.median_filter(z, 3)
    image = ma.masked_where(isnan(z),z)
    title = 'GPSTEC: {} UT'.format(t[i])
    savefn = savefolder +'gpsmap_median_zoom_' + datetime.strftime(t[i],'%H%M%S') + '.png'
    fig = gtec.plotTECmap(xgrid,ygrid,z,title=title,clim=[4,8],
                          figsize=figsize,projection=projection,
                          lonlim=lonlim,latlim=latlim,
                          meridians=meridians,parallels=parallels,
                          tight=False,savefn=savefn,DPI=200)
#    break