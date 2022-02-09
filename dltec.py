#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 14:28:47 2019

@author: smrak
"""
import os, platform, subprocess, yaml, datetime
import numpy as np
import madrigalWeb.madrigalWeb
from dateutil import parser

def dlGPSTEC(t0:str = None, t1:str = None, savedir:str = None,
             fixpath:bool = False,
             los: bool = False):
    
    assert t0 is not None
    assert t1 is not None
    assert savedir is not None
    cfg = os.getcwd() + '.affil.yaml'
    if not os.path.exists(cfg):
        dct = {}
        tmp = input('We need some info for the MardigalWeb interface.\nType you full name: ')
        dct['username'] = tmp
        tmp = input('Type your email: ')
        dct['email'] = tmp
        tmp = input('Type your affiliation: ')
        dct['affiliation'] = tmp
        with open(cfg, 'w') as outfile:
            yaml.dump(dct, outfile)
        del dct, tmp
    stream = yaml.load(open(cfg, 'r'), Loader=yaml.BaseLoader)
    user_fullname = stream.get('username')
    user_email = stream.get('email')
    user_affiliation = stream.get('affiliation')
    
    key = 'los' if los else 'gps'
    
    # Open Madrigal database
    madrigalUrl = 'http://cedar.openmadrigal.org/'
    MD = madrigalWeb.madrigalWeb.MadrigalData(madrigalUrl)
    # ======================== List of instruments ========================== #
    # GPS network ID: 8000
    # IMF and Sw IID: 120
    # Geophysical Indisies ID: 210
    # AE ID: 211
    # DST ID: 212
    
    #instList = MD.getAllInstruments()
    
    #for inst in instList:
    #    if inst.code == 8000:
    #        print((str(inst) + '\n'))
    # ======================= Get/List of Experiments ======================= #
    # GPS Minimum Scallop TEC: 3500
    # GPS LOS: ID: 3505
    # DST ID: 30006
    # Geophysical Ind ID: 30007
    # AE ID: 30008
    
    T0 = parser.parse(t0)
    T1 = parser.parse(t1)
    expList = MD.getExperiments(8000, 
                                T0.year, T0.month, T0.day, 0, 0, 1,
                                T1.year, T1.month, T1.day, 23, 59, 59)
    #for exp in expList:
    #    # should be only one
    #    print((str(exp) + '\n'))
    
    # ==================== Get links-filenames and output paths ============  #
    ids = [n.id for n in expList]
    file_list = np.array([MD.getExperimentFiles(n) for n in ids])
    fnlist = []
    savefnlist = []
    for subarr in file_list:
        for this_file in subarr:
            if this_file.category == 1:
                path = os.path.expanduser(this_file.name)
                parts = path.split('/')
                parts.remove('gps')
                parts[-2] = datetime.datetime.strptime(parts[-2], '%d%b%y').strftime('%m%d')
                path_fn = os.path.split(path)[1]
                if path_fn[:3] == key:
                    if not fixpath:
                        print (parts)
                        p = savedir + os.sep.join(parts[4:])
                        savefn = os.path.join(p)
                        savefnlist.append(savefn)
                    else:
                        savefn = savedir + os.path.split(path)[1]
                        savefnlist.append(savefn)
                    fnlist.append(this_file.name)
    
    # Check for direcotories:
    for ofn in savefnlist:
        print (ofn)
        head = os.path.split(ofn)[0]
        if not os.path.exists(head):
            if platform.system() in ('Linux', 'Darwin'):
                subprocess.call("mkdir -p {}".format(head), timeout=10, shell=True)
            elif platform.system() == 'Windows':
                subprocess.call('mkdir "{}"'.format(head), timeout=10, shell=True)
    
    for i in range(len(savefnlist)):
        if not os.path.exists(savefnlist[i]):
            print ('Downloading {}'.format(os.path.split(savefnlist[i])[1]))
            MD.downloadFile(fnlist[i], savefnlist[i], 
                            user_fullname, user_email, user_affiliation, 
                            format='hdf5')
        else:
            print ("{} already exists".format(savefnlist[i]))

if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument('t0', type=str, help='Time limit 1 YYYY-mm-dd')
    p.add_argument('t1', type=str, help='Time limit 2 YYYY-mm-dd')
    p.add_argument('odir', type=str, help='Output directory root')
    p.add_argument('--fixpath', help='Save to exact directory', action='store_true')
    p.add_argument('--los', help='Get line-of-sight data', action='store_true')
    
    P = p.parse_args()
    
    dlGPSTEC(t0 = P.t0, t1 = P.t1, savedir = P.odir, 
             los = P.los, fixpath=P.fixpath)