#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 14:28:47 2019

@author: smrak
"""
import os
import numpy as np
import madrigalWeb.madrigalWeb
import subprocess
#from dateutil import parser
from datetime import datetime

#fixpath = 1
#savedir = 'C:\\Users\\smrak\\google drive\\bu\\projects\\LWS2019\\'
IX = {'dst': 212, 'aei': 211, 'geo': 210}

def dlIndex(#t0:str = None, #t1:str = None, 
            savedir:str = None,
            index:str = 'dst',
            user_fullname:str = None,
            user_email:str = None,
            user_affiliation:str = None,
            fixpath:bool = False):
    
#    assert t0 is not None
#    assert t1 is not None
    assert savedir is not None
    
    if user_fullname is None: 
        user_fullname = 'Sebastijan Mrak'
    if user_email is None: 
        user_email = 'smrak@bu.edu'
    if user_affiliation is None: 
        user_affiliation = 'BU'
    
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
    #
    #for inst in instList:
    #    if inst.code == 212:
    #        print((str(inst) + '\n'))
    # ======================= Get/List of Experiments ======================= #
    # GPS Minimum Scallop TEC: 3500
    # GPS LOS: ID: 3505
    # DST ID: 30006
    # Geophysical Ind ID: 30007
    # AE ID: 30008

#t0 = '2015-7-22'
#t1 = '2015-9-30'
    T = datetime.now()
#    T0 = parser.parse(t0)
#    T1 = parser.parse(t1)
    expList = MD.getExperiments(IX[index], 
                                1900, 1, 1, 0, 0, 1,
                                T.year, T.month, T.day, 0, 0, 0)
                                #T1.year, T1.month, T1.day, 23, 59, 59)
    #for exp in expList:
        # should be only one
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
                parts = path.split(os.sep)
                path_fn = os.path.split(path)[1]
                print (path_fn)
                if path_fn[:3] == index:
                    if not fixpath:
                        p = savedir + '/'.join(parts[4:])
                        savefn = os.path.join(p)
                        savefnlist.append(os.path.expanduser(savefn))
                    else:
                        savefn = savedir + os.path.split(path)[1]
                        print (savefn)
                        savefnlist.append(os.path.expanduser(savefn))
                    fnlist.append(this_file.name)
    
    # Check for direcotories:
    for ofn in savefnlist:
        head = os.path.split(ofn)[0]
        if not os.path.exists(head):
            subprocess.call("mkdir -p {}".format(head), timeout=10, shell=True)
    
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
#    p.add_argument('t0', type=str, help='Time limit 1 YYYY-mm-dd')
#    p.add_argument('t1', type=str, help='Time limit 2 YYYY-mm-dd')
    p.add_argument('odir', type=str, help='Output directory root')
    p.add_argument('-i', '--index', help='Which index? dst, ae, geo', default='dst')
    p.add_argument('--fixpath', help='Save to exact directory', action='store_true')
    p.add_argument('--name', type=str, help='"Full name"')
    p.add_argument('--email', type=str, help='"email"')
    p.add_argument('--affiliation', type=str, help='"affiliation"')
    
    P = p.parse_args()
    
    dlIndex(#t0 = P.t0, #t1 = P.t1, 
            savedir = P.odir, 
            index=P.index,
             user_fullname = P.name,
             user_email = P.email,
             user_affiliation = P.affiliation,
             fixpath = P.fixpath)
    