#!/usr/bin/env python
import numpy as np
import time
import sys
#def ReadInGro(f = 'conf.gro'):

def ReadInPDBFileCA(f) :
    R = list()
    infile = open(f, 'r')
    lines = [L for L in infile.readlines() if 'ATOM' in L and 'CA' in L ]
    for line in lines:
        resID = int( line[22:26] )
        #print resID
        x  = line[30:38]
        y  = line[38:46]
        z  = line[46:54]
        coords = [float(i)/10. for i in [x,y,z]]
        R.append( np.array(coords, dtype='float64') )
        #print R[-1]
    infile.close()
    return R

def ReadInGroFileCA(f = 'em/em.gro'):
    R = list()
    infile = open(f,'r')
    lines = [x for x in infile.readlines() if 'CA' in x ]
    infile.close()
    for line in lines:
        #if len(line) < 44 or not 'CA' in line: continue
        #print line
        x = line[20:28]
        y = line[28:36]
        z = line[36:44]
        coords = [float(i) for i in [x,y,z] ]
        R.append(np.array(coords,dtype='float64'))
    return R

def ReadInCACoordinates(fileName):
    if fileName[-3:] == 'pdb':
        R = ReadInPDBFileCA(fileName)
    elif fileName[-3:] == 'gro':
        R = ReadInGroFileCA(fileName)
    else:
        sys.exit('File format unknown?')
        R = None
    return R

def ReadInTrajectory(f):
    R = list()
    infile = open(f, 'r')
    lines = [L for L in infile.readlines() if ('ATOM' in L and 'CA' in L) or 'END' in L ]
    infile.close()
    traj = []
    for line in lines:
        if 'END' in line:
            traj.append(list(R))
            R[:] = []
            continue
        #print resID
        resID = int( line[22:26] )
        x  = line[30:38]
        y  = line[38:46]
        z  = line[46:54]
        coords = [float(i)/10.  for i in [x,y,z]]
        R.append( np.array(coords,dtype='float64') )
        #print R[-1]
    return traj
