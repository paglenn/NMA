#!/usr/bin/env python
# take top [numStructures] similarity structures from a pdb trajectory
numStructures = 1
import enm
import numpy as np
import sys
import os
import ptraj
import time

tstart = time.clock()
def WriteReferenceModes(fi = None , ff = None ) :
    if not fi is None  :
        init = enm.enm(fi)
        vi = init.V
        #print vi
        outfile = open('ref_modes.dat','w')
        for i in range(vi.shape[0]-1):
                outfile.write('%.8f '%vi[i])
        outfile.write('%.8f\n'%vi[-1])
    if not ff is None:
        final = enm.enm(ff)
        vf = final.V
        for i in range(vf.shape[0]-1):
                outfile.write('%.8f '%vf[i])
        outfile.write('%.8f\n'%vf[-1])
    outfile.close()

def ReadReferenceModes():
	infile = open('ref_modes.dat','r')
	lines = [x.split() for x in infile.readlines()]
	V0 = [] ; V1 = []
	V0 = [float(x) for x in lines[0] ]
        if len(lines) > 1: V1 = [float(x) for x in lines[1] ]
	V0 = np.array(V0)
	V1 = np.array(V1)
	return V0, V1

if not os.path.isfile('ref_modes.dat'):
    fi = 'em/em3.gro'
    WriteReferenceModes(fi=fi )

# driver part

V0 , V1 = ReadReferenceModes()

trajFile = 'traj.pdb'

traj = ptraj.ReadInTrajectory(trajFile) # inb4AMBER...
init = enm.enm('topol/conf.gro')
X = []
for R in traj :
    current = enm.enm(coords=R )
    sim = current.Similarity(init)
    X.append(sim)
T = [0] + list(np.linspace(1.0,50.0,len(X)-1))

itr = 0
while itr < numStructures :
    M = max(X[1:])
    frame = X.index(M)
    print frame, X[frame]
    X.remove(M)
    itr += 1

outfile = open('inms.dat','w')
for i in range(len(X)):
    outfile.write('%.2f\t%f\n'%(T[i],X[i]))
outfile.close()

tdiff = time.clock() - tstart
print 'time taken: ', tdiff
