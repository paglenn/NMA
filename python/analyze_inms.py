#!/usr/bin/env python
# print top [numStructures] similarity structures from a pdb trajectory
# also pring inm-based similarity as a function of frame
numStructures = 10
import nma
import numpy as np
import sys
import time
prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_prep/'

tstart = time.clock()

#------------------------------------------------
# get distances from minimized structure
em = prefix + 'em/em3.gro'
R0 = nma.ReadInCACoordinates(em)
_Rij = nma.GetDistances(R0)

#------------------------------------------------
# get reference conformation for comparison
ref = prefix + 'topol/conf.gro'
Ri = nma.ReadInGroFileCA(ref)
init = nma.ENM(Ri,_Rij)

#------------------------------------------------
# read in trajectory and calculate similarity
trajFile = prefix + 'traj.pdb'
traj = nma.ReadInPDBTrajectory(trajFile)
X = []
itr = 0
for R in traj :
    curr = nma.ENM(R,_Rij )
    for eig in curr.EigenFreqs: print eig
    sys.exit()
    sim = curr.Similarity(init)
    X.append(sim)
    itr+=1
lenX = len(X)

#------------------------------------------------
# print all similarities
ofn = prefix + 'inms.dat' # outfile name
outfile = open(ofn,'w')
for i in range(len(X)):
    outfile.write('%.2f\t%f\n'%(T[i],X[i]))
outfile.close()

#------------------------------------------------
# print top structures
itr = 0
ofn = prefix + 'top_structures.dat'
outfile = open(ofn, 'w')
while itr < numStructures and itr < lenX :
    M = max(X)
    frame = X.index(M)
    outfile.write('%i %.2f\n'%(frame, X[frame]))
    X.remove(M)
    itr += 1
outfile.close()

tdiff = time.clock() - tstart
print 'time taken: ', tdiff/ 60., 'min'
