#!/usr/bin/env python
# print top [numStructures] similarity structures from a pdb trajectory
# also pring inm-based similarity as a function of frame
numStructures = 10
import cProfile
import re
import nma
import numpy as np
import sys
import time
prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_prep/'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_test/'

em = prefix + 'em/em.gro'
ref = prefix + 'topol/conf.gro'
trajFile = prefix + 'traj.pdb'

tstart = time.clock()
pr = cProfile.Profile()
pr.enable()
X = nma.TrajectoryINMs(em, ref, trajFile)
pr.disable()
pr.print_stats()

tdiff = time.clock() - tstart
print 'time taken: ', tdiff/ 60., 'min'

lenX = len(X)

#------------------------------------------------
# print all similarities
ofn = prefix + 'inms.dat' # outfile name
outfile = open(ofn,'w')
for i in range(len(X)):
    outfile.write('%i %f\n'%(i,X[i]))
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
