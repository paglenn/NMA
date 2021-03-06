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
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/1ptq.pdb'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/1aqb.pdb'
#prefix = '/Users/paulglen/github/NMA/Python/sandbox/1vom.pdb'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_test/em/em.gro'


ref = prefix + 'first.pdb'
em = ref
trajFile = prefix + 'traj.pdb'
#trajFile = prefix + 'tiny.pdb'
#em = prefix
#ref = prefix
#trajFile = prefix

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
ofn = prefix + 'inms_2.dat' # outfile name
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
