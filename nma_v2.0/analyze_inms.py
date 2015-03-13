#!/usr/bin/env python
# usage: ./analyze_inms.py > [output file]
# print top [numStructures] similarity structures from a pdb trajectory
# also bring inm-based similarity as a function of frame
# will add command line arguments...

highPass = 0.1 # percentage of top structures to admit
import cProfile
import re
import nma
import numpy as np
import sys
import time

#prefix = '/home/asmit/cas9/nma_codes/nma_python/'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_prep/'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/1ptq.pdb'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/1aqb.pdb'
#prefix = '/Users/paulglen/github/NMA/Python/sandbox/1vom.pdb'
#prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_test/em/em.gro'
prefix = '/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_test/'


ref = prefix + 'first.pdb'
#trajFile = prefix + 'tiny.pdb'
trajFile = prefix + 'traj.pdb'
#em = prefix
#ref = prefix
#trajFile = prefix

tstart = time.clock()
pr = cProfile.Profile()
pr.enable()
X = nma.TrajectoryINMs(ref, trajFile)
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
numStructures = int( highPass * lenX )
while itr < numStructures and itr < lenX :
    M = max(X)
    frame = X.index(M)
    outfile.write('%i %.2f\n'%(frame, X[frame]))
    X.remove(M)
    itr += 1
outfile.close()

tdiff = time.clock() - tstart
print 'time taken: ', tdiff/ 60., 'min'
