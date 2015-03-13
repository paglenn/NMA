# Paul Glenn
# enm.py
import sys
import numpy as np
import itertools as it
import time
cutoff = 1.3 # [nm] - distance cutoff for hessian calc
tc = 3 # number of chosen NormalModes for comparison

#------------------------------------------------
# file parsing stuff
def ReadInPDBFileCA(fn) :
	R = list()
	fin = open(fn, 'r')
	lines = fin.readlines()
	lines = [L for L in lines if 'CA' in L and 'ATOM' in L]
	fin.close()
	for line in lines:
		resID = int( line[22:26] )
		#print resID
		x  = line[30:38]
		y  = line[38:46]
		z  = line[46:54]
		coords = [float(i)/10. for i in [x,y,z]]
		R.append( np.array(coords) )
	return R

def ReadInGroFileCA(fn = 'em/em.gro'):
	R = list()
	fin = open(fn,'r')
	lines = [x for x in fin.readlines() if 'CA' in x ]
	fin.close()
	for line in lines:
		#if len(line) < 44 or not 'CA' in line: continue
		#print line
		x = line[20:28]
		y = line[28:36]
		z = line[36:44]
		coords = [float(i) for i in [x,y,z] ]
		R.append(np.array(coords,dtype='float64'))
	return R

def ReadInCACoordinates(fn):
	if fn[-3:] == 'pdb':
		R = ReadInPDBFileCA(fn)
	elif fn[-3:] == 'gro':
		R = ReadInGroFileCA(fn)
	else:
		sys.exit('File format unknown?')
		R = None
	return R

def ReadInPDBTrajectory(fn):
	R = list()
	fin = open(fn, 'r')
	lines = [L for L in fin.readlines() if 'CA' == L[13:15] or 'END' == L[:3] ]
	NL = len(lines)
	fin.close()
	traj = []
	for i in range(NL):
		line = lines[i]
		if 'END' in line or i + 1 == NL:
			traj.append(list(R))
			R[:] = []
		else:
			resID = int( line[22:26] )
			x  = line[30:38]
			y  = line[38:46]
			z  = line[46:54]
			coords = [float(i)/10.  for i in [x,y,z]]
			R.append( np.array(coords,dtype='float64') )

	return traj

def GetDistances(R):
	numRes = len(R)
	resList = range(numRes)
	D = np.zeros([ numRes , numRes])
	for i in resList:
		for j in resList[i+1:]:
			rij =  np.linalg.norm(R[i] - R[j])
			D[i,j] = rij
			D[j,i] = rij
	return D

# -----------------------------------------------
# Hessian calculation
# implementation of method from
# Zheng et al, Proteins 2010 78:2469-2481
# -----------------------------------------------
def computeHessian(R):
	NRES = len(R)
	N = 3 * NRES
	H = np.zeros((N, N))
	resList = range(NRES)
	Rij = GetDistances(R)

	PermsOfQ = list(it.product(range(3),range(3)) ) # indices of sub-matrix indicating (x,x), (x,y), (x,z), etc.

	#print 'Computing Hessian Matrix ...'
	for i in resList :

		for j in resList :

			_dij = Rij[i,j]

			if _dij < cutoff:

				if i == j :

					for k in resList:

						_dij = Rij[i,k]

						if _dij < cutoff and i != k:
							if i == 0 and k == 6: print 'Aah!'

							for qi,qj in PermsOfQ:

								dqi = (R[i][qi] - R[k][qi]) / _dij
								dqj = (R[i][qj] - R[k][qj]) / _dij
								h = dqi * dqj
								H[3*i+qi, 3*j + qj] = H[3*i+qi,3*j+qj] + h

				else:
					for qi,qj in PermsOfQ:

						_dij = Rij[i,j]
						if _dij < cutoff:
							dqi = ( R[i][qi] - R[j][qi]) / _dij
							dqj = (R[i][qj] - R[j][qj]) / _dij
							h = - dqi * dqj


						H[3*i+qi, 3*j + qj] = H[3*i+qi,3*j+qj] + h
	# print 'Hessian: '
	# print H
	return H

def DominantMode(evals, evecs):
	imax = len(evals) - 1
	i = 0
	while evals[i] <= 0  and i < imax: i += 1
	return evecs[i]

# --------------------------------
# Elastic Network Model
# --------------------------------
class ENM:

	def __init__(self,R):
		self.NRES = len(R)
		#print 'ENM-R',R
		self.N = 3* self.NRES # Number of modes
		self.R = np.copy(R)   # R=array of x,y,z coord

		# time hessian construction
		start = time.clock()
		self.H = computeHessian(R)
		print 'hessian construction', time.clock() - start, 'sec'

		# time eigendecomposition
		start = time.clock()
		self.EigenFreqs , self.NormalModes = np.linalg.eigh(self.H)
		print 'time for decomposition: ' , time.clock() - start, 'sec'
		#print 'cutoff: ', cutoff
		#print 'modes: ', self.EigenFreqs.shape
		#for i in range(self.N):
		#    print self.EigenFreqs[i], np.linalg.norm(self.NormalModes[i])


	#   self.EigenFreqs = np.fabs(self.EigenFreqs)
		if any (Freq < -1e-6 for Freq in self.EigenFreqs):
			print 'NEGATIVE FREQ',Freq
		idx = self.EigenFreqs.argsort() # idx is an array of indices
		self.EigenFreqs = self.EigenFreqs[idx]
		self.start = list(i > 1e-6 for i in self.EigenFreqs).index(True) # contains index for 1st non-zero mode
		self.NormalModes = self.NormalModes[idx]

	#--------------------------------------------------
	# Structural Similarity and INM overlap calculation
	# is based on Peng et al , Biophys J 2010 ; 98:2356
	# eqs (7) and (8)
	#--------------------------------------------------

	def ComputeOverlap(self,refMode):
		X = []
		firstIndex = self.start # ensure not using zero modes
		lastIndex = firstIndex + min(tc, self.N)
		for itr in xrange(firstIndex,lastIndex):
			vi = self.NormalModes[itr]
			dij =  np.fabs( np.dot(vi,refMode) )
			X.append(dij)
		di = max(X)
		#print 'Degree of Overlap'
		#print di
		return di

	# self needs to be given for python class methods
	def Similarity(self, ref_enm):
		s = 0 # similarity
		SumW = 0 # denominator in eq 4

		firstIndex = ref_enm.start # index of first non-zero mode
		lastIndex = firstIndex + min(tc, self.N)

		for itr in xrange(firstIndex,lastIndex):
			wi2 = ref_enm.EigenFreqs[itr]
			vi = ref_enm.NormalModes[itr]
			di = self.ComputeOverlap(vi)
			weight = 1.0/wi2
			s += weight * di
			SumW  += weight

		s /= SumW # ( normalization)
		return s


def TrajectoryINMs(refPDBFile, trajFile):


	#------------------------------------------------
	# get reference conformation for comparison
	Ri = ReadInCACoordinates(refPDBFile)
	#print Ri
	print 'READ REF STRUCT'
	init = ENM(Ri)

	#------------------------------------------------
	# read in trajectory and calculate similarity
	traj = ReadInPDBTrajectory(trajFile)
	print 'TRAJ PDB READ IN'
	#print traj
	#traj = [ReadInGroFileCA(trajFile)]
	X = []
	for R in traj :
		print 'new frame'
		curr = ENM(R)
# Using class method on curr class
		sim = curr.Similarity(init)
		X.append(sim)


	return X



