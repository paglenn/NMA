# Paul Glenn
# enm.py
import sys
import numpy as np
import itertools as it
import time
cutoff = 1.3 # [nm] - distance cutoff for hessian calc
tc = 1 # number of chosen NormalModes for comparison

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
	fin.close()
	traj = []
	for line in lines:
		if 'END' in line:
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
def computeHessian(R,NRES=0 ):
	if NRES == 0 :
		NRES = len(R)
	N = 3 * NRES
	H = np.zeros((N, N))
	resList = range(NRES)
	Rij = GetDistances(R)

	PermsOfQ = list(it.product(range(3),range(3)) )

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

                                        dqi = ( R[i][qi] - R[k][qi]) / _dij
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
        #print H
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
		self.N = 3* self.NRES
		self.R = np.copy(R)

                start = time.clock()
		self.H = computeHessian(R ,NRES = self.NRES)
                print 'hessian construction', time.clock() - start, 'sec'
		start = time.clock()
		self.EigenFreqs , self.NormalModes = np.linalg.eigh(self.H)
		print 'time for decomposition: ' , time.clock() - start, 'sec'
                print 'cutoff: ', cutoff
		print 'modes: ', self.EigenFreqs.shape
		for i in range(self.N):
                    print self.EigenFreqs[i], np.linalg.norm(self.NormalModes[i])


		self.EigenFreqs = np.fabs(self.EigenFreqs)
		idx = self.EigenFreqs.argsort()
		self.EigenFreqs = self.EigenFreqs[idx]
		self.start = list(i > 1e-6 for i in self.EigenFreqs).index(True)
		self.NormalModes = self.NormalModes[idx]

	#--------------------------------------------------
	# Structural Similarity and INM overlap calculation
	# is based on Peng et al , Biophys J 2010 ; 98:2356
	# eqs (7) and (8)
	#--------------------------------------------------

	def ComputeOverlap(self,refMode):
		X = []
		firstIndex = self.start
		itr = firstIndex
		for wi in self.EigenFreqs[firstIndex:]:
			vi = self.NormalModes[itr]
			# ensure not using zero modes
			dij =  np.fabs( np.dot(vi,refMode) )
			X.append(dij)
			itr += 1
		di = max(X)
                print di
		return di

	def Similarity(self, ref_enm):
		s = 0 # similarity
		SumW = 0
		firstIndex = ref_enm.start
		lastIndex = firstIndex + min(tc, self.N)
		itr = firstIndex
		for wi in ref_enm.EigenFreqs[firstIndex:lastIndex]:
                    vi = ref_enm.NormalModes[itr]
                    weight = 1./wi
                    s += weight * self.ComputeOverlap(vi)
                    SumW  += weight
                    itr += 1
		s /= SumW
		return s


def TrajectoryINMs(em, ref, trajFile):

	#------------------------------------------------
	# get distances from minimized structure
	R0 = ReadInCACoordinates(em)
	_Rij = GetDistances(R0)

	#------------------------------------------------
	# get reference conformation for comparison
	Ri = ReadInCACoordinates(ref)
	init = ENM(Ri)

	#------------------------------------------------
	# read in trajectory and calculate similarity
	traj = ReadInPDBTrajectory(trajFile)
	#traj = [ReadInGroFileCA(trajFile)]
	X = []
	for R in traj :
		print 'new frame'
		curr = ENM(R)

		sim = curr.Similarity(init)
		X.append(sim)


	return X



