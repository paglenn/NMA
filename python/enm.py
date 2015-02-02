import sys
import numpy as np
import ptraj
import scipy.weave
cutoff = 1.5 # [nm]

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

def computeHessian(R, _Rij ):
    #print 'Calculating distances...'
    numRes = len(R)
    H = np.zeros((3*numRes, 3*numRes))
    resList = range(numRes)
    Rij = GetDistances(R)
    print 'Computing Hessian Matrix ...'
    for i in resList :

            for j in resList[i+1:] :

                    Hij = np.zeros([3*numRes,3*numRes])
                    dij = Rij[i,j] ; _dij = _Rij[i,j]
                    ddij = (dij - _dij) / dij

                    if _dij < cutoff:

                            for ip in [i,j]:

                                    for jp in [i,j]:

                                            for qi in range(3):

                                                    for qj in range(3):

                                                            dqi = ( R[ip][qi] - R[jp][qi]) / dij
                                                            dqj = (R[ip][qj] - R[jp][qj]) / dij
                                                            if ip == jp :
                                                                    if qi == qj :
                                                                            h = ddij + (1-ddij) * dqi * dqi
                                                                    else:
                                                                            h = (1-ddij) * dqi * dqj
                                                            else:
                                                                    if qi == qj:
                                                                            h  = - (ddij + (1-ddij) * dqi * dqj)
                                                                    else:
                                                                            h  = - (1-ddij) * dqi * dqj

                                                            H[3*ip+qi, 3*jp + qj] += h
                    #H = H + Hij

            return H

def DominantMode(evals, evecs):
	idx = evals.argsort()
	imax = len(evals) - 1
	i = 0
	while evals[idx[i]] <= 0 and i < imax:
		i += 1
	return evecs[idx[i]]

R0 = ptraj.ReadInCACoordinates('em/em3.gro')
_Rij = GetDistances(R0)


class enm:

	def __init__(self,fileName=None, coords=None ):

		if fileName is None and coords is not None:
			self.R = coords
		elif fileName is not None and coords is None:
			self.R = ptraj.ReadInCACoordinates(fileName)
		else:
			sys.exit('Must supply EITHER path to file with coordinates OR list with coordinates!')
		self.H = computeHessian(self.R,np.copy(_Rij))
		self.freqs , self.modes = np.linalg.eigh(self.H)
		self.V = DominantMode(self.freqs, self.modes )

	def ComputeOverlap(self,Vref):
		X = []
		for i in range(len(self.freqs)):
			wi = self.freqs[i]
			vi = self.modes[i]
			if np.linalg.norm(vi) != 0 and wi != 0:
				dij = np.fabs( np.dot(vi,Vref) )  / (np.linalg.norm(vi) * np.linalg.norm(Vref))
				X.append(dij)
		di = np.amax(X)
		return di


	def Similarity(self, ref_sys):
		s = 0
		SumW = 0
		for i in range(len(ref_sys.freqs)):
			wi = ref_sys.freqs[i]
			vi = ref_sys.modes[i]
			if wi > 0:
				weight = (1./wi) ** 2.
				si = weight * self.ComputeOverlap(vi)
				SumW  += weight
				s += si
		s /= SumW
		return s





