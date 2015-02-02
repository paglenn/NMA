# Paul Glenn
# enm.py
import sys
import numpy as np
import itertools as it
cutoff = 1.5 # [nm] - distance cutoff for hessian calc
tc =  5 # number of chosen NormalModes for comparison

#------------------------------------------------
# file parsing stuff
def ReadInPDBFileCA(fn) :
    R = list()
    fin = open(fn, 'r')
    lines = fin.readlines()
    lines = [L for L in lines if 'CA' == L[13:15]]
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
            continue
        #print resID
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
def computeHessian(R, _Rij ):
    numRes = len(R)
    H = np.zeros((3*numRes, 3*numRes))
    resList = range(numRes)
    Rij = GetDistances(R)

    PermsOfQ = list(it.product(range(3),range(3)) )
    # the it.product is a one-time use object, here we make it permanent

    print 'Computing Hessian Matrix ...'
    for i in resList :

        for j in resList[i+1:] :

            _dij = _Rij[i,j]

            if _dij < cutoff:
                dij = Rij[i,j]
                ddij = (dij - _dij) / dij

                for ip in [i,j]:

                    for jp in [i,j]:

                        for qi,qj in PermsOfQ:

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

    def __init__(self,R,_Rij):
        # eigenvalues: self.EigenFreqs
        # eigenvectors: self.NormalModes

        self.R = R
        self.H = computeHessian(self.R,np.copy(_Rij))
        self.EigenFreqs , self.NormalModes = np.linalg.eigh(self.H)

        self.EigenFreqs = np.fabs(self.EigenFreqs)
        idx = self.EigenFreqs.argsort()
        self.EigenFreqs = self.EigenFreqs[idx]
        self.NormalModes = self.NormalModes[idx]
        self.V = DominantMode(self.EigenFreqs, self.NormalModes )

    #--------------------------------------------------
    # Structural Similarity and INM overlap calculation
    # is based on Peng et al , Biophys J 2010 ; 98:2356
    # eqs (7) and (8)
    #--------------------------------------------------

    def ComputeOverlap(self,Vref):
        X = []
        for i in range(len(self.EigenFreqs)):
            wi = self.EigenFreqs[i]
            vi = self.NormalModes[i]
            if np.linalg.norm(vi) != 0 and wi > 0:
                dij = np.fabs( np.dot(vi,Vref) )  / (np.linalg.norm(vi) * np.linalg.norm(Vref))
                X.append(dij)
        di = np.amax(X)
        return di

    def Similarity(self, ref_enm):
        s = 0 # similarity
        SumW = 0
        for i in range(tc):
            wi = ref_enm.EigenFreqs[i]
            vi = ref_enm.NormalModes[i]
            if wi > 0:
                weight = (1./wi) ** 2.
                si = weight * self.ComputeOverlap(vi)
                SumW  += weight
                s += si
        s /= SumW
        return s





