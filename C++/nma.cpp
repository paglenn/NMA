/* ptraj.cpp 
 * Author : Paul Glenn 
 * Date: Feb 1 2015 
 */ 
#define CUTOFF 1.5
#include<iostream>
#include<fstream>
#include<vector> 
#include<string>
#include<cstdlib>
#include<sstream>
#include<cholmod.h>
#include<Eigen/Dense>
#include<cmath>
using Eigen::MatrixXd;
using std::ifstream; 
using std::cout ; 
using std::endl; 
using std::istringstream; 
using std::vector; 
using std::string;
using std::cerr ; 

double toDouble(string s) {
  double r = 0;
  istringstream ss(s);
  ss >> r;
  return r;
}

MatrixXd ReadInPDBFileCA(string fname )
{
	vector< vector<double> > R ; 
	double x,y,z ; 
	string xs, ys, zs; 
	string line ; 
	ifstream fin; 
	fin.open( fname.c_str() ) ; 

	while(getline(fin,line))
	{ 
		if (line.size() > 30) 
		{
		if (line.substr(13,2) == "CA") 
		{
			x = toDouble( line.substr(30,8) )  ;
			y = toDouble( line.substr(38,8) ) ;
			z = toDouble(  line.substr(46,8) ); 
			vector<double> Ri ; 
			Ri.push_back(x/10.) ; 
			Ri.push_back(y/10.) ; 
			Ri.push_back(z/10.) ; 
			R.push_back(Ri) ; 

		}

		}
	
	}

	int len = R.size() ; 
	MatrixXd Rv(len, 3) ;
	for(int i = 0; i < len; i++) 
	{
		for(int j=0 ; j < 3; j++ ) 
		{
			Rv(i,j) = R[i][j] ; 
		
		}
	
	}
	
	return Rv ; 
}

MatrixXd GetDistances(MatrixXd R) 
{

	int s = R.rows(); 
	double dx, dy, dz; 
	MatrixXd D(s,s);  

	for(int i = 0; i < s ; i++) {
		for(int j = i+1 ; j < s ; j++ ) {
			dx = R(i,0) - R(j,0); 
			dy = R(i,1) - R(j,1);
			dz = R(i,2) - R(j,2);

			D(i,j) = sqrt(dx*dx + dy*dy + dz*dz) ;  
			D(j,i) = D(i,j) ; 
		}
		D(i,i) = 0.0 ; 
		cout << "got all distances for particle "<< i << endl; 
	}
	cout << "Done!" << endl;

	return D ; 
}

/*------------------------------------------
 * Hessian calculation 
 * implementation of method from 
 *	Zheng et al, Proteins 2010 78:2469–2481
 * -----------------------------------------*/
MatrixXd ComputeHessian(MatrixXd R, MatrixXd _Rij) 
{
	cout << "Calculating distances..."<< endl; 
	int numRes, ndf; 
	int i , j , ii, jj, ip, jp, qi, qj ; 
	double dij, _dij, ddij, dqi, dqj,dqij,h ; 
	MatrixXd H, Hij, Rij; 

	numRes= R.rows() ; 
	ndf = 3*numRes; 
	H = MatrixXd::Zero(ndf,ndf) ; 
	Rij = GetDistances(R) ; 

	cout << "computing Hessian" << endl; 
	for(i = 0 ; i < numRes; i++) {
		for(int j = 0; j < numRes; j++) {
			//cout <<"i "<<  i << "j " << j << endl; 

			_dij = _Rij(i,j) ; 
			if(_dij < CUTOFF) {
 				//Hij = MatrixXd::Zero(ndf,ndf); 
				dij = Rij(i,j) ;
				ddij = (dij - _dij)/ dij ; 

				for(ii = 0 ; ii <= 1 ; ii++) {
				for(jj = 0 ; jj <= 1; jj++) {
					if(ii == 0 ) ip = i ; 
					else ip = j ; 
					if(jj == 0 ) jp = i ; 
					else jp = j ; 

					for(qi = 0 ; qi < 3; qi++) {
					for(qj = 0 ; qj < 3; qj++) { 
						
						dqi = (R(ip,qi) - R(jp,qi))/dij; 
						dqj = (R(ip,qj) - R(jp,qi))/dij; 
						dqij = dqi * dqj ; 
						if (ip == jp){ 
							if (qi == qj) {
								h  = ddij + (1-ddij) * dqij; 
							} else { 
								 h = (1 - ddij) * dqij; 
							}
						} else {
							if(qi == qj) {
								h = - (ddij + (1-ddij) * dqij);
							} else {
								h  = - (1-ddij) * dqij ; 
							}
						}
						//Hij(3*ip+qi, 3*jp+qj) = h ; 
						cout << H.rows() << endl; 
						H(3*ip+qi, 3*jp+qj) = H(3*ip+qi, 3*jp+qj) +  h ; 
					}
					}
					
				}
				}
			}
			//H  = H + Hij ; 

		}
	}
	return H; 
}

/*---------------------------------------
 * Elastic Network Model base class [enm] 
 *---------------------------------------*/
class ENM {
	
	public:
	
		MatrixXd R, H , modes; 

		Eigen::EigenSolver<MatrixXd> es; 


		ENM(string fn, MatrixXd Rin, MatrixXd _Rij) {
			if (fn != "") {
				R = ReadInPDBFileCA(fn) ; 
				cout << "a" << endl; 
			} else if (Rin.size() > 0){ 
				R  = Rin ; 
				cout << "b" << endl; 
				//cout << R << endl;
			} else {
				cerr << "Error: " << endl; 
				cerr << "Need a file or coordinates as input to ENM" << endl; 
				cerr << "exiting... " << endl ; 
				exit(1) ; 
			}
			
			H = ComputeHessian(R,_Rij) ;
			//es.compute(H,true)  ; // Compute eigenvectors and eigenvals  
			//cout << es.eigenvectors() << endl;
		} 

		double ComputeOverlap(ENM other) { 
			double O ; 
			return O ; 
		}	

		double Similarity(ENM other) { 
			double S; 
			return S; 
		
		}


} ; 


int main() { 
	
	//ReadInPDBFileCA("/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_prep/init.pdb") ;
	MatrixXd D, R, _D, _R ; 
	R = ReadInPDBFileCA("/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_prep/init.pdb") ;
	_R = ReadInPDBFileCA("/Users/paulglen/Dropbox/CompBiophysics/cas9/inma_prep/em/em.pdb") ; 
	_D = GetDistances(_R ); 
	//cout << _D << endl; 
	ENM * enm = new ENM("" , R , _D);  

	return 0; 


}
