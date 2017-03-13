/* ProgrammedBurn.C

*/

#include "ProgrammedBurn.h"
#include <limits>

ProgrammedBurn::ProgrammedBurn(IoData& ioData, DistSVec<double,3>* _nodeSet) : distInfo(&_nodeSet->info()) , nodeSet0(NULL) {

  nodeSet = _nodeSet;
  std::map<int, SphereData *>::iterator itr;
  std::map<int,int> burnableFluids;
  for (itr = ioData.mf.multiInitialConditions.sphereMap.dataMap.begin(); 
       itr != ioData.mf.multiInitialConditions.sphereMap.dataMap.end();
       ++itr) {
   
    SphereData& S = *(itr->second);
    ProgrammedBurnData& B = S.programmedBurn;
    if (B.unburnedEOS < 0)
      continue;
    
    double x0n[3] = {B.ignitionX0,
		     B.ignitionY0,
		     B.ignitionZ0 };
    
    Burn burn;
    burn.pgData = &B;
    burn.ignited = B.ignited;
    burn.finished = false;
    //burn.finished = B.finished;
    computeNearestNode(x0n, burn.x0,burn.x0subdom,burn.x0id);

    myBurns.push_back(burn);
  }

  std::map<int, PrismData *>::iterator itrp;
  for (itrp = ioData.mf.multiInitialConditions.prismMap.dataMap.begin(); 
       itrp != ioData.mf.multiInitialConditions.prismMap.dataMap.end();
       ++itrp) {
   
    PrismData& S = *(itrp->second);
    ProgrammedBurnData& B = S.programmedBurn;
    if (B.unburnedEOS < 0)
      continue;
    
    double x0n[3] = {B.ignitionX0,
		     B.ignitionY0,
		     B.ignitionZ0 };
    
    Burn burn;


    burn.pgData = &B;
    burn.ignited = B.ignited;
    burn.finished = false;
    //burn.finished = B.finished;
    //std::cout << "burn.ignited = " << burn.ignited << std::endl;
    computeNearestNode(x0n, burn.x0,burn.x0subdom,burn.x0id);

    myBurns.push_back(burn);
  }

  std::map<int, PointData *>::iterator it;
  for (it=ioData.embed.embedIC.pointMap.dataMap.begin();
       it!=ioData.embed.embedIC.pointMap.dataMap.end();
       it++) {

    PointData& S = *(it->second);
    ProgrammedBurnData& B = S.programmedBurn;
    if (B.unburnedEOS < 0)
      continue;

    double x0n[3] = {B.ignitionX0,
		     B.ignitionY0,
		     B.ignitionZ0 };

    Burn burn;
    burn.pgData = &B;
    burn.ignited = B.ignited;
    burn.finished = false;
    //burn.finished = B.finished;
    computeNearestNode(x0n, burn.x0,burn.x0subdom,burn.x0id);

    myBurns.push_back(burn);
  }

  lastTime = 0.0;
}

ProgrammedBurn::ProgrammedBurn(IoData& ioData, SVec<double,1>* _nodeSet) {
   
  nodeSet0 = _nodeSet;
  ProgrammedBurnData& B = ioData.oneDimensionalInfo.programmedBurn;
  if (B.unburnedEOS < 0)
    return;
  
  double x0n[3] = {B.ignitionX0,
		   B.ignitionY0,
		   B.ignitionZ0 };
  
  Burn burn;
  burn.pgData = &B;
  burn.ignited = B.ignited;
  burn.finished = false;
  double min_dist = 1000000000.0;
  int id;
  for (int i = 0; i < nodeSet0->size(); ++i) {
    double dist = fabs((*nodeSet0)[i][0]-x0n[0]);
    if (dist < min_dist) {
      min_dist = dist;
      id = i;
    }
  }
  
  burn.x0[0] = (*nodeSet0)[id][0];
  burn.x0[1] = 0.0;
  burn.x0[2] = 0.0;
  burn.x0id = id;

  //std::cout << "burn location = " << burn.x0[0] << " id = " << id << std::endl;
  
  myBurns.push_back(burn);

  lastTime = 0.0;
}

ProgrammedBurn::~ProgrammedBurn() {

}

int ProgrammedBurn::countBurnableFluids(IoData& ioData) {

  std::map<int, SphereData *>::iterator itr;
  std::map<int,int> burnableFluids;
  for (itr = ioData.mf.multiInitialConditions.sphereMap.dataMap.begin(); 
       itr != ioData.mf.multiInitialConditions.sphereMap.dataMap.end();
       ++itr) {
   
    SphereData& S = *(itr->second);
    if (S.programmedBurn.unburnedEOS < 0)
      continue;

    if (S.programmedBurn.unburnedEOS != S.fluidModelID) {
      fprintf(stderr,"Error: the programmedBurn EOS must match the sphere map EOS\n");
    }
    
    burnableFluids[ S.programmedBurn.unburnedEOS ] = S.programmedBurn.burnedEOS;

  }
  

  std::map<int, PrismData *>::iterator itrp;
  for (itrp = ioData.mf.multiInitialConditions.prismMap.dataMap.begin(); 
       itrp != ioData.mf.multiInitialConditions.prismMap.dataMap.end();
       ++itrp) {
   
    PrismData& S = *(itrp->second);
    if (S.programmedBurn.unburnedEOS < 0)
      continue;

    if (S.programmedBurn.unburnedEOS != S.fluidModelID) {
      fprintf(stderr,"Error: the programmedBurn EOS must match the prism map EOS\n");
    }
    
    burnableFluids[ S.programmedBurn.unburnedEOS ] = S.programmedBurn.burnedEOS;

  }

  
  map<int, PointData *>::iterator it;
  for (it=ioData.embed.embedIC.pointMap.dataMap.begin();
       it!=ioData.embed.embedIC.pointMap.dataMap.end();
       it++) {

    PointData& S = *(it->second);
    if (S.programmedBurn.unburnedEOS < 0)
      continue;

    if (S.programmedBurn.unburnedEOS != S.fluidModelID) {
      fprintf(stderr,"Error: the programmedBurn EOS must match the point map EOS\n");
    }

    burnableFluids[ S.programmedBurn.unburnedEOS ] = S.programmedBurn.burnedEOS;

  }

  //std::cout << "# burnable fluids = " << burnableFluids.size() << std::endl;
 
  return burnableFluids.size();
}

bool ProgrammedBurn::isDetonationInterface(int i,int j,int& tag) const {

  for (int k = 0; k < myBurns.size(); ++k) {
    const Burn& B = myBurns[k];
    if ((B.pgData->unburnedEOS == i && B.pgData->burnedEOS == j) ||
	(B.pgData->unburnedEOS == j && B.pgData->burnedEOS == i)) {
      tag = k;
      return true;
    }
  }

  return false;
}

void ProgrammedBurn::computeNearestNode(const double x0n[3], double x0[3], int& x0subdom, int& x0id) {

  double fmin = std::numeric_limits<double>::max();
  int mini = -1,mins=-1;

  int iSub;

  double *allmin = reinterpret_cast<double *>(alloca(sizeof(double) * distInfo->numLocSub));

#pragma omp parallel for
  for (iSub = 0; iSub < distInfo->numLocSub; ++iSub) {
    double (*x)[3] = nodeSet->subData(iSub);
    double dist = 0.0;
    
    for (int i = 0; i < nodeSet->subSize(iSub); ++i) {
      dist = (x[i][0]-x0n[0])*(x[i][0]-x0n[0])+
	(x[i][1]-x0n[1])*(x[i][1]-x0n[1])+
	(x[i][2]-x0n[2])*(x[i][2]-x0n[2]);
      if (dist < fmin) {
	fmin = dist;
	mini = i;
	mins = iSub;
      }
    }
  }

  double globalmin = fmin;

  distInfo->com->globalMin(1, &globalmin);

  int flag=-1,myflag = -1;
  if (fmin == globalmin) {
    myflag = flag = distInfo->locSubToGlobSub[mins];
  }

  distInfo->com->globalMax(1, &flag);

  if (flag == myflag) {

    // We are the winner.  Set the node to the initial state
    // Snap the initial explosion location to a node location
    x0[0] = x0n[0]; 
    x0[1] = x0n[1];
    x0[2] = x0n[2];
    //std::cout << "Actual ignition location: " << x0[0] << " " << x0[1] << " " << x0[2] << std::endl;
    x0subdom = mins;
    x0id = mini;
  } else {
    x0[0] = x0[1] = x0[2] = 0.0;
  }

  distInfo->com->globalSum(3,x0);
  
}

static double min(double a,double b) {
  if (a < b) return a;
  else return b;
}

void ProgrammedBurn::setFluidIds(double t, DistVec<int>& fluidIds,DistSVec<double,5>& U) {

  int iSub;
  lastTime = t;

  int cnt[5] = {0,0,0,0,0};

#pragma omp parallel for
  for (iSub = 0; iSub < distInfo->numLocSub; ++iSub) {
    int* fid = fluidIds.subData(iSub);
    double (*x)[3] = nodeSet->subData(iSub);
    double r;
    double V[5],v2;

    for (int i = 0; i < fluidIds.subSize(iSub); ++i) {
      
      for (int j = 0; j < myBurns.size(); ++j) {
	Burn& B = myBurns[j];
        if (B.finished)
          continue;

	r = sqrt((x[i][0]-B.x0[0])*(x[i][0]-B.x0[0])+
		 (x[i][1]-B.x0[1])*(x[i][1]-B.x0[1])+
		 (x[i][2]-B.x0[2])*(x[i][2]-B.x0[2]));

	if ((r <= B.pgData->cjDetonationVelocity*(t-B.pgData->ignitionTime) ||
	     (iSub == B.x0subdom && i == B.x0id)) &&
	    fid[i] == B.pgData->unburnedEOS && B.ignited) {
	  //bBurned->subData(iSub)[i] = true;
	  double* v = &(U.subData(iSub)[i][1]);
	  if (B.pgData->limitPeak) {
	    U.subData(iSub)[i][0] = min(U.subData(iSub)[i][0],  B.pgData->cjDensity);
	    U.subData(iSub)[i][4] = min(U.subData(iSub)[i][4],  
					B.pgData->cjDensity*B.pgData->cjEnergy+ 
					0.5*B.pgData->cjDensity*pow(B.pgData->cjDetonationVelocity/B.pgData->factorS,2.0));
	  }

	  //std::cout << "[ " << x[i][0] << " " << x[i][1] << " " << x[i][2] << " ]: ( " << U.subData(iSub)[i][0] <<
	  //  " , " << U.subData(iSub)[i][4] << " )" << std::endl;
	  
	  fid[i] = B.pgData->burnedEOS;
	}
        if (fid[i] == B.pgData->unburnedEOS)
          cnt[j]++;
      }
    }
  }
  
  distInfo->com->globalSum(5,cnt);

//  distInfo->com->fprintf(stderr,"cnt[0] = %d.\n", cnt[0]);

  for (int j = 0; j < myBurns.size(); ++j) {
    if (cnt[j] == 0)
      myBurns[j].finished = true;
  }
}

void ProgrammedBurn::setFluidIds(double t, Vec<int>& fluidIds,SVec<double,5>& U) {

  int iSub;
  double r;
  lastTime = t;
  int cnt[5] = {0,0,0,0,0};

  for (int i = 0; i < fluidIds.size(); ++i) {
    
    for (int j = 0; j < myBurns.size(); ++j) {
      Burn& B = myBurns[j];
      if (B.finished)
        continue;
      r = sqrt(((*nodeSet0)[i][0]-B.x0[0])*((*nodeSet0)[i][0]-B.x0[0])+
	       (B.x0[1])*(B.x0[1])+
	       (B.x0[2])*(B.x0[2]));
      
      if (r <= B.pgData->cjDetonationVelocity*(t-B.pgData->ignitionTime) &&
	  fluidIds[i] == B.pgData->unburnedEOS && B.ignited) {
	//bBurned->subData(iSub)[i] = true;
	/*U[i][0] = min(U[i][0],  B.pgData->cjDensity);
	U[i][4] = min(U[i][4],  
		      B.pgData->cjDensity*B.pgData->cjEnergy+ 
		      0.5*B.pgData->cjDensity*pow(B.pgData->cjDetonationVelocity/B.pgData->factorS,2.0));
	*/
//	std::cout << "[ x = " << (*nodeSet0)[i][0] << " ]: ( " << U[i][0] <<
//	    " , " << U[i][1] <<  " , " << U[i][4] << " )" << std::endl;
	fluidIds[i] = B.pgData->burnedEOS;
      }
      if (fluidIds[i] == B.pgData->unburnedEOS)
        cnt[j]++;
    }
  }
  
  for (int j = 0; j < myBurns.size(); ++j) {
    if (cnt[j] == 0)
      myBurns[j].finished = true;
  } 
}

void ProgrammedBurn::getDetonationNormal(int tag,int i,int j, double xmid[3], double gradphi[3]) {

  // Find the burn for fluid i/j
  Burn* B = &myBurns[tag];
  
  for (int k = 0; k < 3; ++k) {
    gradphi[k] = xmid[k] - B->x0[k];
  }
  
  double r = sqrt( gradphi[0]*gradphi[0] + gradphi[1]*gradphi[1] + gradphi[2]*gradphi[2] );
  for (int k = 0; k < 3; ++k) 
    gradphi[k] /= r;
}


bool ProgrammedBurn::nodeInside(int tag,int iSub, int i) {

  Burn* B = &myBurns[tag];
  double (*x)[3] = nodeSet->subData(iSub);
  double r = sqrt((x[i][0]-B->x0[0])*(x[i][0]-B->x0[0])+
	   (x[i][1]-B->x0[1])*(x[i][1]-B->x0[1])+
	   (x[i][2]-B->x0[2])*(x[i][2]-B->x0[2]));
 
  if ((r <= B->pgData->cjDetonationVelocity*(lastTime-B->pgData->ignitionTime) ||
      (iSub == B->x0subdom && i == B->x0id)) && B->ignited )
    return true;
  else {
    return false;
  }
}

bool ProgrammedBurn::nodeInside(int tag,int i) {

  //std::cout << "last time = " << lastTime << std::endl;
  SVec<double,1>& x = *nodeSet0;
  Burn* B = &myBurns[tag];
  double r = sqrt((x[i][0]-B->x0[0])*(x[i][0]-B->x0[0])+
	   (B->x0[1])*(B->x0[1])+
	   (B->x0[2])*(B->x0[2]));
  
  if ((r <= B->pgData->cjDetonationVelocity*(lastTime-B->pgData->ignitionTime) ||
      (i == B->x0id)) && B->ignited ) {
    return true;
  } else {
    return false;
  }
}

bool ProgrammedBurn::isBurnedEOS(int eos,int& tag) const {

  for (int i = 0; i < myBurns.size(); ++i) {

    if (myBurns[i].pgData->burnedEOS == eos) {
      tag = i;
      return true;
    }
  }
  return false;
}

bool ProgrammedBurn::isUnburnedEOS(int eos,int& tag) const {

  for (int i = 0; i < myBurns.size(); ++i) {

    if (myBurns[i].pgData->unburnedEOS == eos) {
      tag = i;
      return true;
    }
  }
  return false;
}

int ProgrammedBurn::getBurnedEOS(int tag) const {

  return myBurns[tag].pgData->burnedEOS;
}

int ProgrammedBurn::getUnburnedEOS(int tag) const {

  return myBurns[tag].pgData->unburnedEOS;
}

namespace ProgrammedBurn_CJ {

class IdealGasEOS {

public:

	IdealGasEOS(double _g) : gamma(_g) { }

	double computeSoundSpeed(double rho,double e) const {
		return sqrt((gamma)*(gamma-1.0)*e);
	}

	double computePressure(double rho,double e) const {
		return (gamma-1.0)*rho*e;
	}

	double computeDpDrho(double rho,double e) const {
		return (gamma-1.0)*e;
	}

	double computeDpDe(double rho,double e) const {
		return (gamma-1.0)*rho;
	}

	double computeDcDrho(double rho,double e) const {
		return 0.0;
	}

	double computeDcDe(double rho,double e) const {
		return 0.5/computeSoundSpeed(rho,e)*gamma*(gamma-1.0);
	}

private:

	double gamma;
};

class JWLEOS {

public:

	JWLEOS(double _A1,double _A2, double _R1, double _R2, double _omega) : A1(_A1), A2(_A2),
																			R1(_R1), R2(_R2),
																			omega(_omega)
	{ }

	double computeSoundSpeed(double rho,double e) const {
		return sqrt(1.0/rho*( (omega+1.0)*computePressure(rho,e)-F(rho)+rho*Fp(rho)) );
	}

	double computePressure(double rho,double e) const {
		return F(rho)+omega*rho*e;
	}

	double computeDpDrho(double rho,double e) const {
		return Fp(rho)+omega*e;
	}

	double computeDpDe(double rho,double e) const {
		return omega*rho;
	}

	double computeDcDrho(double rho,double e) const {
		double c = computeSoundSpeed(rho,e);
		double q = c*c*rho;
		return 1.0/(2.0*c)*(-1/(rho*rho)*q + 1.0/rho*( (omega+1.0)*computeDpDrho(rho,e)+rho*Fpp(rho) )) ;
	}

	double computeDcDe(double rho,double e) const {
		double c = computeSoundSpeed(rho,e);
		return 1.0/(2.0*c)*(1.0/rho*(omega+1.0)*computeDpDe(rho,e) );
	}

	double F(double rho) const {
		return A1*(1.0-omega*rho/R1)*exp(-R1/rho)+A2*(1.0-omega*rho/R2)*exp(-R2/rho);
	}
	
	double Fp(double rho) const {
		return A1*(-omega/R1+(1.0-omega*rho/R1)*R1/(rho*rho))*exp(-R1/rho)+
			   A2*(-omega/R2+(1.0-omega*rho/R2)*R2/(rho*rho))*exp(-R2/rho);
	}

	double Fpp(double rho) const {
		return A1*( (-omega/R1+(1.0-omega*rho/R1)*R1/(rho*rho))*(R1/(rho*rho))  + (-omega/R1)*R1/(rho*rho) -2.0*(1.0-omega*rho/R1)*R1/(rho*rho*rho) )  *exp(-R1/rho)+
			   A2*( (-omega/R2+(1.0-omega*rho/R2)*R2/(rho*rho))*(R2/(rho*rho))  + (-omega/R2)*R2/(rho*rho) -2.0*(1.0-omega*rho/R2)*R2/(rho*rho*rho) )  *exp(-R2/rho);
	}
private:

	double A1,A2,R1,R2,omega;
};

template <class EOS>
void evaluateF(double rho_0,double e_0,double p_0, double rho_cj, double e_cj, double &f_1, double &f_2, const EOS& theEOS) {

	double c = theEOS.computeSoundSpeed(rho_cj,e_cj);
	double p_cj = theEOS.computePressure(rho_cj,e_cj);
	//s = rho_cj/rho_0*theEOS.computeSoundSpeed(rho_cj, e_cj);

	f_1 = p_0 + c*c*rho_cj*((rho_cj/rho_0)-p_cj/(c*c*rho_cj)-1.0);
	f_2 = (e_0-e_cj) +c*c*(p_0/(rho_0*c*c)+0.5*rho_cj*rho_cj/(rho_0*rho_0)-p_cj/(rho_cj*c*c)-0.5);
}

template <class EOS>
void computeChapmanJouguetState(double p_0,double rho_0,double e_0,
								double &rho_cj, double& e_cj,double&s,
								double &p_cj,const EOS& theEOS) {

	while (1) {
		double c = theEOS.computeSoundSpeed(rho_cj,e_cj);
		p_cj = theEOS.computePressure(rho_cj,e_cj);
		s = rho_cj/rho_0*theEOS.computeSoundSpeed(rho_cj, e_cj);

		double f_1 = p_0 + c*c*rho_cj*((rho_cj/rho_0)-p_cj/(c*c*rho_cj)-1.0);
		double f_2 = (e_0-e_cj) +c*c*(p_0/(rho_0*c*c)+0.5*rho_cj*rho_cj/(rho_0*rho_0)-p_cj/(rho_cj*c*c)-0.5);

		if (fabs(f_1) < 1.0 && fabs(f_2) < 1.0)
			break;

		//std::cout << "Newton res = " << f_1 << " " << f_2 << " rho_cj = " << rho_cj << " e_cj = " << e_cj << std::endl;
		//std::cout << "p = " << p_cj << " s = " << s << " c = " << c << std::endl;

		double dpdrho = theEOS.computeDpDrho(rho_cj,e_cj),
			dpde = theEOS.computeDpDe(rho_cj,e_cj),
			dcdrho = theEOS.computeDcDrho(rho_cj,e_cj),
			dcde = theEOS.computeDcDe(rho_cj,e_cj);

		double J[4] = { (1.0/rho_0*(2.0*rho_cj+2.0*rho_cj*rho_cj*dcdrho/c)-dpdrho/(c*c)-1.0-2.0*rho_cj*dcdrho/c),
						2.0/rho_0*rho_cj*rho_cj*dcde/c-dpde/(c*c)-2.0*rho_cj*dcde/c,

						1.0/(rho_0*rho_0)*(rho_cj+rho_cj*rho_cj*dcdrho/c)-1.0/rho_cj*dpdrho/(c*c)+p_cj/(rho_cj*rho_cj*c*c)-dcdrho/c,
		(1.0/(rho_0*rho_0)*rho_cj*rho_cj*dcde/c-1.0/rho_cj*dpde/(c*c)-dcde/c)-1.0/(c*c) };

		double fp1=0.0,fp2=0.0;
		//evaluateF(rho_0,e_0,p_0, rho_cj+1.0,e_cj, fp1, fp2, theEOS);
		
		//evaluateF(rho_0,e_0,p_0, rho_cj,e_cj+100.0, fp1, fp2, theEOS);

		double det = J[0]*J[3]-J[1]*J[2];
		double drho = 1.0/det*(J[3]*f_1-J[1]*f_2);
		double de = 1.0/det*(-J[2]*f_1+J[0]*f_2);
		rho_cj -= drho/(c*c);
		e_cj -= de/(c*c);
	}

	s = rho_cj/rho_0*theEOS.computeSoundSpeed(rho_cj, e_cj);
}

}

void ProgrammedBurn::computeChapmanJouguetStateJWL(double A1,double A2,double R1,double R2,double omega, // JWL parameters
						   double p0, double rho0, double e0,
						   double& rho_cj,double& p_cj, double& e_cj, double& s) {

  rho_cj=rho0*1.5,p_cj = 0.0, e_cj = e0*1.5,s=0.0;
  ProgrammedBurn_CJ::JWLEOS jwl_eos(A1,A2,R1,R2,omega);
  ProgrammedBurn_CJ::computeChapmanJouguetState( p0,rho0,e0, rho_cj,e_cj,s,p_cj, jwl_eos);

}

void ProgrammedBurn::computeChapmanJouguetStatePG(double gamma, // PG parameters
						  double p0, double rho0, double e0,
						  double& rho_cj,double& p_cj, double& e_cj, double& s) {

  rho_cj=rho0*1.5,p_cj = 0.0, e_cj = e0*1.5,s=0.0;
  ProgrammedBurn_CJ::IdealGasEOS pg_eos(gamma);
  ProgrammedBurn_CJ::computeChapmanJouguetState( p0,rho0,e0, rho_cj,e_cj,s,p_cj, pg_eos);

}
