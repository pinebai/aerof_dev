#include <BCApplier.h>
#include <Domain.h>
#include <SubDomain.h>

#ifdef OLD_STL
#include <list.h>
#include <map.h>
#else
#include <list>
#include <map>
using std::list;
using std::map;
#endif

//#define HB_MESHMOTION_DEBUG

template<int dim>
void
BCApplier::applyP(DistSVec<double,dim> &X)
{

  SubDomain **subDomain = domain->getSubDomain();
  int iSub;
  
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    int *nodeMap = subDomain[iSub]->getNodeMap();
    SVec<double,dim>& Xsub = X(iSub);
    list<ProjData>::iterator it = proj[iSub].begin();
    while(it != proj[iSub].end()){
      it->apply(Xsub[it->node]);
      it++;
    }
  }

}

template<int dim>
void
BCApplier::applyD(DistSVec<double,dim> &X)
{

  SubDomain** subD = domain->getSubDomain();

//Alternative implementation to avoid adding a subdomain level method
//NEED TESTING ...
  if(dofType) {
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; ++iSub) {
      int nn = 0;
      double (*x)[dim] = X.subData(iSub);
      int (*subDofType)[dim] = reinterpret_cast<int (*)[dim]>(dofType[iSub]); 
      for(int i=0;i<X.subSize(iSub); i++)
        for(int l=0; l<dim; l++)
          //PJSA FIX: this projector should be for the constraints from the fluid-structure interface
          // and not the sliding planes which are dealt with in applyP
          if(subDofType[i][l]==BC_MATCHED || subDofType[i][l]==BC_MATCHEDSLIDE || subDofType[i][l]==BC_FIXED) 
          { 
            nn++; 
            x[i][l] = 0.0; 
          } 
          //if(subDofType[i][l]!=BC_FREE) { nn++; x[i][l] = 0.0; }
  
      //fprintf(stderr, "sub %d zeroed %d of %d dofs\n", subD[iSub]->getGlobSubNum(), nn, X.subSize(iSub)*dim);
    }
  }
}

template<int dim>
void
BCApplier::applyD2(DistSVec<double,dim> &X, double dX[dim])
{

  SubDomain** subD = domain->getSubDomain();

  if(dofType) {
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; ++iSub) {
      double (*x)[dim] = X.subData(iSub);
      int (*subDofType)[dim] = reinterpret_cast<int (*)[dim]>(dofType[iSub]); 
      for(int i=0;i<X.subSize(iSub); i++) {
        for(int l=0; l<dim; l++) {
          if(subDofType[i][l]==BC_FIXED) 
            x[i][l] = 0.0; 
	  else if (subDofType[i][l]==BC_MATCHEDSLIDE)
	    x[i][l] = dX[l];
	}
      }
    }
  }
}
