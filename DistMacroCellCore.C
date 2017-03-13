#include <DistMacroCell.h>

#include <MacroCell.h>
#include <Domain.h>
#include <DistVector.h>
#include <SubDomain.h>

//-----------------------------------------------------------------------

DistMacroCellSet::DistMacroCellSet(Domain* dom, double gamma,
				   bool** masterFlag, int scopeWidth,
				   int scopeDepth) : domain(dom)
{

  numLocSub  = domain->getNumLocSub();
  subDomain  = domain->getSubDomain();
  macroCells = new MacroCellSet** [numLocSub];


  #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
       macroCells[iSub] = subDomain[iSub]->findAgglomerateMesh(scopeWidth, scopeDepth, masterFlag[iSub], gamma);
    }
    
  #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) 
      subDomain[iSub]->createMacroCellConnectivities(macroCells[iSub], scopeDepth);

}

//-----------------------------------------------------------------------

DistMacroCellSet::~DistMacroCellSet()
{

  if (macroCells) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      if (macroCells[iSub]) delete macroCells[iSub];
    delete [] macroCells;
  }

}

//-----------------------------------------------------------------------

void DistMacroCellSet::computeDelRatios(DistVec<double> &ctrlVol, int scopeDepth, double *rmin, double *rmax, 
                                        double *rsum, double *ratiosum, int *numNodes)
                                                                                                                                                  
{
  
  double *rrmin, *rrmax, *rrsum, *rtsum, *rratiosum;    
  int *ssumcells;

  rrmin = new double[numLocSub];
  rrmax = new double[numLocSub];
  rrsum = new double[numLocSub];
  rtsum = new double[numLocSub];
  rratiosum = new double[numLocSub];
  ssumcells = new int[numLocSub];

  DistSVec<double,5> *Volume;
  Volume = new DistSVec<double,5>(domain->getNodeDistInfo());

  for (int i = 0; i < scopeDepth; ++i) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      (*macroCells[iSub][i]).computeDelRatios(ctrlVol(iSub), rrmin[iSub], rrmax[iSub], 
                                              rrsum[iSub], ssumcells[iSub], i, (*Volume)(iSub));
    }

    // finding the rmin and rmax of all the subdomains in the same cpu

    rmin[i] = 100000;
    rmax[i] = -100000;

    for (int iSub = 0; iSub < numLocSub; ++iSub) {
       if(rmin[i] > rrmin[iSub]) rmin[i] = rrmin[iSub];
       if(rmax[i] < rrmax[iSub]) rmax[i] = rrmax[iSub];
    }

    rsum[i] = 0.0; numNodes[i] = 0;

    // summing up values of ratios and the number of cells

    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      rsum[i] += rrsum[iSub];
      numNodes[i] += ssumcells[iSub];
    }
  }

  int k = 0;
  for (int i = 0; i < scopeDepth-1; ++i) {
   for (int j =i+1; j<scopeDepth; ++j) {
    ratiosum[k] = 0.0;
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      (*macroCells[iSub][i]).computeDelOverDel((*Volume)(iSub), i, j,  rratiosum[iSub]);
      ratiosum[k] += rratiosum[iSub];
    }
    k = k+1;
   } 
  }

  delete (Volume);
  delete [] rrmin; delete [] rrsum; delete [] rrmax;
  delete [] rtsum; delete [] rratiosum;
  delete [] ssumcells;

}

//-----------------------------------------------------------------------

