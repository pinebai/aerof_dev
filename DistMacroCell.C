#include <DistMacroCell.h>

#include <MacroCell.h>
#include <DistVector.h>
#include <DistGeoState.h>
#include <GeoState.h>


//-----------------------------------------------------------------------

template<int dim>
void DistMacroCellSet::computeDVMS(bool doInitialTasks,
			           DistVec<double> &ctrlVol,
			           DistSVec<double,3> &X,
			           DistSVec<double,dim> &V,
			           DistSVec<double,dim> **VtBar,
			           DistSVec<double,1> **volRatio, 
			           int scopeDepth1, int scopeDepth2)

{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
     (*macroCells[iSub][scopeDepth1-1]).compute(doInitialTasks, ctrlVol(iSub),
                                        X(iSub), V(iSub), (*VtBar[0])(iSub), 
				        (*volRatio[0])(iSub));
     (*macroCells[iSub][scopeDepth2-1]).compute(doInitialTasks, ctrlVol(iSub),
                                        X(iSub), V(iSub), (*VtBar[1])(iSub), 
				        (*volRatio[1])(iSub));
  }

}


//-----------------------------------------------------------------------

template<int dim>
void DistMacroCellSet::computeVMS(bool doInitialTasks,
                               DistVec<double> &ctrlVol,
                               DistSVec<double,3> &X,
                               DistSVec<double,dim> &V,
                               DistSVec<double,dim> &VtBar,
                               DistSVec<double,1> &volRatio, 
			       int scopeDepth)

{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
        (*macroCells[iSub][scopeDepth-1]).compute(doInitialTasks, ctrlVol(iSub),
                                X(iSub), V(iSub), VtBar(iSub), volRatio(iSub));
  }

}

//-----------------------------------------------------------------------

template<int dim>
void DistMacroCellSet::computeVBar(bool doInitialTasks, 
                                   DistGeoState &geoState,
                                   DistSVec<double,dim> &V, 
                                   DistSVec<double,dim> &VBar, 
                                   int scopeDepth, int n)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
     Vec<double>& ctrlVol = geoState(iSub).getCtrlVol_n();
     if (n == 2) ctrlVol = geoState(iSub).getCtrlVol_nm1();
     if (n == 3) ctrlVol = geoState(iSub).getCtrlVol_nm2();
     (*macroCells[iSub][scopeDepth-1]).computeVBar(doInitialTasks, ctrlVol, V(iSub), VBar(iSub));
  }
  
}

//-----------------------------------------------------------------------

