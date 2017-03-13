#include <MacroCell.h>

#include <cstdio>

//---------------------------------------------------------------------------

template<int dim>
void MacroCell::computeBarValues(SVec<double,dim> &V, 
				 Vec<double> &ctrlVol,
				 SVec<double,dim> &VtBar,
				 SVec<double,1> &volRatio)
{

  int i;

  // These will contain the large-scale quantities for the macro-cell
  double rhoBar = 0.0;
  double uBar   = 0.0;
  double vBar   = 0.0;
  double wBar   = 0.0;
  double enBar  = 0.0; 

  // Loop over all cells within the macro-cell
  for (i=0; i<numSubCells; ++i) {
    // Get the node number of the current cell
    int nodeNum = subCells[i];

    // Compute the volume ratio for the current cell
    volRatio[nodeNum][0] = ctrlVol[nodeNum] / volume;

    // These contain the instantaneous pointwise quantities for
    // the current cell contained within the macro-cell
    double nrho = V[nodeNum][0];                      // density
    double nu   = V[nodeNum][1];                      // x-velocity
    double nv   = V[nodeNum][2];                      // y-velocity
    double nw   = V[nodeNum][3];                      // z-velocity
    double ne   = (V[nodeNum][4] * oogamma1) / nrho;  // internal energy

    // Add the contributions to the large-scale quantities
    rhoBar += volRatio[nodeNum][0] * nrho;
    uBar   += volRatio[nodeNum][0] * nu;
    vBar   += volRatio[nodeNum][0] * nv;
    wBar   += volRatio[nodeNum][0] * nw;
    enBar  += volRatio[nodeNum][0] * ne;
  }

  // Now loop over all cells within the macro-cell again to assign
  // the large-scale components of V to each cell within
  for (i=0; i<numSubCells; ++i) {
    int nodeNum = subCells[i];

    VtBar[nodeNum][0] = rhoBar;
    VtBar[nodeNum][1] = uBar;
    VtBar[nodeNum][2] = vBar;
    VtBar[nodeNum][3] = wBar;
    VtBar[nodeNum][4] = enBar*rhoBar/oogamma1;

  }

}

//---------------------------------------------------------------------------

template<int dim>
void MacroCell::computeVBar(SVec<double,dim> &V,
                            Vec<double> &ctrlVol,
                            SVec<double,dim> &VBar)
{

  int i;

  // These will contain the large-scale quantities for the macro-cell
  double rhoBar = 0.0;
  double uBar   = 0.0;
  double vBar   = 0.0;
  double wBar   = 0.0;
  double enBar  = 0.0; 

  // Loop over all cells within the macro-cell
  for (i=0; i<numSubCells; ++i) {
    // Get the node number of the current cell
    int nodeNum = subCells[i];

    // Compute the volume ratio for the current cell
    double volRatio = ctrlVol[nodeNum] / volume;

    // These contain the instantaneous pointwise quantities for
    // the current cell contained within the macro-cell
    double nrho = V[nodeNum][0];                      // density
    double nu   = V[nodeNum][1];                      // x-velocity
    double nv   = V[nodeNum][2];                      // y-velocity
    double nw   = V[nodeNum][3];                      // z-velocity
    double ne   = (V[nodeNum][4] * oogamma1) / nrho;  // internal energy

    // Add the contributions to the large-scale quantities
    rhoBar += volRatio * nrho;
    uBar   += volRatio * nu;
    vBar   += volRatio * nv;
    wBar   += volRatio * nw;
    enBar  += volRatio * ne;
  }

  // Now loop over all cells within the macro-cell again to assign
  // the large-scale components of V to each cell within
  for (i=0; i<numSubCells; ++i) {
    int nodeNum = subCells[i];

    VBar[nodeNum][0]  = rhoBar;
    VBar[nodeNum][1]  = uBar;
    VBar[nodeNum][2]  = vBar;
    VBar[nodeNum][3]  = wBar;
    VBar[nodeNum][4]  = enBar*rhoBar/oogamma1;

  }

}

//---------------------------------------------------------------------------

template<int dim>
void MacroCellSet::compute(bool doInitialTasks,
			   Vec<double> &ctrlVol,
			   SVec<double,3> &X,
			   SVec<double,dim> &V,
			   SVec<double,dim> &VtBar,
			   SVec<double,1> &volRatio)
{

  VtBar    = 0.0;
  volRatio = 0.0;

  for (int k=0; k<numMacroCells; ++k) {
    // Only compute macro-cell volumes if the mesh has changed
    // or if performing a re-start
    if (doInitialTasks)
      macroCells[k]->computeVolume(ctrlVol);
    // Compute the large-scale components of V and the volume ratios
    macroCells[k]->computeBarValues(V, ctrlVol, VtBar, volRatio);
  }

}

//---------------------------------------------------------------------------

template<int dim>
void MacroCellSet::computeVBar(bool doInitialTasks,
                               Vec<double> &ctrlVol,
                               SVec<double,dim> &V,
                               SVec<double,dim> &VBar)
{

  VBar    = 0.0;

  for (int k=0; k<numMacroCells; ++k) {
    if (doInitialTasks)
      macroCells[k]->computeVolume(ctrlVol);
    macroCells[k]->computeVBar(V, ctrlVol, VBar);
  }

}

//---------------------------------------------------------------------------
