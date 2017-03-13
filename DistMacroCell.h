#ifndef _DIST_MACRO_CELL__H_
#define _DIST_MACRO_CELL__H_

class Domain;
class MacroCellSet;
class Domain;
class SubDomain;
class GeoState;
class DistGeoState;

#include <DistVector.h>


//-----------------------------------------------------------------------

class DistMacroCellSet {

  int numLocSub;
  Domain* domain;
  SubDomain** subDomain;
  MacroCellSet*** macroCells;
	 
public:

  DistMacroCellSet(Domain*, double, bool **, int, int);
  ~DistMacroCellSet();

  MacroCellSet* obtainMacroCell (int iSub, int scopeDepth) { return macroCells[iSub][scopeDepth];}

  void computeDelRatios(DistVec<double> &, int, double *, double *, double *, double *, int *);

  template<int dim>
  void computeDVMS(bool, DistVec<double> &, DistSVec<double,3> &,
	       DistSVec<double,dim> &, DistSVec<double,dim> **,
	       DistSVec<double,1> **, int, int);

  template<int dim>
  void computeVMS(bool, DistVec<double> &, DistSVec<double,3> &,
               DistSVec<double,dim> &, DistSVec<double,dim> &,
               DistSVec<double,1> &, int);

  template<int dim>
  void computeVBar(bool, DistGeoState &, DistSVec<double,dim> &, DistSVec<double,dim> &, int, int);

};

//-----------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistMacroCell.C>
#endif

#endif
