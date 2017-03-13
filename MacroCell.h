#ifndef _MACRO_CELL_H_
#define _MACRO_CELL_H_

#include <Vector.h>
#include <Connectivity.h>

//---------------------------------------------------------------------------

class MacroCell {

  int constructionProgress;

  int IDtag;
  int numSubCells;
  int* subCells;

  double volume;
  double gamma;
  double oogamma1;

public:

  int& operator[](int i) { return subCells[i]; }

  MacroCell(int, int, double);
  ~MacroCell();

  int whereis(int);
  double computeVolume(Vec<double> &);
  bool addCell(int);

  void computeDelRatios(Vec<double> &, double &, double &, double &, int &, int, SVec<double,5> &);

  void computeDelOverDel(SVec<double,5> &, int, int, double &);

  template<int dim>
  void computeBarValues(SVec<double,dim> &, Vec<double> &,
			SVec<double,dim> &, SVec<double,1> &);

  template<int dim>
  void computeVBar(SVec<double,dim> &, Vec<double> &, SVec<double,dim> &);

  int whoami() { return IDtag; }

  int size() { return numSubCells; }

  double getVolume() { return volume; }

};

//---------------------------------------------------------------------------

class MacroCellSet {

  int numMacroCells;
  int numTotalCells;
  int *cellsPerMacroCell;
  int *nodeToMacroCellMap;
  MacroCell **macroCells;

public:

  MacroCell &operator[](int i) const { return *(macroCells[i]); }

  MacroCellSet(int, int *, int, double);
  MacroCellSet(int nCells, Connectivity *cellToNode, int nNodes, double gam);
  ~MacroCellSet();

  void computeDelRatios(Vec<double> &, double &, double &, double &, int &, int, SVec<double,5> &);

  void computeDelOverDel(SVec<double,5> &, int, int, double &);

  template<int dim>
  void compute(bool, Vec<double> &, SVec<double,3> &,
	       SVec<double,dim> &, SVec<double,dim> &, SVec<double,1> &);

  template<int dim>
  void computeVBar(bool, Vec<double> &, SVec<double,dim> &, SVec<double,dim> &);

  int size() const { return numMacroCells; }

  int containing(int i) { return nodeToMacroCellMap[i]; }

  double getVolume(int i) { return macroCells[i]->getVolume(); }

};

//---------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <MacroCell.C>
#endif

#endif
