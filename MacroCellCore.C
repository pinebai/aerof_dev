#include <MacroCell.h>

#include <cstdio>

//---------------------------------------------------------------------------

MacroCell::MacroCell(int val, int macroCellID, double gam)
{

  IDtag                = macroCellID;
  gamma                = gam;
  oogamma1             = 1.0/(gamma - 1.0);
  constructionProgress = 0;
  volume               = 0.0;
  numSubCells          = val;
  subCells             = new int [val];

  for (int k=0; k<val; ++k) {
    subCells[k] = -1;
  }
}

//---------------------------------------------------------------------------

MacroCell::~MacroCell()
{

  if (subCells) delete [] subCells;

}

//---------------------------------------------------------------------------

int MacroCell::whereis(int nodeNum)
{

  for (int i=0; i<numSubCells; ++i)
    if (subCells[i]==nodeNum)
      return i;

  return -1;

}

//---------------------------------------------------------------------------

double MacroCell::computeVolume(Vec<double> &ctrlVol)
{

  double vol = 0.0;

  for (int i=0; i<numSubCells; ++i)
    vol += ctrlVol[ subCells[i] ];

  volume = vol;
  return vol;

}

//---------------------------------------------------------------------------

bool MacroCell::addCell(int subCellID)
{

  if (constructionProgress < numSubCells) {
    subCells[constructionProgress] = subCellID;
    constructionProgress++;
    return true;
  }
  else
    return false;

}

//---------------------------------------------------------------------------

MacroCellSet::MacroCellSet(int numCells,
			   int *cellMap,
			   int numNodes,
			   double gam)

{

  numMacroCells      = numCells;
  numTotalCells      = numNodes;
  cellsPerMacroCell  = new int [numMacroCells];
  nodeToMacroCellMap = new int [numNodes];
  int m, k;

  // Initialize cellsPerMacroCell
  for (k=0; k < numMacroCells; ++k)
    cellsPerMacroCell[k] = 0;

  // Make local copy of node-to-macro-cell map
  // and count number of cells in each macro-cell
  for (m=0; m<numNodes; ++m) {
    nodeToMacroCellMap[m] = cellMap[m];

    // A cell only counts if it is included
    if (nodeToMacroCellMap[m] != -1)
      cellsPerMacroCell[ nodeToMacroCellMap[m] ]++;
  }

  // Create a new array of pointers to MacroCell,
  // and create a new MacroCell in each element of the array
  macroCells = new MacroCell *[numMacroCells];
  for (k=0; k<numMacroCells; ++k) {
    macroCells[k] = new MacroCell(cellsPerMacroCell[k], k, gam);
  }

  bool stat;

  // Fill macro-cells
  for (m=0; m<numNodes; ++m) {
    if (nodeToMacroCellMap[m] != -1) {
      stat = macroCells[ nodeToMacroCellMap[m] ]->addCell(m);
      if (!stat)
	printf("Warning: Could not add cell defined by node %d to macrocell %d\n",
	       m,nodeToMacroCellMap[m]);
    }
  }

}

//---------------------------------------------------------------------------

MacroCellSet::MacroCellSet(int nCells, Connectivity *cellToNode, int nNodes, double gam)
{
  numMacroCells = nCells;
  numTotalCells = nNodes;
  cellsPerMacroCell  = new int [numMacroCells];
  nodeToMacroCellMap = new int [nNodes];
  int m, k;

  for(m = 0 ; m < nNodes; ++m)
    nodeToMacroCellMap[m] = -1;
    
  // Create a new array of pointers to MacroCell,
  // and create a new MacroCell in each element of the array
  macroCells = new MacroCell *[numMacroCells];
  for (k=0; k < numMacroCells; ++k) {
    cellsPerMacroCell[k] = cellToNode->num(k);
    macroCells[k] = new MacroCell(cellsPerMacroCell[k], k, gam);
    for(m = 0; m < cellToNode->num(k); ++m) {
       nodeToMacroCellMap[(*cellToNode)[k][m]] = k;
       bool stat = macroCells[k]->addCell((*cellToNode)[k][m]);
       if (!stat)
	printf("Warning: Could not add cell defined by node %d to macrocell %d\n",
	       m,nodeToMacroCellMap[m]);
    }
  }

}
//---------------------------------------------------------------------------

MacroCellSet::~MacroCellSet()
{

  if (macroCells) {
    for (int k=0; k<numMacroCells; ++k)
      if (macroCells[k]) delete macroCells[k];
    delete [] macroCells;
  }

  if (cellsPerMacroCell) delete [] cellsPerMacroCell;
  if (nodeToMacroCellMap) delete [] nodeToMacroCellMap;

}

//---------------------------------------------------------------------------
                                                                                                                                                  
void MacroCellSet::computeDelRatios(Vec<double> &ctrlVol, double &rmin, double &rmax, double &rsum, int &sumcells, 
                                    int scd, SVec<double,5> &Volume)
{

  double *rrmax, *rrmin, *rrsum;
  int *ssumcells;

  rrmax = new double[numMacroCells];
  rrmin = new double[numMacroCells];
  rrsum = new double[numMacroCells];
  ssumcells = new int[numMacroCells];

  rmin = 10000; rmax = -10000; rsum = 0.0; sumcells = 0;

  for (int k=0; k<numMacroCells; ++k) {
    macroCells[k]->computeVolume(ctrlVol);
    macroCells[k]->computeDelRatios(ctrlVol, rrmax[k], rrmin[k], rrsum[k], ssumcells[k], scd, Volume);
    rsum = rsum + rrsum[k]; sumcells = sumcells + ssumcells[k];
  }

  for (int k=0; k<numMacroCells; ++k) {
    if(rrmin[k]<rmin) rmin = rrmin[k];
    if(rrmax[k]>rmax) rmax = rrmax[k];
  }

  delete rrmax; delete rrmin; delete rrsum;
  delete ssumcells;

}

//---------------------------------------------------------------------------

void MacroCell::computeDelRatios(Vec<double> &ctrlVol, double &rmax, double &rmin, double &rsum, int &sumcells, 
                                 int scd, SVec<double,5> &Volume)
{

  int i;
  double *value;

  value = new double[numSubCells];

  rsum = 0.0; rmin = 10000; rmax = -10000;
  sumcells = numSubCells;

  // Loop over all cells within the macro-cell
  for (i=0; i<numSubCells; ++i) {
    int nodeNum = subCells[i];
    value[i] = (ctrlVol[nodeNum]*100) / volume; // percentage value
    rsum = rsum + value[i];
    Volume[nodeNum][scd] = value[i];
  }

  for (i=0; i<numSubCells; ++i) {
    if (value[i]<rmin) rmin = value[i];
    if (value[i]>rmax) rmax = value[i];
  }

  delete value;

}

//---------------------------------------------------------------------------

void MacroCell::computeDelOverDel(SVec<double,5> &Volume, int scd1, int scd2, double &rsum)

{

  int i;
  double *value;

  value = new double[numSubCells];

  rsum = 0.0; 

  // Loop over all cells within the macro-cell
  for (i=0; i<numSubCells; ++i) {
    int nodeNum = subCells[i];
    rsum = rsum + (Volume[nodeNum][scd2] / Volume[nodeNum][scd1])*100;
  }

  delete value;

}

//---------------------------------------------------------------------------

void MacroCellSet::computeDelOverDel(SVec<double,5> &Volume, int scd1, int scd2, double &rsum)
                                    
{

  double *rrsum;
  rrsum = new double[numMacroCells];

  rsum = 0.0;

  for (int k=0; k<numMacroCells; ++k) {
    macroCells[k]->computeDelOverDel(Volume, scd1, scd2, rrsum[k]);
    rsum = rsum + rrsum[k]; 
  }

  delete rrsum;

}

//---------------------------------------------------------------------------
