#ifndef _MATCH_NODE_H_
#define _MATCH_NODE_H_

class BinFileHandler;

//------------------------------------------------------------------------------

class MatchNodeSet {

  int numNodes, totalSize;

  int (*index)[3];
  double (*gap)[3];
  double (*xi)[3];
  double *normgap;

public:

  MatchNodeSet(const char *fileName);
  MatchNodeSet() { numNodes = 0; totalSize = 0; index = 0; gap = 0; xi = 0; normgap = 0; }
  MatchNodeSet(int);
  ~MatchNodeSet();

  void read(BinFileHandler &, int, int (*)[2]);
  void autoInit(int);
  void updateNumStNodes(int nn) {numNodes = nn;}

  template<class NodeMap>
  void renumberNodes(NodeMap &);

  void exportInfo(int, int (*)[3]);
  void setBufferPosition(int, int);
  void setBufferPosition(int, int, double[2], int (*)[3], double (*)[3]);
  void getDisplacement(int, double, double, double, bool *, double (*)[2][3], double (*)[3], 
		       double (*)[3], double (*)[3], double (*)[3], double *,
		       bool isEmbedded); 
  void getDisplacementSensitivity(bool *, double, double (*)[2][3], double (*)[3], double *);
                       
  void getDisplacement(double (*)[3], int (*)[3], double (*)[3], double (*)[3], double (*)[3], double);

  double getTemperature(int, double, double, bool*, double*, double*);

  template<int dim>
  void send(double, double (*)[dim], double (*)[dim]);

  template<int dim>
  void sendWithMasterFlag(double, double (*)[dim], double (*)[dim], bool *, int);

  double (*getGap(int, int *))[3];

  int size() const { return numNodes; }
  int totSize() const { return totalSize; }
  int subMatchNode(int i) { return(index[i][0]); } //HB: returns the subdomain node of the ith matched node 	
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <MatchNode.C>
#endif

#endif
