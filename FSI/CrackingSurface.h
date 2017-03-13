/*
 * CrackingSurface.h
 *
 *  Created on: Feb 2, 2011 (the night before the Year of Rabbit)
 *  Author: Kevin Wang
 */

#ifndef CRACKINGSURFACE_H_
#define CRACKINGSURFACE_H_
#include<LOCAL_LEVELSET.h>
#include<map>
#include<set>
#include<fstream>

//------------------------------------------------------------------------------

struct PhantomElement {
  int nNodes;
  double *phi;
  int *nodes;

  // constructors
  PhantomElement(): nNodes(-1), phi(0), nodes(0) {}
  PhantomElement(int n, int* nod, double* ph);
  PhantomElement(int a, int b, int c, int d, 
                 double phia, double phib, double phic, double phid);
  // destructor
  ~PhantomElement() {if(phi) delete[] phi;  if(nodes) delete[] nodes;}

  // update nodes
  void update(int* nod, double* ph); //KW: Doesn't have to update both nodes and phi. Use "NULL".

  void writeCrackingData(std::ofstream& restart_file) const ;
  static PhantomElement* readCrackingData(std::ifstream& restart_file) ;
};

//------------------------------------------------------------------------------

struct LatestCracking {
  std::set<int> phantomQuads;
  std::map<int,int> phantomNodes; //Note: "phantomNodes" are NOT equivalent to "nodes of phantomQuads"!
  LatestCracking() {/*nothing*/}
};

//------------------------------------------------------------------------------

class CrackingSurface : public LocalLevelSet {
  const int elemType; //currently only support quadrangles.
  int nTotalQuads, nUsedQuads;
  int nTotalTrias, nUsedTrias;
  int nTotalNodes, nUsedNodes;

  std::map<int,PhantomElement*> phantoms; //size: number of cracked (quad) elements
  LatestCracking latest;
  bool gotNewCracking;
  int (*tria2quad)[2]; //size: nTotalTrias
  int (*quad2tria)[2]; //size: nTotalQuads
  bool *cracked; //size: nTotalQuads
  bool *deleted; //size: nTotalQuads, in the case of Element Deletion, a "cracked" element is "deleted".

  /// For cracking simulations, this contains a map from a set of purely undeleted
  /// triangles to the full list of triangles.
  /// 
  int* triangle_id_map;
  int numRealTriangles;

  void constructTriangleMap();

public:
  CrackingSurface(int eType, int nUsed, int nTotal, int nUsedNd, int nTotNodes);
  ~CrackingSurface();

  //called by EmbeddedStructure only!
  int splitQuads(int* quadTopo, int nQuads, int(*triaTopo)[3]);
  int updateCracking(int numConnUpdate, int numLSUpdate, int* connUpdate, double* phi,
                     int* phiIndex, int(*triaTopo)[3], int nUsedNd, int* new2old, int numNewNodes);

  int numCrackedElements() {return phantoms.size();}
  bool hasCracked(int trId);
  double getPhi(int trId, double xi1, double xi2, bool* hasCracked=0, bool debug=false);
  
  //
  double getPhiPhysBAM(int trId, double xi1, double xi2, bool* hasCracked=0, bool debug=false);

  bool purelyPhantom(int trId);
  bool purelyPhantomPhysBAM(int trId);

  bool getNewCrackingFlag() const {return gotNewCracking;}
  void setNewCrackingFlag(bool flag) {gotNewCracking = flag;}

  int totNodes()  const {return nTotalNodes;}
  int usedNodes() const {return nUsedNodes;}
  int totTrias()  const {return nTotalTrias;}
  int usedTrias() const {return nUsedTrias;}
  std::set<int> getLatestPhantomQuads() const {return latest.phantomQuads;}
  std::map<int,int> getLatestPhantomNodes() const {return latest.phantomNodes;}
  void getQuad2Tria(int quad, int &trId1, int &trId2) {trId1=quad2tria[quad][0]; trId2=quad2tria[quad][1];}

  // for debug
  void printInfo(char* filename);

  void writeCrackingData(std::ofstream& restart_file) const;
  void readCrackingData(std::ifstream& restart_file);

  // see triangle_id_map
  int mapTriangleID(int);

  int numberRealTriangles() { 
    if (triangle_id_map)
      return numRealTriangles;
    else
      return nUsedTrias;
  }
};

//------------------------------------------------------------------------------
#endif /* CRACKINGSURFACE_H_ */

