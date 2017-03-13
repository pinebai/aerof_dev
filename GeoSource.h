#ifndef _GEO_SOURCE_H_
#define _GEO_SOURCE_H_

class IoData;
class MatchNodeSet;
class SubDomain;
class Connectivity;
class Communicator;
class BinFileHandler;
class BCondSet;
class ElemSet;
class FaceSet;

//------------------------------------------------------------------------------

class GeoSource {

  char *conName;
  char *geoName;
  char *decName;
  char *mapName;
  char *matchName;

  int numGlobSub, numCPU, numClusters;
  int numLocSub, cpuNum, numLocThreads;

  int (*subToCluster)[2];
  Connectivity *subToSub;
  Connectivity *cpuToSub;
  int *cpuToThreads;

  MatchNodeSet **matchNodes;

  double oolscale;
  IoData *iod;

  BCondSet *clus_bconds;

public:

  GeoSource(IoData &);
  ~GeoSource();

  void readConnectivityInfo(Communicator *);
  void readCpuToSub();
  int getRangeInfo(BinFileHandler &, int &, int (*&)[2]);
  SubDomain *getSubDomain(int);

  int getNumGlobSub() const { return numGlobSub; }
  int getNumLocSub() const { return numLocSub; }
  int getNumLocThreads() const { return numLocThreads; }
  Connectivity *getSubToSub() const { return subToSub; }
  Connectivity *getCpuToSub() const { return cpuToSub; }
  MatchNodeSet **getMatchNodes() const { return matchNodes; }
  template<class MapType> 
    void distributeBCs(SubDomain *sub, MapType &cl2LocNodeMap);
  template<class MapType>
    void getBC(BCondSet *&subBC, MapType &cl2LocNodeMap);

};

//------------------------------------------------------------------------------

#endif
