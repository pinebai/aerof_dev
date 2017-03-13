#include <cstdio>

#ifdef OLD_STL
#include <map.h>
#else
#include <map>
using std::map;
#endif

#include <GeoSource.h>

#include <IoData.h>
#include <MatchNode.h>
#include <SubDomain.h>
#include <Connectivity.h>
#include <Communicator.h>
#include <BinFileHandler.h>
#include <BCond.h>

//------------------------------------------------------------------------------

GeoSource::GeoSource(IoData &ioData)
{

  iod = &ioData;

  int sp = strlen(ioData.input.prefix) + 1;

  conName = new char[sp + strlen(ioData.input.connectivity)];
  sprintf(conName, "%s%s", ioData.input.prefix, ioData.input.connectivity);

  geoName = new char[sp + strlen(ioData.input.geometry)];
  sprintf(geoName, "%s%s", ioData.input.prefix, ioData.input.geometry);

  decName = new char[sp + strlen(ioData.input.decomposition)];
  sprintf(decName, "%s%s", ioData.input.prefix, ioData.input.decomposition);

  mapName = new char[sp + strlen(ioData.input.cpumap)];
  sprintf(mapName, "%s%s", ioData.input.prefix, ioData.input.cpumap);

  if (ioData.problem.framework==ProblemData::BODYFITTED && //KW: For Embedded, matcher file is in ASCII
      (ioData.problem.type[ProblemData::AERO] ||        //      and it's loaded in FSI/Dynam... 
       ioData.problem.type[ProblemData::THERMO])) {
    matchName = new char[sp + strlen(ioData.input.match)];
    sprintf(matchName, "%s%s", ioData.input.prefix, ioData.input.match);
  }
  else if (ioData.problem.framework==ProblemData::EMBEDDEDALE && ioData.dmesh.type == DefoMeshMotionData::COROTATIONAL && strcmp(ioData.input.embmeshmatch,"") != 0) {
    matchName = new char[sp + strlen(ioData.input.embmeshmatch)];
    sprintf(matchName, "%s%s", ioData.input.prefix, ioData.input.embmeshmatch);
  }
  else {
    matchName = new char[1];
    sprintf(matchName, "");
  }

  subToCluster = 0;
  subToSub = 0;
  cpuToSub = 0;
  cpuToThreads = 0;

  matchNodes = 0;

  oolscale = 1.0 / ioData.ref.rv.tlength;

  clus_bconds = 0;
}

//------------------------------------------------------------------------------

GeoSource::~GeoSource()
{

  if (conName) delete [] conName;
  if (geoName) delete [] geoName;
  if (decName) delete [] decName;
  if (mapName) delete [] mapName;
  if (matchName) delete [] matchName;

  if (subToCluster) delete [] subToCluster;
  if (subToSub) delete subToSub;
  if (cpuToSub) delete cpuToSub;
  if (cpuToThreads) delete [] cpuToThreads;
  if (matchNodes) {
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub)
      if (matchNodes[iSub]) delete matchNodes[iSub];
    delete [] matchNodes;
  }
  if(clus_bconds) delete clus_bconds;
}

//------------------------------------------------------------------------------

void GeoSource::readConnectivityInfo(Communicator *com)
{

  numCPU = com->size();
  cpuNum = com->cpuNum();

  BinFileHandler conFile(conName, "rb");

  Connectivity clusterToSub(conFile);

  subToSub = new Connectivity(conFile);

  numGlobSub = clusterToSub.numConnect();

  subToCluster = new int[numGlobSub][2];

  numClusters = clusterToSub.csize();

  for (int i=0; i<numClusters; ++i) {
    for (int iSub=0; iSub<clusterToSub.num(i); ++iSub) {
      subToCluster[ clusterToSub[i][iSub] ][0] = i;
      subToCluster[ clusterToSub[i][iSub] ][1] = iSub;
    }
  }

  if (cpuNum == -1) {
    fprintf(stdout, "Global decomposition with %d subdomain%s\n", 
	    numGlobSub, numGlobSub>1? "s":"");
    fprintf(stdout, "Cluster to subdomain:\n");
    clusterToSub.print(stdout);
    fprintf(stdout, "Subdomain to subdomain:\n");
    subToSub->print(stdout);
  }

  readCpuToSub();

  numLocSub = cpuToSub->num(cpuNum);
  numLocThreads = cpuToThreads[cpuNum];

  if (matchName[0] != 0){
    matchNodes = new MatchNodeSet *[numLocSub];
    for (int iSub=0; iSub<numLocSub; ++iSub) matchNodes[iSub] = 0;
  }
}

//------------------------------------------------------------------------------

void GeoSource::readCpuToSub()
{

  int *ptr = new int[numCPU+1];
  int *target = new int[numGlobSub];
  cpuToThreads = new int[numCPU];

#if defined(_OPENMP) || defined(USE_MPI)
  FILE *fp = fopen(mapName, "r");
  if (!fp) {
    fprintf(stderr, "*** Error: could not open CPU Map file \'%s\'\n", mapName);
    exit(1);
  }

  char line[MAXLINE];
  char* toto = fgets(line, MAXLINE, fp);
  int nCPU;
  sscanf(line, "%d", &nCPU);

  if (nCPU != numCPU) {
    fprintf(stderr, "*** Error: incorrect number of CPUs (%d) in the map file (It should be %d)\n", nCPU, numCPU);
    exit(1);
  }

  int curSub = 0;
  for (int iCPU=0; iCPU<numCPU; ++iCPU) {

    int numSub;
    toto = fgets(line, MAXLINE, fp);
    sscanf(line, "%d %d", &numSub, &cpuToThreads[iCPU]);
    ptr[iCPU] = curSub;

    if (curSub + numSub > numGlobSub) {
      if (cpuNum == 0)
	fprintf(stderr, "*** Error: this decomposition contains more elements"
		" than the original mesh\n");
      exit(1);
    }

    for (int iSub=0; iSub<numSub; ++iSub) {
      toto = fgets(line, MAXLINE, fp);
      sscanf(line, "%d", target+curSub);
      target[curSub] -= 1;
      curSub++;
    }

  }

  ptr[numCPU] = curSub;
#else
  ptr[0] = 0;
  ptr[1] = numGlobSub;
  for (int iSub=0; iSub<numGlobSub; ++iSub)
    target[iSub] = iSub;
  cpuToThreads[0] = 1;
#endif

  cpuToSub = new Connectivity(numCPU, ptr, target);

}

//------------------------------------------------------------------------------

int GeoSource::getRangeInfo(BinFileHandler &file, int &numRanges, int (*&ranges)[2])
{

  // read the numbe of ranges and singles

  int nr, ns;
  file.read(&nr, 1);
  file.read(&ns, 1);

  numRanges = nr + ns;
  ranges = new int[numRanges][2];

  // read the ranges for this subdomain
  
  file.read(reinterpret_cast<int *>(ranges), 2*nr);

  // read the singles for this subdomain

  int *s = new int[ns];
  file.read(s, ns);

  // merge ranges and singles

  int i;
  for (i=0; i<ns; ++i) {
    ranges[nr + i][0] = s[i];
    ranges[nr + i][1] = s[i];
  }

  delete [] s;

  // compute the number of values for this subdomain
  int numValues = 0;
  for (i=0; i<numRanges; ++i)
    numValues += ranges[i][1] - ranges[i][0] + 1;

  return numValues;

}

//------------------------------------------------------------------------------

SubDomain *GeoSource::getSubDomain(int iSub)
{

  int globSubNum = (*cpuToSub)[cpuNum][iSub];
  int clusterNum = subToCluster[globSubNum][0];
  int clusterSubNum = subToCluster[globSubNum][1];

  char fullgeoName[MAXLINE], fulldecName[MAXLINE], fullmatchName[MAXLINE];
  char *suffix = computeClusterSuffix(clusterNum+1, numClusters);
  sprintf(fullgeoName, "%s%s", geoName, suffix);
  sprintf(fulldecName, "%s%s", decName, suffix);
  sprintf(fullmatchName, "%s%s", matchName, suffix);

  // get the decomposition

  BinFileHandler decFile(fulldecName, "rb");

  // get number of subdomains in cluster (not used)
  int numSubCluster;
  decFile.read(&numSubCluster, 1);

  // get starting location of the subdomain
  decFile.seek(sizeof(int) + clusterSubNum*sizeof(BinFileHandler::OffType));
  BinFileHandler::OffType offSet;
  decFile.read(&offSet, 1);
  decFile.seek(offSet);

  // get node ranges
  int numNodeRanges, (*nodeRanges)[2] = 0;
  int numNodes = getRangeInfo(decFile, numNodeRanges, nodeRanges);

  // get elem ranges
  int numElemRanges, (*elemRanges)[2] = 0;
  int numElems = getRangeInfo(decFile, numElemRanges, elemRanges);
  
  // get face ranges
  int numFaceRanges, (*faceRanges)[2] = 0;
  int numFaces = getRangeInfo(decFile, numFaceRanges, faceRanges);

  // get the neighboring subs and their interface nodes
  int numNeighb;
  decFile.read(&numNeighb, 1);
  int *neighb;
  Connectivity *sharedNodes;
  if (numNeighb > 0) {
    neighb = new int[numNeighb];
    decFile.read(neighb, numNeighb);
    sharedNodes = new Connectivity(decFile);
  } 
  else {
    neighb = 0;
    int *ptr = new int[2];
    ptr[0] = 0;
    ptr[1] = 0;
    int *target = new int[1];
    sharedNodes = new Connectivity(0, ptr, target);
  }

  // get match ranges
  int numMatchedNodes, numMatchRanges, (*matchRanges)[2] = 0;
  if (matchNodes)
    numMatchedNodes = getRangeInfo(decFile, numMatchRanges, matchRanges);

  // allocate memory for the nodes, elements and faces
  NodeSet *nodes = new NodeSet(numNodes);
  ElemSet *elems = new ElemSet(numElems);
  FaceSet *faces = new FaceSet(numFaces);

  int *locToGlobNodeMap = new int[numNodes];
  int *locToGlobElemMap = new int[numElems];
  int *locToGlobFaceMap = new int[numFaces];

  // get the geometry file
  BinFileHandler geoFile(fullgeoName, "rb");

  // read the offset for the TOC
  BinFileHandler::OffType tocOffset[3];
  geoFile.read(tocOffset, 3);
  // read start of BCs
  BinFileHandler::OffType startOfBCsLoc;
  geoFile.read(&startOfBCsLoc, 1);

  int *locToClusNodeMap = new int[numNodes];

  // read the nodes
  geoFile.seek(tocOffset[0]);
  int numClusNodes = nodes->read(geoFile, numNodeRanges, nodeRanges, 
				 locToGlobNodeMap, locToClusNodeMap);
  (*nodes) *= oolscale;

  // read the elems
  geoFile.seek(tocOffset[1]);
  elems->read(geoFile, numElemRanges, elemRanges, locToGlobElemMap, iod->volumes.volumeMap.dataMap);

  // read the faces
  geoFile.seek(tocOffset[2]);
  faces->read(geoFile, numFaceRanges, faceRanges, locToGlobFaceMap);

  // read the boundary conditions for the entire cluster if necessary
  geoFile.seek(startOfBCsLoc);
  if(!clus_bconds) {
    clus_bconds = new BCondSet();
    clus_bconds->read(geoFile);
  }

  // read the match points
  if (matchNodes) {
    matchNodes[iSub] = new MatchNodeSet(numMatchedNodes);
    BinFileHandler matchFile(fullmatchName, "rb");
    matchFile.read(tocOffset, 1);
    matchFile.seek(tocOffset[0]);
    matchNodes[iSub]->read(matchFile, numMatchRanges, matchRanges);
  }

#ifdef OLD_STL
  map<int, int, less<int> > clusToLocNodeMap;
#else
  map<int, int> clusToLocNodeMap;
#endif

  int i;

  // make a map from cluster to local
  for (i=0; i<numNodes; ++i) 
    clusToLocNodeMap[ locToClusNodeMap[i] ] = i;

  // remap the elem nodes
  for (i=0; i<numElems; ++i)
    (*elems)[i].renumberNodes(clusToLocNodeMap);

  // remap the face nodes
  for (i=0; i<numFaces; ++i)
    (*faces)[i].renumberNodes(clusToLocNodeMap);

  // remap the shared nodes
  sharedNodes->renumberTargets(clusToLocNodeMap);

  // remap the matched nodes
  if (matchNodes)
    matchNodes[iSub]->renumberNodes(clusToLocNodeMap);

  // create the ranges for reading and writing the nodal solutions

  Vec<int> flag(numNodes);

  flag = 1;

  for (int jSub = 0; jSub < numNeighb; ++jSub)
    if ((subToCluster[neighb[jSub]][0] == clusterNum) && (globSubNum > neighb[jSub]))
      for (i=0; i<sharedNodes->num(jSub); ++i)
	flag[ (*sharedNodes)[jSub][i] ] = 0;

  int numRanges = 1;

  for (i=1; i<numNodes; ++i)
    if ((locToClusNodeMap[i] != locToClusNodeMap[i-1]+1) || (flag[i] != flag[i-1]))
      ++numRanges;

  int (*ranges)[3] = new int[numRanges][3];

  numRanges = 1;
  ranges[0][0] = 1;
  ranges[0][1] = locToClusNodeMap[0];
  ranges[0][2] = flag[0];
  for (i=1; i<numNodes; ++i) {
    if ((locToClusNodeMap[i] != locToClusNodeMap[i-1]+1) || (flag[i] != flag[i-1])) {
      ++numRanges;
      ranges[numRanges-1][0] = 1;
      ranges[numRanges-1][1] = locToClusNodeMap[i];
      ranges[numRanges-1][2] = flag[i];
    }
    else
      ranges[numRanges-1][0] += 1;
  }

  // create the subdomain
  SubDomain *sub = new SubDomain(iSub, clusterSubNum, globSubNum, numClusNodes, suffix, 
				 nodes, faces, elems,
                                 numNeighb, neighb, sharedNodes, 
				 locToGlobNodeMap, locToGlobFaceMap, locToGlobElemMap, 
				 numRanges, ranges);

  // get the local subdomain boundary conditions from the cluster boundary conditions
  distributeBCs(sub, clusToLocNodeMap);

  if (nodeRanges) delete [] nodeRanges;
  if (elemRanges) delete [] elemRanges;
  if (faceRanges) delete [] faceRanges;
  if (matchRanges) delete [] matchRanges;
  if (locToClusNodeMap) delete [] locToClusNodeMap;
  delete[] suffix;
  return sub;

}

//------------------------------------------------------------------------------

template<class MapType>
void GeoSource::distributeBCs(SubDomain *sub, MapType &cl2LocNodeMap)
{
  BCondSet *subBC;
  getBC(subBC, cl2LocNodeMap);
  if(subBC->size() > 0) sub->setBCond(subBC);
  else delete subBC;
}

template<class MapType>
void GeoSource::getBC(BCondSet *&subBC, MapType &cl2LocNodeMap)
{
  int numLocBC = 0;
  for(int iBC = 0; iBC < clus_bconds->size(); iBC++)  {
    if(cl2LocNodeMap[(*clus_bconds)[iBC].nnum] >= 0)
      numLocBC++;
  }

  // set BC
  subBC = new BCondSet(numLocBC);
  int count = 0;
  for(int iBC = 0; iBC < clus_bconds->size(); iBC++)  {
    if(cl2LocNodeMap[(*clus_bconds)[iBC].nnum] >= 0)  {
      (*subBC)[count] = (*clus_bconds)[iBC];
      // renumber to local node number
      (*subBC)[count++].nnum = cl2LocNodeMap[(*clus_bconds)[iBC].nnum];
    }
  }
}

//------------------------------------------------------------------------------
