#ifndef _FLOODFILL_H_
#define _FLOODFILL_H_

#include <set>
#include <map>
#include <Vector.h>
#include <Mpi_Utilities.h>

#include "Communicator.h"

class Domain;
class EdgeSet;
class SubDomain;

using std::pair;
using PhysBAM::GLOBAL_SUBD_ID;

// TODO(jontg): const correctness
class FloodFill {
private:
// Data useful for distributing the color map -- populated by generateSubDToProcessorMap(...)
GLOBAL_SUBD_ID nGlobalSubDomains;
int* subDomainToProcessorMap;
std::set<pair<pair<GLOBAL_SUBD_ID,int>,pair<GLOBAL_SUBD_ID,int> > > connections;

public:
FloodFill() : nGlobalSubDomains(-1),subDomainToProcessorMap(0)
{}

~FloodFill()
{delete [] subDomainToProcessorMap;}

void generateSubDToProcessorMap(Domain& domain,Communicator& com);

void generateConnectionsSet(Domain& domain,Communicator& com,DistVec<int>& nodeColors);
void unionColors(Domain& domain,Communicator& com,const int nLocalSubDomains,const int* localColorCount,
     std::map<pair<PhysBAM::GLOBAL_SUBD_ID,int>,int>& localToGlobalColorMap);

static int floodFillSubDomain(SubDomain& subDomain, 
										EdgeSet& edges, 
										const Vec<bool>& edgeIntersections, 
										const Vec<bool>& occluded_node, 
										Vec<int>& nodeColors);

static void dumpColorsToFile(const std::string& prefix,SubDomain& subDomain,const SVec<double,3>& X,Vec<int>& color);
};

#endif
