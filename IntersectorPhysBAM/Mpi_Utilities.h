#ifndef __PHYSBAM_MPI_UTILITIES_H__
#define __PHYSBAM_MPI_UTILITIES_H__

#include <cstdio>
#include <iostream>
#include <set>
#include <map>
#include <queue>

#include <Connectivity.h>
#include <Domain.h>

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>

using std::set;
using std::map;
using std::pair;

//----------------------------------------------------------------------------
// ELEMENT_ID to help keep track of global sub-domains.
//----------------------------------------------------------------------------
namespace PhysBAM{
    PHYSBAM_DECLARE_ELEMENT_ID(GLOBAL_SUBD_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);

    GLOBAL_SUBD_ID getGlobSubNum(SubDomain& subD);
};
//----------------------------------------------------------------------------
// local helper functions...
//----------------------------------------------------------------------------
namespace PHYSBAM_MPI_UTILITIES{
static const int masterProcessor=0;
#define GATHER_NUMBERS_TAG 10050
#define DISTRIBUTE_NUMBERS_TAG 20050
#define GATHER_SET_TAG 30050
#define DISTRIBUTE_SET_TAG 40050
#define DISTRIBUTE_COLOR_MAP_TAG 50050
#define SYNC_MAP_SEND_TAG 60050
#define SYNC_MAP_RECV_TAG 70050
#define UNION_COLORS_TAG 80050

//----------------------------------------------------------------------------
template<class T>
void gatherNumbers(Communicator& com,T localValue,T* toGather)
{
  toGather[com.cpuNum()]=localValue;
  if(com.cpuNum()==masterProcessor){
      MPI_Status *status = reinterpret_cast<MPI_Status *>(alloca(sizeof(MPI_Status)*(com.size()-1)));
      MPI_Request *rcvId = reinterpret_cast<MPI_Request *>(alloca(sizeof(MPI_Request)*(com.size()-1)));
      for(int i=1;i<com.size();++i)
          MPI_Irecv(toGather+i, CommTrace<T>::multiplicity, CommTrace<T>::MPIType, i, GATHER_NUMBERS_TAG, com.comm, rcvId+i-1);
      MPI_Waitall(com.size()-1,rcvId,status);}
  else MPI_Send(toGather+com.cpuNum(), CommTrace<T>::multiplicity, CommTrace<T>::MPIType, 0, GATHER_NUMBERS_TAG, com.comm);
}

//----------------------------------------------------------------------------
template<class T>
void distributeValue(Communicator& com,T& localValue,T toDistribute)
{
    localValue=toDistribute;
    com.broadcast(1,&localValue);
}

//----------------------------------------------------------------------------
template<class T>
void synchronizeMaxNumber(Communicator& com, T& numberToMax)
{
    com.globalMax(1,(T*)&numberToMax);
}

//----------------------------------------------------------------------------
void gatherSet(Communicator& com, set<pair<pair<PhysBAM::GLOBAL_SUBD_ID,int>,pair<PhysBAM::GLOBAL_SUBD_ID,int> > >& set);
void distributeColorMap(Domain& domain, Communicator& com,const int* colorCount,const PhysBAM::GLOBAL_SUBD_ID nGlobalSubDomains,
    const int* subDomainToProcessorMap,map<pair<PhysBAM::GLOBAL_SUBD_ID,int>,int>& localToGlobalColorMap);
void syncMap(Domain& domain, Communicator& com,map<int,int>& localMap);
void syncNodeColors(Domain& domain, Communicator& com,const DistVec<int>& nodeColors_input,
    set<pair<pair<PhysBAM::GLOBAL_SUBD_ID,int>,pair<PhysBAM::GLOBAL_SUBD_ID,int> > >& connections_set);

void dump_stack_trace();
}

#endif
