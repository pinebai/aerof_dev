#include <cstdio>
#include <execinfo.h>
#include <iostream>
#include <set>
#include <map>
#include <queue>

#include <Communicator.h>
#include <Connectivity.h>
#include <Domain.h>
#include <SubDomain.h>

#include "Mpi_Utilities.h"

using std::set;
using std::map;
using std::pair;

using PhysBAM::GLOBAL_SUBD_ID;

//----------------------------------------------------------------------------
// ELEMENT_ID to help keep track of global sub-domains.
//----------------------------------------------------------------------------
namespace PhysBAM{
    GLOBAL_SUBD_ID getGlobSubNum(SubDomain& subD)
    {return GLOBAL_SUBD_ID(subD.getGlobSubNum());}
};
//----------------------------------------------------------------------------
// local helper functions...
//----------------------------------------------------------------------------
namespace PHYSBAM_MPI_UTILITIES{
//----------------------------------------------------------------------------
void gatherSet(Communicator& com, set<pair<pair<GLOBAL_SUBD_ID,int>,pair<GLOBAL_SUBD_ID,int> > >& set)
{
    typedef std::set<pair<pair<GLOBAL_SUBD_ID,int>,pair<GLOBAL_SUBD_ID,int> > > LIST;
    typedef LIST::const_iterator set_iterator;

    int *setLengths = reinterpret_cast<int*>(alloca(com.size()*sizeof(int)));
    for(int i=0; i<com.size(); ++i) setLengths[i]=0;
    gatherNumbers(com,(int)set.size(),(int*)setLengths);

    int bufferSize=0;
    if(com.cpuNum()==masterProcessor) for(int i=0;i<com.size();++i) bufferSize=std::max(bufferSize,4*setLengths[i]);++bufferSize;
    distributeValue(com,bufferSize,bufferSize);

    int *buffer = reinterpret_cast<int*>(alloca((com.cpuNum()==masterProcessor?com.size()*bufferSize:bufferSize)*sizeof(int)));
    int i=0;
    for(set_iterator it=set.begin();it != set.end();it++){const pair<pair<GLOBAL_SUBD_ID,int>,pair<GLOBAL_SUBD_ID,int> >& p=*it;
        buffer[i++]=p.first.first.Value();buffer[i++]=p.first.second; buffer[i++]=p.second.first.Value();buffer[i++]=p.second.second;}

    if(com.cpuNum()==masterProcessor){
      MPI_Status *status = reinterpret_cast<MPI_Status *>(alloca(sizeof(MPI_Status)*(com.size()-1)));
      MPI_Request *rcvId = reinterpret_cast<MPI_Request *>(alloca(sizeof(MPI_Request)*(com.size()-1)));
      for(int i=1;i<com.size();++i)
          MPI_Irecv((int*)(buffer+i*bufferSize), 4*setLengths[i], MPI_INT, i, GATHER_SET_TAG, com.comm, rcvId+i-1);
      MPI_Waitall(com.size()-1,rcvId,status);}
    else MPI_Send((int*)buffer, 4*setLengths[com.cpuNum()], MPI_INT, 0, GATHER_SET_TAG, com.comm);

    if(masterProcessor == com.cpuNum())
        for(int iCPU=1;iCPU<com.size();++iCPU){ pair<pair<GLOBAL_SUBD_ID,int>,pair<GLOBAL_SUBD_ID,int> > buffer_input;
            for(int i=0;i<setLengths[iCPU];++i){
                buffer_input.first.first=GLOBAL_SUBD_ID(buffer[iCPU*bufferSize+4*i]);buffer_input.first.second=buffer[iCPU*bufferSize+4*i+1];
                buffer_input.second.first=GLOBAL_SUBD_ID(buffer[iCPU*bufferSize+4*i+2]);buffer_input.second.second=buffer[iCPU*bufferSize+4*i+3];
                set.insert(buffer_input);}}
}

//----------------------------------------------------------------------------
void distributeColorMap(Domain& domain, Communicator& com,const int* colorCount,const GLOBAL_SUBD_ID nGlobalSubDomains,const int* subDomainToProcessorMap,
    map<pair<GLOBAL_SUBD_ID,int>,int>& localToGlobalColorMap)
{
    int bufferSize=-1;
    if(com.cpuNum()==masterProcessor) for(GLOBAL_SUBD_ID gSub(0);gSub<nGlobalSubDomains;++gSub)
        bufferSize=std::max(bufferSize,colorCount[gSub.Value()]);
    distributeValue(com,bufferSize,bufferSize);

    int *buffer = reinterpret_cast<int*>(alloca(2*nGlobalSubDomains.Value()*bufferSize*sizeof(int)));
    if(com.cpuNum()==masterProcessor)
	for(GLOBAL_SUBD_ID gSub(0);gSub<nGlobalSubDomains;++gSub) for(int i=0;i<colorCount[gSub.Value()];++i){
            buffer[gSub.Value()*bufferSize+i]=localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(gSub,i+1)];}

    int index=0;
    MPI_Status *status = reinterpret_cast<MPI_Status *>(alloca(sizeof(MPI_Status)*nGlobalSubDomains.Value()));
    MPI_Request *mpiId = reinterpret_cast<MPI_Request *>(alloca(sizeof(MPI_Request)*nGlobalSubDomains.Value()));
    for(int i=0;i<nGlobalSubDomains.Value();++i){
        if(subDomainToProcessorMap[i]==masterProcessor) continue;
        if(com.cpuNum()==masterProcessor)
            MPI_Isend(buffer+bufferSize*i, colorCount[i], MPI_INT, subDomainToProcessorMap[i], DISTRIBUTE_COLOR_MAP_TAG, com.comm, mpiId+(index++));
        else if(com.cpuNum()==subDomainToProcessorMap[i]){
            MPI_Status status;MPI_Recv(buffer+bufferSize*i, colorCount[i], MPI_INT, 0, DISTRIBUTE_COLOR_MAP_TAG, com.comm, &status);}}
    if(com.cpuNum()==masterProcessor) MPI_Waitall(index,mpiId,status);

    for(GLOBAL_SUBD_ID gSub(0);gSub<nGlobalSubDomains;++gSub)
        if(com.cpuNum()==subDomainToProcessorMap[gSub.Value()]) for(int i=0;i<colorCount[gSub.Value()];++i)
            localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(gSub,i+1)]=buffer[gSub.Value()*bufferSize+i];

#if 0 // Debug output
    for(map<pair<GLOBAL_SUBD_ID,int>,int>::const_iterator iter=localToGlobalColorMap.begin();iter!=localToGlobalColorMap.end();iter++)
        fprintf(stderr,"%02d: Mapping (Global SubD %d, Local Color %d) To Global Color %d\n",com.cpuNum(),iter->first.first.Value(),iter->first.second,iter->second);
#endif
}

//----------------------------------------------------------------------------
void syncMap(Domain& domain, Communicator& com,map<int,int>& localMap){
    map<int,int> globalMap;
    { // Gather all the local maps to the master processor.
        int bufferSize=0;
        int *mapLengths = reinterpret_cast<int*>(alloca(com.size()*sizeof(int)));
        gatherNumbers(com,(int)localMap.size(),(int*)mapLengths);
        if(com.cpuNum()==masterProcessor) for(int i=0;i<com.size();++i) bufferSize=std::max(bufferSize,2*mapLengths[i]);
        distributeValue(com,bufferSize,bufferSize);

        int offset=0;
        int *buffer = reinterpret_cast<int*>(alloca((com.cpuNum()==masterProcessor?com.size()*bufferSize:bufferSize)*sizeof(int)));
        for(map<int,int>::const_iterator iter=localMap.begin();iter!=localMap.end();offset+=2,iter++){
            // fprintf(stderr,"%d REPORTING globalMap[%d] = %d\n", com.cpuNum(), iter->first, iter->second);
            buffer[offset]=iter->first;buffer[offset+1]=iter->second;}

        if(com.cpuNum()==masterProcessor){
            MPI_Status *status = reinterpret_cast<MPI_Status *>(alloca(sizeof(MPI_Status)*(com.size()-1)));
            MPI_Request *rcvId = reinterpret_cast<MPI_Request *>(alloca(sizeof(MPI_Request)*(com.size()-1)));
            for(int i=1;i<com.size();++i)
                MPI_Irecv((int*)(buffer+i*bufferSize), 2*mapLengths[i], MPI_INT, i, GATHER_SET_TAG, com.comm, rcvId+i-1);
            MPI_Waitall(com.size()-1,rcvId,status);}
        else MPI_Send((int*)buffer, 2*mapLengths[com.cpuNum()], MPI_INT, 0, GATHER_SET_TAG, com.comm);

        if(com.cpuNum()==masterProcessor)
            for(int iCPU=0;iCPU<com.size();++iCPU)
                for(int i=0;i<mapLengths[iCPU];++i){
                    std::pair<map<int,int>::iterator,bool> result=globalMap.insert(std::make_pair(buffer[iCPU*bufferSize+2*i],buffer[iCPU*bufferSize+2*i+1]));
                    if(!result.second && globalMap[buffer[iCPU*bufferSize+2*i]] != buffer[iCPU*bufferSize+2*i+1]){
                        if(result.first->second == 0){
			    fprintf(stderr,"Tried to assign globalMap[%d] = %d, but globalMap[%d] = %d.  OVERWRITING VALUE!\n",
					    buffer[iCPU*bufferSize+2*i], buffer[iCPU*bufferSize+2*i+1], buffer[iCPU*bufferSize+2*i], result.first->second);
			    globalMap.erase(result.first);
			    globalMap.insert(std::make_pair(buffer[iCPU*bufferSize+2*i],buffer[iCPU*bufferSize+2*i+1]));
			} else if(buffer[iCPU*bufferSize+2*i+1] != 0){
			    fprintf(stderr,"CONFLICT: Tried to set globalMap[color %d]=%d, but %d was already there!\n",buffer[iCPU*bufferSize+2*i],buffer[iCPU*bufferSize+2*i+1],globalMap[buffer[iCPU*bufferSize+2*i]]);exit(-1);
			}
#if 0 // Debug output
			else fprintf(stderr,"Tried to insert FF color; ignoring...\n");
#endif
		    }
		}
    }

#if 0 // Debug output
    if(com.cpuNum()==masterProcessor) for(map<int,int>::const_iterator iter=globalMap.begin();iter!=globalMap.end();iter++)
        fprintf(stderr,"ON MASTER PROC, Mapping %d to %d\n",iter->first,iter->second);
#endif

    { // Distribute the global map to every processor.
        int bufferSize=0;
        if(com.cpuNum()==masterProcessor) bufferSize=2*globalMap.size();
        distributeValue(com,bufferSize,bufferSize);
        int *buffer = reinterpret_cast<int*>(alloca(bufferSize*sizeof(int)));

        if(com.cpuNum()==masterProcessor){int offset=0;for(map<int,int>::const_iterator iter=globalMap.begin();iter!=globalMap.end();offset+=2,iter++){
            buffer[offset]=iter->first;buffer[offset+1]=iter->second;}}

        if(com.cpuNum()==masterProcessor){
            MPI_Status *status = reinterpret_cast<MPI_Status *>(alloca(sizeof(MPI_Status)*(com.size()-1)));
            MPI_Request *rcvId = reinterpret_cast<MPI_Request *>(alloca(sizeof(MPI_Request)*(com.size()-1)));
            for(int i=1;i<com.size();++i) MPI_Isend(buffer, bufferSize, MPI_INT, i, DISTRIBUTE_SET_TAG,com.comm,rcvId+i-1);
            MPI_Waitall(com.size()-1,rcvId,status);}
        else{
            MPI_Status status;
            MPI_Recv(buffer, bufferSize, MPI_INT, 0, DISTRIBUTE_SET_TAG, com.comm,&status);
            for(int i=0;i<bufferSize;i+=2) globalMap[buffer[i]]=buffer[i+1];}
    }

    localMap.swap(globalMap);
}

//----------------------------------------------------------------------------
static CommPattern<pair<int,int> >* colorPattern=0;
//the following two lines cause "cannot be specialized in the current scope" error using intel compiler version 11.1
//template <> MPI_Datatype CommTrace<pair<int,int> >::MPIType=MPI_INT;
//template <> int CommTrace<pair<int,int> >::multiplicity=2;

void assembleIntoSets(Domain& domain, Communicator& com,DistVec<pair<int,int> >& localColor,DistVec<set<pair<int,int> > >& nodeLocalColors){
    int numLocSub=domain.getNumLocSub();
    if(!colorPattern){
        colorPattern=new CommPattern<pair<int,int> >(domain.getSubTopo(),&com,CommPattern<pair<int,int> >::CopyOnSend);
        for(int iSub=0; iSub<numLocSub; ++iSub) (domain.getSubDomain()[iSub])->setComLenNodes(2,*colorPattern);
        colorPattern->finalize();}

#pragma omp parallel for
    for(int iSub=0; iSub < numLocSub; ++iSub)
        for(int i=0; i<localColor(iSub).size(); ++i)
            nodeLocalColors(iSub)[i].insert(localColor(iSub)[i]);

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
        domain.getSubDomain()[iSub]->sndData(*colorPattern, reinterpret_cast<pair<int,int> (*)[1]>(localColor.subData(iSub)));
    colorPattern->exchange();

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub){
        SubDomain& sub=*(domain.getSubDomain()[iSub]);
        set<pair<int,int> > (*rcvBuffer)[1]=reinterpret_cast<set<pair<int,int> > (*)[1]>(nodeLocalColors.subData(iSub));
        for(int nSub=0; nSub<sub.getNumNeighb(); ++nSub){
            SubRecInfo<pair<int,int> > sInfo=colorPattern->recData(sub.getRcvChannel()[nSub]);
            pair<int,int> (*buffer)[1]=reinterpret_cast<pair<int,int> (*)[1]>(sInfo.data);
            for(int iNode=0; iNode<sub.getSharedNodes()->num(nSub); ++iNode){
                rcvBuffer[(*sub.getSharedNodes())[nSub][iNode]][0].insert(buffer[iNode][0]);}}}
}

//----------------------------------------------------------------------------
void syncNodeColors(Domain& domain, Communicator& com,const DistVec<int>& nodeColors_input,set<pair<pair<GLOBAL_SUBD_ID,int>,pair<GLOBAL_SUBD_ID,int> > >& connections){
    int numLocSub=domain.getNumLocSub();
    DistVec<pair<int,int> > localColor(domain.getNodeDistInfo());
    DistVec<set<pair<int,int> > > nodeLocalColor(domain.getNodeDistInfo());
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
        SubDomain& sub=*(domain.getSubDomain()[iSub]);
        for(int i=0; i<nodeLocalColor(iSub).size(); ++i)
            localColor(iSub)[i]=pair<int,int>(PhysBAM::getGlobSubNum(sub).Value(),nodeColors_input(iSub)[i]);}

    assembleIntoSets(domain,com,localColor,nodeLocalColor);

    for(int iSub=0; iSub<numLocSub; ++iSub)
        for(int i=0; i<nodeLocalColor(iSub).size(); ++i)
            if(nodeLocalColor(iSub)[i].size() > 1){
                set<pair<int,int> >::const_iterator first_element=nodeLocalColor(iSub)[i].begin();
                if(first_element->second <= 0) continue;
                for(set<pair<int,int> >::const_iterator iter=first_element; iter!=nodeLocalColor(iSub)[i].end();iter++){
                    if(iter==first_element || iter->second <= 0) continue; // connections to uncolored, undecided and occluded nodes should be neglected
                    else connections.insert(pair<pair<GLOBAL_SUBD_ID,int>,pair<GLOBAL_SUBD_ID,int> >(
                            pair<GLOBAL_SUBD_ID,int>(GLOBAL_SUBD_ID(first_element->first),first_element->second),
                            pair<GLOBAL_SUBD_ID,int>(GLOBAL_SUBD_ID(iter->first),iter->second)));}}
}

//----------------------------------------------------------------------------
void dump_stack_trace()
{
       void *array[100];
       size_t size;
       char **strings;
       size_t i;
     
       size = backtrace (array, 100);
       strings = backtrace_symbols (array, size);
     
       printf ("Obtained %zd stack frames.\n", size);
     
       for (i = 0; i < size; i++)
          printf ("%s\n", strings[i]);
     
       free (strings);

}
}
template <> MPI_Datatype CommTrace<pair<int,int> >::MPIType=MPI_INT;
template <> int CommTrace<pair<int,int> >::multiplicity=2;
