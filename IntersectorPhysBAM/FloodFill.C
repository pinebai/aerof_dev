#include <cstdio>
#include <iostream>
#include <set>
#include <map>
#include <queue>

#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>

#include <Connectivity.h>
#include <Domain.h>
#include <Edge.h>
#include <Vector.h>
#include <DistInfo.h>
#include <DistVector.h>
#include "FloodFill.h"
#include "Mpi_Utilities.h"

static double total_time=0;

using std::set;
using std::map;
using std::pair;
using PHYSBAM_MPI_UTILITIES::masterProcessor;
using PhysBAM::GLOBAL_SUBD_ID;
//----------------------------------------------------------------------------
int FloodFill::floodFillSubDomain(SubDomain& sub, 
											 EdgeSet& edges, 
											 const Vec<bool>& edgeIntersections, 
											 const Vec<bool>& occluded_node, 
											 Vec<int>&nodeColors) 
{

   // A color of -2 means that this node has no neighbors

    int color=0,nColoredNodes=0,seedNode=0;

    const Connectivity &nToN = *sub.getNodeToNode(); 
	
	for(int i=0; i<sub.numNodes(); ++i){

		nodeColors[i]=-1;

		if(occluded_node[i]) {
			nodeColors[i]=-2;
			++nColoredNodes;
		}
	}

    std::queue<int> nodeQueue;
	while(nColoredNodes < sub.numNodes()){

		++color;

        while(seedNode < sub.numNodes() && nodeColors[seedNode] != -1) ++seedNode;
        nodeQueue.push(seedNode);

		while(!nodeQueue.empty()){

			++nColoredNodes;

			int currentNode = nodeQueue.front();
			nodeQueue.pop();
			nodeColors[currentNode] = color;

			for(int i=0; i<nToN.num(currentNode); ++i){

				int neighborNode = nToN[currentNode][i];

               if(currentNode == neighborNode) continue;
               if(nodeColors[neighborNode] == -1 && !edgeIntersections[edges.find(currentNode,neighborNode)]){
					nodeColors[neighborNode] = 0;
					nodeQueue.push(neighborNode);
				}
			}

		}

	}

#if 0 // Debug output
    bool has_bad_nodes=false;
    int colorCount[color+3];
    for(int i=0;i<color+3;++i) colorCount[i]=0;
    fprintf(stderr,"%02d Colored nodes = %d/%d, colors = %d\n",sub.getGlobSubNum(),nColoredNodes,sub.numNodes(),color);

	for(int i=0;i<sub.numNodes();++i){
		if(nodeColors[i] == 0 || nodeColors[i] == -1){
			fprintf(stderr,"There's still a bug!\n");has_bad_nodes=true;continue;
		}
		colorCount[nodeColors[i]+2]++;
	}

    for(int i=0;i<color+3;i++){
        fprintf(stderr,"%02d color %d has %d nodes\n",sub.getGlobSubNum(),i-2,colorCount[i]);
        if(colorCount[i]==1)
			for(int j=0;j<sub.numNodes();++j) 
				if(nodeColors[j]==i-2)
                fprintf(stderr,"\t(%02d) color %d is node %d\n",sub.getGlobSubNum(),i-2,j);}
    if(has_bad_nodes) exit(-1);
#endif
    return color; 
}

//----------------------------------------------------------------------------
void FloodFill::dumpColorsToFile(const std::string& prefix,SubDomain& subDomain,const SVec<double,3>& X,Vec<int>& color){
    GLOBAL_SUBD_ID gSub = PhysBAM::getGlobSubNum(subDomain);
    int max_color_number=color.max(),initial_offset=color.min();
    int *mapping = reinterpret_cast<int*>(alloca(color.size()*sizeof(int)));
    for(int current_color=initial_offset;current_color<=max_color_number;++current_color){
        int number_of_nodes=0,number_of_elements=0;
        for(int i=0;i<color.size();++i) if(color[i] == current_color) mapping[i]=++number_of_nodes;
        for(int i=0;i<subDomain.numElems();++i){int* nodes=subDomain.getElemNodeNum(i);
            if(color[nodes[0]] == current_color && color[nodes[1]] == current_color &&
               color[nodes[2]] == current_color && color[nodes[3]] == current_color) ++number_of_elements;}
        if(!number_of_nodes || !number_of_elements) continue;

        std::string fileName = PhysBAM::STRING_UTILITIES::string_sprintf("%s/subD-%d_color-%d",prefix.c_str(),gSub.Value(),current_color);
        FILE* tet_volume = fopen(fileName.c_str(),"w"); if(!tet_volume){fprintf(stderr,"SubD %d unable to write to '%s', exiting early\n",subDomain.getGlobSubNum(),fileName.c_str());return;}
        fprintf(tet_volume, "Nodes Nodes%d\n",number_of_nodes);
        for(int i=0;i<color.size();++i) if(color[i] == current_color)
            fprintf(tet_volume,"%d %e %e %e\n",mapping[i],X[i][0],X[i][1],X[i][2]);

        fprintf(tet_volume, "Elements Tets%d Using Nodes%d\n",number_of_elements,number_of_nodes);int element=0;
        for(int i=0;i<subDomain.numElems();++i){int* nodes=subDomain.getElemNodeNum(i);
            if(color[nodes[0]] == current_color && color[nodes[1]] == current_color &&
               color[nodes[2]] == current_color && color[nodes[3]] == current_color)
                fprintf(tet_volume,"%d 5 %d %d %d %d\n",++element,mapping[nodes[0]],mapping[nodes[1]],mapping[nodes[2]],mapping[nodes[3]]);}
        fclose(tet_volume);
    }
}

//----------------------------------------------------------------------------
// local helper functions...
//----------------------------------------------------------------------------
namespace {
void unionize(const set<pair<pair<GLOBAL_SUBD_ID,int>,pair<GLOBAL_SUBD_ID,int> > >& connections, map<pair<GLOBAL_SUBD_ID,int>,int>& localToGlobalColorMap)
{
    PhysBAM::UNION_FIND<int> union_find(localToGlobalColorMap.size()+1);
    typedef std::set<pair<pair<GLOBAL_SUBD_ID,int>,pair<GLOBAL_SUBD_ID,int> > > LIST;

    for(LIST::const_iterator iter=connections.begin();iter!=connections.end();++iter){
        union_find.Union(localToGlobalColorMap[iter->first],localToGlobalColorMap[iter->second]);}
    for(map<pair<GLOBAL_SUBD_ID,int>,int>::iterator iter=localToGlobalColorMap.begin();iter!=localToGlobalColorMap.end();++iter)
        localToGlobalColorMap[iter->first] = union_find.Find(iter->second);

#if 0 // Debug output
    for(LIST::const_iterator iter=connections.begin();iter!=connections.end();++iter){
        int firstColor=localToGlobalColorMap[iter->first],secondColor=localToGlobalColorMap[iter->second];
        if(firstColor != secondColor){fprintf(stderr,"Flood Fill Unionize failed to union all the colors...\n");exit(-1);}}

    set<int> globalColors;
    for(map<pair<GLOBAL_SUBD_ID,int>,int>::iterator iter=localToGlobalColorMap.begin();iter!=localToGlobalColorMap.end();++iter){
#if 0 // Detailed debug output
        fprintf(stderr,"(%d,%d) maps to global color %d\n",iter->first.first.Value(),iter->first.second,iter->second);
#endif
        globalColors.insert(iter->second);}
    for(set<int>::const_iterator iter=globalColors.begin();iter!=globalColors.end();iter++)
        fprintf(stderr,"Unionize determined %d was a global color\n",*iter);
#endif
}

}

//----------------------------------------------------------------------------
void FloodFill::generateConnectionsSet(Domain& domain,Communicator& com,DistVec<int>& nodeColors)
{
    connections.clear();
    PHYSBAM_MPI_UTILITIES::syncNodeColors(domain,com,nodeColors,connections);
    PHYSBAM_MPI_UTILITIES::gatherSet(com,connections);

#if 0 // Detailed debug output
    com.barrier();
    if(com.cpuNum()==masterProcessor)
		for(set<pair<pair<GLOBAL_SUBD_ID,int>,pair<GLOBAL_SUBD_ID,int> > >::const_iterator iter=connections.begin();
			 iter!=connections.end(); iter++)
			com.fprintf(stderr,"CONNECTED: (%d,%d) to (%d,%d)\n", iter->first.first.Value(),iter->first.second,
							iter->second.first.Value(),iter->second.second);
    com.barrier();
#endif
}

//----------------------------------------------------------------------------
void FloodFill::generateSubDToProcessorMap(Domain& domain,Communicator& com)
{
    int number_global_sub=-1;
#pragma omp parallel for
    for(int iSub=0;iSub<domain.getNumLocSub();++iSub)
	number_global_sub=std::max(number_global_sub,(domain.getSubDomain())[iSub]->getGlobSubNum()+1);

    PHYSBAM_MPI_UTILITIES::synchronizeMaxNumber(com,number_global_sub);
    nGlobalSubDomains=GLOBAL_SUBD_ID(number_global_sub);

    SubDomain **subD = domain.getSubDomain();
    subDomainToProcessorMap = new int[number_global_sub];

    for(GLOBAL_SUBD_ID gSub(0);gSub<nGlobalSubDomains;++gSub) subDomainToProcessorMap[gSub.Value()]=-1;
    // TODO(jontg): Probably better to get rid of this window communication too, but it's only ever called once...
    Communication::Window<int> windowSub2Proc(com,number_global_sub*sizeof(int),(int*)subDomainToProcessorMap);

    windowSub2Proc.fence(true);
#pragma omp parallel for
    for(int iSub=0; iSub<domain.getNumLocSub(); iSub++){
        subDomainToProcessorMap[PhysBAM::getGlobSubNum(*subD[iSub]).Value()]=com.cpuNum(); // This could be filled only once.
        int myCPU = com.cpuNum();
        windowSub2Proc.put((int*)&myCPU, 0, 1, masterProcessor, subD[iSub]->getGlobSubNum());}

    windowSub2Proc.fence(false);
#if 0 // Debug output
    if(com.cpuNum()==masterProcessor){
        fprintf(stderr,"subDomainToProcessorMap(%d): ",number_global_sub);
        for(GLOBAL_SUBD_ID gSub(0); gSub<nGlobalSubDomains; ++gSub)
            fprintf(stderr,"%d ",subDomainToProcessorMap[gSub.Value()]);
        fprintf(stderr,"\n");}
    com.barrier();
#endif
}

//----------------------------------------------------------------------------
void FloodFill::unionColors(Domain& domain, Communicator& com,
									 const int nLocalSubDomains, const int* localColorCount,
         map<pair<GLOBAL_SUBD_ID,int>,int>& localToGlobalColorMap)
{
    if(nGlobalSubDomains.Value()<0) generateSubDToProcessorMap(domain,com);

    // Only filled on the master processor, and valid on local processors.
    int *colorCount = reinterpret_cast<int*>(alloca(nGlobalSubDomains.Value()*sizeof(int))); 
    { // Send color counts to the master processor.
        // TODO(jontg): There should be a more efficient way to do this...
        for(int iSub=0;iSub<nLocalSubDomains;++iSub) colorCount[(domain.getSubDomain())[iSub]->getGlobSubNum()]=localColorCount[iSub];
        if(com.cpuNum()==masterProcessor){int rcvNum=0;
            MPI_Request *rcvId = reinterpret_cast<MPI_Request *>(alloca(sizeof(MPI_Request)*(nGlobalSubDomains.Value()-nLocalSubDomains)));
            MPI_Status *status = reinterpret_cast<MPI_Status *>(alloca(sizeof(MPI_Status)*(nGlobalSubDomains.Value()-nLocalSubDomains)));
            for(GLOBAL_SUBD_ID gSub(0);gSub<nGlobalSubDomains;++gSub) if(subDomainToProcessorMap[gSub.Value()]!=masterProcessor){
                MPI_Irecv(colorCount+gSub.Value(),1,MPI_INT,subDomainToProcessorMap[gSub.Value()],UNION_COLORS_TAG,com.comm,rcvId+(rcvNum++));}
            MPI_Waitall(rcvNum,rcvId,status);}
        else{int sndNum=0;
            MPI_Request *sndId = reinterpret_cast<MPI_Request *>(alloca(sizeof(MPI_Request)*(nLocalSubDomains)));
            MPI_Status *status = reinterpret_cast<MPI_Status *>(alloca(sizeof(MPI_Status)*(nLocalSubDomains)));
            for(GLOBAL_SUBD_ID gSub(0);gSub<nGlobalSubDomains;++gSub) if(subDomainToProcessorMap[gSub.Value()]==com.cpuNum()){
                MPI_Isend(colorCount+gSub.Value(),1,MPI_INT,masterProcessor,UNION_COLORS_TAG,com.comm,sndId+(sndNum++));}
            MPI_Waitall(sndNum,sndId,status);}
    }

#if 0 // Debug output
    com.barrier();
    if(com.cpuNum()==masterProcessor){
        fprintf(stderr,"colorCount: ");for(GLOBAL_SUBD_ID gSub(0);gSub<nGlobalSubDomains;++gSub) fprintf(stderr,"%d ",colorCount[gSub.Value()]); fprintf(stderr,"\n");}
    com.barrier();
#endif

    localToGlobalColorMap.clear();
    if(com.cpuNum()==masterProcessor){
        int globalColorCount=0;
        for(GLOBAL_SUBD_ID gSub(0);gSub<nGlobalSubDomains;++gSub) for(int localColor=1;localColor<=colorCount[gSub.Value()];++localColor)
            localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(gSub,localColor)]=++globalColorCount;
        unionize(connections,localToGlobalColorMap);}

    PHYSBAM_MPI_UTILITIES::distributeColorMap(domain,com,colorCount,nGlobalSubDomains,subDomainToProcessorMap,localToGlobalColorMap);
}
