#include <BCApplier.h>
#include <Domain.h>
#include <SubDomain.h>
#include <IoData.h>
#include <MatchNode.h>
#include <DistVector.h>

#define HB_MESHMOTION_WARNINGS

BCApplier::BCApplier(int nLocSub, Domain* dom, IoData& ioData):domain(dom),surfaceMap(ioData.surfaces.surfaceMap.dataMap),
                                                               slidingSurfaceTreatment(ioData.dmesh.slidingSurfaceTreatment)
{
 numLocSub = nLocSub;
 dofType   = 0;
 proj      = new list<ProjData>[numLocSub];
}

void
BCApplier::addProj(int locSubNum, int nd, double nrm[3]) 
{
  proj[locSubNum].push_back(ProjData(nd, nrm));
}

void
BCApplier::setDofType(MatchNodeSet** matchNodes)
{
  SubDomain** subD  = domain->getSubDomain();
  SubDTopo* subTopo = domain->getSubTopo();
  Communicator* com = domain->getCommunicator();

#ifdef HB_MESHMOTION_WARNINGS
  if(matchNodes) {
    int** ndType = domain->getNodeType();
    int nMoving = 0; // counts (global) number of nodes "labelled" as moving
    int nMatched= 0; // counts (global) number of "matched" nodes
    int i;
#pragma omp parallel for private(i) reduction(+: nMoving,nMatched)
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      bool* ndFlag = domain->getNodeDistInfo().getMasterFlag(iSub); // nodes master flag array -> to count each node ONLY once
      // Count number of nodes "labelled" as moving in iSub 
      for (i=0; i<domain->getNodeDistInfo().subSize(iSub); i++) // Loop over sub's nodes
        if(ndFlag[i] && ndType[iSub][i]<BC_INTERNAL) nMoving++; // iSub is the "master" for this node which is labelled as moving 
      // Count number of "matched" nodes iSub 
      for (i=0; i<matchNodes[iSub]->size(); i++)  // Loop over sub's matched nodes
        if(ndFlag[matchNodes[iSub]->subMatchNode(i)]) nMatched++; // iSub is the "master" for this node which is "matched" 
    }
    // sum across the CPUs
    int ndCount[2] = {nMoving,nMatched};
    com->globalSum(2, ndCount); // ndCount[0] <- global number of "moving" nodes
    if(ndCount[0]!=ndCount[1])  // ndCount[1] <- global number of "matched" nodes 
      com->fprintf(stderr," ... WARNING: there are %d moving nodes that are not matched.\n",ndCount[0]-ndCount[1]);
  }
#endif
 
  dofType = new int*[numLocSub];
  dofTypeStep1 = 0;
  dofTypeStep2 = 0;
  CommPattern<int> ndC(subTopo, com, CommPattern<int>::CopyOnSend);

#pragma omp parallel for
  for (int iSub = 0; iSub<numLocSub; ++iSub)
    subD[iSub]->setComLenNodes(3, ndC);

  ndC.finalize();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    MatchNodeSet* subMatchNodes = (matchNodes) ? matchNodes[iSub] : 0;
    dofType[iSub] = subD[iSub]->getMeshMotionDofType(surfaceMap, ndC, subMatchNodes);
  }
  
  ndC.exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub<numLocSub; ++iSub)
    subD[iSub]->completeMeshMotionDofType(dofType[iSub], ndC);
}


//---------------------------------------------------------------------------------------------------------------
void
BCApplier::setEmbeddedALEDofType(MatchNodeSet** matchNodes)
{
  SubDomain** subD  = domain->getSubDomain();
  SubDTopo* subTopo = domain->getSubTopo();
  Communicator* com = domain->getCommunicator();

  dofType = new int*[numLocSub];
  CommPattern<int> ndC(subTopo, com, CommPattern<int>::CopyOnSend);

#pragma omp parallel for
  for (int iSub = 0; iSub<numLocSub; ++iSub)
    subD[iSub]->setComLenNodes(3, ndC);

  ndC.finalize();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    MatchNodeSet* subMatchNodes = (matchNodes) ? matchNodes[iSub] : 0;
    dofType[iSub] = subD[iSub]->getEmbeddedALEMeshMotionDofType(surfaceMap, ndC, subMatchNodes, slidingSurfaceTreatment);
  }
  
  ndC.exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub<numLocSub; ++iSub)
    subD[iSub]->completeEmbeddedALEMeshMotionDofType(dofType[iSub], ndC);
}


//---------------------------------------------------------------------------------------------------------------

void 
BCApplier::print()
{
  SubDomain** subD = domain->getSubDomain();
  Communicator* com = domain->getCommunicator();
//#ifdef HB_MESHMOTION_DEBUG
  com->sync();
  for(int iCPU=0; iCPU<com->size(); iCPU++) {
    if(iCPU==com->cpuNum()) {
      for (int iSub = 0; iSub < numLocSub; ++iSub) {
        if(proj[iSub].size()) { fprintf(stderr," ... BCApplier of subd %4d contains %5d projections\n.",
                                subD[iSub]->getGlobSubNum(),int(proj[iSub].size())); fflush(stderr); }
      }
    }
    com->sync();
  }
//#endif
  //int count = 0;
  //if(!dofType[iSub]) fprintf(stderr," ERROR: Subd %d, dofType = 0\n");
  //for(int i=0;i<subD[iSub]->numNodes(); i++) {
  //  for(int l=0; l<3; l++) 
  //    if(dofType[iSub][3*i+l]!=BC_FREE) count++;
  //}
  //fprintf(stderr," Subd %d has %d constrained dofs\n",subD[iSub]->getGlobSubNum(),count);
}

BCApplier::~BCApplier()
{
  if(dofType) {
#pragma omp parallel for
    for (int iSub = 0; iSub<numLocSub; ++iSub)
      if(dofType[iSub]) { delete [] dofType[iSub]; dofType[iSub] = 0; }
    delete [] dofType; dofType = 0;
  }
  if(proj) { delete[] proj; proj = 0; }
  numLocSub = 0;
  domain    = 0;
}

// Rotate the projection normal according to the given (incremental) rotation matrix dRot
// Note that due to numerical round-off, we potentially lost the orthogonality of the 
// (at most 3) normal define at a node. If this is an issue, we may re-orthogonalize them
// periodically. Also as this rotation (& re-orthogonalization if any) is done independantly 
// (i.e. no communication) for shared nodes, those nodes may potentially end up having 
// slightly different projections on different subdomains. 
void
BCApplier::rotateProjNormal(double dRot[3][3])
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    list<ProjData>::iterator it = proj[iSub].begin();
    while(it != proj[iSub].end()){
      it->rotateNormal(dRot); 
      it++;
    }
  }
}
