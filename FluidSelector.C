#include "LevelSet/LevelSetStructure.h"
//#include <Domain.h>
#include <DebugTools.h>
#include <cassert>

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::initializeFluidIds(DistSVec<double,dim> &Phin, DistSVec<double,dim> &Phinm1, DistSVec<double,dim> &Phinm2){
  getFluidId(Phin);
  *fluidIdn = *fluidId;
  if(fluidIdnm1) getFluidId(*fluidIdnm1, Phinm1);
  if(fluidIdnm2) getFluidId(*fluidIdnm2, Phinm2);
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::reinitializeFluidIds(DistVec<int> &fsId, DistSVec<double,dim> &Phin)
{
  getFluidId(*fluidId, Phin, &fsId);
  //if (fluidIdnm1) getFluidId(*fluidIdm1, Phinm1, &fsId);
  //if(fluidIdnm2) getFluidId(*fluidIdm2, Phinm1, &fsId);
  *fluidIdn = *fluidId;
  if(fluidIdnm1) *fluidIdnm1 = *fluidId;
  if(fluidIdnm2) *fluidIdnm2 = *fluidId;
}

template<int dim>
void FluidSelector::reinitializeFluidIdsWithCracking(DistVec<int> &fsId, DistSVec<double,dim> &Phin)
{
  getFluidId(*fluidId, Phin, NULL);
  //if (fluidIdnm1) getFluidId(*fluidIdm1, Phinm1, &fsId);
  //if(fluidIdnm2) getFluidId(*fluidIdm2, Phinm1, &fsId);
  *fluidIdn = *fluidId;
  if(fluidIdnm1) *fluidIdnm1 = *fluidId;
  if(fluidIdnm2) *fluidIdnm2 = *fluidId;
}


//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(DistSVec<double,dim> &Phi){
  assert(dim<=numPhases-1);
  int numLocSub = Phi.numLocSub();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int     *tag       = fluidId->subData(iSub);
    int burnTag;
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      //if (programmedBurn && programmedBurn->isBurnedEOS(tag[iNode],burnTag))
      //	continue;
      tag[iNode] = 0;
      for(int i=0; i<dim; i++) {
        if(phi[iNode][i]>0.0) {
	  if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
				 programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
	    if (programmedBurn->nodeInside(burnTag,iSub,iNode) || 
                programmedBurn->isFinished(burnTag))
	      tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
	    else
	      tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
	    break;
	  }
	  else 
	    {  tag[iNode] = i+1; break; }
	}
      }
    }
  }
}
//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(int &tag, double *phi){
  assert(dim<=numPhases-1);
  tag = 0;
  for(int i=0; i<dim; i++)
    if(phi[i]>0.0) {tag = i+1; return; }
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(Vec<int> &tag, SVec<double,dim> &phi){
  assert(dim<=numPhases-1);
  int burnTag;
  for(int iNode=0; iNode<phi.size(); iNode++){
    tag[iNode] = 0;
    for(int i=0; i<dim; i++) {
      if(phi[iNode][i]>0.0) {
	if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
			       programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
	  if (programmedBurn->nodeInside(burnTag,iNode) ||
              programmedBurn->isFinished(burnTag))
	    tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
	  else
	    tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
	  break;
	}
	else
	  { tag[iNode] = i+1; break; }
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void FluidSelector::getFluidId(DistVec<int> &Tag, DistSVec<double,dim> &Phi, DistVec<int>* fsId){
  assert(dim<=numPhases-1);
  //std::cout << "Dim = " << dim << std::endl;
  int numLocSub = Phi.numLocSub();
  int oldtag;
  int iSub;
  //std::cout << programmedBurn << std::endl;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int     *tag       = Tag.subData(iSub);
    int     *fsid      = fsId ? fsId->subData(iSub) : 0;
    int burnTag;
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      //if (programmedBurn && programmedBurn->isBurnedEOS(tag[iNode],burnTag))
      //	continue;

      tag[iNode] = 0;
      if(fsid && fsid[iNode]!=0) //isolated by structure. use fsId as fluidId.
        tag[iNode] = fsid[iNode];
      else for(int i=0; i<dim; i++) {
	if(phi[iNode][i]>0.0) { 
	  if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
				 programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
	    if (programmedBurn->nodeInside(burnTag,iSub,iNode) ||
                programmedBurn->isFinished(burnTag))
	      tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
	    else
	      tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
	    break;
	  }
	  else{
	    tag[iNode] = i+1; break; 
	  }
	}
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim> /*this dim is actually dimLS*/
void FluidSelector::updateFluidIdFS(DistLevelSetStructure *distLSS, DistSVec<double,dim> &PhiV)
{
  assert(dim<=numPhases-1);
  DistVec<int> &fsId(distLSS->getStatus());

  int burnTag;
  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<PhiV.numLocSub(); ++iSub) {
    Vec<int> &subfsId(fsId(iSub));
    Vec<int> &subId((*fluidId)(iSub));
    SVec<double,dim> &subPhiV(PhiV(iSub));
    LevelSetStructure &LSS((*distLSS)(iSub));

    for(int i=0; i<subPhiV.size(); i++) {
      int Id = LSS.fluidModel(0.0,i);
      bool swept = LSS.isSwept(0.0,i);
      //if(subId[i]==2) fprintf(stderr,"my subId = %d, swept = %d, Id = %d.\n", subId[i], swept, Id);

      if(Id==0) { // not isolated by structure. need to consider level-set
        if(swept) {
          subId[i] = 0;
          for(int k=0; k<dim; k++) {
            if(subPhiV[i][k]>0.0) {
	      if (programmedBurn && (programmedBurn->isUnburnedEOS(k+1,burnTag) ||
				     programmedBurn->isBurnedEOS(k+1,burnTag)) ) {
		if (programmedBurn->nodeInside(burnTag,iSub,i) ||
                    programmedBurn->isFinished(burnTag)) {
                  //fprintf(stderr,"Inside updateFluidIdFS, I am here!\n");
		  subId[i] = programmedBurn->getBurnedEOS(burnTag);}
		else
		  subId[i] = programmedBurn->getUnburnedEOS(burnTag);
		break;
	      }
	      else {
		subId[i] = k+1;
		break;
	      }
            }
	  }
	} 
      } else // isolated by structure. Id determined by intersector
        subId[i] = Id;
    }


    //for(int i=0; i<subPhiV.size(); i++) {
    //  if(subId[i]==2) fprintf(stderr,"In FluidSelector::updateFluidIdFS(...), found Id = 2!\n");}
  }
}

//------------------------------------------------------------------------------

template<int dim> /*this dim is actually dimLS*/
void FluidSelector::updateFluidIdFS2(DistLevelSetStructure *distLSS, DistSVec<double,dim> &PhiV, DistSVec<bool,4> &pollp)
{
  // ------- Determine status for grid-points swept by FS interface -------
  // Rule No.1: If this grid point is "occluded", set its status to "numPhases".
  // Rule No.2: If its "visible && !occluded && !swept" neighbors have the same status, use this one. 
  // Rule No.3: Otherwise, consider the sign of "PhiV". (PhiV should have been "blurred".)
  //            KW(TODO)(03/06/2014): Really? Is it correct to use "PhiV" here?

  int burnTag, rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#pragma omp parallel for
  for (int iSub=0; iSub<PhiV.numLocSub(); ++iSub) {
    Vec<int> &subId((*fluidId)(iSub));
    SVec<double,dim> &subPhiV(PhiV(iSub));
    LevelSetStructure &LSS((*distLSS)(iSub));

    SVec<bool,4> &poll(pollp(iSub));

    for(int i=0; i<subPhiV.size(); i++) {
      bool swept = LSS.isSwept(0.0,i);
      bool occluded = LSS.isOccluded(0.0,i);
/*
      if(rank==31 && i==2402)
        fprintf(stderr,"Rank 31, i = %d, swept = %d, occluded = %d, fluidId = %d, dimLS = %d, poll = %d %d %d %d. LSS.numOfFluids = %d, subPhiV = %e.\n", i, swept, occluded, subId[i], dim, poll[i][0], poll[i][1], poll[i][2], poll[i][3], LSS.numOfFluids(), subPhiV[i][0]);
*/
      if(!swept) {//nothing to be done
        if(!occluded && subId[i]!=LSS.numOfFluids())
          continue;
      }

      //if (i == 4566)
      //  fprintf(stderr,"Rank %d, i = %d, swept = %d, occluded = %d, fluidId = %d, dimLS = %d, poll = %d %d %d %d. LSS.numOfFluids = %d, subPhiV = %e.\n", rank,i, swept, occluded, subId[i], dim, poll[i][0], poll[i][1], poll[i][2], poll[i][3], LSS.numOfFluids(), subPhiV[i][0]);

      if(occluded) { // Rule No.1
        subId[i] = LSS.numOfFluids();
        if(!poll[i][3]) {
          fprintf(stderr,"Rank %d, i = %d, swept = %d, occluded = %d, fluidId = %d, dimLS = %d, poll = %d %d %d %d. LSS.numOfFluids = %d, subPhiV = %e.\n", rank,i, swept, occluded, subId[i], dim, poll[i][0], poll[i][1], poll[i][2], poll[i][3], LSS.numOfFluids(), subPhiV[i][0]);
          fprintf(stderr,"TOO BAD TOO! %i %i %i %i\n",poll[i][0], poll[i][1],poll[i][2],poll[i][3]);
          DebugTools::PrintBacktrace();
        }
        continue;
      }

      int count = (int)poll[i][0] + (int)poll[i][1] + (int)poll[i][2] + (int)poll[i][3];
      switch (count) {
        case 0: //no info
          fprintf(stderr,"More than one layer of nodes are swept in one step: Rank %d, i = %d.\n", rank, i);
  //        LSS.forceOccluded(0.0,i);
  //        poll[i][3] = 1;
  //        subId[i] = LSS.numOfFluids(); 
          DebugTools::SpitRank();
          break;
        case 1: // Rule No.2
          for(int j=0; j<3; j++)
            if(poll[i][j]) subId[i] = j;
          if(poll[i][3]) fprintf(stderr,"WARNING: poll[occluded] = true. Could be a software bug.\n");
          break;
      }
    
      if(count==1) //already applied Rule No.2
        continue;

      // Rule No.3
//-------------------------
      //KW: I think this part is kind of arbitrary...
/*
      bool settled = false;
      int burnedEOS, unburnedEOS;
      if(programmedBurn) {
        for(int burnTag=0; burnTag<programmedBurn->numberOfBurns(); burnTag++) {
          burnedEOS = programmedBurn->getBurnedEOS(burnTag); 
          unburnedEOS = programmedBurn->getUnburnedEOS(burnTag);
          if(poll[i][burnedEOS] && poll[i][unburnedEOS]) {
            if(programmedBurn->nodeInside(burnTag,iSub,i) || programmedBurn->isFinished(burnTag))
              subId[i] = burnedEOS;
            else
              subId[i] = unburnedEOS;
            settled = true;
            break;
          } else if(poll[i][burnedEOS]) {
            if(programmedBurn->nodeInside(burnTag,iSub,i) || programmedBurn->isFinished(burnTag)) {
              subId[i] = burnedEOS;
              settled = true;
              break;
            } 
          }
        }
      } 
      if(settled)
        continue; 
*/
      if(dim==1) { // consider detonation inside structure
        if(subPhiV[i][0]>0.0 || (subPhiV[i][0]<=0.0 && poll[i][0]==0)) {
          if(programmedBurn && (programmedBurn->isUnburnedEOS(1,burnTag) || programmedBurn->isBurnedEOS(1,burnTag)) ) {
            if(programmedBurn->nodeInside(burnTag,iSub,i) || programmedBurn->isFinished(burnTag)) {
              subId[i] = programmedBurn->getBurnedEOS(burnTag);
              if(!poll[i][subId[i]]) {
                if(poll[i][0]) {
                  subId[i] = 0;
                  fprintf(stderr,"WARNING: Setting node tag to 0 based on poll.\n");
                } else {
                  fprintf(stderr,"WARNING: poll = %d %d %d %d, fluidId = %d.\n", 
                          poll[i][0], poll[i][1], poll[i][2], poll[i][3], subId[i]);
                }
              }
            } else {
              subId[i] = programmedBurn->getUnburnedEOS(burnTag);
              if(!poll[i][subId[i]]) {
                if(poll[i][0]) {
                  subId[i] = 0;
                  fprintf(stderr,"WARNING2: Setting node tag to 0 based on poll.\n");
                } else {
                  fprintf(stderr,"WARNING2: poll = %d %d %d %d, fluidId = %d.\n", 
                          poll[i][0], poll[i][1], poll[i][2], poll[i][3], subId[i]);
                }
              }
            }
          } else {
            subId[i] = 1;
          }
        } else
          subId[i] = 0;
      } else {
        bool done = false;
        for(int k=0; k<dim; k++)
          if(subPhiV[i][k]>0.0) {
            subId[i] = k+1;
            done = true;
            break;
          }
        if(!done)
          subId[i] = 0;
      }
//--------------------------------
    }
  } 
}

//------------------------------------------------------------------------------

template<int dim> /*this dim is actually dimLS*/
void FluidSelector::updateFluidIdFF(DistLevelSetStructure *distLSS, DistSVec<double,dim> &Phi)
{
  assert(dim<=numPhases-1);
  int numLocSub = Phi.numLocSub();
  int burnTag;
  DistVec<int> &fsId(distLSS->getStatus());
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int     *tag       = fluidId->subData(iSub);
    int     *fsid      = fsId.subData(iSub);
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      if(fsid[iNode]!=0) {
        if(fsid[iNode]!=tag[iNode]){
          fprintf(stderr,"This must be a bug!\n"); exit(-1);}
        continue;
      }
      tag[iNode] = 0;
      for(int i=0; i<dim; i++) {
	if(phi[iNode][i]>0.0) {
	  if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
				 programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
	    if (programmedBurn->nodeInside(burnTag,iSub,iNode) ||
                programmedBurn->isFinished(burnTag)) {
	      tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
              //fprintf(stderr,"I am here... \n");
            }
	    else {
	      tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
            }
	    break;
	  }
	  else{ 
	    tag[iNode] = i+1; break; 
	  }
	}
      }
    }
  }
}

//------------------------------------------------------------------------------
// This function allows fracture. It does not look at fsid.
template<int dim> /*this dim is actually dimLS*/
void FluidSelector::updateFluidIdFF2(DistLevelSetStructure *distLSS, DistSVec<double,dim> &Phi)
{
  int numLocSub = Phi.numLocSub();
  int burnTag;
  int iSub;
 
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim]     = Phi.subData(iSub);
    int     *tag           = fluidId->subData(iSub);
    LevelSetStructure &LSS = (*distLSS)(iSub);

    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++){
      if(LSS.isOccluded(0.0,iNode)) {
        //phi[iNode][dim-1] = 0.0;
        tag[iNode] = LSS.numOfFluids();
        continue;
      }
      tag[iNode] = 0;
      for(int i=0; i<dim; i++) {
        if(phi[iNode][i]>0.0) {
          if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
                                 programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
            if (/*programmedBurn->isIgnited(burnTag) &&*/
                (programmedBurn->nodeInside(burnTag,iSub,iNode) || programmedBurn->isFinished(burnTag)))
              tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
            else {
              tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
            }
            break;
          }
          else{
            tag[iNode] = i+1; break;
          }
        }
      }
    }

  }
}

//------------------------------------------------------------------------------

template<int dim> /*this dim is actually dimLS*/
void FluidSelector::checkLSConsistency(DistSVec<double,dim> &Phi)
{
  int numLocSub = Phi.numLocSub();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double (*phi)[dim] = Phi.subData(iSub);
    int *tag           = fluidId->subData(iSub);
    for(int i=0; i<Phi.subSize(iSub); i++) {
      if(tag[i]==0) {
        if(phi[i][dim-1]>0.0) {
          fprintf(stderr,"BUG: Inconsistency between fluidId (%d) and phi (%e). numPhases = %d.\n", tag[i], phi[i][dim-1], numPhases);
          exit(-1);}}
      else if(tag[i]==numPhases) {
        if(fabs(phi[i][dim-1])>1.0e-10) {
          fprintf(stderr,"BUG: Inconsistency between fluidId (%d) and phi (%e). numPhases = %d.\n", tag[i], phi[i][dim-1], numPhases);
          exit(-1);}}
      else if (tag[i] == numPhases-1) {
        if(fabs(phi[i][dim-1]<=0.0)) {
          fprintf(stderr,"BUG: Inconsistency between fluidId (%d) and phi (%e). numPhases = %d.\n", tag[i], phi[i][dim-1], numPhases);
          exit(-1);}}
    }
  }
}

//------------------------------------------------------------------------------

