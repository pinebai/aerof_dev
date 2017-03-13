#ifdef OLD_STL
#include <defalloc.h>
#include <algo.h>
#else
#include <algorithm>
using std::stable_sort;
using std::lower_bound;
#endif

#include <Communicator.h>
#include <Connectivity.h>

SubDTopo::SubDTopo(int myCPU, Connectivity *subToSub, Connectivity *CPUToSub)
{
 int iCPU, iSub;

 numCPU = CPUToSub->csize();
 cpuNum = myCPU;
 int glNumSub  = subToSub->csize();
 localNumSub = CPUToSub->num(myCPU);

 // fprintf(stderr, "We got %d CPUs\n", numCPU);
 // list of local subdomains in global numbering
 localSubToGl = (*CPUToSub)[myCPU];

 glSubToLocal = new int[glNumSub];

 // build the map from global to local.
 for(iSub = 0; iSub < glNumSub; ++iSub)
   glSubToLocal[iSub] = -1;
 
 // fprintf(stderr, "Got %d subs and %d globals\n", localNumSub, glNumSub);
 for(iSub = 0; iSub < localNumSub; ++iSub)
   glSubToLocal[localSubToGl[iSub]] = iSub;

 // Map from global number to local number
 glSubToCPU = new int[glNumSub];

 for(iSub = 0; iSub < glNumSub; ++iSub)
   glSubToCPU[iSub] = -1;

 for(iCPU = 0; iCPU < numCPU; ++iCPU)
   for(iSub = 0; iSub < CPUToSub->num(iCPU); ++iSub)
     glSubToCPU[(*CPUToSub)[iCPU][iSub]] = iCPU;

 //fprintf(stderr, "Going to make pairs\n");
 makeLocalPairIndex(subToSub);
 //fprintf(stderr, "Done\n");
 // create the cross-CPU send and receive connectivities if necessary.
 if(numCPU > 1) {
   makeCrossConnect();
 } else{
   cpuSC = cpuRC = 0;
 }
}

SubDTopo::~SubDTopo()
{
  delete[] glSubToLocal;
  delete[] glSubToCPU;
  delete[] allPairs;
}
/******************************************************************************
  makeConnect creates the send and receive connectivities for neighboring
  CPUs 
 ******************************************************************************/
void
SubDTopo::makeCrossConnect()
{
 int iCPU, iPair;
 // First step, find all the neighboring CPUs
 int *glToLocCPU = new int[numCPU];
 for(iCPU = 0; iCPU < numCPU; ++iCPU)
   glToLocCPU[iCPU] = -1;

 int numNeighbCPU = 0;
 for(iPair = 0; iPair < numPairs; ++iPair) {
    int cpu;
    if((cpu = allPairs[iPair].cpuID) >= 0) {
      if(glToLocCPU[cpu] < 0) glToLocCPU[cpu] = numNeighbCPU++;
    }
 }
 neighbCPU =  new int[numNeighbCPU];
 for(iCPU = 0; iCPU < numCPU; ++iCPU)
   if(glToLocCPU[iCPU] >= 0)
     neighbCPU[glToLocCPU[iCPU]] = iCPU;

 // Now count the send and receives (actually they should be the same for
 // a symmetric communication
 int *sndPtr = new int[numNeighbCPU+1];
 int *rcvPtr = new int[numNeighbCPU+1];
 for(iCPU = 0; iCPU < numNeighbCPU+1; ++iCPU)
   sndPtr[iCPU] = rcvPtr[iCPU] = 0;

 for(iPair = 0; iPair < numPairs; ++iPair) {
    int cpu;
    if((cpu = allPairs[iPair].cpuID) >= 0) {
     if(glSubToCPU[allPairs[iPair].from] == cpuNum) // send pair
       sndPtr[glToLocCPU[cpu]]++;
     else // receive
       rcvPtr[glToLocCPU[cpu]]++;
    }
 }
 // Make the actual pointers (shifted for easy fill in)
 for(iCPU = 0; iCPU < numNeighbCPU; ++iCPU) {
   rcvPtr[iCPU+1] += rcvPtr[iCPU];
   sndPtr[iCPU+1] += rcvPtr[iCPU];
 }

 int *sndTrgt = new int[sndPtr[numNeighbCPU]];
 int *rcvTrgt = new int[rcvPtr[numNeighbCPU]];
 // Final fill in
 for(iPair = 0; iPair < numPairs; ++iPair) {
    int cpu;
    if((cpu = allPairs[iPair].cpuID) >= 0) {
     if(glSubToCPU[allPairs[iPair].from] == cpuNum) { // send pair
       sndPtr[glToLocCPU[cpu]]--;
       sndTrgt[sndPtr[glToLocCPU[cpu]]] = iPair;
     } else { // receive
       rcvPtr[glToLocCPU[cpu]]--;
       rcvTrgt[rcvPtr[glToLocCPU[cpu]]] = iPair;
     }
    }
 }
 cpuSC = new Connectivity(numNeighbCPU, sndPtr, sndTrgt);
 cpuRC = new Connectivity(numNeighbCPU, rcvPtr, rcvTrgt);
}

/******************************************************************************
   makeLocalPairIndex associate with every communication pair local to this
   subdomain:
      - A unique tag in the communication table 
   this index is the most convenient way to get access to data in send or
   receive.
   There is a different tag for the send (from a to b) and for receive
   (from b to a)
 *****************************************************************************/
void
SubDTopo::makeLocalPairIndex(Connectivity *subToSub)
{
  int iSub, jSub;
  // we just build a list of tripples from, to, cpu and sort them;
  // We create the tripples by
  // looping over the local subdomains subI and then loop over its neighbors
  // subJ. We add the tripple (subI, subJ, -1) if subJ is local.
  // (subI, subJ, cpu[subJ) if not. We also
  // add its symmetric (subJ, subI, [subJ]) if subJ is not in the current CPU

  //Step one count the number of tripples.
  numPairs = 0;
  for(iSub = 0; iSub < localNumSub; ++iSub) {
     int subI = localSubToGl[iSub];
     for(jSub = 0; jSub < subToSub->num(subI); ++jSub) {
       int subJ = (*subToSub)[subI][jSub];
       if(subI == subJ) continue;
       numPairs += 1; // internal pair
       if(glSubToLocal[subJ] < 0) numPairs +=2;  // external pair
     }
  }

  // create the (unsorted) pair list
  allPairs = new CPair[numPairs];
  numPairs = 0;
  for(iSub = 0; iSub < localNumSub; ++iSub) {
     int subI = localSubToGl[iSub];
     for(jSub = 0; jSub < subToSub->num(subI); ++jSub) {
       int subJ = (*subToSub)[subI][jSub];
       if(subI == subJ) continue;
       if(glSubToLocal[subJ] < 0) {   // external pair
          allPairs[numPairs++] = CPair(subI, subJ, glSubToCPU[subJ]);
          allPairs[numPairs++] = CPair(subJ, subI, glSubToCPU[subJ]);
       } else
          allPairs[numPairs++] = CPair(subI,subJ, -1); // internal pair
     }
  }
   
  // for easy finding, we now sort the pairs (a,b) < (c,d) iif a<c or
  // a == c &&  b < d
#ifdef OLD_STL
  sort(allPairs, allPairs+numPairs);
#else
  stable_sort(allPairs, allPairs+numPairs);
#endif

}

int
SubDTopo::getChannelID(int glFrom, int glTo)
{
 int cpu = -1;
 if(glSubToLocal[glFrom] < 0) cpu = glSubToCPU[glFrom];
 if(glSubToLocal[glTo  ] < 0) cpu = glSubToCPU[glTo  ];
 CPair cp(glFrom, glTo, cpu);
 CPair *found = lower_bound(allPairs, allPairs+numPairs, cp);
 int ID = found-allPairs;
 if(ID >= numPairs) return -1; // Did not find this pair
 if(*found == cp) return ID;
 return -1; // Did not find this pair
}

int
SubDTopo::reverseChannel(int channel)
{
 return getChannelID( allPairs[channel].to, allPairs[channel].from );
}
/*
Communicator::Communicator(int CPU, Connectivity *subToSub,
                                        Connectivity *CPUToSub)
{
 refTopo = new SubDTopo(CPU, subToSub, CPUToSub);
}
*/

