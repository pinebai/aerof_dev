#include <InletNode.h>

#include <Node.h>
#include <Vector.h>
#include <Vector3D.h>
#include <Extrapolation.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>

//------------------------------------------------------------------------------

InletNodeSet::InletNodeSet()
{
  locSubNum = 0;
  numInletNodes = 0;
  inletNodes = 0;
}

//------------------------------------------------------------------------------

InletNodeSet::~InletNodeSet()
{
  if(inletNodes) delete [] inletNodes;
}

//------------------------------------------------------------------------------

void InletNodeSet::setup(int loc, int num, IoData& ioData)
{
  locSubNum = loc;
  numInletNodes = num;
  inletNodes = new InletNode[num];
}

//------------------------------------------------------------------------------

void InletNodeSet::checkInletNodes(int* locToGlobNodeMap)
{
  fprintf(stderr, "#### total # of inletnodes is %d\n", numInletNodes);	
  for (int i = 0; i<numInletNodes; ++i)
    inletNodes[i].checkInletNodes(i, locToGlobNodeMap);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

InletNode::InletNode()
{
	node = -1;
	numTets = -1;
	numFaces = -1;
	tets = 0;
	tets2 = 0;
	faces = 0;
	master = false;
}

//------------------------------------------------------------------------------

InletNode::~InletNode()
{
	if (tets) delete [] tets;
	if (tets2) delete [] tets2;
	if (faces) delete [] faces;	
}

//------------------------------------------------------------------------------

void InletNode::setInletNode(int i, int fnum, int tnum, int* flist, int* tlist)
{
	node = i;
	numTets = tnum;
	numFaces = fnum;
	tets = tlist;
	faces = flist;
}

//------------------------------------------------------------------------------

void InletNode::addSecondTets( int* listtets)
{
        tets2 = listtets;
}
//------------------------------------------------------------------------------

void InletNode::chooseExtrapolation(int dim, double u, int* tet, int* mask, double* V1, double* V2, double* V){

  if (tet[0]==-1) {
    for(int i=0; i<dim; i++)
      V[i] = V1[i];
    return;
  }

  int type;
  if(u==0.0) type = 3;
  else if(u>0) type = 2;
  else type = 1;

  bool linear = true;

  for (int j=1; j<3; j++)
    if (mask[tet[j]]!=type && mask[tet[j]]!=0)
      linear = false;

  if(linear){
    for (int i=0; i<dim; i++)
      V[i]=V2[i];
    for (int j=0; j<3; j++)
      mask[tet[j]] = type;
  }
  else{
    for (int i=0; i<dim; i++)
      V[i]=V1[i];
  }
}

//------------------------------------------------------------------------------
 
void InletNode::checkInletNodes(int i, int* locToGlobNodeMap)
{
  fprintf(stderr, "------------------------------------------\n");
  fprintf(stderr, "-       node %d %d (master = %d)\n", locToGlobNodeMap[node]+1, node, master);
  fprintf(stderr, "# of tets = %d --- # of faces = %d\n", numTets, numFaces);
  for (int i=0; i<numTets; i++)
    fprintf(stderr, "%d  ", tets[i]);
//  fprintf(stderr, "\n----\n");
//  for (int j=0; j<numTets; j++)
//    fprintf(stderr, "%d  ", tets2[j]);
//  fprintf(stderr, "\n----\n");
//  for (int k=0; k<numFaces; k++)
//    fprintf(stderr, "%d  ", faces[k]);
  fprintf(stderr, "\n-----------------------------------------\n");

}

//------------------------------------------------------------------------------
