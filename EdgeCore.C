#include <Edge.h>
#include "LevelSet/LevelSetStructure.h"

#include <Vector.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::swap;
#endif
#define EDGE_LENGTH
#ifdef EDGE_LENGTH
#include <Vector3D.h>
#endif

//------------------------------------------------------------------------------

EdgeSet::EdgeSet()
{
  numEdges  = 0;
  numSampledEdges = 0;
  ptr       = 0;
  masterFlag= 0;
  mp        = new MapPair;
  sampleMesh = false;

  programmedBurn = 0;
  higherOrderMF = 0;
  higherOrderFSI = 0;

#ifdef EDGE_LENGTH  //HB
  edgeLength= 0;
#endif

  mfRiemannNormal = MF_RIEMANN_NORMAL_REAL;

  triangulatedLSS = NULL;
}

//------------------------------------------------------------------------------

EdgeSet::~EdgeSet()
{
  if (ptr)        delete [] ptr;
  if (masterFlag) delete [] masterFlag;
  if (mp)         delete mp;
#ifdef EDGE_LENGTH
  if (edgeLength) delete [] edgeLength;
#endif
}

//------------------------------------------------------------------------------

int EdgeSet::find(int first, int second)
{

  if (first > second) swap(first, second);

  MapPair::iterator v = mp->find(Pair(first, second));

  if (v == mp->end()) { 
    (*mp)[ Pair(first, second) ] = numEdges;
		++numSampledEdges;
    return numEdges++;
  }

  return (*v).second;

}

//------------------------------------------------------------------------------

int EdgeSet::findOnly(int first, int second) const
{

  if (first > second) swap(first, second);

  MapPair::iterator v = mp->find(Pair(first, second));

  if (v == mp->end()) {
    fprintf(stderr,"ERROR: Failed in finding an edge connecting %d and %d. (local index).\n", first, second);
    exit(-1);
  }

  return (*v).second;

}

//-----------------------------------------------------------

void EdgeSet::createPointers(Vec<int> &newNum)
{

  ptr = new int[numEdges][2];
  edgeLength = new double[numEdges];

  MapPair::iterator it = mp->begin();
  MapPair::iterator last = mp->end();

  int l = 0;

  while (it != last) {

    // XML I am not sure why one would need to renumber...
    ptr[l][0] = (*it).first.first;
    ptr[l][1] = (*it).first.second;
    edgeLength[l] = 0.0;

    newNum[(*it).second] = l;

    (*it).second = l;

    ++l;
    ++it;

  }

}
//------------------------------------------------------------------------------

#ifdef EDGE_LENGTH //HB: compute the edges'length for given nodes'coordinates
void EdgeSet::updateLength(SVec<double,3>& X)
{
  if(!edgeLength) edgeLength = new double[numEdges];
 
  for(int iedge=0; iedge<numEdges; iedge++) {
    Vec3D d(X[ptr[iedge][0]]);
    Vec3D e(X[ptr[iedge][1]]);
    Vec3D ed = d-e;
    edgeLength[iedge] = ed.norm();
  }
}
#endif

//------------------------------------------------------------------------------
int EdgeSet::checkReconstructedValues(int i, int j, double *Vi, double *Vj, VarFcn *vf, 
				int *locToGlobNodeMap, int failsafe, SVec<int,2>& tag,
                                double *originalVi, double *originalVj,
												  int IDi, int IDj, bool iAct, bool jAct)
{
//proceed to checking positivity of certain quantities required in the computation of the fluxes for both nodes of an edge.
//these quantities are the ones involved in the computation of the sound speed (most often pressure and density but not
//always, see Tait for example).

  int ierr = 0;

  if(iAct)
  {
	  if(vf->checkReconstructedValues(Vi, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1, IDi, IDj, failsafe, IDi))
	  {
    ++ierr;
		  if(failsafe) tag[i][0] = 1;
	  }
  }

  if(jAct)
  {
	  if(vf->checkReconstructedValues(Vj, locToGlobNodeMap[j]+1, locToGlobNodeMap[i]+1, IDj, IDi, failsafe, IDj))
	  {
    ++ierr;
		  if(failsafe) tag[j][0] = 1;
	  }
  }

  return ierr;

}

//------------------------------------------------------------------------------

#ifdef EDGE_LENGTH
void EdgeSet::computeCharacteristicEdgeLength(SVec<double,3> &X, double &minLength, double &aveLength, double &maxLength, int &numInsideEdges, const double xmin, const double xmax, const double ymin, const double ymax, const double zmin, const double zmax)
{
  if (!this->edgeLength) this->updateLength(X);
  numInsideEdges = 0;
  bool inside = true, start = true;
  for (int iEdge=0; iEdge<numEdges; iEdge++) {
    double* position1 = X[ptr[iEdge][0]];
    double* position2 = X[ptr[iEdge][1]];
    if ( (position1[0] < xmin || position1[0] > xmax || position1[1] < ymin || position1[1] > ymax || position1[2] < zmin || position1[2] > zmax ) && (position2[0] < xmin || position2[0] > xmax || position2[1] < ymin || position2[1] > ymax || position2[2] < zmin || position2[2] > zmax ) )
      inside = false; else inside = true;
    if (inside == true) {
      if (start == true) {
        minLength = length(iEdge);  aveLength = length(iEdge);  maxLength = length(iEdge);
        numInsideEdges = 1; start = false;
      }
      else {
        aveLength += length(iEdge);  
        if (length(iEdge)<minLength) minLength = length(iEdge);
        if (length(iEdge)>maxLength) maxLength = length(iEdge);
        numInsideEdges++;
      }
    } 
  } 
  aveLength /= numInsideEdges;
}
#endif
//---------------------------------------------------------------------------------

void EdgeSet::attachProgrammedBurn(ProgrammedBurn* p) {

  programmedBurn = p;
}

void EdgeSet::attachHigherOrderMultiFluid(HigherOrderMultiFluid* mf) {

  higherOrderMF = mf;
}

void EdgeSet::attachHigherOrderFSI(HigherOrderFSI* fsi) {

  higherOrderFSI = fsi;
}

//------------------------------------------------------------------------------

void EdgeSet::computeConnectedEdges(const std::vector<int> &locSampleNodes) 
{

	sampleMesh = true;
  edgesConnectedToSampleNode.clear();
  for(int l=0; l<numEdges; l++){
		//compute which nodes are attached
    int i = ptr[l][0];
    int j = ptr[l][1];
		// check if either node is a sample node
		for (int iSampleNode = 0; iSampleNode < locSampleNodes.size(); ++iSampleNode) {
			if (i == locSampleNodes[iSampleNode] || j == locSampleNodes[iSampleNode]) {
				edgesConnectedToSampleNode.push_back(l);
				break;
			}
		}
	}
	numSampledEdges = edgesConnectedToSampleNode.size();
}

//------------------------------------------------------------------------------

void EdgeSet::computeGlobalConnectedEdges(const std::vector<int> &globalNeighborNodes, const int *locToGlobNodeMap) 
{

	sampleMesh = true;
  edgesTwoLayersSampleNode.clear();
  for(int l=0; l<numEdges; l++){
		//compute which nodes are attached
    int i = ptr[l][0];
    int j = ptr[l][1];
		// check if either node is a sample node
		for (int iSampleNode = 0; iSampleNode < globalNeighborNodes.size(); ++iSampleNode) {
			if (locToGlobNodeMap[i] == globalNeighborNodes[iSampleNode]  ||
				 	locToGlobNodeMap[j] == globalNeighborNodes[iSampleNode] ) {
				edgesTwoLayersSampleNode.push_back(l);
				break;
			}
		}
	}
	numTwoLayerEdges = edgesTwoLayersSampleNode.size();
}

void EdgeSet::setMultiFluidRiemannNormal(MultifluidRiemannNormal m) {

  mfRiemannNormal = m;
}

void EdgeSet::attachTriangulatedInterfaceLSS(LevelSetStructure* LSS) {

  triangulatedLSS = LSS;
}
