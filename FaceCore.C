#include <Face.h>

#include <RefVal.h>
#include <BcDef.h>
#include <Edge.h>
#include <Vector3D.h>
#include <Vector.h>
#include <BinFileHandler.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>

#ifdef OLD_STL
#include <defalloc.h>
#include <algo.h>
#else
#include <algorithm>
using std::stable_sort;
using std::min;
using std::max;
using std::swap;
#endif

//------------------------------------------------------------------------------

Face::Face()
{
  higherOrderMF = NULL;
  //ffWeight = 1.0;  //for nonlinear ROM simulations requiring additional weighting for far field BCs
}

//------------------------------------------------------------------------------

void Face::setup(int fc, int *nn, int nnum, int sid)
{

  if (fc == BC_OUTLET_MOVING || fc == BC_OUTLET_FIXED || 
      fc == BC_INLET_MOVING || fc == BC_INLET_FIXED ||
      fc == BC_DIRECTSTATE_OUTLET_MOVING || fc == BC_DIRECTSTATE_OUTLET_FIXED || 
      fc == BC_DIRECTSTATE_INLET_MOVING || fc == BC_DIRECTSTATE_INLET_FIXED ||
      fc == BC_MASSFLOW_OUTLET_MOVING || fc == BC_MASSFLOW_OUTLET_FIXED || 
      fc == BC_MASSFLOW_INLET_MOVING || fc == BC_MASSFLOW_INLET_FIXED ||
      fc == BC_ADIABATIC_WALL_MOVING || fc == BC_ADIABATIC_WALL_FIXED ||
      fc == BC_SLIP_WALL_MOVING || fc == BC_SLIP_WALL_FIXED ||
      fc == BC_POROUS_WALL_MOVING || fc == BC_POROUS_WALL_FIXED ||
      fc == BC_ISOTHERMAL_WALL_MOVING || fc == BC_ISOTHERMAL_WALL_FIXED ||
      fc == BC_SYMMETRY || fc == BC_KIRCHHOFF_SURFACE)
  {
    code = fc;
  }
  else {
    fprintf(stderr, "*** Error: incorrect boundary face code (%d)\n", fc);
    exit(1);
  }

  for (int l=0; l<numNodes(); ++l)
    nodeNum(l) = nn[l];
  
  surface_id = sid;
  normNum = nnum;

  // for HH Farfield BC
  faceCenter[0] = faceCenter[1] = faceCenter[2] = 0.0;
}

//------------------------------------------------------------------------------

void Face::setType(int *facemap)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)

  if ((code >= BC_MIN_CODE) && (code <= BC_MAX_CODE))
    code = facemap[code];

}

//------------------------------------------------------------------------------

void Face::setNodeType(int* priority, int* nodeType)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)

  if ((code >= BC_MIN_CODE) && (code <= BC_MAX_CODE))
  {
    for (int j=0; j<numNodes(); ++j)
      if (priority[code] > priority[ nodeType[ nodeNum(j) ] ])
        nodeType[ nodeNum(j) ] = code;
  }

}

//------------------------------------------------------------------------------

void Face::setNodeFaceType(int* nodeFaceType)
{
  int nft;

  if( code==BC_INLET_FIXED  || code==BC_INLET_MOVING  ||
      code==BC_OUTLET_FIXED || code==BC_OUTLET_MOVING ||
      code==BC_DIRECTSTATE_INLET_FIXED  || code==BC_DIRECTSTATE_INLET_MOVING ||
      code==BC_DIRECTSTATE_OUTLET_FIXED || code==BC_DIRECTSTATE_OUTLET_MOVING ||
      code==BC_MASSFLOW_INLET_FIXED  || code==BC_MASSFLOW_INLET_MOVING ||
      code==BC_MASSFLOW_OUTLET_FIXED || code==BC_MASSFLOW_OUTLET_MOVING ) {
    for (int j=0; j<numNodes(); j++){
      nft = nodeFaceType[nodeNum(j)];
      if(nft == 0)
        nodeFaceType[nodeNum(j)] = 1;
      if(nft == -1)
        nodeFaceType[nodeNum(j)] = 2;
    }
  } 
  else if ((code >= BC_MIN_CODE) && (code <= BC_MAX_CODE))
  {
    //only possibilities left should be wall and symmetry
    for (int j=0; j<numNodes(); j++)
    {
      nft = nodeFaceType[nodeNum(j)];
      if(nft == 0)
        nodeFaceType[nodeNum(j)] = -1;
      if(nft == 1)
        nodeFaceType[nodeNum(j)] = 2;
    } // for (int j=0; j<numNodes(); j++)
  }

}

//------------------------------------------------------------------------------

void Face::setElementNumber(int num, int rotDir)
{ 
  elemNum = num;

  if (rotDir < 0) {
    int i, j, nndo2 = numNodes()/2;
    for (i=0, j=numNodes()-1; i<nndo2; ++i, --j) {
      int nodeNum_j = nodeNum(j);
      nodeNum(j) = nodeNum(i);
      nodeNum(i) = nodeNum_j ;
    }
  }

}

//------------------------------------------------------------------------------

void Face::tagNodesOnBoundaries(Vec<bool> &tagNodes)
{
  for (int j=0; j<numNodes(); j++)
    tagNodes[ nodeNum(j) ] = true;
}

//------------------------------------------------------------------------------

void Face::tagEdgesOnBoundaries(Vec<bool> &tagEdges)
{
  for (int j=0; j<numNodes(); j++)
    tagEdges[ edgeNum(j) ] = true;
}

//------------------------------------------------------------------------------

void Face::reorder()
{ 
  int *nN = nodes();
#ifdef OLD_STL
  sort(&nN[0], &nN[numNodes()]); 
#else
  stable_sort(&nN[0], &nN[numNodes()]); 
#endif
}

//------------------------------------------------------------------------------

void Face::numberEdges(EdgeSet &edges)
{
  
  int numEdges = edges.size();
  
  for (int j=0; j<numNodes(); j++)
    edgeNum(j) = edges.find(nodeNum( edgeEnd(j,0) ), nodeNum( edgeEnd(j,1) ));

  if (edges.size() != numEdges) {
    fprintf(stderr, "*** Error: could not find face edges\n");
    exit(1);
  }

}

//------------------------------------------------------------------------------
// WARNING: THIS IS A FUNCTION FOR TFACES ONLY
double Face::computeVolume(Vec3D &xa_n, Vec3D &xb_n, Vec3D &xc_n, 
			   Vec3D &xa_np1, Vec3D &xb_np1, Vec3D &xc_np1)
{
  const static double sixth = 1.0/6.0;
  const static double third = 1.0/3.0;

  Vec3D dxa = xa_np1 - xa_n;
  Vec3D dxb = xb_np1 - xb_n;
  Vec3D dxc = xc_np1 - xc_n;
  Vec3D xac_n = xa_n - xc_n;
  Vec3D xbc_n = xb_n - xc_n;
  Vec3D xac_np1 = xa_np1 - xc_np1;
  Vec3D xbc_np1 = xb_np1 - xc_np1;
  
  Vec3D eta = sixth * ((xac_n ^ xbc_n) + (xac_np1 ^ xbc_np1) + 
		       0.5 * ((xac_n ^ xbc_np1) + (xac_np1 ^ xbc_n)));
  
  return third * (dxa + dxb + dxc) * eta;

}

//------------------------------------------------------------------------------

void Face::computeEdgeNormals(SVec<double,3>& X, int* l2gl, SVec<double,6>& normals)
{

  Vec3D faceNorm;
  computeNormal(X, faceNorm);
  faceNorm *= -1.0 / sqrt(faceNorm*faceNorm);

  for (int l=0; l<numNodes(); ++l) {
    bool ori;
    if (l2gl[ nodeNum( edgeEnd(l,0) ) ] < l2gl[ nodeNum( edgeEnd(l,1) ) ])
      ori = true;
    else
      ori = false;

    if (ori) {
      normals[ edgeNum(l) ][0] = faceNorm[0];
      normals[ edgeNum(l) ][1] = faceNorm[1];
      normals[ edgeNum(l) ][2] = faceNorm[2];
    }
    else {
      normals[ edgeNum(l) ][3] = faceNorm[0];
      normals[ edgeNum(l) ][4] = faceNorm[1];
      normals[ edgeNum(l) ][5] = faceNorm[2];
    }
  }

}

//------------------------------------------------------------------------------
FaceSet::FaceSet(int value)
{
  numFaces = value;
  faces = new Face*[value];
	sampleMesh = false;

  // Set total number of face normals to 0: 
  // it will be incremented when reading faces
  numFaceNorms = 0;

  higherOrderMF = NULL;
}

//------------------------------------------------------------------------------
FaceSet::~FaceSet()
{
  delete [] faces;
}

//------------------------------------------------------------------------------
int FaceSet::read(BinFileHandler &file, int numRanges, int (*ranges)[2], int *map)
{
  // read in number of faces in cluster (not used)
  int numClusFaces;
  file.read(&numClusFaces, 1);

  // read in the offset for the first face
  BinFileHandler::OffType start, tocStart;
  file.read(&start, 1);

  // set the location of first face in the TOC
  tocStart = file.tell();

  // intialize counters
  int count = 0;
  int countNorm = 0;

  // read in ranges
  int nSymm = 0;
  for (int iRange = 0; iRange < numRanges; ++iRange) {

    // compute number of faces in range
    int nFaces = ranges[iRange][1] - ranges[iRange][0] + 1;

    // seek to correct position in file thanks to the table of contents
    file.seek(tocStart + sizeof(int) * ranges[iRange][0]);
    int toc;
    file.read(&toc, 1);
    file.seek(start + sizeof(int) * toc);
    for (int i = 0; i < nFaces; i++)  {
      int type;
      int nodeNum[Face::MaxNumNd];
      int code;
      int surface_id;

      // read in the face type
      file.read(&type, 1);	

      // create face in faceset depending on type
      switch (type) {
      case Face::TRIA:
	addFace(count, new(memFaces) FaceTria);
	break;

      default:
	fprintf(stderr, "Error: Could not find face type %d\n", type);
	exit(1);
      }

      // read in global face number
      file.read(&map[count], 1);

      // read in face nodes
      file.read( nodeNum, faces[count]->numNodes());

      // read in the face code
      file.read(&code, 1);
      if (code == 6) {
        code = BC_SYMMETRY;
        nSymm++;
      }

      // read in the surface id
      file.read(&surface_id, 1);
      
      // setup face
      faces[count]->setup(code, nodeNum, countNorm, surface_id);

      // count total number of normals to be stored
      countNorm += faces[count]->numNorms();

      // count number of faces
      count++;
    }
  }

  // total number of normals to be stored
  numFaceNorms = countNorm;

  if (count != numFaces) {
    fprintf(stderr, "*** Error: wrong number of faces read (%d instead of %d)\n",
	    count, numFaces);
    exit(1);
  }

  return numClusFaces;
}

//------------------------------------------------------------------------------

void FaceSet::computeConnectedFaces(const std::vector<int> &locSampleNodes) 
{

	sampleMesh = true;
  facesConnectedToSampleNode.clear();
	int nSampleNode = locSampleNodes.size();
  for(int l=0; l<numFaces; ++l) {
		bool connectedFace = false;
		for (int iNode = 0; iNode < (*faces[l]).numNodes(); ++iNode) {
			for (int iSampleNode = 0; iSampleNode < nSampleNode; ++iSampleNode) {
				if ((*faces[l])[iNode] == locSampleNodes[iSampleNode]) {
					facesConnectedToSampleNode.push_back(l);
					connectedFace = true;
					break;
				}
			}
			if (connectedFace)
				break;
		}
	}

	numSampledFaces = facesConnectedToSampleNode.size();
	int tmp;
}

void FaceSet::attachHigherOrderMF(class HigherOrderMultiFluid* mf) {

  higherOrderMF = mf;
  for(int l=0; l<numFaces; ++l) {
    faces[l]->attachHigherOrderMF(mf);
  }
}
