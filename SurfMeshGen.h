#ifndef _SURF_MESH_GEN_H_
#define _SURF_MESH_GEN_H_

#include <GappyPreprocessing.h>
template <int dim>
class SurfMeshGen : public GappyPreprocessing<dim> {

	void setUpGreedy() { ;}

	void determineSampleNodes() { ;}

	void computePseudoInverse() { ;}

	void assembleOnlineMatrices() { ;}

	void outputOnlineMatrices() { ;}

	void outputSampleNodes() { ;}

//	void addSampleNodesAndNeighbors() { ;}

	public:
	SurfMeshGen(Communicator *, IoData &, Domain &, DistGeoState *);
};
#include "SurfMeshGen.C"
#endif
