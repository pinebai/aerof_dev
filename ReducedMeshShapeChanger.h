#ifndef _REDUCED_MESH_SHAPE_CHANGER_H_
#define _REDUCED_MESH_SHAPE_CHANGER_H_

#include <SubDomain.h>
#include <GappyPreprocessing.h>
#include <iostream>
#include <string>
template <int dim>
class ReducedMeshShapeChanger : public GappyPreprocessing<dim> {

	double **xyz;	// xyz[iReducedNodes][iXYZ] is the iXYZ coordinte of the iReducedNodes node

	void readReducedNodes(const char *);
	void readWriteTopFile(const char *);
	void fillXYZ();

	public:
	void buildReducedModel();	// build all offline info (do everything)
	ReducedMeshShapeChanger(Communicator *, IoData &, Domain &, DistGeoState *);
	~ReducedMeshShapeChanger();
};
#include "ReducedMeshShapeChanger.C"
#endif
