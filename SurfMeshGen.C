#include <SurfMeshGen.h>

template<int dim>
SurfMeshGen<dim>::SurfMeshGen(Communicator *_com, IoData &_ioData, Domain &dom, DistGeoState *_geoState) : 
GappyPreprocessing<dim>(_com, _ioData, dom, _geoState) {

  _com->fprintf(stdout, "\n\nBuilding surface mesh quantities\n");

  this->surfaceMeshConstruction = true;
	this->nSampleNodes = 1;
	// by default, include all lift and drag faces in the surface mesh
	if (this->includeLiftFaces == 0)
	 	this->includeLiftFaces = 2;
}
