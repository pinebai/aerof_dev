/* Multigrid.h

 */

#pragma once

#include <Vector.h>

#include <IoData.h>

class Multigrid {

 public:

  Multigrid();

  ~Multigrid();

  static void createMappingFromMeshes(IoData& ioData);
  
  static void constructMapping(SVec<double,3>& Xfine,
			       SVec<double,3>& Xcoarse,
			       int** decFine, int nSubDFine,int* subDsizeFine,
			       int** decCoarse,int nSubDCoarse,int* subDsizeCoarse,
			       const char* packageFile,
			       const char* collectionFile,
			       double,int, double);


};
  
