#ifndef __REINITIALIZE_DISTANCE_TO_WALL_H__
#define __REINITIALIZE_DISTANCE_TO_WALL_H__

#include <DistVector.h>
#include <LevelSet/LevelSetStructure.h>
#include <Vector.h>

class Domain;
class SubDomain;
class DistGeoState;
class DistLevelSetStructure;
class GeoState;
class LevelSetStructure;

template<int dimLS>
class ReinitializeDistanceToWall
{
  IoData& iod;
  Domain& dom;
  DistVec<bool> done;
  DistSVec<double,1> d2wall;
  DistVec<int> sortedNodes;
  int* nSortedNodes,*firstCheckedNode;

  DistVec<int> tag;
  DistSVec<double,dimLS> dummyPhi;

public:
  ReinitializeDistanceToWall(IoData &ioData, Domain& domain);
  ~ReinitializeDistanceToWall();

  void ComputeWallFunction(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState);
  void DistanceToClosestPointOnMovingStructure(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState);
  void PrescribedValues(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState);
  void GetLevelsFromInterfaceAndMarchForward(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState);
  void PseudoFastMarchingMethod(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState,int iterativeLevel);
  void computeExactErrors(DistLevelSetStructure& LSS,DistSVec<double,3>& X,DistGeoState& distGeoState);

private:
  void InitializeWallFunction(SubDomain& subD,LevelSetStructure& LSS,Vec<bool>& done,SVec<double,3>& X,SVec<double,1>& d2w,Vec<int>& tag);
};

#endif
