#include <GeoState.h>

#include <GeoData.h>
#include <Vector3D.h>
#include <Vector.h>

//------------------------------------------------------------------------------

// Included (MB)
GeoState::GeoState(const GeoData &gd, Vec<double> &ctrlvoln, Vec<double> &ctrlvolnm1,
                   Vec<double> &ctrlvolnm2,
                   Vec<double> &d2w,
                   Vec<Vec3D> &edgenorm, Vec<Vec3D> &facenorm,
                   Vec<double> &edgenormvel, Vec<double> &facenormvel,
                   Vec<Vec3D> &inletnodenorm, Vec<int>& numfaceneighb,
                   Vec<Vec3D> &dEdgenorm, Vec<Vec3D> &dFacenorm,
                   Vec<double> &dEdgenormvel, Vec<double> &dFacenormvel,
                   SVec<double,3> &x, SVec<double,3> &dx) :
                   data(gd), ctrlVol_n(ctrlvoln), ctrlVol_nm1(ctrlvolnm1), ctrlVol_nm2(ctrlvolnm2),
                   d2wall(d2w), edgeNorm(edgenorm), faceNorm(facenorm),
                   edgeNormVel(edgenormvel), faceNormVel(facenormvel), inletNodeNorm(inletnodenorm),
                   numFaceNeighb(numfaceneighb)
{
                   dEdgeNorm=&dEdgenorm;
                   dFaceNorm=&dFacenorm;
                   dEdgeNormVel=&dEdgenormvel;
                   dFaceNormVel=&dFacenormvel;
                   X=&x;
                   dX=&dx;
}

GeoState::GeoState(const GeoData &gd, Vec<double> &ctrlvoln, Vec<double> &ctrlvolnm1,
		   Vec<double> &ctrlvolnm2,
                   Vec<double> &d2w,
		   Vec<Vec3D> &edgenorm, Vec<Vec3D> &facenorm,
		   Vec<double> &edgenormvel, Vec<double> &facenormvel,
		   Vec<Vec3D> &inletnodenorm, Vec<int>& numfaceneighb) :
  data(gd), ctrlVol_n(ctrlvoln), ctrlVol_nm1(ctrlvolnm1), ctrlVol_nm2(ctrlvolnm2),
  d2wall(d2w), edgeNorm(edgenorm), faceNorm(facenorm),
  edgeNormVel(edgenormvel), faceNormVel(facenormvel), inletNodeNorm(inletnodenorm),
  numFaceNeighb(numfaceneighb)
{

}

//------------------------------------------------------------------------------
