#include <CurvatureDetection.h>
#include <Domain.h>
#include <SubDomain.h>
#include <Vector3D.h>
#include <DistVector.h>

//------------------------------------------------------------------------------

CurvatureDetection::CurvatureDetection(Domain* domain)
{

  numLocSub = domain->getNumLocSub();
  subDomain = domain->getSubDomain();
  tag = new DistVec<double>(domain->getEdgeDistInfo());
  normals = new DistSVec<double,6>(domain->getEdgeDistInfo());
  vec1 = domain->getVolPat();
  vec6 = new CommPattern<double>(domain->getSubTopo(), domain->getCommunicator(),
				 CommPattern<double>::CopyOnSend);
#pragma omp parallel for
  for (int iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setComLenEdges(6, *vec6);
  vec6->finalize();
  
}

//------------------------------------------------------------------------------

CurvatureDetection::~CurvatureDetection()
{

  if (tag) delete tag;
  if (normals) delete normals;
  if (vec6) delete vec6;

}

//------------------------------------------------------------------------------

void CurvatureDetection::compute(double threshold, DistSVec<double,3>& X, DistVec<bool>& t)
{

  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->computeFaceEdgeNormals(X(iSub), (*normals)(iSub));
    subDomain[iSub]->sndEdgeData(*vec6, normals->subData(iSub));
  }
  vec6->exchange();
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->addRcvEdgeData(*vec6, normals->subData(iSub));
    subDomain[iSub]->computeEdgeDihedralAngle(threshold, (*normals)(iSub), (*tag)(iSub));
    subDomain[iSub]->sndData(*vec1, reinterpret_cast<double (*)[1]>(tag->subData(iSub)));
  }
  vec1->exchange();
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*vec1, reinterpret_cast<double (*)[1]>(tag->subData(iSub)));
    subDomain[iSub]->propagateInfoAlongEdges((*tag)(iSub));
    subDomain[iSub]->sndData(*vec1, reinterpret_cast<double (*)[1]>(tag->subData(iSub)));
  }
  vec1->exchange();
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*vec1, reinterpret_cast<double (*)[1]>(tag->subData(iSub)));
    double* _tag = tag->subData(iSub);
    bool* _t = t.subData(iSub);
    for (int i=0; i<t.subSize(iSub); ++i) {
      if (_tag[i] > 0.0)
	_t[i] = true;
    }
  }

}

//------------------------------------------------------------------------------
