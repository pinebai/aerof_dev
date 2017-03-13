#include <IoData.h>
#include <EdgeGrad.h>
#include <Domain.h>
#include <DistVector.h>
#include <MacroCell.h>
#include <CurvatureDetection.h>

#include <cstdio>

//------------------------------------------------------------------------------

template<int dim>
DistEdgeGrad<dim>::DistEdgeGrad(IoData& iod, Domain* domain)
{

  tag = new DistVec<bool>(domain->getNodeDistInfo());
  *tag = false;

  failSafeNewton = iod.ts.implicit.newton.failsafe;

  numLocSub  = domain->getNumLocSub();
  subDomain  = domain->getSubDomain();
  it0        = iod.restart.iteration;
  lastIt     = it0;
  lastConfig = -1;

  double spheres[SchemeFixData::num * 2][4];
  double boxes[SchemeFixData::num * 2][2][3];
  double cones[SchemeFixData::num * 2][2][4];
  int j, nspheres = 0, nboxes = 0, ncones = 0;

  for (j=0; j<iod.schemes.fixes.num; ++j) {
    if (iod.schemes.fixes.spheres[j]->r > 0.0) {
      spheres[nspheres][0] = iod.schemes.fixes.spheres[j]->x0;
      spheres[nspheres][1] = iod.schemes.fixes.spheres[j]->y0;
      spheres[nspheres][2] = iod.schemes.fixes.spheres[j]->z0;
      spheres[nspheres][3] = iod.schemes.fixes.spheres[j]->r;
      ++nspheres;
      if (iod.schemes.fixes.symmetry == SchemeFixData::X) {
        spheres[nspheres][0] = - iod.schemes.fixes.spheres[j]->x0;
        spheres[nspheres][1] = iod.schemes.fixes.spheres[j]->y0;
        spheres[nspheres][2] = iod.schemes.fixes.spheres[j]->z0;
        spheres[nspheres][3] = iod.schemes.fixes.spheres[j]->r;
        ++nspheres;
      }
      else if (iod.schemes.fixes.symmetry == SchemeFixData::Y) {
        spheres[nspheres][0] = iod.schemes.fixes.spheres[j]->x0;
        spheres[nspheres][1] = - iod.schemes.fixes.spheres[j]->y0;
        spheres[nspheres][2] = iod.schemes.fixes.spheres[j]->z0;
        spheres[nspheres][3] = iod.schemes.fixes.spheres[j]->r;
        ++nspheres;
      }
      else if (iod.schemes.fixes.symmetry == SchemeFixData::Z) {
        spheres[nspheres][0] = iod.schemes.fixes.spheres[j]->x0;
        spheres[nspheres][1] = iod.schemes.fixes.spheres[j]->y0;
        spheres[nspheres][2] = - iod.schemes.fixes.spheres[j]->z0;
        spheres[nspheres][3] = iod.schemes.fixes.spheres[j]->r;
        ++nspheres;
      }
    }
    if (iod.schemes.fixes.boxes[j]->x0 < iod.schemes.fixes.boxes[j]->x1) {
      boxes[nboxes][0][0] = iod.schemes.fixes.boxes[j]->x0;
      boxes[nboxes][0][1] = iod.schemes.fixes.boxes[j]->y0;
      boxes[nboxes][0][2] = iod.schemes.fixes.boxes[j]->z0;
      boxes[nboxes][1][0] = iod.schemes.fixes.boxes[j]->x1;
      boxes[nboxes][1][1] = iod.schemes.fixes.boxes[j]->y1;
      boxes[nboxes][1][2] = iod.schemes.fixes.boxes[j]->z1;
      ++nboxes;
      if (iod.schemes.fixes.symmetry == SchemeFixData::X) {
        boxes[nboxes][0][0] = -iod.schemes.fixes.boxes[j]->x1;
        boxes[nboxes][0][1] = iod.schemes.fixes.boxes[j]->y0;
        boxes[nboxes][0][2] = iod.schemes.fixes.boxes[j]->z0;
        boxes[nboxes][1][0] = -iod.schemes.fixes.boxes[j]->x0;
        boxes[nboxes][1][1] = iod.schemes.fixes.boxes[j]->y1;
        boxes[nboxes][1][2] = iod.schemes.fixes.boxes[j]->z1;
        ++nboxes;
      }
      if (iod.schemes.fixes.symmetry == SchemeFixData::Y) {
        boxes[nboxes][0][0] = iod.schemes.fixes.boxes[j]->x0;
        boxes[nboxes][0][1] = -iod.schemes.fixes.boxes[j]->y1;
        boxes[nboxes][0][2] = iod.schemes.fixes.boxes[j]->z0;
        boxes[nboxes][1][0] = iod.schemes.fixes.boxes[j]->x1;
        boxes[nboxes][1][1] = -iod.schemes.fixes.boxes[j]->y0;
        boxes[nboxes][1][2] = iod.schemes.fixes.boxes[j]->z1;
        ++nboxes;
      }
     if (iod.schemes.fixes.symmetry == SchemeFixData::Z) {
        boxes[nboxes][0][0] = iod.schemes.fixes.boxes[j]->x0;
        boxes[nboxes][0][1] = iod.schemes.fixes.boxes[j]->y0;
        boxes[nboxes][0][2] = -iod.schemes.fixes.boxes[j]->z1;
        boxes[nboxes][1][0] = iod.schemes.fixes.boxes[j]->x1;
        boxes[nboxes][1][1] = iod.schemes.fixes.boxes[j]->y1;
        boxes[nboxes][1][2] = -iod.schemes.fixes.boxes[j]->z0;
        ++nboxes;
      }
    }
    if (iod.schemes.fixes.cones[j]->r0 >= 0.0 && iod.schemes.fixes.cones[j]->r1 >= 0.0) {
      cones[ncones][0][0] = iod.schemes.fixes.cones[j]->x0;
      cones[ncones][0][1] = iod.schemes.fixes.cones[j]->y0;
      cones[ncones][0][2] = iod.schemes.fixes.cones[j]->z0;
      cones[ncones][0][3] = iod.schemes.fixes.cones[j]->r0;
      cones[ncones][1][0] = iod.schemes.fixes.cones[j]->x1;
      cones[ncones][1][1] = iod.schemes.fixes.cones[j]->y1;
      cones[ncones][1][2] = iod.schemes.fixes.cones[j]->z1;
      cones[ncones][1][3] = iod.schemes.fixes.cones[j]->r1;
      ++ncones;
      if (iod.schemes.fixes.symmetry == SchemeFixData::X) {
        cones[ncones][0][0] = -iod.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = iod.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = iod.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = iod.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = -iod.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = iod.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = iod.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = iod.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
      if (iod.schemes.fixes.symmetry == SchemeFixData::Y) {
        cones[ncones][0][0] = iod.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = -iod.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = iod.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = iod.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = iod.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = -iod.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = iod.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = iod.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
      if (iod.schemes.fixes.symmetry == SchemeFixData::Z) {
        cones[ncones][0][0] = iod.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = iod.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = -iod.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = iod.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = iod.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = iod.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = -iod.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = iod.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
    }
  }

 if (nspheres > 0 || nboxes > 0 || ncones > 0) {
  Communicator* com = domain->getCommunicator();
    for (j=0; j<nspheres; ++j)
      com->printf(1, "*** Warning: set the gradients to zero in [(%g, %g, %g), %g]\n",
                  spheres[j][0], spheres[j][1], spheres[j][2], spheres[j][3]);
    for (j=0; j<nboxes; ++j)
      com->printf(1, "*** Warning: set the gradients to zero in [(%g, %g, %g), (%g, %g, %g)]\n",
                  boxes[j][0][0], boxes[j][0][1], boxes[j][0][2],
                  boxes[j][1][0], boxes[j][1][1], boxes[j][1][2]);

    for (j=0; j<ncones; ++j)
      com->printf(1, "*** Warning: set the gradients to zero in cone [(%g, %g, %g), %g; (%g, %g, %g), %g]\n",
                  cones[j][0][0], cones[j][0][1], cones[j][0][2], cones[j][0][3],
                  cones[j][1][0], cones[j][1][1], cones[j][1][2], cones[j][1][3]);


    DistSVec<double,3> X0(domain->getNodeDistInfo());
    domain->getReferenceMeshPosition(X0);

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      bool* loctag = tag->subData(iSub);
      double (*x0)[3] = X0.subData(iSub);
      for (int i=0; i<tag->subSize(iSub); ++i) {
        loctag[i] = false;
        for (j=0; j<nspheres; ++j) {
          double r = sqrt( (x0[i][0] - spheres[j][0])*(x0[i][0] - spheres[j][0]) +
                           (x0[i][1] - spheres[j][1])*(x0[i][1] - spheres[j][1]) +
                           (x0[i][2] - spheres[j][2])*(x0[i][2] - spheres[j][2]) );
          if (r <= spheres[j][3])
            loctag[i] = true;
        }
        for (j=0; j<nboxes; ++j) {
          if ((x0[i][0] >= boxes[j][0][0]) && (x0[i][0] <= boxes[j][1][0]) &&
              (x0[i][1] >= boxes[j][0][1]) && (x0[i][1] <= boxes[j][1][1]) &&
              (x0[i][2] >= boxes[j][0][2]) && (x0[i][2] <= boxes[j][1][2]))
            loctag[i] = true;
        }
        for (j=0; j<ncones; ++j)  {
          Vec3D dr(cones[j][1][0]-cones[j][0][0], cones[j][1][1]-cones[j][0][1], cones[j][1][2]-cones[j][0][2]);
          double height = dr.norm();
          dr /= height;
          Vec3D xp;
          Vec3D pr0(x0[i][0]-cones[j][0][0], x0[i][1]-cones[j][0][1], x0[i][2]-cones[j][0][2]);
          double h = pr0*dr;
          if (h >= 0.0 && h <= height)  {
            xp = pr0 - (h*dr);
            double r = cones[j][0][3] + (cones[j][1][3]-cones[j][0][3]) * h / height;
            if (xp.norm() < r)
              loctag[i] = true;
          }
        }
      }
    }
  }

  if (iod.schemes.fixes.dihedralAngle > 0.0) {
    DistSVec<double,3> X0(domain->getNodeDistInfo());
    domain->getReferenceMeshPosition(X0);
    CurvatureDetection crvdet(domain);
    crvdet.compute(iod.schemes.fixes.dihedralAngle, X0, *tag);
  }

  subEdgeGrad = new EdgeGrad<dim>*[numLocSub];

  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subEdgeGrad[iSub] = new EdgeGrad<dim>(iod);
    subEdgeGrad[iSub]->settag(tag->subData(iSub));
  }

  if (tag && failSafeNewton==1) {
    backuptag = new DistVec<bool>(domain->getNodeDistInfo());
    *backuptag = *tag;
  } else {
    backuptag = 0;
  }

}

//------------------------------------------------------------------------------

template<int dim>
DistEdgeGrad<dim>::~DistEdgeGrad()
{

  if (subEdgeGrad) {

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub){

      if (subEdgeGrad[iSub])
	delete subEdgeGrad[iSub];

      delete [] subEdgeGrad;
    }
  }

  if (tag) delete(tag);
  if (backuptag) delete(backuptag);
}

//------------------------------------------------------------------------------
extern void exactinit();

template<int dim>
void DistEdgeGrad<dim>::compute(int config, DistSVec<double,3> &X)
{

  if ((config != lastConfig) || (lastIt == it0)) {
    exactinit(); //HB: to initialize all the constants used in the estimation of the rounding
                 // truncature error in floating point operations in orient3d (see utils/Predicate.C)
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->findEdgeTetrahedra(X(iSub), subEdgeGrad[iSub]->getV6NodeData());
      //subEdgeGrad[iSub]->findEdgeTetrahedra(subDomain[iSub], X(iSub));
    if (lastConfig != config)
      lastConfig = config;
    if (lastIt == it0)
      lastIt = -1;
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistEdgeGrad<dim>::computeDerivative(int config, DistSVec<double,3> &X, DistSVec<double,3> &dX)
{

  fprintf(stderr, "***** DistEdgeGrad<dim>::computeDerivative is not implemented!\n");
  exit(1);

}

//------------------------------------------------------------------------------
