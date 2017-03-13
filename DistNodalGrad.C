#include <RecFcnDesc.h>
#include <NodalGrad.h>
#include <CurvatureDetection.h>
#include <Domain.h>
#include <DistVector.h>
#include <Vector3D.h>
#include <MatVecProd.h>

#include <cstdio>
#include <cmath>
#include <alloca.h>
#include <cassert>

//------------------------------------------------------------------------------

template<int dim, class Scalar>
DistNodalGrad<dim, Scalar>::DistNodalGrad(IoData &ioData, Domain *dom) : domain(dom)
{

  int iSub;

  myIoData = &ioData;

  typeGradient = ioData.schemes.ns.gradient;

  failSafeNewton = ioData.ts.implicit.newton.failsafe;

  numLocSub = domain->getNumLocSub();

  tag = 0;

  if (ioData.schemes.ns.limiter == SchemeData::BARTH ||
      ioData.schemes.ns.limiter == SchemeData::VENKAT) 
  {
    Vmin = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
    Vmax = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
    phi = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());

// Included (MB)
    if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
        ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
		ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
		ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_)
	 {
      dVmin = new DistSVec<double,dim>(domain->getNodeDistInfo());
      dVmax = new DistSVec<double,dim>(domain->getNodeDistInfo());
      dphi = new DistSVec<double,dim>(domain->getNodeDistInfo());
    }
    else {
      dVmin = 0;
      dVmax = 0;
      dphi = 0;
    }

  }
  else {
    Vmin = 0;
    Vmax = 0;
    phi = 0;

// Included (MB)
    dVmin = 0;
    dVmax = 0;
    dphi = 0;

  }

  if (ioData.schemes.ns.limiter == SchemeData::P_SENSOR) {
    sensor = new DistSVec<Scalar,3>(domain->getNodeDistInfo());
    sigma = new DistVec<Scalar>(domain->getNodeDistInfo());
  }
  else {
    sensor = 0;
    sigma = 0;
  }

  ddx = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
  ddy = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
  ddz = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());

  dTdx = new DistVec<Scalar>(domain->getNodeDistInfo());
  dTdy = new DistVec<Scalar>(domain->getNodeDistInfo());
  dTdz = new DistVec<Scalar>(domain->getNodeDistInfo());

// Included (MB)
  if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
      ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
	  ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
	  ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
    dddx = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
    dddy = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
    dddz = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
    *dddx = 0.0;
    *dddy = 0.0;
    *dddz = 0.0;
  }
  else {
    dddx = 0;
    dddy = 0;
    dddz = 0;
  }

  *ddx = 0.0;
  *ddy = 0.0;
  *ddz = 0.0;

  *dTdx = 0.0;
  *dTdy = 0.0;
  *dTdz = 0.0;

  if (typeGradient == SchemeData::LEAST_SQUARES) {
    R = new DistSVec<double,6>(domain->getNodeDistInfo());
    wii = wij = wji = 0;

// Included (MB)
    dwii = 0;
    dwij = 0;
    dwji = 0;
    if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
        ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
		ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
		ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
      dR = new DistSVec<double,6>(domain->getNodeDistInfo());
    }
    else {
      dR = 0;
    }

  }
  else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL) {
    R = 0;
    wii = new DistSVec<double,3>(domain->getNodeDistInfo());
    wij = new DistSVec<double,3>(domain->getEdgeDistInfo());
    wji = new DistSVec<double,3>(domain->getEdgeDistInfo());

// Included (MB)
    dR = 0;
    if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
        ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
		ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
		ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
      dwii = new DistSVec<double,3>(domain->getNodeDistInfo());
      dwij = new DistSVec<double,3>(domain->getEdgeDistInfo());
      dwji = new DistSVec<double,3>(domain->getEdgeDistInfo());
    }
    else {
      dwii = 0;
      dwij = 0;
      dwji = 0;
    }

  }

  subNodalGrad = new NodalGrad<dim, Scalar>*[numLocSub];

// Included (MB)
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
        ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
		ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
		ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
      subNodalGrad[iSub] = new NodalGrad<dim, Scalar>((*ddx)(iSub), (*ddy)(iSub), (*ddz)(iSub), (*dTdx)(iSub), (*dTdy)(iSub), (*dTdz)(iSub),
                                             (*dddx)(iSub), (*dddy)(iSub), (*dddz)(iSub));
      lastConfigSA = -1;
    }
    else
      subNodalGrad[iSub] = new NodalGrad<dim, Scalar>((*ddx)(iSub), (*ddy)(iSub), (*ddz)(iSub), (*dTdx)(iSub), (*dTdy)(iSub), (*dTdz)(iSub));

  lastConfig = -1;

  int j, nspheres = 0, nboxes = 0, ncones = 0;
  for (j=0; j<ioData.schemes.fixes.num; ++j) {

    if (ioData.schemes.fixes.spheres[j]->r > 0.0) {
      if (ioData.schemes.fixes.spheres[j]->failsafe == 
          SFixData::OFF)
        continue;
      spheres[nspheres][0] = ioData.schemes.fixes.spheres[j]->x0;
      spheres[nspheres][1] = ioData.schemes.fixes.spheres[j]->y0;
      spheres[nspheres][2] = ioData.schemes.fixes.spheres[j]->z0;
      spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
      ++nspheres;
      if (ioData.schemes.fixes.symmetry == SchemeFixData::X) {
	spheres[nspheres][0] = - ioData.schemes.fixes.spheres[j]->x0;
	spheres[nspheres][1] = ioData.schemes.fixes.spheres[j]->y0;
	spheres[nspheres][2] = ioData.schemes.fixes.spheres[j]->z0;
	spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
	++nspheres;
      }
      else if (ioData.schemes.fixes.symmetry == SchemeFixData::Y) {
	spheres[nspheres][0] = ioData.schemes.fixes.spheres[j]->x0;
	spheres[nspheres][1] = - ioData.schemes.fixes.spheres[j]->y0;
	spheres[nspheres][2] = ioData.schemes.fixes.spheres[j]->z0;
	spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
	++nspheres;
      }
      else if (ioData.schemes.fixes.symmetry == SchemeFixData::Z) {
	spheres[nspheres][0] = ioData.schemes.fixes.spheres[j]->x0;
	spheres[nspheres][1] = ioData.schemes.fixes.spheres[j]->y0;
	spheres[nspheres][2] = - ioData.schemes.fixes.spheres[j]->z0;
	spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
	++nspheres;
      }
    }
    if (ioData.schemes.fixes.boxes[j]->x0 < ioData.schemes.fixes.boxes[j]->x1) {
      if (ioData.schemes.fixes.boxes[j]->failsafe == 
          BFixData::OFF)
        continue;
      boxes[nboxes][0][0] = ioData.schemes.fixes.boxes[j]->x0;
      boxes[nboxes][0][1] = ioData.schemes.fixes.boxes[j]->y0;
      boxes[nboxes][0][2] = ioData.schemes.fixes.boxes[j]->z0;
      boxes[nboxes][1][0] = ioData.schemes.fixes.boxes[j]->x1;
      boxes[nboxes][1][1] = ioData.schemes.fixes.boxes[j]->y1;
      boxes[nboxes][1][2] = ioData.schemes.fixes.boxes[j]->z1;
      ++nboxes;
      if (ioData.schemes.fixes.symmetry == SchemeFixData::X) {
	boxes[nboxes][0][0] = -ioData.schemes.fixes.boxes[j]->x1;
	boxes[nboxes][0][1] = ioData.schemes.fixes.boxes[j]->y0;
	boxes[nboxes][0][2] = ioData.schemes.fixes.boxes[j]->z0;
	boxes[nboxes][1][0] = -ioData.schemes.fixes.boxes[j]->x0;
	boxes[nboxes][1][1] = ioData.schemes.fixes.boxes[j]->y1;
	boxes[nboxes][1][2] = ioData.schemes.fixes.boxes[j]->z1;
	++nboxes;
      }
      if (ioData.schemes.fixes.symmetry == SchemeFixData::Y) {
	boxes[nboxes][0][0] = ioData.schemes.fixes.boxes[j]->x0;
	boxes[nboxes][0][1] = -ioData.schemes.fixes.boxes[j]->y1;
	boxes[nboxes][0][2] = ioData.schemes.fixes.boxes[j]->z0;
	boxes[nboxes][1][0] = ioData.schemes.fixes.boxes[j]->x1;
	boxes[nboxes][1][1] = -ioData.schemes.fixes.boxes[j]->y0;
	boxes[nboxes][1][2] = ioData.schemes.fixes.boxes[j]->z1;
	++nboxes;
      }
     if (ioData.schemes.fixes.symmetry == SchemeFixData::Z) {
	boxes[nboxes][0][0] = ioData.schemes.fixes.boxes[j]->x0;
	boxes[nboxes][0][1] = ioData.schemes.fixes.boxes[j]->y0;
	boxes[nboxes][0][2] = -ioData.schemes.fixes.boxes[j]->z1;
	boxes[nboxes][1][0] = ioData.schemes.fixes.boxes[j]->x1;
	boxes[nboxes][1][1] = ioData.schemes.fixes.boxes[j]->y1;
	boxes[nboxes][1][2] = -ioData.schemes.fixes.boxes[j]->z0;
	++nboxes;
      }
    }
    if (ioData.schemes.fixes.cones[j]->r0 >= 0.0 && ioData.schemes.fixes.cones[j]->r1 >= 0.0) {
      if (ioData.schemes.fixes.cones[j]->failsafe == 
          CFixData::OFF)
        continue;
      cones[ncones][0][0] = ioData.schemes.fixes.cones[j]->x0;
      cones[ncones][0][1] = ioData.schemes.fixes.cones[j]->y0;
      cones[ncones][0][2] = ioData.schemes.fixes.cones[j]->z0;
      cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
      cones[ncones][1][0] = ioData.schemes.fixes.cones[j]->x1;
      cones[ncones][1][1] = ioData.schemes.fixes.cones[j]->y1;
      cones[ncones][1][2] = ioData.schemes.fixes.cones[j]->z1;
      cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
      ++ncones;
      if (ioData.schemes.fixes.symmetry == SchemeFixData::X) {
        cones[ncones][0][0] = -ioData.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = ioData.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = ioData.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = -ioData.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = ioData.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = ioData.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
      if (ioData.schemes.fixes.symmetry == SchemeFixData::Y) {
        cones[ncones][0][0] = ioData.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = -ioData.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = ioData.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = ioData.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = -ioData.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = ioData.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
      if (ioData.schemes.fixes.symmetry == SchemeFixData::Z) {
        cones[ncones][0][0] = ioData.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = ioData.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = -ioData.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = ioData.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = ioData.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = -ioData.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
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


    tag = new DistVec<bool>(domain->getNodeDistInfo());
    DistSVec<double,3> X0(domain->getNodeDistInfo());
    domain->getReferenceMeshPosition(X0);

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
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

  if (ioData.schemes.fixes.dihedralAngle > 0.0) {
    if (!tag) {
      tag = new DistVec<bool>(domain->getNodeDistInfo());
      *tag = false;
    }
    DistSVec<double,3> X0(domain->getNodeDistInfo());
    domain->getReferenceMeshPosition(X0);
    CurvatureDetection crvdet(domain);
    crvdet.compute(ioData.schemes.fixes.dihedralAngle, X0, *tag);
  }

  if (tag && failSafeNewton==1) {
    backuptag = new DistVec<bool>(domain->getNodeDistInfo());
    *backuptag = *tag;
  }

/*
  if (tag) {
    DistSVec<double,1> nt(domain->getNodeDistInfo());
#pragma omp parallel for
    for (iSub = 0; iSub<numLocSub; ++iSub) {
      double (*_nt)[1] = nt.subData(iSub);
      bool* _tag = tag->subData(iSub);
      for (int i=0; i<nt.subSize(iSub); ++i) {
	if (_tag[i])
	  _nt[i][0] = 1.0;
	else
	  _nt[i][0] = 0.0;
      }
    }
    domain->writeVectorToFile("nodetag", 0, 0.0, nt);
  }
*/

  iteration = 0;
}

//------------------------------------------------------------------------------
// constructor for levelset variable (ie phi or rhophi)
// no need for fixes since phi can be both positive or negative
template<int dim, class Scalar>
DistNodalGrad<dim, Scalar>::DistNodalGrad(IoData &ioData, Domain *dom, int whichone) : domain(dom)
{

  int iSub;

  typeGradient = ioData.schemes.ls.gradient;

  numLocSub = domain->getNumLocSub();

  tag = 0;

  if (ioData.schemes.ls.limiter == SchemeData::BARTH ||
      ioData.schemes.ls.limiter == SchemeData::VENKAT) {
    Vmin = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
    Vmax = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
    phi = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
  }
  else {
    Vmin = 0;
    Vmax = 0;
    phi = 0;
  }

  if (ioData.schemes.ls.limiter == SchemeData::P_SENSOR) {
    sensor = new DistSVec<Scalar,3>(domain->getNodeDistInfo());
    sigma = new DistVec<Scalar>(domain->getNodeDistInfo());
  }
  else {
    sensor = 0;
    sigma = 0;
  }

  ddx = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
  ddy = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());
  ddz = new DistSVec<Scalar,dim>(domain->getNodeDistInfo());

  dTdx = new DistVec<Scalar>(domain->getNodeDistInfo());
  dTdy = new DistVec<Scalar>(domain->getNodeDistInfo());
  dTdz = new DistVec<Scalar>(domain->getNodeDistInfo());

  *ddx = 0.0;
  *ddy = 0.0;
  *ddz = 0.0;

  *dTdx = 0.0;
  *dTdy = 0.0;
  *dTdz = 0.0;

  if(whichone==1 || whichone==2){
    //select galerkin gradient
    R = 0;
    wii = new DistSVec<double,3>(domain->getNodeDistInfo());
    wij = new DistSVec<double,3>(domain->getEdgeDistInfo());
    wji = new DistSVec<double,3>(domain->getEdgeDistInfo());
    typeGradient = SchemeData::GALERKIN;
  }else if(whichone==0){
    //select least square gradient
    R = new DistSVec<double,6>(domain->getNodeDistInfo());
    wii = wij = wji = 0;
    typeGradient = SchemeData::LEAST_SQUARES;
  } else {
    if (typeGradient == SchemeData::LEAST_SQUARES) {
      R = new DistSVec<double,6>(domain->getNodeDistInfo());
      wii = wij = wji = 0;
    }
    else if (typeGradient == SchemeData::GALERKIN) {
      R = 0;
      wii = new DistSVec<double,3>(domain->getNodeDistInfo());
      wij = new DistSVec<double,3>(domain->getEdgeDistInfo());
      wji = new DistSVec<double,3>(domain->getEdgeDistInfo());
    }
  }

  subNodalGrad = new NodalGrad<dim, Scalar>*[numLocSub];
  lastConfig = -1;

  lastConfigSA = -1;
  dVmin = 0;
  dVmax = 0;
  dphi  = 0;
  dR    = 0;
  dwii  = 0;
  dwij  = 0;
  dwji  = 0;
  dddx  = 0;
  dddy  = 0;
  dddz  = 0;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subNodalGrad[iSub] = new NodalGrad<dim, Scalar>((*ddx)(iSub), (*ddy)(iSub), (*ddz)(iSub), (*dTdx)(iSub), (*dTdy)(iSub), (*dTdz)(iSub));
}
//------------------------------------------------------------------------------

template<int dim, class Scalar>
DistNodalGrad<dim, Scalar>::~DistNodalGrad()
{

  if (tag) delete tag;

  if (Vmin) delete Vmin;
  if (Vmax) delete Vmax;
  if (phi) delete phi;

  if (sensor) delete sensor;
  if (sigma) delete sigma;

  if (R) delete R;

  if (wii) delete wii;
  if (wij) delete wij;
  if (wji) delete wji;

  if (ddx) delete ddx;
  if (ddy) delete ddy;
  if (ddz) delete ddz;

  if (dTdx) delete dTdx;
  if (dTdy) delete dTdy;
  if (dTdz) delete dTdz;

  if (subNodalGrad) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      if (subNodalGrad[iSub]) delete subNodalGrad[iSub];

    delete [] subNodalGrad;
  }

// Included (MB)
  if (dR) delete dR;
  if (dwii) delete dwii;
  if (dwij) delete dwij;
  if (dwji) delete dwji;
  if (dddx) delete dddx;
  if (dddy) delete dddy;
  if (dddz) delete dddz;

}

template<int dim, class Scalar>
void DistNodalGrad<dim,Scalar>::updateFixes() {

  int j,iSub;

  for (j=0; j<myIoData->schemes.fixes.num; ++j) {

    if (myIoData->schemes.fixes.spheres[j]->r > 0.0) {
      if (myIoData->schemes.fixes.spheres[j]->failsafe ==
          SFixData::ALWAYSON)
        continue;

      if (iteration == myIoData->schemes.fixes.spheres[j]->failsafeN) {
        DistSVec<double,3> X0(domain->getNodeDistInfo());
        domain->getReferenceMeshPosition(X0);

        domain->getCommunicator()->fprintf(stdout, "Sphere fix at [%lf %lf %lf], r = %lf expired\n",
                        spheres[j][0],spheres[j][1],spheres[j][2],spheres[j][3]);

#pragma omp parallel for
        for (iSub = 0; iSub < numLocSub; ++iSub) {
          bool* loctag = tag->subData(iSub);
          double (*x0)[3] = X0.subData(iSub);
          for (int i=0; i<tag->subSize(iSub); ++i) {
	    double r = sqrt((x0[i][0] - spheres[j][0])*(x0[i][0] - spheres[j][0]) +
	                    (x0[i][1] - spheres[j][1])*(x0[i][1] - spheres[j][1]) +
			    (x0[i][2] - spheres[j][2])*(x0[i][2] - spheres[j][2]) );
	    if (r <= spheres[j][3])
	        loctag[i] = false;
	  }
        }
      }
    }
    
    if (myIoData->schemes.fixes.boxes[j]->x0 < myIoData->schemes.fixes.boxes[j]->x1) {
      if (myIoData->schemes.fixes.boxes[j]->failsafe == 
          BFixData::ALWAYSON)
        continue;

      if (iteration == myIoData->schemes.fixes.boxes[j]->failsafeN) {
        DistSVec<double,3> X0(domain->getNodeDistInfo());
        domain->getReferenceMeshPosition(X0);
        domain->getCommunicator()->fprintf(stdout, "Box fix at [[%lf %lf %lf], [%lf %lf %lf]] expired\n",
                        boxes[j][0],boxes[j][1],boxes[j][2],
                        boxes[j][3],boxes[j][4],boxes[j][5]);

#pragma omp parallel for
        for (iSub = 0; iSub < numLocSub; ++iSub) {
          bool* loctag = tag->subData(iSub);
          double (*x0)[3] = X0.subData(iSub);
          for (int i=0; i<tag->subSize(iSub); ++i) {
            if ((x0[i][0] >= boxes[j][0][0]) && (x0[i][0] <= boxes[j][1][0]) &&
                (x0[i][1] >= boxes[j][0][1]) && (x0[i][1] <= boxes[j][1][1]) &&
                (x0[i][2] >= boxes[j][0][2]) && (x0[i][2] <= boxes[j][1][2]))
              loctag[i] = false;
          }
        }
      }
    }

    if (myIoData->schemes.fixes.cones[j]->r0 >= 0.0 && myIoData->schemes.fixes.cones[j]->r1 >= 0.0) {
      if (myIoData->schemes.fixes.cones[j]->failsafe == 
          CFixData::ALWAYSON)
        continue;
 
      if (iteration == myIoData->schemes.fixes.cones[j]->failsafeN) {
        DistSVec<double,3> X0(domain->getNodeDistInfo());
        domain->getReferenceMeshPosition(X0);

        domain->getCommunicator()->fprintf(stdout,"Cone fix expired\n");

#pragma omp parallel for
        for (iSub = 0; iSub < numLocSub; ++iSub) {
          bool* loctag = tag->subData(iSub);
          double (*x0)[3] = X0.subData(iSub);
          for (int i=0; i<tag->subSize(iSub); ++i) {
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
                loctag[i] = false;
            }
          }
        }
      }
    }
  }    

  ++iteration;
}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::computeWeights(DistSVec<double,3> &X)
{

  if (typeGradient == SchemeData::LEAST_SQUARES)
    domain->computeWeightsLeastSquares(X, *R);
  else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL)
    domain->computeWeightsGalerkin(X, *wii, *wij, *wji);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::computeDerivativeOfWeights(DistSVec<double,3> &X, DistSVec<double,3> &dX)
{

  if (typeGradient == SchemeData::LEAST_SQUARES) {
//Remark: Error mesage for pointers
    if (dR == 0) {
      fprintf(stderr, "*** Error: Variable dR does not exist!\n");
      exit(1);
    }
    domain->computeDerivativeOfWeightsLeastSquares(X, dX, *dR);
  }
  else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL) {
//Remark: Error mesage for pointers
    if (dwii == 0) {
      fprintf(stderr, "*** Error: Variable dwii does not exist!\n");
      exit(1);
    }
    if (dwij == 0) {
      fprintf(stderr, "*** Error: Variable dwij does not exist!\n");
      exit(1);
    }
    if (dwji == 0) {
      fprintf(stderr, "*** Error: Variable dwji does not exist!\n");
      exit(1);
    }

    domain->computeDerivativeOfWeightsGalerkin(X, dX, *dwii, *dwij, *dwji);

  }

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::computeDerivativeOfWeights(dRdXoperators<dim> &dRdXop, DistSVec<double,3> &dX, DistSVec<double,6> &dR2)
{

  if (typeGradient == SchemeData::LEAST_SQUARES) {
//Remark: Error mesage for pointers
    if (dR == 0) {
      fprintf(stderr, "*** Error: Variable dR does not exist!\n");
      exit(1);
    }
    domain->computeDerivativeOfWeightsLeastSquares(dRdXop.dRdX,dRdXop.dRdR, dX, dR2);
   
  }
  else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL) {
//Remark: Error mesage for pointers
    if (dwii == 0) {
      fprintf(stderr, "*** Error: Variable dwii does not exist!\n");
      exit(1);
    }
    if (dwij == 0) {
      fprintf(stderr, "*** Error: Variable dwij does not exist!\n");
      exit(1);
    }
    if (dwji == 0) {
      fprintf(stderr, "*** Error: Variable dwji does not exist!\n");
      exit(1);
    }

//    domain->computeDerivativeOfWeightsGalerkin(X, dX, *dwii, *dwij, *dwji);
    fprintf(stderr, "*** Error: DistNodalGrad<dim, Scalar>::computeDerivativeOfWeights is not implemented yet.\n");
    fprintf(stderr, "           Turn off SparseFlag under SensitivityAnalysis.\n");  exit(-1);

  }

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::computeTransposeDerivativeOfWeights(dRdXoperators<dim> &dRdXop, DistSVec<double,6> &dR2, DistSVec<double,3> &dX)
{

  if (typeGradient == SchemeData::LEAST_SQUARES) {
//Remark: Error mesage for pointers
    if (dR == 0) {
      fprintf(stderr, "*** Error: Variable dR does not exist!\n");
      exit(1);
    }
    domain->computeTransposeDerivativeOfWeightsLeastSquares(dRdXop.dRdX, dRdXop.dRdR, dR2, dX);
  }
  else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL) {
//Remark: Error mesage for pointers
    if (dwii == 0) {
      fprintf(stderr, "*** Error: Variable dwii does not exist!\n");
      exit(1);
    }
    if (dwij == 0) {
      fprintf(stderr, "*** Error: Variable dwij does not exist!\n");
      exit(1);
    }
    if (dwji == 0) {
      fprintf(stderr, "*** Error: Variable dwji does not exist!\n");
      exit(1);
    }

//    domain->computeDerivativeOfWeightsGalerkin(X, dX, *dwii, *dwij, *dwji);
    fprintf(stderr, "*** Error: DistNodalGrad<dim, Scalar>::computeTransposeDerivativeOfWeights is not implemented yet.\n");
    fprintf(stderr, "           Turn off SparseFlag under SensitivityAnalysis.\n");  exit(-1);

  }

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::computeDerivativeOfWeightsOperators(DistSVec<double,3> &X, dRdXoperators<dim> &dRdXop)
{

  if (typeGradient == SchemeData::LEAST_SQUARES) 
    domain->computeDerivativeOperatorsOfWeightsLeastSquares(X, dRdXop.dRdX, dRdXop.dRdR);
  else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL) {
//Remark: Error mesage for pointers
    if (dwii == 0) {
      fprintf(stderr, "*** Error: Variable dwii does not exist!\n");
      exit(1);
    }
    if (dwij == 0) {
      fprintf(stderr, "*** Error: Variable dwij does not exist!\n");
      exit(1);
    }
    if (dwji == 0) {
      fprintf(stderr, "*** Error: Variable dwji does not exist!\n");
      exit(1);
    }

    fprintf(stderr, "*** Error: DistNodalGrad<dim, Scalar>::computeDerivativeOfWeightsOperators is not implemented yet.\n");
    fprintf(stderr, "           Turn off SparseFlag under SensitivityAnalysis.\n");  exit(-1);

  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
template<class Scalar2>
void DistNodalGrad<dim, Scalar>::compute(int config, DistSVec<double,3> &X,
													  DistVec<double> &ctrlVol, 
													  DistSVec<Scalar2, dim> &V)
{

	if(config != lastConfig) 
	{
    computeWeights(X);
    lastConfig = config;
  }

  if (typeGradient == SchemeData::LEAST_SQUARES)
    domain->computeGradientsLeastSquares(X, *R, V, *ddx, *ddy, *ddz);	// KTC fix
  else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL)
    domain->computeGradientsGalerkin(ctrlVol, *wii, *wij, *wji, V, *ddx, *ddy, *ddz);	// KTC fix

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
template<class Scalar2>
void DistNodalGrad<dim, Scalar>::computeDerivative(int configSA, DistSVec<double,3> &X, DistSVec<double,3> &dX,
		DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol, DistSVec<Scalar2,dim> &V, DistSVec<Scalar2,dim> &dV)
{

//Remark: Error mesage for pointers
  if (dddx == 0) {
    fprintf(stderr, "*** Error: Variable dddx does not exist!\n");
    exit(1);
  }
  if (dddy == 0) {
    fprintf(stderr, "*** Error: Variable dddy does not exist!\n");
    exit(1);
  }
  if (dddz == 0) {
    fprintf(stderr, "*** Error: Variable dddz does not exist!\n");
    exit(1);
  }

  if (configSA != lastConfigSA) {
    computeDerivativeOfWeights(X, dX);
    lastConfigSA = configSA;
  }

  if (typeGradient == SchemeData::LEAST_SQUARES) {
    domain->computeDerivativeOfGradientsLeastSquares(X, dX, *R, *dR, V, dV,*dddx, *dddy, *dddz);
  } else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL) {
    domain->computeDerivativeOfGradientsGalerkin(ctrlVol, dCtrlVol, *wii, *wij, *wji, *dwii, *dwij, *dwji, V, dV, *dddx, *dddy, *dddz);
  }
}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::computeDerivative(dRdXoperators<dim> *dRdXop, DistSVec<double,3> &dX, 
                                                   DistVec<double> &dCtrlVol, DistSVec<double,dim> &dV_r, 
                                                   DistSVec<double,6>& dR_r, DistSVec<double,dim>& dddx_r, 
                                                   DistSVec<double,dim>& dddy_r, DistSVec<double,dim>& dddz_r)
{

  computeDerivativeOfWeights(*dRdXop, dX, dR_r);

  *dR = dR_r;

  if (typeGradient == SchemeData::LEAST_SQUARES) {
    domain->computeDerivativeOfGradientsLeastSquares(*dRdXop, dX, dR_r, dV_r, dddx_r, dddy_r, dddz_r);
    *dddx = dddx_r;    *dddy = dddy_r;    *dddz = dddz_r;
  } else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL) {
	fprintf(stderr, " *** ERROR : DistNodalGrad<dim, Scalar>::computeDerivative is not implemented\n");
	exit(-1);
  }
}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::computeTransposeDerivative(dRdXoperators<dim> *dRdXop,  
                                                            DistSVec<double,dim> &dddx2, 
                                                            DistSVec<double,dim> &dddy2, 
                                                            DistSVec<double,dim> &dddz2, 
                                                            DistSVec<double,6> &dR2,
				                                                    DistVec<double> &dCtrlVol2,
                                                            DistSVec<double,dim> &dV2,
                                                            DistSVec<double,3> &dX2)
{

  if (typeGradient == SchemeData::LEAST_SQUARES) {
    Communicator* com = domain->getCommunicator();
    domain->computeTransposeDerivativeOfGradientsLeastSquares(*dRdXop, dddx2, dddy2, dddz2, dX2, dR2, dV2);
  } else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL) {
    fprintf(stderr, " *** ERROR : DistNodalGrad<dim, Scalar>::computeTransposeDerivative is not implemented\n");
    exit(-1);
  }

  computeTransposeDerivativeOfWeights(*dRdXop, dR2, dX2);
}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim, class Scalar>
template<class Scalar2>
void DistNodalGrad<dim, Scalar>::computeDerivativeOperators(DistSVec<double,3> &X, DistVec<double> &ctrlVol, 
                                                            DistSVec<Scalar2,dim> &V, dRdXoperators<dim> &dRdXop)
{

//Remark: Error mesage for pointers
  if (dddx == 0) {
    fprintf(stderr, "*** Error: Variable dddx does not exist!\n");
    exit(1);
  }
  if (dddy == 0) {
    fprintf(stderr, "*** Error: Variable dddy does not exist!\n");
    exit(1);
  }
  if (dddz == 0) {
    fprintf(stderr, "*** Error: Variable dddz does not exist!\n");
    exit(1);
  }

  computeDerivativeOfWeightsOperators(X,dRdXop);

  if (typeGradient == SchemeData::LEAST_SQUARES)
    domain->computeDerivativeOperatorsOfGradientsLeastSquares(X, *R, V, dRdXop);
  else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL) {
    fprintf(stderr, "*** ERROR: DistNodalGrad<dim, Scalar>::computeDerivativeOperators is not implemented yet\n");
    exit(-1);
  }
}

//------------------------------------------------------------------------------
// least square gradient involving only nodes of same fluid (multiphase flow)
// $d2d
template<int dim, class Scalar>
template<class Scalar2>
void DistNodalGrad<dim, Scalar>::compute(int config, DistSVec<double,3> &X,
					 DistVec<double> &ctrlVol, DistVec<int> &fluidId,
					 DistSVec<Scalar2, dim> &V, bool linFSI,
					 DistLevelSetStructure *distLSS,
					 bool includeSweptNodes)
{

	if(typeGradient == SchemeData::LEAST_SQUARES)
	{    
		domain->computeWeightsLeastSquares(X, fluidId, *R, distLSS, includeSweptNodes);
    
    domain->computeGradientsLeastSquares(X, fluidId, *R, V, *ddx, *ddy, *ddz, linFSI, 
					 distLSS, includeSweptNodes);
	}
	else if(typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL)
	{    
		domain->computeWeightsGalerkin(X, fluidId, *wii, *wij, *wji, distLSS, includeSweptNodes);

    domain->computeGradientsGalerkin(ctrlVol, *wii, *wij, *wji, V, *ddx, *ddy, *ddz);       

  }

}

//------------------------------------------------------------------------------
// least square gradient of temperature involving only nodes of same fluid (multiphase flow)
// $dd
template<int dim, class Scalar>
template<class Scalar2>
void DistNodalGrad<dim, Scalar>::computeTemperatureGradient(int config, DistSVec<double,3> &X,
                                 DistVec<double> &ctrlVol, DistVec<int> &fluidId,
                                 DistVec<Scalar2> &T,
                                 DistLevelSetStructure *distLSS)
{
  if (typeGradient == SchemeData::LEAST_SQUARES){
    
    //domain->computeWeightsLeastSquares(X, fluidId, *R, distLSS);
    domain->computeGradientsLeastSquares(X, fluidId, *R, T, *dTdx, *dTdy, *dTdz, distLSS);

  }else if(typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL){
    
    domain->computeWeightsGalerkin(X, fluidId, *wii, *wij, *wji, distLSS);
    domain->computeGradientsGalerkin(ctrlVol, *wii, *wij, *wji, T, *dTdx, *dTdy, *dTdz);

  }
  
}

//------------------------------------------------------------------------------
// least square gradient involving only nodes of fluid (FSI)
// Wstar is involved in gradient computation
// $dd
template<int dim, class Scalar>
template<class Scalar2>
void DistNodalGrad<dim, Scalar>::compute(int config, DistSVec<double,3> &X,
					 DistVec<double> &ctrlVol, DistVec<int> &fluidId,
					 DistSVec<Scalar2, dim> &V, 
					 DistSVec<Scalar2, dim> &Wstarij, DistSVec<Scalar2, dim> &Wstarji, 
					 DistVec<int> &countWstarij, DistVec<int> &countWstarji,
					 bool linFSI, DistLevelSetStructure *distLSS)

{
  if (typeGradient == SchemeData::LEAST_SQUARES){
    
    domain->computeWeightsLeastSquares(X, fluidId, *R, countWstarij, countWstarji, distLSS);
    domain->computeGradientsLeastSquares(X, fluidId, *R, V, Wstarij, Wstarji, 
					 countWstarij, countWstarji, *ddx, *ddy, *ddz, linFSI, distLSS);

  }else if(typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL){
    
    //domain->computeWeightsGalerkin(X, fluidId, *wii, *wij, *wji, distLSS);
    //domain->computeGradientsGalerkin(ctrlVol, *wii, *wij, *wji, V, *ddx, *ddy, *ddz);

  }

}

//------------------------------------------------------------------------------
template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::compute(int config, DistSVec<double,3> &X,
				         DistSVec<double,dim> &Psi)
{

  //$dd ???
  assert(typeGradient == SchemeData::LEAST_SQUARES);

  if (config != lastConfig)
    computeWeights(X);

  domain->computeGradientsLeastSquares(X, *R, Psi, *ddx, *ddy, *ddz);

}
//------------------------------------------------------------------------------

template<int dim, class Scalar>
template<class Scalar2>
void DistNodalGrad<dim, Scalar>::computeT(int config, DistSVec<double,3> &X,
                DistVec<double> &ctrlVol, DistSVec<Scalar2,dim> &Vx,
                DistSVec<Scalar2,dim> &Vy, DistSVec<Scalar2,dim> &Vz)
{

  if (config != lastConfig) {
    computeWeights(X);
    lastConfig = config;
  }

  if (typeGradient == SchemeData::LEAST_SQUARES){
    Communicator* com = domain->getCommunicator();
    com->fprintf(stderr," *** Error : The transpose of the Jacobian is not supported for Least Squares gradients types\n");
    exit(-1);
    //domain->computeGradientsLeastSquares(X, *R, Vx, *ddx, *ddy, *ddz); buggish routine
   }
  else if (typeGradient == SchemeData::GALERKIN || typeGradient == SchemeData::NON_NODAL)
    domain->computeGradientsGalerkinT(ctrlVol, *wii, *wij, *wji, Vx, Vy, Vz, *ddx, *ddy, *ddz);

  if (tag) {

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      bool *loctag = tag->subData(iSub);
      Scalar2 (*locddx)[dim] = ddx->subData(iSub);
      Scalar2 (*locddy)[dim] = ddy->subData(iSub);
      Scalar2 (*locddz)[dim] = ddz->subData(iSub);
      for (int i=0; i<tag->subSize(iSub); ++i) {
        if (loctag[i]) {
          for (int j=0; j<dim; ++j) {
            locddx[i][j] = 0.0;
            locddy[i][j] = 0.0;
            locddz[i][j] = 0.0;
          }
        }
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
template<class Scalar2>
void DistNodalGrad<dim, Scalar>::limit(RecFcn *recFcn, DistSVec<double,3> &X,
			       DistVec<double> &ctrlVol, DistSVec<Scalar2,dim> &V,
                               DistVec<int>* additionalFirstOrderNodes)
{

  RecFcnLtdMultiDim<dim>* ltdmd = dynamic_cast<RecFcnLtdMultiDim<dim>*>(recFcn);
  if (ltdmd)
    domain->computeMultiDimLimiter(ltdmd, X, ctrlVol, V, *ddx, *ddy, *ddz, *Vmin, *Vmax, *phi);

  RecFcnLtdSensor* ltdsensor = dynamic_cast<RecFcnLtdSensor*>(recFcn);
  if (ltdsensor)
    domain->computePressureSensor(ltdsensor->getThreshold(), X, V, *ddx, *ddy, *ddz, *sensor, *sigma);

  if (tag) 
  {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) 
	 {
      bool *loctag = tag->subData(iSub);
      Scalar (*locddx)[dim] = ddx->subData(iSub);
      Scalar (*locddy)[dim] = ddy->subData(iSub);
      Scalar (*locddz)[dim] = ddz->subData(iSub);

		 for (int i=0; i<tag->subSize(iSub); ++i) 
		 {
			 if (loctag[i]) 
			 {
				 for (int j=0; j<dim; ++j) 
				 {
	    locddx[i][j] = 0.0;
	    locddy[i][j] = 0.0;
	    locddz[i][j] = 0.0;
	  }
	}
      }
    }
  }

  if (additionalFirstOrderNodes) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      int *loctag = additionalFirstOrderNodes->subData(iSub);
      Scalar (*locddx)[dim] = ddx->subData(iSub);
      Scalar (*locddy)[dim] = ddy->subData(iSub);
      Scalar (*locddz)[dim] = ddz->subData(iSub);
      for (int i=0; i<additionalFirstOrderNodes->subSize(iSub); ++i) {
 	if (loctag[i]) {
  //        std::cout << "Gradient set to zero additionally at x = [" << X(iSub)[i][0] << ", " <<
  //                   X(iSub)[i][1] << ", " << X(iSub)[i][2] << "] " << std::endl;
 	  for (int j=0; j<dim; ++j) {
	    locddx[i][j] = 0.0;
	    locddy[i][j] = 0.0;
	    locddz[i][j] = 0.0;
	  }
	}
      }
    }
  }


}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
template<class Scalar2>
void DistNodalGrad<dim, Scalar>::limitDerivative(RecFcn *recFcn, DistSVec<double,3> &X, DistSVec<double,3> &dX,
			       DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol, DistSVec<Scalar2,dim> &V, DistSVec<Scalar2,dim> &dV)
{
  RecFcnLtdMultiDim<dim>* ltdmd = dynamic_cast<RecFcnLtdMultiDim<dim> *>(recFcn);
  Communicator* com = domain->getCommunicator();

  if (ltdmd) {
    //Remark: Error mesage for pointers
     if (dVmin == 0) {
       fprintf(stderr, "*** Error: Variable dVmin does not exist!\n");
       exit(1);
     }
     if (dVmax == 0) {
       fprintf(stderr, "*** Error: Variable dVmax does not exist!\n");
       exit(1);
     }
     if (dphi == 0) {
       fprintf(stderr, "*** Error: Variable dphi does not exist!\n");
      exit(1);
    }

    domain->computeDerivativeOfMultiDimLimiter(ltdmd, X, dX, ctrlVol, dCtrlVol, V, dV, *ddx, *ddy, *ddz,
                                               *dddx, *dddy, *dddz, *Vmin, *dVmin, *Vmax, *dVmax, *phi, *dphi);
  }

  RecFcnLtdSensor* ltdsensor = dynamic_cast<RecFcnLtdSensor*>(recFcn);
  if (ltdsensor) {
    fprintf(stderr, "*** Error: The derivative of the function computePressure does not exist!\n");
    exit(1);
  }

  if (tag) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      bool *loctag = tag->subData(iSub);
      Scalar (*locddx)[dim] = ddx->subData(iSub);
      Scalar (*locddy)[dim] = ddy->subData(iSub);
      Scalar (*locddz)[dim] = ddz->subData(iSub);
      Scalar (*dlocddx)[dim] = dddx->subData(iSub);
      Scalar (*dlocddy)[dim] = dddy->subData(iSub);
      Scalar (*dlocddz)[dim] = dddz->subData(iSub);
      for (int i=0; i<tag->subSize(iSub); ++i) {
        if (loctag[i]) {
          for (int j=0; j<dim; ++j) {
            locddx[i][j] = 0.0;
            locddy[i][j] = 0.0;
            locddz[i][j] = 0.0;
            dlocddx[i][j] = 0.0;
            dlocddy[i][j] = 0.0;
            dlocddz[i][j] = 0.0;
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
template<class Scalar2>
void DistNodalGrad<dim, Scalar>::limit(RecFcn *recFcn, DistSVec<double,3> &X,
													DistVec<double> &ctrlVol, 
													DistLevelSetStructure *distLSS,
													DistSVec<Scalar2,dim> &V)
{

	RecFcnLtdMultiDim<dim>* ltdmd = dynamic_cast<RecFcnLtdMultiDim<dim>*>(recFcn);
	if (ltdmd)
		domain->computeMultiDimLimiter(X, ctrlVol, V, *ddx, *ddy, *ddz, distLSS);
		
	RecFcnLtdSensor* ltdsensor = dynamic_cast<RecFcnLtdSensor*>(recFcn);
	if (ltdsensor)
		domain->computePressureSensor(ltdsensor->getThreshold(), X, V, *ddx, *ddy, *ddz, *sensor, *sigma);

	if (tag) 
	{
#pragma omp parallel for
		for (int iSub = 0; iSub < numLocSub; ++iSub) 
		{
			bool *loctag = tag->subData(iSub);

			Scalar (*locddx)[dim] = ddx->subData(iSub);
			Scalar (*locddy)[dim] = ddy->subData(iSub);
			Scalar (*locddz)[dim] = ddz->subData(iSub);

			for (int i=0; i<tag->subSize(iSub); ++i) 
			{
				if (loctag[i]) 
				{
					for (int j=0; j<dim; ++j) 
					{
						locddx[i][j] = 0.0;
						locddy[i][j] = 0.0;
						locddz[i][j] = 0.0;
					}
				}
			}
		}
	}

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::fix(DistSVec<bool,2>& fstag)
{

  if (!tag) {
    tag = new DistVec<bool>(domain->getNodeDistInfo());
    *tag = false;
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    bool* t = tag->subData(iSub);
    bool (*fst)[2] = fstag.subData(iSub);
    for (int i=0; i<tag->subSize(iSub); ++i)
      t[i] = t[i] || fst[i][1];
  }

}

//------------------------------------------------------------------------------


template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::fix(DistSVec<int,2>& fstag)
{

  if (!tag) {
    tag = new DistVec<bool>(domain->getNodeDistInfo());
    *tag = false;
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    bool* t = tag->subData(iSub);
    int (*fst)[2] = fstag.subData(iSub);
    for (int i=0; i<tag->subSize(iSub); ++i) {
      t[i] = t[i] || bool(fst[i][1]);
    }
  }

}


//------------------------------------------------------------------------------


template<int dim, class Scalar>
void DistNodalGrad<dim, Scalar>::resetTag()
{


#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    bool* t = tag->subData(iSub);
    bool* backupt = backuptag->subData(iSub);
    for (int i=0; i<tag->subSize(iSub); ++i)
      t[i] = backupt[i];
  }

}

//------------------------------------------------------------------------------


