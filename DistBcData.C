#include <DistBcData.h>

#include <IoData.h>
#include <VarFcn.h>
#include <BcData.h>
#include <SubDomain.h>
#include <Domain.h>

#include <cmath>
#include <cassert>

//------------------------------------------------------------------------------

template<int dim>
DistBcData<dim>::DistBcData(IoData &ioData, VarFcn *varFcn, Domain *domain,
                            DistInfo* __nodeDistInfo,DistInfo* __inletNodeDistInfo,
                            DistInfo* __faceDistInfo ) : 
  nodeDistInfo(__nodeDistInfo?*__nodeDistInfo:domain->getNodeDistInfo()),
  inletNodeDistInfo(__inletNodeDistInfo?*__inletNodeDistInfo:domain->getInletNodeDistInfo()),
  faceDistInfo(__faceDistInfo?*__faceDistInfo:domain->getFaceDistInfo()),
  Xdot(nodeDistInfo), Temp(nodeDistInfo),vf(varFcn),
  Ufarin(nodeDistInfo), Ufarout(nodeDistInfo), Uporouswall(nodeDistInfo),
  Uface(faceDistInfo), Unode(nodeDistInfo),
  Uinletnode(inletNodeDistInfo), rotInfo(ioData.rotations.rotationMap.dataMap)
{
  this->boundaryStateHH = 0;

// Included (MB)
  if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
      ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
	  ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
	  ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
    this->dXdot = new DistSVec<double,3>(nodeDistInfo);
    this->dTemp = new DistVec<double>(nodeDistInfo);
    this->dUface = new DistSVec<double,dim>(faceDistInfo);
    this->dUnode = new DistSVec<double,dim>(nodeDistInfo);
    this->dUinletnode = new DistSVec<double,dim>(faceDistInfo);
    this->dUfarin = new DistSVec<double,dim>(nodeDistInfo);
    this->dUfarout = new DistSVec<double,dim>(nodeDistInfo);
    this->dUporouswall = new DistSVec<double,dim>(nodeDistInfo);
  }
  else {
    this->dXdot = 0;
    this->dTemp = 0;
    this->dUface = 0;
    this->dUnode = 0;
    this->dUinletnode = 0;
    this->dUfarin = 0;
    this->dUfarout = 0;
    this->dUporouswall = 0;
  }
  if ((ioData.eqs.type == EquationsData::NAVIER_STOKES) && (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY)) {
    if ((ioData.bc.wall.integration == BcsWallData::WALL_FUNCTION) && (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS)) {
      this->dUfaceSA = new DistSVec<double,dim>(faceDistInfo);
      this->dUnodeSA = new DistSVec<double,dim>(nodeDistInfo);
    } 
    else {
      this->dUfaceSA = 0;
      this->dUnodeSA = 0;
    }
  }
  else {
    this->dUfaceSA = 0;
    this->dUnodeSA = 0;
  }

  if (ioData.schemes.bc.type == BoundarySchemeData::MODIFIED_GHIDAGLIA) {

    this->boundaryStateHH = new DistVec<double>(faceDistInfo);
  }

  //KW: "insensitive" Euler starts here.
  this->numLocSub = domain->getNumLocSub();
  this->subDomain = domain->getSubDomain();
  this->com = domain->getCommunicator();

  subBcData = new BcData<dim>*[this->numLocSub];
#pragma omp parallel for
  for (int iSub=0; iSub<this->numLocSub; ++iSub)
// Included (MB)
    if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
        ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
		ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
		ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
      if ((ioData.eqs.type == EquationsData::NAVIER_STOKES) && (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY)) {
        if ((ioData.bc.wall.integration == BcsWallData::WALL_FUNCTION) && (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS)) {
          subBcData[iSub] = new BcData<dim>(this->Uface(iSub), this->Unode(iSub), this->Uinletnode(iSub), this->Ufarin(iSub), this->Ufarout(iSub), this->Uporouswall(iSub), (*dUface)(iSub), (*dUnode)(iSub), (*dUinletnode)(iSub), (*dUfarin)(iSub), (*dUfarout)(iSub), (*dUporouswall)(iSub), (*dUfaceSA)(iSub), (*dUnodeSA)(iSub));
        }
        else {
          subBcData[iSub] = new BcData<dim>(this->Uface(iSub), this->Unode(iSub), this->Uinletnode(iSub), this->Ufarin(iSub), this->Ufarout(iSub), this->Uporouswall(iSub), (*dUface)(iSub), (*dUnode)(iSub), (*dUinletnode)(iSub), (*dUfarin)(iSub), (*dUfarout)(iSub), (*dUporouswall)(iSub));
	}
      }
      else {
        subBcData[iSub] = new BcData<dim>(this->Uface(iSub), this->Unode(iSub), this->Uinletnode(iSub), this->Ufarin(iSub), this->Ufarout(iSub), this->Uporouswall(iSub), (*dUface)(iSub), (*dUnode)(iSub), (*dUinletnode)(iSub), (*dUfarin)(iSub), (*dUfarout)(iSub), (*dUporouswall)(iSub));
      }
    }
    else {
      if ((ioData.eqs.type == EquationsData::NAVIER_STOKES) && (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY)) {
        if ((ioData.bc.wall.integration == BcsWallData::WALL_FUNCTION) && (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS)) {
          subBcData[iSub] = new BcData<dim>(this->Uface(iSub), this->Unode(iSub), this->Uinletnode(iSub), this->Ufarin(iSub), this->Ufarout(iSub), this->Uporouswall(iSub), (*dUfaceSA)(iSub), (*dUnodeSA)(iSub));
        }
        else {
          subBcData[iSub] = new BcData<dim>(this->Uface(iSub), this->Unode(iSub), this->Uinletnode(iSub), this->Ufarin(iSub), this->Ufarout(iSub), this->Uporouswall(iSub));
	}
      }
      else {
        subBcData[iSub] = new BcData<dim>(this->Uface(iSub),this->Unode(iSub), this->Uinletnode(iSub), this->Ufarin(iSub), this->Ufarout(iSub), this->Uporouswall(iSub));
      }
    }

  if (ioData.schemes.bc.type == BoundarySchemeData::MODIFIED_GHIDAGLIA) {

#pragma omp parallel for
    for (int iSub=0; iSub<this->numLocSub; ++iSub)
      subBcData[iSub]->setBoundaryStateVectorHH(&(*(this->boundaryStateHH))(iSub));
  }

  // it is assumed that at only fluid can be found at the far-field boundary of an external flow
  // and it is always the fluid with id=0 (by default)
  if(varFcn->getType()==VarFcnBase::STIFFENEDGAS || varFcn->getType()==VarFcnBase::PERFECTGAS)
    boundaryFluid = GAS;
  else if (varFcn->getType()==VarFcnBase::JWL)
    boundaryFluid = JWL;
  else
    boundaryFluid = TAIT;

  angles[0] = ioData.bc.inlet.alpha;
  angles[1] = ioData.bc.inlet.beta;

  depth     = ioData.bc.hydro.depth;
  ngravity[0] = ioData.eqs.gravity_x;
  ngravity[1] = ioData.eqs.gravity_y;
  ngravity[2] = ioData.eqs.gravity_z;
  gravity   = sqrt(ngravity[0]*ngravity[0]+ngravity[1]*ngravity[1]+ngravity[2]*ngravity[2]);


  Xdot      = 0.0;
  Temp      = ioData.bc.wall.temperature;
  this->Ufarin     = 0.0;
  this->Ufarout    = 0.0;
  this->Uporouswall= 0.0;
  this->Unode      = 0.0;
  this->Uface      = 0.0;
  this->Uinletnode = 0.0;
  tref = ioData.ref.rv.time;
  vref = ioData.ref.rv.velocity;

  /** Update the temperature for walls for which a different temperature was specified. */

  // The input __nodeDistInfo is only specified for multigrid problems, 
  // in which we specify our own DistInfo.  For these problems, this block is turned
  // off.
  if (!__nodeDistInfo) {

    DistVec<int> faceFlag(nodeDistInfo);
    faceFlag = 0;
    CommPattern<int> ndC(domain->getSubTopo(), this->com, CommPattern<int>::CopyOnSend);

#pragma omp parallel for
   for (int iSub = 0; iSub<numLocSub; ++iSub)
      subDomain[iSub]->setComLenNodes(1, ndC);

    ndC.finalize();

#pragma omp parallel for
    for (int iSub = 0; iSub<numLocSub; ++iSub)
      subDomain[iSub]->markFaceBelongsToSurface(faceFlag(iSub), ndC);

    ndC.exchange();

#pragma omp parallel for
    for (int iSub = 0; iSub<numLocSub; ++iSub)
      subDomain[iSub]->completeFaceBelongsToSurface(faceFlag(iSub), Temp(iSub),
                          ioData.surfaces.surfaceMap.dataMap, ndC);


   this->com->sync();
  }


// Included (MB)
  if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
      ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
	  ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
	  ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
    (*dXdot) = 0.0;
    (*dTemp) = 0.0;
    (*dUface) = 0.0;
    (*dUnode) = 0.0;
    (*dUinletnode) = 0.0;
    (*dUfarin) = 0.0;
    (*dUfarout) = 0.0;
    (*dUporouswall) = 0.0;
  }
  dtrefdMach = ioData.ref.rv.dtimedMach;
  dvrefdMach = ioData.ref.rv.dvelocitydMach;
  if ((ioData.eqs.type == EquationsData::NAVIER_STOKES) && (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY))
    if ((ioData.bc.wall.integration == BcsWallData::WALL_FUNCTION) && (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS)) {
      (*dUfaceSA) = 0.0;
      (*dUnodeSA) = 0.0;
    } 

}

//------------------------------------------------------------------------------

template<int dim>
DistBcData<dim>::~DistBcData()
{
 
  if (subBcData) {
#pragma omp parallel for
    for (int iSub=0; iSub<this->numLocSub; ++iSub)
      if (subBcData[iSub]) 
	delete subBcData[iSub];
    delete [] subBcData;
  }

// Included (MB)
  if (this->dXdot) delete this->dXdot;
  if (this->dTemp) delete this->dTemp;
  if (this->dUface) delete this->dUface;
  if (this->dUnode) delete this->dUnode;
  if (this->dUfaceSA) delete this->dUfaceSA;
  if (this->dUnodeSA) delete this->dUnodeSA;
  if (this->dUinletnode) delete this->dUinletnode;
  if (this->dUfarin) delete this->dUfarin;
  if (this->dUfarout) delete this->dUfarout;
  if (this->dUporouswall) delete this->dUporouswall;

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcData<dim>::finalize(DistSVec<double,3> &X)
{
  //assumption : only one phase is found at the far-field boundary
//  com->printf(2, "Conservative Inlet : %e %e %e\n", Uin[0],Uin[1],Uin[4]);
  this->vf->conservativeToPrimitive(this->Uin, Vin);
  this->vf->conservativeToPrimitive(this->Uout, Vout);

  double Pressure = this->vf->getPressure(Vin);
/*  com->printf(2, "Non-dimensionalized primitive state vector:\n");
  com->printf(2, "Inlet: ");
  for (int k=0; k<dim; ++k)
    com->printf(2, " %g", Vin[k]);
  com->printf(2, "\nOutlet:");
  for (int k=0; k<dim; ++k)
    com->printf(2, " %g", Vout[k]);
  com->printf(2,"\n");
*/
#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->assignFreeStreamValues2(this->Ufarin(iSub),
                                 this->Ufarout(iSub), this->Uface(iSub),
                                 this->Uinletnode(iSub));
    this->subDomain[iSub]->assignPorousWallValues(this->Uporouswall(iSub), this->Uface(iSub));
  }
  update(X);
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistBcData<dim>::finalizeSA(DistSVec<double,3> &X, DistSVec<double,3> &dX, double &dMach)
{

//Remark: Error mesage for pointers
  if (dUface == 0) {
    fprintf(stderr, "*** Error: Variable dUface does not exist!\n");
    exit(1);
  }

  this->vf->conservativeToPrimitiveDerivative(this->Uin, this->dUin, Vin, dVin);
  this->vf->conservativeToPrimitiveDerivative(this->Uout, this->dUout, Vout, dVout);

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->assignFreeStreamValues2((*dUfarin)(iSub),
                                 (*dUfarout)(iSub), (*dUface)(iSub),
                                 (*dUinletnode)(iSub));
    this->subDomain[iSub]->assignPorousWallValues((*dUporouswall)(iSub), (*dUface)(iSub));
  }

  updateSA(X, dX, dMach);

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcData<dim>::update(DistSVec<double,3> &X)  {

#pragma omp parallel for
  for (int iSub=0; iSub<this->numLocSub; ++iSub) {
    double (*xdot)[3] = Xdot.subData(iSub);
    double *temp = Temp.subData(iSub);
    double (*unode)[dim] = this->Unode.subData(iSub);
    int *rotOwn = this->subDomain[iSub]->getRotOwn();
    NodeSet &nodes = this->subDomain[iSub]->getNodes();
    //int count = 0;
    if (rotOwn)  { 
      for (int i=0; i<this->Unode.subSize(iSub); ++i) {
        if (rotOwn[i]>=0) {    // node belongs to a (potential) "rotating" surface
          map<int,RotationData *>::iterator it = rotInfo.find(rotOwn[i]);
          if(it != rotInfo.end()) { // the rotation data have been defined
	    if(it->second->infRadius == RotationData::TRUE) {
              double vel = it->second->omega / vref;
	      unode[i][1] = vel*it->second->nx;
	      unode[i][2] = vel*it->second->ny;
	      unode[i][3] = vel*it->second->nz;
	    } 
            else {
	      double ox = tref*it->second->omega*it->second->nx;
	      double oy = tref*it->second->omega*it->second->ny;
	      double oz = tref*it->second->omega*it->second->nz;
	      double xd = nodes[i][0] - it->second->x0;
	      double yd = nodes[i][1] - it->second->y0;
	      double zd = nodes[i][2] - it->second->z0;
	      unode[i][1] = oy*zd-oz*yd + xdot[i][0];
	      unode[i][2] = oz*xd-ox*zd + xdot[i][1];
	      unode[i][3] = ox*yd-oy*xd + xdot[i][2];
	    }
	  } 
          else  { // no rotation data -> use velocity from mesh motion if any
            unode[i][1] = xdot[i][0];
            unode[i][2] = xdot[i][1];
            unode[i][3] = xdot[i][2];
	  }
        }
        else  { // no rotation data -> use velocity from mesh motion if any
          unode[i][1] = xdot[i][0];
          unode[i][2] = xdot[i][1];
          unode[i][3] = xdot[i][2];
        }
        unode[i][4] = temp[i];
      }
    }
    else { // node does not belong to a "rotating" surface
      for (int i=0; i<this->Unode.subSize(iSub); ++i) {
        unode[i][1] = xdot[i][0];
        unode[i][2] = xdot[i][1];
        unode[i][3] = xdot[i][2];
        unode[i][4] = temp[i];
      }
    }
    //if(count) { fprintf(stderr," In DistBcData<dim>::update(): subd %3d has %6d 'rotating' nodes\n",this->subDomain[iSub]->getGlobSubNum(),count); fflush(stderr); }
#if defined(STRONG_INLET_BC)
    this->subDomain[iSub]->setNodeBcValue(Vin, this->Unode(iSub));
#endif
#if defined(STRONG_FARFIELD_BC)
    this->subDomain[iSub]->setNodeBcValue2(this->Uin, this->Unode(iSub));
#endif
    this->subDomain[iSub]->computeFaceBcValue(this->Unode(iSub), this->Uface(iSub));
  }
  if ( this->gravity > 0.0 )
    updateFarField(X);
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistBcData<dim>::updateSA(DistSVec<double,3> &X, DistSVec<double,3> &dX, double &dMach)
{

//Remark: Error mesage for pointers
  if (dUface == 0) {
    fprintf(stderr, "*** Error: Variable dUface does not exist!\n");
    exit(1);
  }
  if (dUnode == 0) {
    fprintf(stderr, "*** Error: Variable dUnode does not exist!\n");
    exit(1);
  }

#pragma omp parallel for
  for (int iSub=0; iSub<this->numLocSub; ++iSub) {
    double (*dxdot)[3] = dXdot->subData(iSub);
    double *dtemp = dTemp->subData(iSub);
    double (*dunode)[dim] = dUnode->subData(iSub);
    int *rotOwn = this->subDomain[iSub]->getRotOwn();
    NodeSet &nodes = this->subDomain[iSub]->getNodes();
    //int count = 0;
    if (rotOwn)  { 
      for (int i=0; i<dUnode->subSize(iSub); ++i) {
        if(rotOwn[i]>=0) { // node belongs to a (potential) "rotating" surface
          map<int,RotationData *>::iterator it = rotInfo.find(rotOwn[i]);
          if(it != rotInfo.end()) { // the rotation data have been defined
	    if(it->second->infRadius == RotationData::TRUE) {
              double vel = it->second->omega / vref;
              double dvel = -it->second->omega / (vref*vref)*dvrefdMach*dMach;
	      dunode[i][1] = dvel*it->second->nx;
	      dunode[i][2] = dvel*it->second->ny;
	      dunode[i][3] = dvel*it->second->nz;
	    }
	    else {
	      double ox = tref*it->second->omega*it->second->nx;
	      double oy = tref*it->second->omega*it->second->ny;
	      double oz = tref*it->second->omega*it->second->nz;
	      double dox = dtrefdMach*dMach*it->second->omega*it->second->nx;
	      double doy = dtrefdMach*dMach*it->second->omega*it->second->ny;
	      double doz = dtrefdMach*dMach*it->second->omega*it->second->nz;
	      double xd = nodes[i][0] - it->second->x0;
	      double yd = nodes[i][1] - it->second->y0;
	      double zd = nodes[i][2] - it->second->z0;
	      dunode[i][1] = doy*zd-doz*yd + dxdot[i][0];
	      dunode[i][2] = doz*xd-dox*zd + dxdot[i][1];
	      dunode[i][3] = dox*yd-doy*xd + dxdot[i][2];
	    }
	  }
	  else {
            dunode[i][1] = dxdot[i][0];
            dunode[i][2] = dxdot[i][1];
            dunode[i][3] = dxdot[i][2];
	  }
        }
        dunode[i][4] = dtemp[i];
      }
    }
    else { // node does not belong to a "rotating" surface
      for (int i=0; i<dUnode->subSize(iSub); ++i) {
        dunode[i][1] = dxdot[i][0];
        dunode[i][2] = dxdot[i][1];
        dunode[i][3] = dxdot[i][2];
        dunode[i][4] = dtemp[i];
      }
    }
    //if(count) { fprintf(stderr," In DistBcData<dim>::update(): subd %3d has %6d 'rotating' nodes\n",subDomain[iSub]->getGlobSubNum(),count); fflush(stderr); }
#if defined(STRONG_INLET_BC)
    this->subDomain[iSub]->setNodeBcValue(dVin, (*dUnode)(iSub));
#endif
#if defined(STRONG_FARFIELD_BC)
    this->subDomain[iSub]->setNodeBcValue2(this->dUin, (*dUnode)(iSub));
#endif
    this->subDomain[iSub]->computeFaceBcValue((*dUnode)(iSub), (*dUface)(iSub));
  }

  if ( this->gravity > 0.0 )
    updateFarFieldSA(X, dX, dMach);
}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::updateFarField(DistSVec<double,3> &X)
{
// it does not matter what kind of simulation is done.
// we only need to know which fluid is assumed to be at the
//    boundary.
  if (this->boundaryFluid == this->GAS)
    updateFarFieldGas(X);

  else if(this->boundaryFluid == this->TAIT)
    updateFarFieldLiquid(X);

  else if(this->boundaryFluid == this->JWL)
    updateFarFieldJWL(X);

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->assignFreeStreamValues2(this->Ufarin(iSub),
                                 this->Ufarout(iSub), this->Uface(iSub),
                                 this->Uinletnode(iSub));
    this->subDomain[iSub]->assignPorousWallValues(this->Uporouswall(iSub), this->Uface(iSub));
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistBcDataEuler<dim>::updateFarFieldSA(DistSVec<double,3> &X, DistSVec<double,3> &dX, double &dMach)
{
// it does not matter what kind of simulation is done.
// we only need to know which fluid is assumed to be at the
//    boundary.
  if (this->boundaryFluid == this->GAS) {
    fprintf(stderr, " ... boundary fluid is GAS!\n");
    updateFarFieldGasSA(X, dX, dMach);
  } else if(this->boundaryFluid == this->TAIT || this->boundaryFluid == this->JWL) {
    fprintf(stderr, "*** Error: Function updateFarFieldSA (at DistBcData.C) does not support this option!\n");
    exit(1);
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->assignFreeStreamValues2((*this->dUfarin)(iSub),
                                 (*this->dUfarout)(iSub), (*this->dUface)(iSub),
                                 (*this->dUinletnode)(iSub));
    this->subDomain[iSub]->assignPorousWallValues((*this->dUporouswall)(iSub), (*this->dUface)(iSub));
  }
}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::updateFarFieldGas(DistSVec<double,3> &X)
{

  // flow properties
  double gam = this->vf->getGamma();
  double Pstiff = this->vf->getPressureConstant();

#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double ptempin, ptempout, un, velin2, velout2;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = this->Vin[4] + this->Vin[0] *un;
      ptempout = this->Vout[4]+ this->Vout[0]*un;
      velin2  = this->Vin[1]*this->Vin[1]+this->Vin[2]*this->Vin[2]+this->Vin[3]*this->Vin[3];
      velout2 = this->Vout[1]*this->Vout[1]+this->Vout[2]*this->Vout[2]+this->Vout[3]*this->Vout[3];

      uin[inode][4] = 0.5*this->Vin[0]*velin2 + (ptempin+gam*Pstiff)/(gam-1.0);
      uout[inode][4] = 0.5*this->Vout[0]*velout2 + (ptempout+gam*Pstiff)/(gam-1.0);

    }

  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistBcDataEuler<dim>::updateFarFieldGasSA(DistSVec<double,3> &X, DistSVec<double,3> &dX, double &dMach)
{

  // flow properties
  double gam = this->vf->getGamma();
  double Pstiff = this->vf->getPressureConstant();
  double dPstiff = this->vf->getDerivativeOfPressureConstant()*dMach;

//  assert(this->Vin[0]>0.0);
//  assert(this->Vout[0]>0.0);

#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]       = X.subData(iSub);
    double (*dx)[3]      = dX.subData(iSub);
    double (*duin)[dim]  = this->dUfarin->subData(iSub);
    double (*duout)[dim] = this->dUfarout->subData(iSub);
    double ptempin, dptempin, ptempout, dptempout, un, dun, velin2, dvelin2, velout2, dvelout2;

    for(int inode = 0; inode<this->dUnode->subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      dun = (dx[inode][0]*this->ngravity[0]+dx[inode][1]*this->ngravity[1]+dx[inode][2]*this->ngravity[2]);
      ptempin  = this->Vin[4] + this->Vin[0] *un;
      ptempout = this->Vout[4]+ this->Vout[0]*un;
      velin2  = this->Vin[1]*this->Vin[1]+this->Vin[2]*this->Vin[2]+this->Vin[3]*this->Vin[3];
      velout2 = this->Vout[1]*this->Vout[1]+this->Vout[2]*this->Vout[2]+this->Vout[3]*this->Vout[3];
      dptempin  = this->dVin[4] + this->dVin[0] *un + this->Vin[0] *dun;
      dptempout = this->dVout[4] + this->dVout[0]*un + this->Vout[0]*dun;
      dvelin2  = 2.0*(this->Vin[1]*this->dVin[1]+this->Vin[2]*this->dVin[2]+this->Vin[3]*this->dVin[3]);
      dvelout2 = 2.0*(this->Vout[1]*this->dVout[1]+this->Vout[2]*this->dVout[2]+this->Vout[3]*this->dVout[3]);

      duin[inode][4] = 0.5*this->dVin[0]*velin2 + 0.5*this->Vin[0]*dvelin2 + (dptempin+gam*dPstiff)/(gam-1.0);
      duout[inode][4] = 0.5*this->dVout[0]*velout2 + 0.5*this->Vout[0]*dvelout2 + (dptempout+gam*dPstiff)/(gam-1.0);

    }

  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::updateFarFieldLiquid(DistSVec<double,3> &X)
{

  //flow properties
  double a = this->vf->getAlphaWater();
  double b = this->vf->getBetaWater();
  double P = this->vf->getPrefWater();
  double c = this->vf->getCv();
  double coeff = (b-1.0)/(a*b);

#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double rtempin, rtempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      rtempin  = pow(this->Vin[0], b-1.0)+coeff*(un-this->gravity*this->depth);
      rtempout = pow(this->Vout[0],b-1.0)+coeff*(un-this->gravity*this->depth);
      rtempin  = pow(rtempin, 1.0/(b-1.0));
      rtempout = pow(rtempout,1.0/(b-1.0));

      uin[inode][0] = rtempin;
      uin[inode][1] = rtempin*this->Vin[1];
      uin[inode][2] = rtempin*this->Vin[2];
      uin[inode][3] = rtempin*this->Vin[3];
      uin[inode][4] = rtempin*this->Uin[4]/this->Vin[0];

      uout[inode][0] = rtempout;
      uout[inode][1] = rtempout*this->Vout[1];
      uout[inode][2] = rtempout*this->Vout[2];
      uout[inode][3] = rtempout*this->Vout[3];
      uout[inode][4] = rtempout*this->Uout[4]/this->Vout[0];

    }

  }
}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::updateFarFieldJWL(DistSVec<double,3> &X)
{

  // flow properties
  double omega  = this->vf->getOmega();

#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double ptempin, ptempout, un, velin2, velout2;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = this->Vin[4] + this->Vin[0] *un;
      ptempout = this->Vout[4]+ this->Vout[0]*un;
      velin2  = this->Vin[1]*this->Vin[1]+this->Vin[2]*this->Vin[2]+this->Vin[3]*this->Vin[3];
      velout2 = this->Vout[1]*this->Vout[1]+this->Vout[2]*this->Vout[2]+this->Vout[3]*this->Vout[3];

      uin[inode][4] = 0.5*this->Vin[0]*velin2 + (ptempin-this->vf->computeFrho(this->Vin[0]))/omega;
      uout[inode][4] = 0.5*this->Vout[0]*velout2 + (ptempout-this->vf->computeFrho(this->Vout[0]))/omega;

    }

  }

}

//------------------------------------------------------------------------------
template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsGas(IoData &iod,
                                DistSVec<double,3> &X)
{

/* case 1: no gravity (typically aero simulation)
 *   pressure and density are constants and imposed as are
 *   velocity is constant but imposed as Mach*speed of sound
 *
 * case 2: gravity and depth (typically hydro simulation)
 *   density is constant and imposed as is
 *   pressure is a function of the depth, 
 *          ie P = P_surface + rho*g*(h-z)
 *   velocity is constant and imposed as Mach*speed of sound at depth "depth"
 *
 */

  // flow properties
  double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  double Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant;
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;
  double pressurein  = iod.bc.inlet.pressure + rhoin*this->gravity*this->depth;
  double pressureout = iod.bc.outlet.pressure + rhoout*this->gravity*this->depth;

  double velin2 = 0.0;
  double velout2 = 0.0;

  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = gam * (pressurein + Pstiff)*
      iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin;
    velout2 = gam * (pressureout + Pstiff)*
      iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout;
  }else if(iod.bc.inlet.velocity >= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
           iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin = sqrt(velin2);
  double velout= sqrt(velout2);
  
// computation of boundary values "on average", ie at node of coordinates (0,0,0)
  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0] * velin * sin(iod.bc.inlet.alpha);
  this->Uin[4] = (pressurein+gam*Pstiff)/(gam - 1.0) + 0.5 * this->Uin[0] * velin2;

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0] * velout * sin(iod.bc.outlet.alpha);
  this->Uout[4] = (pressureout+gam*Pstiff)/(gam - 1.0) + 0.5 * this->Uout[0] * velout2;

// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation 
#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double (*uporouswall)[dim] = this->Uporouswall.subData(iSub);
    double ptempin, ptempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = pressurein  + rhoin *un;
      ptempout = pressureout + rhoout*un;

      uin[inode][0] = this->Uin[0];
      uin[inode][1] = this->Uin[1];
      uin[inode][2] = this->Uin[2];
      uin[inode][3] = this->Uin[3];
      uin[inode][4] = 0.5*rhoin*velin2 + (ptempin+gam*Pstiff)/(gam-1.0); 

      uout[inode][0] = this->Uout[0];
      uout[inode][1] = this->Uout[1];
      uout[inode][2] = this->Uout[2];
      uout[inode][3] = this->Uout[3];
      uout[inode][4] = 0.5*rhoout*velout2 + (ptempout+gam*Pstiff)/(gam-1.0);


// pressure pulse for channel flow along x
      
//      double amplitude = 5.0;
//      double dx = 0.075*0.075;
//     double radius2 = (x[inode][2]-0.125)*(x[inode][2]-0.125);// +
//                       (x[inode][1]-0.)*(x[inode][1]-0.) +
//                       (x[inode][2]-0.125)*(x[inode][2]-0.125);
//      uin[inode][0] = rhoin*(1.0 + amplitude*exp(-radius2/dx));
//      uin[inode][4] = 0.5*uin[inode][0]*velin2 +(ptempin+gam*Pstiff)/(gam-1.0);
//      /*double amplitude = 0.2;
//      double dx = 0.1*0.1;
//      double radius2 = (x[inode][0]-0.)*(x[inode][0]-0.) +
//                       (x[inode][1]-0.)*(x[inode][1]-0.) +
//                       (x[inode][2]-0.)*(x[inode][2]-0.);
//      ptempin = pressurein*(1.0 + amplitude*exp( -pow((x[inode][0]-0.5),2)/dx));
//      ptempin = pressurein*(1.0 + amplitude*exp(-radius2/dx));
//      uin[inode][4] = 0.5*rhoin*velin2 +(ptempin+gam*Pstiff)/(gam-1.0);
//      */

    }

    FaceSet& faces = this->subDomain[iSub]->getFaces();
    for (int i=0;i<faces.size(); i++) { //loop over faces 
      map<int,SurfaceData *> &surfaceMap = iod.surfaces.surfaceMap.dataMap;
      map<int,SurfaceData*>::iterator it = surfaceMap.find(faces[i].getSurfaceID());
      if(it!=surfaceMap.end()) { // surface has attribut in the input file
        map<int,BoundaryData *> &bcMap = iod.bc.bcMap.dataMap;
        map<int,BoundaryData *>::iterator it2 = bcMap.find(it->second->bcID);
        if(it2 != bcMap.end()) { // the bc data have been defined
          if(it2->second->type == BoundaryData::DIRECTSTATE || 
             it2->second->type == BoundaryData::MASSFLOW) {
            if (faces[i].getCode() ==  BC_DIRECTSTATE_INLET_MOVING ||
                faces[i].getCode() ==  BC_DIRECTSTATE_INLET_FIXED ) {
              for (int l = 0; l<faces[i].numNodes();++l) {
                int k = faces[i][l];
                uin[k][0] = it2->second->totalTemperature;
                uin[k][1] = it2->second->velocityX;
                uin[k][2] = it2->second->velocityY;
                uin[k][3] = it2->second->velocityZ;
                uin[k][4] = it2->second->totalPressure;
              }
            }
            if (faces[i].getCode() ==  BC_MASSFLOW_INLET_MOVING ||
                faces[i].getCode() ==  BC_MASSFLOW_INLET_FIXED ) {
            }

            if (faces[i].getCode() ==  BC_DIRECTSTATE_OUTLET_MOVING ||
                faces[i].getCode() ==  BC_DIRECTSTATE_OUTLET_FIXED ) {
              for (int l = 0; l<faces[i].numNodes();++l) {
                int k = faces[i][l];
                uout[k][4] = it2->second->pressure;
              }
            }
            if (faces[i].getCode() ==  BC_MASSFLOW_OUTLET_MOVING ||
                faces[i].getCode() ==  BC_MASSFLOW_OUTLET_FIXED ) {
            }
          }

          else if(it2->second->type == BoundaryData::POROUSWALL ) {
            for (int l = 0; l<faces[i].numNodes();++l) {
              int k = faces[i][l];
              uporouswall[k][0] = it2->second->mdot;
              uporouswall[k][1] = it2->second->porosity;
              uporouswall[k][4] = it2->second->temperature;
            }
          }
        }
      }
    }

  }
  for (int idim=0; idim<dim; idim++)
    this->Ub[idim] = 0.0;

}

//------------------------------------------------------------------------------


// Included (MB)
template<int dim>
void DistBcDataEuler<dim>::setDerivativeOfBoundaryConditionsGas(IoData &iod,
                                DistSVec<double,3> &X, DistSVec<double,3> &dX, double dM, double dA, double dB)
{

 // flow properties
  double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  double Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant;

  double dPstiff = -iod.eqs.fluidModel.gasModel.pressureConstant*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach)) * dM;
//double dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant*(-2.0 / (iod.bc.inlet.mach)) * dM;//TODO WRONG

  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;
  double pressurein  = iod.bc.inlet.pressure + rhoin*this->gravity*this->depth;

//double dpressurein = iod.bc.inlet.pressure * (-2.0 / (iod.bc.inlet.mach)) * dM;//TODO WRONG
  double dpressurein = -2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach) * dM;

  double pressureout = iod.bc.outlet.pressure + rhoout*this->gravity*this->depth;

//double dpressureout = iod.bc.outlet.pressure * (-2.0 / (iod.bc.outlet.mach)) * dM;//TODO WRONG
  double dpressureout = -2.0 / (gam * iod.bc.outlet.mach * iod.bc.outlet.mach * iod.bc.outlet.mach) * dM;

  double velin2 = gam * (pressurein + Pstiff) * iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin;

  double dvelin2 = (gam * dPstiff * iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin + 2.0 * gam * pressurein * iod.bc.inlet.mach / rhoin - 2.0 / (rhoin*iod.bc.inlet.mach)) * dM;
//double dvelin2 = (gam * (dpressurein + dPstiff) * iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin + 2.0 * gam * (pressurein+Pstiff) * iod.bc.inlet.mach / rhoin) * dM;//TODO WRONG

  double velout2 = gam * (pressureout + Pstiff) * iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout;

  double dvelout2 = (gam * dPstiff * iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout + 2.0 * gam * pressureout * iod.bc.outlet.mach / rhoout  - 2.0 / (rhoout*iod.bc.outlet.mach)) * dM;
//  double dvelout2 = (gam * (dpressureout + dPstiff) * iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout + 2.0 * gam * (pressureout+Pstiff) * iod.bc.outlet.mach / rhoout) * dM;//TODO wrong

  double velin = sqrt(velin2);
  double dvelin = dvelin2/(2*velin);
  double velout = sqrt(velout2);
  double dvelout = dvelout2/(2.0*velout);

  double alpha_in  = iod.bc.inlet.alpha;
  double alpha_out = iod.bc.outlet.alpha;
  double beta_in   = iod.bc.inlet.beta;
  double beta_out  = iod.bc.outlet.beta;

  this->dUin[0] = 0.0;
  this->dUin[1] = rhoin * dvelin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->dUin[2] = rhoin * dvelin * cos(iod.bc.inlet.alpha)  * sin(iod.bc.inlet.beta);
  this->dUin[3] = rhoin * dvelin * sin(iod.bc.inlet.alpha);
  this->dUin[4] = (dpressurein + gam*dPstiff)/(gam - 1.0) + 0.5 * rhoin * dvelin2;

  this->dUout[0] = 0.0;
  this->dUout[1] = rhoout * dvelout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->dUout[2] = rhoout * dvelout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->dUout[3] = rhoout * dvelout * sin(iod.bc.outlet.alpha);
  this->dUout[4] = (dpressureout + gam*dPstiff)/(gam - 1.0)  + 0.5 * rhoout * dvelout2;

  this->dUin[0] += 0.0;
  this->dUin[1] += rhoin * velin * (-sin(iod.bc.inlet.alpha)*dA) * cos(iod.bc.inlet.beta);
  this->dUin[2] += rhoin * velin * (-sin(iod.bc.inlet.alpha)*dA) * sin(iod.bc.inlet.beta);
  this->dUin[3] += rhoin * velin * (cos(iod.bc.inlet.alpha)*dA);
  this->dUin[4] += 0.0;

  this->dUout[0] += 0.0;
  this->dUout[1] += rhoout * velout * (-sin(iod.bc.outlet.alpha)*dA) * cos(iod.bc.outlet.beta);
  this->dUout[2] += rhoout * velout * (-sin(iod.bc.outlet.alpha)*dA) * sin(iod.bc.outlet.beta);
  this->dUout[3] += rhoout * velout * (cos(iod.bc.outlet.alpha)*dA);
  this->dUout[4] += 0.0;

  this->dUin[0] += 0.0;
  this->dUin[1] += rhoin * velin * cos(iod.bc.inlet.alpha) * (-sin(iod.bc.inlet.beta)*dB);
  this->dUin[2] += rhoin * velin * cos(iod.bc.inlet.alpha) * (cos(iod.bc.inlet.beta)*dB);
  this->dUin[3] += 0.0;
  this->dUin[4] += 0.0;

  this->dUout[0] += 0.0;
  this->dUout[1] += rhoout * velout * cos(iod.bc.outlet.alpha) * (-sin(iod.bc.outlet.beta)*dB);
  this->dUout[2] += rhoout * velout * cos(iod.bc.outlet.alpha) * (cos(iod.bc.outlet.beta)*dB);
  this->dUout[3] += 0.0;
  this->dUout[4] += 0.0;

// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation
#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]       = X.subData(iSub);
    double (*dx)[3]      = dX.subData(iSub);
    double (*duin)[dim]  = this->dUfarin->subData(iSub);
    double (*duout)[dim] = this->dUfarout->subData(iSub);
    double ptempin, ptempout, dptempin, dptempout, un, dun;

    for(int inode = 0; inode<this->dUnode->subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      dun = (dx[inode][0]*this->ngravity[0]+dx[inode][1]*this->ngravity[1]+dx[inode][2]*this->ngravity[2]);
      ptempin  = pressurein  + rhoin *un;
      ptempout = pressureout + rhoout*un;
      dptempin  = dpressurein + rhoin*dun;
      dptempout = dpressureout + rhoout*dun;

      duin[inode][0] = this->dUin[0];
      duin[inode][1] = this->dUin[1];
      duin[inode][2] = this->dUin[2];
      duin[inode][3] = this->dUin[3];
      duin[inode][4] = 0.5*rhoin*dvelin2 + (dptempin+gam*dPstiff)/(gam-1.0);

      duout[inode][0] = this->dUout[0];
      duout[inode][1] = this->dUout[1];
      duout[inode][2] = this->dUout[2];
      duout[inode][3] = this->dUout[3];
      duout[inode][4] = 0.5*rhoout*dvelout2 + (dptempout+gam*dPstiff)/(gam-1.0);

    }
  }

}

////TODO this is the working version so far
//// Included (MB)
//template<int dim>
//void DistBcDataEuler<dim>::setDerivativeOfBoundaryConditionsGas(IoData &iod,
//                                DistSVec<double,3> &X, DistSVec<double,3> &dX, double dM, double dA, double dB)
//{
//
// // reading flow properties from input
//  double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;// \f$(x_1,y_1)\f$ //
//  double Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant;
//  double rhoin = iod.bc.inlet.density;
//  double rhoout = iod.bc.outlet.density;
//  double pressurein  = iod.bc.inlet.pressure + rhoin*this->gravity*this->depth;
//  double pressureout = iod.bc.outlet.pressure + rhoout*this->gravity*this->depth;
//
//  //squared values of
//  double velin2 = gam * (pressurein + Pstiff) * iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin;
//  double velout2 = gam * (pressureout + Pstiff) * iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout;
//
//  double alpha_in  = iod.bc.inlet.alpha;
//  double alpha_out = iod.bc.outlet.alpha;
//  double beta_in   = iod.bc.inlet.beta;
//  double beta_out  = iod.bc.outlet.beta;
//
//  double dvelin2=0.0;
//  double dvelout2=0.0;
//  double dpressurein=0.0;
//  double dpressureout=0.0;
//  double dPstiff=0.0;
//  if(true){//new derivatives
//    double dPstiff = -iod.eqs.fluidModel.gasModel.pressureConstant*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach)) * dM;
//    double dpressurein = -2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach) * dM;
//    double dpressureout = -2.0 / (gam * iod.bc.outlet.mach * iod.bc.outlet.mach * iod.bc.outlet.mach) * dM;
//    double dvelin2 = (gam * dPstiff * iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin + 2.0 * gam * pressurein * iod.bc.inlet.mach / rhoin - 2.0 / (rhoin*iod.bc.inlet.mach)) * dM;
//    double dvelout2 = (gam * dPstiff * iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout + 2.0 * gam * pressureout * iod.bc.outlet.mach / rhoout  - 2.0 / (rhoout*iod.bc.outlet.mach)) * dM;
//  }else{//old derivatives
//    double dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant*(-2.0 / (iod.bc.inlet.mach)) * dM;
//    double dpressurein = iod.bc.inlet.pressure * (-2.0 / (iod.bc.inlet.mach)) * dM;
//    double dpressureout = iod.bc.outlet.pressure * (-2.0 / (iod.bc.outlet.mach)) * dM;
//    double dvelin2 = (gam * (dpressurein + dPstiff) * iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin + 2.0 * gam * (pressurein+Pstiff) * iod.bc.inlet.mach / rhoin) * dM;
//    double dvelout2 = (gam * (dpressureout + dPstiff) * iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout + 2.0 * gam * (pressureout+Pstiff) * iod.bc.outlet.mach / rhoout) * dM;
//  }
//
//  double velin = sqrt(velin2);
//  double dvelin = dvelin2/(2*velin);
//  double velout = sqrt(velout2);
//  double dvelout = dvelout2/(2.0*velout);
//
//  //Initialization
//  this->dUin[0] = 0.0;
//  this->dUin[1] = rhoin * dvelin * cos(alpha_in) * cos(beta_in);
//  this->dUin[2] = rhoin * dvelin * cos(alpha_in)  * sin(beta_in);
//  this->dUin[3] = rhoin * dvelin * sin(alpha_in);
//  this->dUin[4] = (dpressurein + gam*dPstiff)/(gam - 1.0) + 0.5 * rhoin * dvelin2;
//
//  //TODO why is there an independent derivative?
//  this->dUout[0] = 0.0;
//  this->dUout[1] = rhoout * dvelout * cos(alpha_out) * cos(beta_out);
//  this->dUout[2] = rhoout * dvelout * cos(alpha_out) * sin(beta_out);
//  this->dUout[3] = rhoout * dvelout * sin(alpha_out);
//  this->dUout[4] = (dpressureout + gam*dPstiff)/(gam - 1.0)  + 0.5 * rhoout * dvelout2;
//
//  //derivative contribution in case of alpha sensitivity
//  if(dA!=0.0){
//	this->dUin[0] += 0.0;
//	this->dUin[1] += rhoin * velin * (-sin(alpha_in)*dA) * cos(beta_in);
//	this->dUin[2] += rhoin * velin * (-sin(alpha_in)*dA) * sin(beta_in);
//	this->dUin[3] += rhoin * velin * (cos(alpha_in)*dA);
//	this->dUin[4] += 0.0;
//
//	this->dUout[0] += 0.0;
//	this->dUout[1] += rhoout * velout * (-sin(alpha_out)*dA) * cos(beta_out);
//	this->dUout[2] += rhoout * velout * (-sin(alpha_out)*dA) * sin(beta_out);
//	this->dUout[3] += rhoout * velout * (cos(alpha_out)*dA);
//	this->dUout[4] += 0.0;
//  }
//
//  //derivative contribution in case of beta sensitivity
//  if(dB!=0.0){
//    this->dUin[0] += 0.0;
//    this->dUin[1] += rhoin * velin * cos(alpha_in) * (-sin(beta_in)*dB);
//    this->dUin[2] += rhoin * velin * cos(alpha_in) * (cos(beta_in)*dB);
//    this->dUin[3] += 0.0;
//    this->dUin[4] += 0.0;
//
//    this->dUout[0] += 0.0;
//    this->dUout[1] += rhoout * velout * cos(alpha_out) * (-sin(beta_out)*dB);
//    this->dUout[2] += rhoout * velout * cos(alpha_out) * (cos(beta_out)*dB);
//    this->dUout[3] += 0.0;
//    this->dUout[4] += 0.0;
//  }
//
//// computation for each node according to its depth
//// this will be passed in DistTimeState to initialize simulation
//#pragma omp parallel for
//  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
//    double (*x)[3]       = X.subData(iSub);
//    double (*dx)[3]      = dX.subData(iSub);
//    double (*duin)[dim]  = this->dUfarin->subData(iSub);
//    double (*duout)[dim] = this->dUfarout->subData(iSub);
//    double ptempin, ptempout, dptempin, dptempout, un, dun;
//
//    for(int inode = 0; inode<this->dUnode->subSize(iSub); inode++){
//      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
//      dun = (dx[inode][0]*this->ngravity[0]+dx[inode][1]*this->ngravity[1]+dx[inode][2]*this->ngravity[2]);
//      ptempin  = pressurein  + rhoin *un;
//      ptempout = pressureout + rhoout*un;
//      dptempin  = dpressurein + rhoin*dun;
//      dptempout = dpressureout + rhoout*dun;
//
//      duin[inode][0] = this->dUin[0];
//      duin[inode][1] = this->dUin[1];
//      duin[inode][2] = this->dUin[2];
//      duin[inode][3] = this->dUin[3];
//      duin[inode][4] = 0.5*rhoin*dvelin2 + (dptempin+gam*dPstiff)/(gam-1.0);
//
//      duout[inode][0] = this->dUout[0];
//      duout[inode][1] = this->dUout[1];
//      duout[inode][2] = this->dUout[2];
//      duout[inode][3] = this->dUout[3];
//      duout[inode][4] = 0.5*rhoout*dvelout2 + (dptempout+gam*dPstiff)/(gam-1.0);
//
//    }
//  }
//
//}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsLiquid(IoData &iod,
                                DistSVec<double,3> &X)
{
  // flow properties
  double a = this->vf->getAlphaWater();
  double b = this->vf->getBetaWater();
  double c = this->vf->getCv();
  double P = this->vf->getPrefWater();
  double coeff = -(b-1.0)/(a*b);
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = iod.bc.inlet.mach*iod.bc.inlet.mach * a*b*pow(rhoin, b-1.0);
    velout2 = iod.bc.outlet.mach*iod.bc.outlet.mach * a*b*pow(rhoout, b-1.0);
  }else if (iod.bc.inlet.velocity >= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);
  double pin    = P + a*pow(rhoin, b);
  double pout   = P + a*pow(rhoout, b);

// computation of boundary values "on average", ie at node of coordinates (0,0,0)
  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0] * velin * sin(iod.bc.inlet.alpha);
  this->Uin[4] = this->Uin[0]*(c*iod.bc.inlet.temperature + 0.5 * velin2) - pin;

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0] * velout * sin(iod.bc.outlet.alpha);
  this->Uout[4] = this->Uout[0]*(c*iod.bc.outlet.temperature + 0.5 * velout2) - pout;


// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation 
#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double rtempin, rtempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      rtempin  = pow(rhoin ,b-1.0)+coeff*(un-this->gravity*this->depth);
      rtempout = pow(rhoout,b-1.0)+coeff*(un-this->gravity*this->depth);
      rtempin  = pow(rtempin, 1.0/(b-1.0));
      rtempout = pow(rtempout,1.0/(b-1.0));

      uin[inode][0] = rtempin;
      uin[inode][1] = rtempin*this->Uin[1]/this->Uin[0];
      uin[inode][2] = rtempin*this->Uin[2]/this->Uin[0];
      uin[inode][3] = rtempin*this->Uin[3]/this->Uin[0];
      uin[inode][4] = rtempin*this->Uin[4]/this->Uin[0];

      uout[inode][0] = rtempout;
      uout[inode][1] = rtempout*this->Uout[1]/this->Uout[0];
      uout[inode][2] = rtempout*this->Uout[2]/this->Uout[0];
      uout[inode][3] = rtempout*this->Uout[3]/this->Uout[0];
      uout[inode][4] = rtempout*this->Uout[4]/this->Uout[0];

    }

  }
  for (int idim=0; idim<dim; idim++)
    this->Ub[idim] = 0.0;

}

//------------------------------------------------------------------------------
template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsJWL(IoData &iod,
                                DistSVec<double,3> &X)
{

/* case 1: no gravity (typically aero simulation)
 *   pressure and density are constants and imposed as are
 *   velocity is constant but imposed as Mach*speed of sound
 *
 * case 2: gravity and depth (typically hydro simulation)
 *   density is constant and imposed as is
 *   pressure is a function of the depth, 
 *          ie P = P_surface + rho*g*(h-z)
 *   velocity is constant and imposed as Mach*speed of sound at depth "depth"
 *
 */

  // flow properties
  double omega  = iod.eqs.fluidModel.jwlModel.omega;
  double A1     = iod.eqs.fluidModel.jwlModel.A1;
  double A2     = iod.eqs.fluidModel.jwlModel.A2;
  double R1     = iod.eqs.fluidModel.jwlModel.R1;
  double R2     = iod.eqs.fluidModel.jwlModel.R2;
  double rhoref = iod.eqs.fluidModel.jwlModel.rhoref;
  double R1r = R1*rhoref;
  double R2r = R2*rhoref;

  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;
  double pressurein  = iod.bc.inlet.pressure + rhoin*this->gravity*this->depth;
  double pressureout = iod.bc.outlet.pressure + rhoout*this->gravity*this->depth;

  double velin2  = 0.0;
  double velout2 = 0.0;
  double frhoin, frhoout, frhop;

  //compute the speed of sounds based on inlet and outlet conditions
  frhoin   = A1*(1.0-omega*rhoin/R1r)*exp(-R1r/rhoin) + A2*(1.0-omega*rhoin/R2r)*exp(-R2r/rhoin);
  frhop    = A1*(-omega/R1r+(1.0-omega*rhoin/R1r)*R1r/(rhoin*rhoin))*exp(-R1r/rhoin) 
           + A2*(-omega/R2r+(1.0-omega*rhoin/R2r)*R2r/(rhoin*rhoin))*exp(-R2r/rhoin); 
  double cin2  = (omega+1.0)*pressurein/rhoin - frhoin/rhoin + frhop;
  frhoout  = A1*(1.0-omega*rhoout/R1r)*exp(-R1r/rhoout) + A2*(1.0-omega*rhoout/R2r)*exp(-R2r/rhoout);
  frhop    = A1*(-omega/R1r+(1.0-omega*rhoout/R1r)*R1r/(rhoout*rhoout))*exp(-R1r/rhoout) 
           + A2*(-omega/R2r+(1.0-omega*rhoout/R2r)*R2r/(rhoout*rhoout))*exp(-R2r/rhoout); 
  double cout2 = (omega+1.0)*pressureout/rhoout - frhoout/rhoout + frhop;

  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2  = iod.bc.inlet.mach * iod.bc.inlet.mach  * cin2;
    velout2 = iod.bc.outlet.mach* iod.bc.outlet.mach * cout2;
  }else if(iod.bc.inlet.velocity >= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
           iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2  = iod.bc.inlet.velocity * iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity* iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin = sqrt(velin2);
  double velout= sqrt(velout2);
  
// computation of boundary values "on average", ie at node of coordinates (0,0,0)
  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0] * velin * sin(iod.bc.inlet.alpha);
  this->Uin[4] = (pressurein-frhoin)/omega + 0.5 * this->Uin[0] * velin2;

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0] * velout * sin(iod.bc.outlet.alpha);
  this->Uout[4] = (pressureout-frhoout)/omega + 0.5 * this->Uout[0] * velout2;

// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation 
#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double ptempin, ptempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = pressurein  + rhoin *un;
      ptempout = pressureout + rhoout*un;

      uin[inode][0] = this->Uin[0];
      uin[inode][1] = this->Uin[1];
      uin[inode][2] = this->Uin[2];
      uin[inode][3] = this->Uin[3];
      uin[inode][4] = 0.5*rhoin*velin2 + (ptempin-frhoin)/omega;

      uout[inode][0] = this->Uout[0];
      uout[inode][1] = this->Uout[1];
      uout[inode][2] = this->Uout[2];
      uout[inode][3] = this->Uout[3];
      uout[inode][4] = 0.5*rhoout*velout2 + (ptempout-frhoout)/omega;

    }
  }
  for (int idim=0; idim<dim; idim++)
    this->Ub[idim] = 0.0;

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsGasGas(IoData &iod,
                                DistSVec<double,3> &X)
{

// this->Uin set up with inlet values and properties of fluidModel1 = GasModel1
  double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  double Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant;
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;
  double pressurein = iod.bc.inlet.pressure + rhoin*this->gravity*this->depth;
  double pressureout = iod.bc.outlet.pressure+ rhoout*this->gravity*this->depth;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = gam * (pressurein+Pstiff) * 
      iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin;
    velout2 = gam * (pressureout+ Pstiff)*
      iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout;
  }else if (iod.bc.inlet.velocity>= 0.0 && iod.bc.outlet.velocity >= 0.0  &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

  //std::cout << "velout = " << iod.bc.outlet.velocity << std::endl;

  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0]*velin*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0]*velin*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0]*velin*sin(iod.bc.inlet.alpha);
  this->Uin[4] = (pressurein+gam*Pstiff)/(gam-1.0) + 0.5 * this->Uin[0] * velin2;

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0]*velout*cos(iod.bc.outlet.alpha)*cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0]*velout*cos(iod.bc.outlet.alpha)*sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0]*velout*sin(iod.bc.outlet.alpha);
  this->Uout[4] = (pressureout+gam*Pstiff)/(gam-1.0) + 0.5 * this->Uout[0] * velout2;

#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double ptempin, ptempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = pressurein + rhoin *un;
      ptempout = pressureout+ rhoout*un;

      uin[inode][0] = this->Uin[0];
      uin[inode][1] = this->Uin[1];
      uin[inode][2] = this->Uin[2];
      uin[inode][3] = this->Uin[3];
      uin[inode][4] = 0.5*this->Uin[0]*velin2 + (ptempin+gam*Pstiff)/(gam-1.0); 

      uout[inode][0] = this->Uout[0];
      uout[inode][1] = this->Uout[1];
      uout[inode][2] = this->Uout[2];
      uout[inode][3] = this->Uout[3];
      uout[inode][4] = 0.5*this->Uout[0]*velout2 + (ptempout+gam*Pstiff)/(gam-1.0);

    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsLiquidLiquid(IoData &iod,
                                DistSVec<double,3> &X)
{

  // flow properties for fluid at boundary
  double a = this->vf->getAlphaWater();
  double b = this->vf->getBetaWater();
  double c = this->vf->getCv();
  double P = this->vf->getPrefWater();
  double coeff = (b-1.0)/(a*b);
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = iod.bc.inlet.mach*iod.bc.inlet.mach * a*b*pow(rhoin, b-1.0);
    velout2 = iod.bc.outlet.mach*iod.bc.outlet.mach * a*b*pow(rhoout, b-1.0);
  }else if (iod.bc.inlet.velocity>= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

// Uin and Uout set up with iod.bc.let.values and properties of fluidModel1 = LiquidModel1
// computation of boundary values "on average", ie at node of coordinates (0,0,0)
  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0] * velin * sin(iod.bc.inlet.alpha);
  this->Uin[4] = this->Uin[0]*(c*iod.bc.inlet.temperature + 0.5 * velin2);

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0] * velout * sin(iod.bc.outlet.alpha);
  this->Uout[4] = this->Uout[0]*(c*iod.bc.outlet.temperature + 0.5 * velout2);


// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation
#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double rtempin, rtempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      rtempin  = pow(rhoin ,b-1.0)-coeff*(un-this->gravity*this->depth);
      rtempout = pow(rhoout,b-1.0)-coeff*(un-this->gravity*this->depth);
      rtempin  = pow(rtempin, 1.0/(b-1.0));
      rtempout = pow(rtempout,1.0/(b-1.0));

      uin[inode][0] = rtempin;
      uin[inode][1] = rtempin*this->Uin[1]/this->Uin[0];
      uin[inode][2] = rtempin*this->Uin[2]/this->Uin[0];
      uin[inode][3] = rtempin*this->Uin[3]/this->Uin[0];
      uin[inode][4] = rtempin*this->Uin[4]/this->Uin[0];

      uout[inode][0] = rtempout;
      uout[inode][1] = rtempout*this->Uout[1]/this->Uout[0];
      uout[inode][2] = rtempout*this->Uout[2]/this->Uout[0];
      uout[inode][3] = rtempout*this->Uout[3]/this->Uout[0];
      uout[inode][4] = rtempout*this->Uout[4]/this->Uout[0];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsGasLiquid(IoData &iod,
                                                          DistSVec<double,3> &X)
{
  int gas = 1, liq = 0;

  // flow properties for fluid at boundary (ie liquid)
  double a = this->vf->getAlphaWater(liq);
  double b = this->vf->getBetaWater(liq);
  double c = this->vf->getCv(liq);
  double P = this->vf->getPrefWater(liq);
  double coeff = (b-1.0)/(a*b);
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = iod.bc.inlet.mach*iod.bc.inlet.mach * a*b*pow(rhoin, b-1.0);
    velout2 = iod.bc.outlet.mach*iod.bc.outlet.mach * a*b*pow(rhoout, b-1.0);
  }else if (iod.bc.inlet.velocity >= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

// Uin and Uout set up with iod.bc.let.values and properties of fluidModel1 = LiquidModel1
// computation of boundary values "on average", ie at node of coordinates (0,0,0)
  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0] * velin * sin(iod.bc.inlet.alpha);
  this->Uin[4] = this->Uin[0]*(c*iod.bc.inlet.temperature + 0.5 * velin2);

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0] * velout * sin(iod.bc.outlet.alpha);
  this->Uout[4] = this->Uout[0]*(c*iod.bc.outlet.temperature + 0.5 * velout2);


// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation
#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double rtempin, rtempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      rtempin  = pow(rhoin ,b-1.0)+coeff*(un-this->gravity*this->depth);
      rtempout = pow(rhoout,b-1.0)+coeff*(un-this->gravity*this->depth);
      rtempin  = pow(rtempin, 1.0/(b-1.0));
      rtempout = pow(rtempout,1.0/(b-1.0));

      uin[inode][0] = rtempin;
      uin[inode][1] = rtempin*this->Uin[1]/this->Uin[0];
      uin[inode][2] = rtempin*this->Uin[2]/this->Uin[0];
      uin[inode][3] = rtempin*this->Uin[3]/this->Uin[0];
      uin[inode][4] = rtempin*this->Uin[4]/this->Uin[0];

      uout[inode][0] = rtempout;
      uout[inode][1] = rtempout*this->Uout[1]/this->Uout[0];
      uout[inode][2] = rtempout*this->Uout[2]/this->Uout[0];
      uout[inode][3] = rtempout*this->Uout[3]/this->Uout[0];
      uout[inode][4] = rtempout*this->Uout[4]/this->Uout[0];
    }
  }

}

//------------------------------------------------------------------------------


template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsLiquidGas(IoData &iod,
				DistSVec<double,3> &X)
{

// this->Uin set up with inlet values and properties of fluidModel1 = GasModel1
  double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  double Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant;
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;
  double pressurein = iod.bc.inlet.pressure + rhoin*this->gravity*this->depth;
  double pressureout = iod.bc.outlet.pressure+ rhoout*this->gravity*this->depth;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = gam * (pressurein+Pstiff) * 
      iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin;
    velout2 = gam * (pressureout+ Pstiff)*
      iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout;
  }else if (iod.bc.inlet.velocity >= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0]*velin*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0]*velin*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0]*velin*sin(iod.bc.inlet.alpha);
  this->Uin[4] = (pressurein+gam*Pstiff)/(gam-1.0) + 0.5 * this->Uin[0] * velin2;

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0]*velout*cos(iod.bc.outlet.alpha)*cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0]*velout*cos(iod.bc.outlet.alpha)*sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0]*velout*sin(iod.bc.outlet.alpha);

#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double ptempin, ptempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = pressurein + rhoin *un;
      ptempout = pressureout+ rhoout*un;

      uin[inode][0] = this->Uin[0];
      uin[inode][1] = this->Uin[1];
      uin[inode][2] = this->Uin[2];
      uin[inode][3] = this->Uin[3];
      uin[inode][4] = 0.5*this->Uin[0]*velin2 + (ptempin+gam*Pstiff)/(gam-1.0); 

      uout[inode][0] = this->Uout[0];
      uout[inode][1] = this->Uout[1];
      uout[inode][2] = this->Uout[2];
      uout[inode][3] = this->Uout[3];
      uout[inode][4] = 0.5*this->Uout[0]*velout2 + (ptempout+gam*Pstiff)/(gam-1.0);

    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsJWLGas(IoData &iod,
                                DistSVec<double,3> &X)
{

// this->Uin set up with inlet values and properties of fluidModel1 = GasModel1
  double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  double Pstiff = iod.eqs.fluidModel.gasModel.pressureConstant;
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;
  double pressurein = iod.bc.inlet.pressure + rhoin*this->gravity*this->depth;
  double pressureout = iod.bc.outlet.pressure+ rhoout*this->gravity*this->depth;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = gam * (pressurein+Pstiff) * 
      iod.bc.inlet.mach*iod.bc.inlet.mach / rhoin;
    velout2 = gam * (pressureout+ Pstiff)*
      iod.bc.outlet.mach*iod.bc.outlet.mach / rhoout;
  }else if (iod.bc.inlet.velocity>= 0.0 && iod.bc.outlet.velocity >= 0.0  &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0]*velin*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0]*velin*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0]*velin*sin(iod.bc.inlet.alpha);
  this->Uin[4] = (pressurein+gam*Pstiff)/(gam-1.0) + 0.5 * this->Uin[0] * velin2;
//  this->com->printf(2, "Conservative Inlet : %e %e %e %e %e\n", pressurein, gam, Pstiff, velin2, this->Uin[4]);

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0]*velout*cos(iod.bc.outlet.alpha)*cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0]*velout*cos(iod.bc.outlet.alpha)*sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0]*velout*sin(iod.bc.outlet.alpha);
  this->Uout[4] = (pressureout+gam*Pstiff)/(gam-1.0) + 0.5 * this->Uout[0] * velout2;

#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double ptempin, ptempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      ptempin  = pressurein + rhoin *un;
      ptempout = pressureout+ rhoout*un;

      uin[inode][0] = this->Uin[0];
      uin[inode][1] = this->Uin[1];
      uin[inode][2] = this->Uin[2];
      uin[inode][3] = this->Uin[3];
      uin[inode][4] = 0.5*this->Uin[0]*velin2 + (ptempin+gam*Pstiff)/(gam-1.0); 

      uout[inode][0] = this->Uout[0];
      uout[inode][1] = this->Uout[1];
      uout[inode][2] = this->Uout[2];
      uout[inode][3] = this->Uout[3];
      uout[inode][4] = 0.5*this->Uout[0]*velout2 + (ptempout+gam*Pstiff)/(gam-1.0);

    }
  }

}
//------------------------------------------------------------------------------
template<int dim>
void DistBcDataEuler<dim>::setBoundaryConditionsJWLLiquid(IoData &iod,
                                DistSVec<double,3> &X)
{
 // flow properties for fluid at boundary
  double a = this->vf->getAlphaWater();
  double b = this->vf->getBetaWater();
  double c = this->vf->getCv();
  double P = this->vf->getPrefWater();
  double coeff = (b-1.0)/(a*b);
  double rhoin = iod.bc.inlet.density;
  double rhoout = iod.bc.outlet.density;

  double velin2 = 0.0;
  double velout2 = 0.0;
  if(iod.bc.inlet.mach >= 0.0 && iod.bc.outlet.mach >= 0.0){
    velin2 = iod.bc.inlet.mach*iod.bc.inlet.mach * a*b*pow(rhoin, b-1.0);
    velout2 = iod.bc.outlet.mach*iod.bc.outlet.mach * a*b*pow(rhoout, b-1.0);
  }else if (iod.bc.inlet.velocity>= 0.0 && iod.bc.outlet.velocity >= 0.0 &&
            iod.bc.inlet.velocity+iod.bc.outlet.velocity>0.0){
    velin2 = iod.bc.inlet.velocity*iod.bc.inlet.velocity;
    velout2 = iod.bc.outlet.velocity*iod.bc.outlet.velocity;
  }else{
    this->com->fprintf(stdout, " no proper velocity or mach number specified\n");
    exit(1);
  }
  double velin  = sqrt(velin2);
  double velout = sqrt(velout2);

// Uin and Uout set up with iod.bc.let.values and properties of fluidModel1 = LiquidModel1
// computation of boundary values "on average", ie at node of coordinates (0,0,0)
  this->Uin[0] = rhoin;
  this->Uin[1] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * cos(iod.bc.inlet.beta);
  this->Uin[2] = this->Uin[0] * velin * cos(iod.bc.inlet.alpha) * sin(iod.bc.inlet.beta);
  this->Uin[3] = this->Uin[0] * velin * sin(iod.bc.inlet.alpha);
  this->Uin[4] = this->Uin[0]*(c*iod.bc.inlet.temperature + 0.5 * velin2);

  this->Uout[0] = rhoout;
  this->Uout[1] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * cos(iod.bc.outlet.beta);
  this->Uout[2] = this->Uout[0] * velout * cos(iod.bc.outlet.alpha) * sin(iod.bc.outlet.beta);
  this->Uout[3] = this->Uout[0] * velout * sin(iod.bc.outlet.alpha);
  this->Uout[4] = this->Uout[0]*(c*iod.bc.outlet.temperature + 0.5 * velout2);


// computation for each node according to its depth
// this will be passed in DistTimeState to initialize simulation
#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*x)[3]      = X.subData(iSub);
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    double rtempin, rtempout, un;

    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      un = (x[inode][0]*this->ngravity[0]+x[inode][1]*this->ngravity[1]+x[inode][2]*this->ngravity[2]);
      rtempin  = pow(rhoin ,b-1.0)-coeff*(un-this->gravity*this->depth);
      rtempout = pow(rhoout,b-1.0)-coeff*(un-this->gravity*this->depth);
      rtempin  = pow(rtempin, 1.0/(b-1.0));
      rtempout = pow(rtempout,1.0/(b-1.0));

      uin[inode][0] = rtempin;
      uin[inode][1] = rtempin*this->Uin[1]/this->Uin[0];
      uin[inode][2] = rtempin*this->Uin[2]/this->Uin[0];
      uin[inode][3] = rtempin*this->Uin[3]/this->Uin[0];
      uin[inode][4] = rtempin*this->Uin[4]/this->Uin[0];

      uout[inode][0] = rtempout;
      uout[inode][1] = rtempout*this->Uout[1]/this->Uout[0];
      uout[inode][2] = rtempout*this->Uout[2]/this->Uout[0];
      uout[inode][3] = rtempout*this->Uout[3]/this->Uout[0];
      uout[inode][4] = rtempout*this->Uout[4]/this->Uout[0];
    }
  }

}
//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistBcDataEuler<dim>::initialize(IoData &iod, DistSVec<double,3> &X)
{

  if (iod.eqs.numPhase == 1){
    if (iod.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
        iod.eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS)
      setBoundaryConditionsGas(iod, X);
    else if(iod.eqs.fluidModel.fluid == FluidModelData::LIQUID)
      setBoundaryConditionsLiquid(iod, X);
    else if(iod.eqs.fluidModel.fluid == FluidModelData::JWL)
      setBoundaryConditionsJWL(iod, X);
    
  }else if (iod.eqs.numPhase > 1){
    map<int, FluidModelData *>::iterator it = iod.eqs.fluidModelMap.dataMap.find(1);
    if(it == iod.eqs.fluidModelMap.dataMap.end()){
      fprintf(stderr, "*** Error: no FluidModel[1] was specified\n");
      exit(1);
    }
    if ((iod.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
        iod.eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS) &&
        (it->second->fluid == FluidModelData::PERFECT_GAS ||
        it->second->fluid == FluidModelData::STIFFENED_GAS))
      setBoundaryConditionsGasGas(iod, X);
    
    else if (iod.eqs.fluidModel.fluid == FluidModelData::LIQUID &&
        it->second->fluid == FluidModelData::LIQUID)
      setBoundaryConditionsLiquidLiquid(iod, X);
    
    else if (iod.eqs.fluidModel.fluid == FluidModelData::LIQUID &&
        (it->second->fluid == FluidModelData::PERFECT_GAS ||
        it->second->fluid == FluidModelData::STIFFENED_GAS))
      setBoundaryConditionsGasLiquid(iod, X);

    else if ((iod.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
	    iod.eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS) &&
        it->second->fluid == FluidModelData::LIQUID)
      setBoundaryConditionsLiquidGas(iod, X);

    else if ((iod.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
	    iod.eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS) &&
        it->second->fluid == FluidModelData::JWL)
      setBoundaryConditionsJWLGas(iod,X);

    else if (iod.eqs.fluidModel.fluid == FluidModelData::LIQUID &&
        it->second->fluid == FluidModelData::JWL)
      setBoundaryConditionsJWLLiquid(iod,X);

    else{
      fprintf(stderr, "*** Error : no such two-phase simulation is supported by AERO-F\n");
      exit(1);
    }

  }

  if (dim == 5)
    this->finalize(X);

  if (dim == 6) {
    this->Uin[5] = this->Uin[0] * iod.bc.inlet.nutilde;
    this->Uout[5] = this->Uout[0] * iod.bc.outlet.nutilde;
    this->finalize(X);
  }

  if (dim == 7) {
    this->Uin[5] = this->Uin[0] * iod.bc.inlet.kenergy;
    this->Uin[6] = this->Uin[0] * iod.bc.inlet.eps;
    this->Uout[5] = this->Uout[0] * iod.bc.outlet.kenergy;
    this->Uout[6] = this->Uout[0] * iod.bc.outlet.eps;
    this->finalize(X);
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistBcDataEuler<dim>::initializeSA(IoData &iod, DistSVec<double,3> &X,
                                        DistSVec<double,3> &dX,
                                        double & dM, double & dA, double & dB)//indicators for mach, alpha and beta sensitivity
{

  if (iod.eqs.numPhase == 1){
    if (iod.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS)
      setDerivativeOfBoundaryConditionsGas(iod, X, dX, dM, dA, dB);
    else if(iod.eqs.fluidModel.fluid == FluidModelData::LIQUID || iod.eqs.fluidModel.fluid == FluidModelData::JWL) {
      fprintf(stderr, "*** Error: Function initializeSA (at DistBcData.C) does not support this option!\n");
      exit(1);
    }
    
  }else {
    fprintf(stderr, "*** Error: Function initializeSA (at DistBcData.C) does not support this option!\n");
    exit(1);
  }

  if (dim == 5)
    this->finalizeSA(X,dX,dM);

  if (dim == 6) {
    this->dUin[5] = 0.0;
    this->dUout[5] = 0.0;
    this->finalizeSA(X,dX,dM);
  }

  if (dim == 7) {
// Remark: it has to be derived correctly see function "optRestartBcFluxs", kenergy and eps are dependent on Mach
//    Uin[5] = Uin[0] * iod.bc.inlet.kenergy;
//    Uin[6] = Uin[0] * iod.bc.inlet.eps;
//    Uout[5] = Uout[0] * iod.bc.outlet.kenergy;
//    Uout[6] = Uout[0] * iod.bc.outlet.eps;

    this->dUin[5] = 0.0;
    this->dUin[6] = 0.0;
    this->dUout[5] = 0.0;
    this->dUout[6] = 0.0;
    this->finalizeSA(X,dX,dM);
  }
//  double dUfacenorm = this->dUface->norm();
//  if(dUfacenorm != 0) fprintf(stderr, " *********** norm of dUface is %e\n", dUfacenorm);

}


//------------------------------------------------------------------------------
// unode[i][5] contains k and unode[i][6] contains eps
// TODO this function is empty becaus there is no nutilde term here.
// Also the name EULER is misleading, since it is actually a laminar simulation
template<int dim>
void DistBcDataEuler<dim>::computeDerivativeOfNodeValue(DistSVec<double,3> &X, DistSVec<double,3> &dX)
{ }

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
DistBcDataEuler<dim>::DistBcDataEuler(IoData &iod, VarFcn *vf, Domain *dom, DistSVec<double,3> &X) :
    DistBcData<dim>(iod, vf, dom)
{

  initialize(iod, X);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistBcData<dim>::rstVar(IoData &ioData)
{

  angles[0] = ioData.bc.inlet.alpha;
  angles[1] = ioData.bc.inlet.beta;

  tref = ioData.ref.rv.time;
  dtrefdMach = ioData.ref.rv.dtimedMach;

}

//------------------------------------------------------------------------------

template<int dim>
DistBcDataSA<dim>::DistBcDataSA(IoData &iod, VarFcn *vf, Domain *dom, DistSVec<double,3> &X) : 
  DistBcDataEuler<dim>(iod, vf, dom, X)
{

  tmp = 0;
  vec2Pat = 0;

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION) {
    tmp = new DistSVec<double,2>(dom->getNodeDistInfo());
    vec2Pat = new CommPattern<double>(dom->getSubTopo(), this->com, CommPattern<double>::CopyOnSend);

// Included (MB)
    if (iod.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
        iod.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
		iod.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
		iod.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
      dtmp = new DistSVec<double,2>(dom->getNodeDistInfo());
    }
    else {
      dtmp = 0;
    }

#pragma omp parallel for
    for (int iSub = 0; iSub<this->numLocSub; ++iSub)
      this->subDomain[iSub]->setComLenNodes(2, *vec2Pat);
    vec2Pat->finalize();
  }
// Included (MB)
  else {
    dtmp = 0;
  }  

  this->Uin[5] = this->Uin[0] * iod.bc.inlet.nutilde;
  this->Uout[5] = this->Uout[0] * iod.bc.outlet.nutilde;

#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      uin[inode][5] = this->Uin[5];
      uout[inode][5] = this->Uout[5];
    }

    FaceSet& faces = this->subDomain[iSub]->getFaces();
    for (int i=0;i<faces.size(); i++) { //loop over faces 
      map<int,SurfaceData *> &surfaceMap = iod.surfaces.surfaceMap.dataMap;
      map<int,SurfaceData*>::iterator it = surfaceMap.find(faces[i].getSurfaceID());
      if(it!=surfaceMap.end()) { // surface has attribut in the input file
        map<int,BoundaryData *> &bcMap = iod.bc.bcMap.dataMap;
        map<int,BoundaryData *>::iterator it2 = bcMap.find(it->second->bcID);
        if(it2 != bcMap.end()) { // the bc data have been defined
          if(it2->second->type == BoundaryData::DIRECTSTATE || 
             it2->second->type == BoundaryData::MASSFLOW) {
            if (faces[i].getCode() ==  BC_DIRECTSTATE_INLET_MOVING ||
                faces[i].getCode() ==  BC_DIRECTSTATE_INLET_FIXED ) {
              for (int l = 0; l<faces[i].numNodes();++l) {
                int k = faces[i][l];
                uin[k][5] = it2->second->nutilde;
              }
            }
            if (faces[i].getCode() ==  BC_MASSFLOW_INLET_MOVING ||
                faces[i].getCode() ==  BC_MASSFLOW_INLET_FIXED ) {
            }
          }
        }
      }
    }

  }

// Included (MB)
  if ((iod.eqs.type == EquationsData::NAVIER_STOKES) && (iod.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY)) {
    if ((iod.bc.wall.integration == BcsWallData::WALL_FUNCTION) && (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS)) {
      dnormsa = new DistSVec<double,1>(dom->getNodeDistInfo());
      dtmpsa = new DistSVec<double,dim>(dom->getNodeDistInfo());
      vec1PatSA = new CommPattern<double>(dom->getSubTopo(), this->com, CommPattern<double>::CopyOnSend);
      vec2PatSA = new CommPattern<double>(dom->getSubTopo(), this->com, CommPattern<double>::CopyOnSend);

#pragma omp parallel for
      for (int iSub = 0; iSub<this->numLocSub; ++iSub) {
        this->subDomain[iSub]->setComLenNodes(dim, *vec1PatSA);
        this->subDomain[iSub]->setComLenNodes(1, *vec2PatSA);
      }
      vec1PatSA->finalize();
      vec2PatSA->finalize();
    }
    else {
      dnormsa = 0;
      dtmpsa = 0;
      vec1PatSA = 0;
      vec2PatSA = 0;
    }
  }
  else {
    dnormsa = 0;
    dtmpsa = 0;
    vec1PatSA = 0;
    vec2PatSA = 0;
  }

  this->finalize(X);

}

//------------------------------------------------------------------------------

template<int dim>
DistBcDataSA<dim>::~DistBcDataSA()
{

  if (tmp) delete tmp;
  if (vec2Pat) delete vec2Pat;

// Included (MB)
  if (dtmp) delete dtmp;
  delete dtmpsa;
  delete dnormsa;
  delete vec1PatSA;
  delete vec2PatSA;

}

//------------------------------------------------------------------------------
// unode[i][5] contains mutilde = rho*nutilde (=0.0 by default)

template<int dim>
void DistBcDataSA<dim>::computeNodeValue(DistSVec<double,3> &X)
{

  if (tmp) {
    int iSub;
#pragma omp parallel for
    for (iSub=0; iSub<this->numLocSub; ++iSub) {
      this->subDomain[iSub]->computeNodeBcValue(X(iSub), this->Uface(iSub), (*tmp)(iSub));
      this->subDomain[iSub]->sndData(*vec2Pat, tmp->subData(iSub));
    }

    vec2Pat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < this->numLocSub; ++iSub) {
      this->subDomain[iSub]->addRcvData(*vec2Pat, tmp->subData(iSub));
      double (*t)[2] = tmp->subData(iSub);
      double (*unode)[dim] = this->Unode.subData(iSub);
      for (int i=0; i<this->Unode.subSize(iSub); ++i) {
	    if (t[i][0] != 0.0) {
	      double w = 1.0 / t[i][0];
	      unode[i][5] = w * t[i][1];
	    }
      }
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistBcDataSA<dim>::computeDerivativeOfNodeValue(DistSVec<double,3> &X, DistSVec<double,3> &dX)
{

//Remark: Error mesage for pointers
  if (dtmp == 0) {
    fprintf(stderr, "*** Warning: Variable dtmp does not exist!\n");
    //exit(1);
  }

  if (dtmp) {
    int iSub;
#pragma omp parallel for
    for (iSub=0; iSub<this->numLocSub; ++iSub) {
      this->subDomain[iSub]->computeDerivativeOfNodeBcValue(X(iSub), dX(iSub), this->Uface(iSub), (*this->dUface)(iSub), (*dtmp)(iSub));
      this->subDomain[iSub]->sndData(*vec2Pat, dtmp->subData(iSub));
    }

    vec2Pat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < this->numLocSub; ++iSub) {
      this->subDomain[iSub]->addRcvData(*vec2Pat, dtmp->subData(iSub));
      double (*t)[2] = tmp->subData(iSub);
      double (*dt)[2] = dtmp->subData(iSub);
      double (*dunode)[dim] = this->dUnode->subData(iSub);
      for (int i=0; i<this->dUnode->subSize(iSub); ++i) {
  	    if (t[i][0] != 0.0) {
	      double w = 1.0 / t[i][0];
	      double dw = -1.0 / ( t[i][0] * t[i][0] ) * dt[i][0];
	      dunode[i][5] = dw * t[i][1] + w * dt[i][1];
	    }
      }
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistBcDataSA<dim>::computeNodeWallValues(DistSVec<double,3> &X)
{

  if (dtmpsa) {
    int iSub;
#pragma omp parallel for
    for (iSub=0; iSub<this->numLocSub; ++iSub) {
      this->subDomain[iSub]->computeNodeBCsWallValues(X(iSub), (*dnormsa)(iSub), (*this->dUfaceSA)(iSub), (*dtmpsa)(iSub));
      this->subDomain[iSub]->sndData(*vec1PatSA, dtmpsa->subData(iSub));
      this->subDomain[iSub]->sndData(*vec2PatSA, dnormsa->subData(iSub));
    }

    vec1PatSA->exchange();
    vec2PatSA->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < this->numLocSub; ++iSub) {
      this->subDomain[iSub]->addRcvData(*vec1PatSA, dtmpsa->subData(iSub));
      this->subDomain[iSub]->addRcvData(*vec2PatSA, dnormsa->subData(iSub));
      double (*dnsa)[1] = dnormsa->subData(iSub);
      double (*dtsa)[dim] = dtmpsa->subData(iSub);
      double (*dunodesa)[dim] = this->dUnodeSA->subData(iSub);
      for (int i=0; i<this->dUnodeSA->subSize(iSub); ++i) {
	if (dnsa[i][0] != 0.0) {
	  for (int k=0; k<dim; ++k)
	    dunodesa[i][k] = dtsa[i][k] / dnsa[i][0];
	}
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
DistBcDataKE<dim>::DistBcDataKE(IoData &iod, VarFcn *vf, Domain *dom, DistSVec<double,3> &X) : 
  DistBcDataEuler<dim>(iod, vf, dom, X)
{

  tmp = new DistSVec<double,3>(dom->getNodeDistInfo());
  vec3Pat = dom->getVec3DPat();

// Included (MB)
  if (iod.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
      iod.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
	  iod.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
	  iod.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
    dtmp = new DistSVec<double,3>(dom->getNodeDistInfo());
  }
  else {
    dtmp = 0;
  }

  this->Uin[5] = this->Uin[0] * iod.bc.inlet.kenergy;
  this->Uin[6] = this->Uin[0] * iod.bc.inlet.eps;
  this->Uout[5] = this->Uout[0] * iod.bc.outlet.kenergy;
  this->Uout[6] = this->Uout[0] * iod.bc.outlet.eps;

#pragma omp parallel for
  for(int iSub = 0; iSub<this->numLocSub; ++iSub) {
    double (*uin)[dim]  = this->Ufarin.subData(iSub);
    double (*uout)[dim] = this->Ufarout.subData(iSub);
    for(int inode = 0; inode<this->Unode.subSize(iSub); inode++){
      uin[inode][5] = this->Uin[5];
      uin[inode][6] = this->Uin[6];
      uout[inode][5] = this->Uout[5];
      uout[inode][6] = this->Uout[6];
    }
  }

  this->finalize(X);

}

//------------------------------------------------------------------------------

template<int dim>
DistBcDataKE<dim>::~DistBcDataKE()
{

  if (tmp) delete tmp;

// Included (MB)
  if (dtmp) delete dtmp;

}
//------------------------------------------------------------------------------
// unode[i][5] contains k and unode[i][6] contains eps

template<int dim>
void DistBcDataKE<dim>::computeNodeValue(DistSVec<double,3> &X)
{

  int iSub;

#pragma omp parallel for
  for (iSub=0; iSub<this->numLocSub; ++iSub) {
    this->subDomain[iSub]->computeNodeBcValue(X(iSub), this->Uface(iSub), (*tmp)(iSub));
    this->subDomain[iSub]->sndData(*vec3Pat, tmp->subData(iSub));
  }

  vec3Pat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->addRcvData(*vec3Pat, tmp->subData(iSub));
    double (*t)[3] = tmp->subData(iSub);
    double (*unode)[dim] = this->Unode.subData(iSub);
    for (int i=0; i<this->Unode.subSize(iSub); ++i) {
      if (t[i][0] != 0.0) {
	double w = 1.0 / t[i][0];
	unode[i][5] = w * t[i][1];
	unode[i][6] = w * t[i][2];
      }
    }
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistBcDataKE<dim>::computeDerivativeOfNodeValue(DistSVec<double,3> &X, DistSVec<double,3> &dX)
{

//Remark: Error mesage for pointers
  if (dtmp == 0) {
    fprintf(stderr, "*** Warning: Variable dtmp does not exist!\n");
    //fprintf(stderr, "*** Error: Variable dtmp does not exist!\n");
    //exit(1);
  }

  int iSub;

#pragma omp parallel for
  for (iSub=0; iSub<this->numLocSub; ++iSub) {
    this->subDomain[iSub]->computeDerivativeOfNodeBcValue(X(iSub), dX(iSub), this->Uface(iSub), (*this->dUface)(iSub), (*dtmp)(iSub));
    this->subDomain[iSub]->sndData(*vec3Pat, dtmp->subData(iSub));
  }

  vec3Pat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->addRcvData(*vec3Pat, dtmp->subData(iSub));
    double (*t)[3] = tmp->subData(iSub);
    double (*dt)[3] = dtmp->subData(iSub);
    double (*dunode)[dim] = this->dUnode->subData(iSub);
    for (int i=0; i<this->dUnode->subSize(iSub); ++i) {
      if (t[i][0] != 0.0) {
	double w = 1.0 / t[i][0];
	double dw = -1.0 / ( t[i][0] * t[i][0] ) * dt[i][0];
	dunode[i][5] = dw * t[i][1] + w * dt[i][1];
	dunode[i][6] = dw * t[i][2] + w * dt[i][2];
      }
    }
  }

}

//------------------------------------------------------------------------------
