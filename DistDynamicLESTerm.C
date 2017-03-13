#include <DistDynamicLESTerm.h>

#include <DistVector.h>
#include <Domain.h>

#include <cstdio>

//------------------------------------------------------------------------

template<int dim>
DistDynamicLESTerm<dim>::DistDynamicLESTerm(VarFcn *vf, IoData &iod, Domain *dom) : domain(dom), varFcn(vf)    

{

  numLocSub  = domain->getNumLocSub();

  gam = iod.eqs.fluidModel.gasModel.specificHeatRatio; // ratio of cp/cv
  Rideal = iod.eqs.fluidModel.gasModel.idealGasConstant; // ideal gas cons

  VCap = new DistSVec<double,dim>(domain->getNodeDistInfo()); 
  Mom_Test = new DistSVec<double,16>(domain->getNodeDistInfo());
  Sij_Test = new DistSVec<double,6>(domain->getNodeDistInfo());
  modS_Test = new DistVec<double>(domain->getNodeDistInfo());
  Eng_Test = new DistSVec<double,8>(domain->getNodeDistInfo());
  Ni = new DistVec<int>(domain->getNodeDistInfo());
  Cs = new DistSVec<double,2>(domain->getNodeDistInfo());

  dlest = new DynamicLESTerm(iod, varFcn);
  
  *Ni = 0;
  domain->computeTetsConnectedToNode(*Ni); //computing no of tets conneted to node i

}

//------------------------------------------------------------------------

template<int dim>
DistDynamicLESTerm<dim>::~DistDynamicLESTerm()

{

  if (VCap) delete VCap;
  if (Mom_Test) delete Mom_Test;
  if (Sij_Test) delete Sij_Test;
  if (modS_Test) delete modS_Test;
  if (Eng_Test) delete Eng_Test;
  if (Cs) delete Cs;
  if (Ni) delete Ni;
  if (dlest) delete dlest;

}

//------------------------------------------------------------------------

template<int dim>
void DistDynamicLESTerm<dim>::compute(DistVec<double> &ctrlVol, 
                                      DistBcData<dim> &bcData,
				      DistSVec<double,3> &X, 
				      DistSVec<double,dim> &V, 
                                      DistSVec<double,dim> &R,
                                      DistVec<GhostPoint<dim>*> *ghostPoints, 
                                      DistLevelSetStructure *LSS,
	                                   bool externalSI)
{

  *VCap = 0.0;       // contains test filtered values of the primitive variables //
  *Mom_Test = 0.0;   // contains test filtered values for computing cs from momentum equation //
  *Sij_Test = 0.0;   // contains test filtered value of Sij //
  *modS_Test = 0.0;  // contains test filtered modulus of Sij //
  *Eng_Test = 0.0;   // contains test filtered values for computing pt from energy equation //
  *Cs = 0.0;         // contains Cs and Pt values //

  domain->computeTestFilterValues(ctrlVol, *VCap, *Mom_Test, *Sij_Test, *modS_Test, 
                                  *Eng_Test, *Cs, *Ni, bcData, X, V, gam, Rideal, 
                                  ghostPoints, LSS, externalSI);

  domain->computeDynamicLESTerm(dlest, *Cs, X, V, R, ghostPoints, LSS, externalSI);

}

//------------------------------------------------------------------------

template<int dim>
void DistDynamicLESTerm<dim>::computeMutOMu(DistVec<double> &ctrlVol, 
                                            DistBcData<dim> &bcData,
                                            DistSVec<double,3> &X, 
                                            DistSVec<double,dim> &V, 
                                            DistVec<double> &mutOmu)
{

  *VCap = 0.0;       // contains test filtered values of the primitive variables //
  *Mom_Test = 0.0;   // contains test filtered values for computing cs from momentum equation //
  *Sij_Test = 0.0;   // contains test filtered value of Sij //
  *modS_Test = 0.0;  // contains test filtered modulus of Sij //
  *Eng_Test = 0.0;   // contains test filtered values for computing pt from energy equation //

  *Cs = 0.0; 

  domain->computeTestFilterValues(ctrlVol, *VCap, *Mom_Test, *Sij_Test, *modS_Test,
                                  *Eng_Test, *Cs, *Ni, bcData, X, V, gam, Rideal);
  domain->computeMutOMuDynamicLES(dlest, ctrlVol, *Cs, X, V, mutOmu);

}

//------------------------------------------------------------------------

template<int dim>
void DistDynamicLESTerm<dim>::computeCsValue(DistVec<double> &ctrlVol,
                                             DistBcData<dim> &bcData,
                                             DistSVec<double,3> &X,
                                             DistSVec<double,dim> &V,
                                             DistVec<double> &CsVal)
{

  *VCap = 0.0;       // contains test filtered values of the primitive variables //
  *Mom_Test = 0.0;   // contains test filtered values for computing cs from momentum equation //
  *Sij_Test = 0.0;   // contains test filtered value of Sij //
  *modS_Test = 0.0;  // contains test filtered modulus of Sij //
  *Eng_Test = 0.0;   // contains test filtered values for computing pt from energy equation //

  *Cs = 0.0;

  domain->computeTestFilterValues(ctrlVol, *VCap, *Mom_Test, *Sij_Test, *modS_Test,
                                  *Eng_Test, *Cs, *Ni, bcData, X, V, gam, Rideal);

  domain->outputCsDynamicLES(dlest, ctrlVol, *Cs, X, CsVal); 

}

//------------------------------------------------------------------------
