#ifndef _IMPLICIT_ROM_POSTPRO_TS_DESC_H_
#define _IMPLICIT_ROM_POSTPRO_TS_DESC_H_

#include <ImplicitRomTsDesc.h>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitRomPostproTsDesc : public ImplicitRomTsDesc<dim> {

protected:

	FILE *reducedCoordsFile;	// file of reduced coordinates

  void saveNewtonSystemVectors(const int totalTimeSteps) {this->saveNewtonSystemVectorsAction(totalTimeSteps);}
	void solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim>&, const int& totalTimeSteps = 0);
	virtual void computeFullResidual(int it, DistSVec<double, dim> &Q, bool applyWeighting,  DistSVec<double, dim> *R, bool includeHomotopy);
	virtual void computeAJ(int it, DistSVec<double, dim> &Q, bool applyWeighting,  DistSVec<double, dim> *R);
  DistSVec<double, dim> Uinitial;	// solution increment at EACH NEWTON ITERATION in full coordinates
	virtual void postProStep(DistSVec<double,dim> &, int);	// by default, do not do post processing
  
  double dt;

  double meritFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double) {return 0.0;} // not applicable 

public:

  void checkLocalRomStatus(DistSVec<double, dim> &, const int);
  double computeTimeStep(int, double*, DistSVec<double,dim> &, double);
  bool monitorConvergence(int, DistSVec<double,dim> &);
 
  ImplicitRomPostproTsDesc(IoData &, GeoSource &, Domain *);
  void formInterpolatedInitialCondition(DistSVec<double,dim> *U, IoData &iod)  {
    if (this->ioData->romDatabase.files.surfacePrefix[0]==0 && this->ioData->romDatabase.files.surfaceRefStateName[0]==0) {
      dynamic_cast<TsDesc<dim>*>(this)->TsDesc<dim>::formInterpolatedInitialCondition(U, iod);
    } else {
      this->rom->formInterpolatedInitialCondition(U, iod);
    }
  }


};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitRomPostproTsDesc.C>
#endif

#endif
