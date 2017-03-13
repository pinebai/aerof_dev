#include <IoData.h>

#include <parser/Assigner.h>

#include<limits>

//------------------------------------------------------------------------------

template<class GenericKrylov>
NewtonData<GenericKrylov>::NewtonData()
{

  failsafe = NO;
  maxIts = 1;
  eps = 1.e-2;
  JacSkip = 1;
  epsAbsRes = std::numeric_limits<double>::epsilon();
  epsAbsInc = std::numeric_limits<double>::epsilon();
  output = "";
}

//------------------------------------------------------------------------------

template<class GenericKrylov>
void NewtonData<GenericKrylov>::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassToken<NewtonData>
    (ca, "FailSafe", this, reinterpret_cast<int NewtonData::*>(&NewtonData::failsafe), 3,
     "Off", 0, "On", 1, "AlwaysOn", 2);
  new ClassInt<NewtonData>(ca, "MaxIts", this, &NewtonData::maxIts);
  new ClassDouble<NewtonData>(ca, "Eps", this, &NewtonData::eps);
  new ClassInt<NewtonData>(ca, "JacobianFrequency", this, &NewtonData::JacSkip);
  new ClassDouble<NewtonData>(ca, "EpsAbsRes", this, &NewtonData::epsAbsRes);
  new ClassDouble<NewtonData>(ca, "EpsAbsInc", this, &NewtonData::epsAbsInc);
  new ClassStr<NewtonData>(ca, "Output", this, &NewtonData::output);

  ksp.setup("LinearSolver", ca);
  lineSearch.setup("LineSearch",ca);
}


//------------------------------------------------------------------------------

template<class DataType>
void ObjectMap<DataType>::setup(const char *name, ClassAssigner *p)  {

  SysMapObj<DataType> *smo = new SysMapObj<DataType>(name, &dataMap);

  if (p) p->addSmb(name, smo);
  else addSysSymbol(name, smo);

}

