#include <RefVal.h>
#include <IoData.h>
#include <cmath>

//------------------------------------------------------------------------------
//CHANGES_FOR_WATER
//	new way to define reference values for viscosity, velocity and 
//		temperature 
//------------------------------------------------------------------------------

RefVal::RefVal()
{

  mode = NON_DIMENSIONAL;
  length = 1.0;
  density = 1.0;
  velocity = 1.0;
  pressure = 1.0;
  temperature = 1.0;
  viscosity_mu = 1.0;
  nutilde = 1.0;
  kenergy = 1.0;
  epsilon = 1.0;
  time = 1.0;
  force = 1.0;
  energy = 1.0;
  power = 1.0;
  entropy = 1.0;

  tlength = 1.0;
  tvelocity = 1.0;
  tforce = 1.0;
  tpower = 1.0;

// Included (MB)
  dvelocitydMach = 1.0;
  dtimedMach = 1.0;

}

//------------------------------------------------------------------------------
/*RefVal::RefVal(IoData &ioData)
{


  if (ioData.problem.mode == ProblemData::NON_DIMENSIONAL) {
    mode = NON_DIMENSIONAL;
    length = 1.0;
    density = 1.0;
    velocity = 1.0;
    pressure = 1.0;
    temperature = 1.0;
    viscosity_mu = 1.0;
    nutilde = 1.0;
    kenergy = 1.0;
    epsilon = 1.0;
    time = 1.0;
    force = 1.0;
    energy = 1.0;
    power = 1.0;
    entropy = 1.0;

    tlength = 1.0;
    tvelocity = 1.0;
    tforce = 1.0;
    tpower = 1.0;

// Included (MB)
    dvelocitydMach = 1.0;
    dtimedMach = 1.0;

  }
  else if (ioData.problem.mode  == ProblemData::DIMENSIONAL) {
    mode = DIMENSIONAL;
    if(ioData.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
       ioData.eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS){
      double gam = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
      double R = ioData.eqs.fluidModel.gasModel.idealGasConstant;
      double Pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;
      mach = ioData.ref.mach;
      density = ioData.ref.density;
      velocity = mach * sqrt(gam * (ioData.ref.pressure + Pstiff) / density);
      pressure = density * velocity*velocity;
      temperature = gam*(gam - 1.0) * mach*mach * (ioData.ref.pressure + Pstiff)/(density*R);
      entropy = pow(density,1.0-gam)*velocity*velocity;
// Included (MB)
      dvelocitydMach = sqrt(gam * ioData.ref.pressure / density);
      dtimedMach = - length / (velocity * velocity) * dvelocitydMach;
    }
    else if(ioData.eqs.fluidModel.fluid == FluidModelData::JWL){
      double omegajwl  = ioData.eqs.fluidModel.jwlModel.omega;
      double A1jwl     = ioData.eqs.fluidModel.jwlModel.A1;
      double A2jwl     = ioData.eqs.fluidModel.jwlModel.A2;
      double R1jwl     = ioData.eqs.fluidModel.jwlModel.R1;
      double R2jwl     = ioData.eqs.fluidModel.jwlModel.R2;
      double rhorefjwl = ioData.eqs.fluidModel.jwlModel.rhoref;
      mach = ioData.ref.mach;
      density = ioData.ref.density;
      double frhoref   = A1jwl*(1-omegajwl*density/(R1jwl*rhorefjwl))*exp(-R1jwl*rhorefjwl/density) +
                         A2jwl*(1-omegajwl*density/(R2jwl*rhorefjwl))*exp(-R2jwl*rhorefjwl/density);
      double frhorefp  = A1jwl*(-omegajwl/(R1jwl*rhorefjwl) + (1.0-omegajwl*density/(R1jwl*rhorefjwl))*R1jwl*rhorefjwl/(density*density))
                          *exp(-R1jwl*rhorefjwl/density)
                       + A2jwl*(-omegajwl/(R2jwl*rhorefjwl) + (1.0-omegajwl*density/(R2jwl*rhorefjwl))*R2jwl*rhorefjwl/(density*density))
                          *exp(-R2jwl*rhorefjwl/density);
      velocity = mach * sqrt(((omegajwl+1.0)*ioData.ref.pressure - frhoref + density*frhorefp)/density);
      pressure = density * velocity*velocity;
      temperature = omegajwl*(omegajwl + 1.0) * mach*mach * ((omegajwl+1.0)*ioData.ref.pressure - frhoref + density*frhorefp);
      entropy = pow(density,-omegajwl)*velocity*velocity;
// Included (MB)
      dvelocitydMach = 1.0;
      dtimedMach = 1.0;
    }
    else if(ioData.eqs.fluidModel.fluid == FluidModelData::LIQUID){
      double C = ioData.eqs.fluidModel.liquidModel.specificHeat;
      double awater = ioData.eqs.fluidModel.liquidModel.alpha;
      double bwater = ioData.eqs.fluidModel.liquidModel.beta;
      mach = ioData.ref.mach;
      density = ioData.ref.density;
      velocity = mach * sqrt(awater*bwater*pow(density, bwater - 1.0));
      pressure = density * velocity*velocity;
      temperature = velocity*velocity/C;
      entropy = pow(density,1.0-bwater)*velocity*velocity;
    }

    length = ioData.ref.length;
    //surface = ioData.ref.surface;

    force = density * velocity*velocity * length*length;
    moment = force * length;
    energy = force * length;

    //cf2force = force * surface / (2.0 * length*length);
    //cm2moment = moment * surface / (2.0 * length*length);
    cf2force = force  / (2.0 * length*length);
    cm2moment = moment / (2.0 * length*length);

    time = length / velocity;

  }

}
*/
//------------------------------------------------------------------------------

// Included (MB)
void RefVal::rstVar(IoData &ioData)
{

  if (ioData.problem.mode == ProblemData::NON_DIMENSIONAL) {
    mode = NON_DIMENSIONAL;
    length = 1.0;
    density = 1.0;
    velocity = 1.0;
    pressure = 1.0;
    temperature = 1.0;
    viscosity_mu = 1.0;
    nutilde = 1.0;
    kenergy = 1.0;
    epsilon = 1.0;
    time = 1.0;
    force = 1.0;
    energy = 1.0;
    power = 1.0;
    entropy = 1.0;

    tlength = 1.0;
    tvelocity = 1.0;
    tforce = 1.0;
    tpower = 1.0;

    dvelocitydMach = 1.0;
    dtimedMach = 1.0;

  }
  else if (ioData.problem.mode  == ProblemData::DIMENSIONAL) {
    mode = DIMENSIONAL;
    if(ioData.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
       ioData.eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS){
      double gam = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
      double R = ioData.eqs.fluidModel.gasModel.idealGasConstant;
      double Pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;
      mach = ioData.ref.mach;
      density = ioData.ref.density;
      velocity = mach * sqrt(gam * (ioData.ref.pressure + Pstiff) / density);
      pressure = density * velocity*velocity;
      temperature = gam*(gam - 1.0) * mach*mach * (ioData.ref.pressure + Pstiff)/(density*R);
      dvelocitydMach = sqrt(gam * ioData.ref.pressure / density);
      dtimedMach = - length / (velocity * velocity) * dvelocitydMach;
    }

    length = ioData.ref.length;
    //surface = ioData.ref.surface;

    force = density * velocity*velocity * length*length;
    moment = force * length;
    energy = force * length;

    //cf2force = force * surface / (2.0 * length*length);
    //cm2moment = moment * surface / (2.0 * length*length);
    cf2force = force  / (2.0 * length*length);
    cm2moment = moment / (2.0 * length*length);

    time = length / velocity;

  }

}

//------------------------------------------------------------------------------
