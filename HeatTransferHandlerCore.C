#include <HeatTransferHandler.h>

#include <IoData.h>
#include <MatchNode.h>
#include <StructExc.h>
#include <Domain.h>

//------------------------------------------------------------------------------

HeatTransferHandler::HeatTransferHandler(IoData& iod, MatchNodeSet** matchNodes, Domain* dom) : 
  domain(dom), P(dom->getNodeDistInfo())
{

  steady = !iod.problem.type[ProblemData::UNSTEADY];
  it0 = iod.restart.iteration;
  com = dom->getCommunicator();

  strExc = new StructExc(iod, matchNodes, 1, domain->getHeatCommunicator(), domain->getCommunicator(),
           domain->getNumLocSub());
  strExc->negotiate();
  strExc->getInfo();

}

//------------------------------------------------------------------------------

void HeatTransferHandler::setup(int* rstrt, double* maxTime)
{
  *rstrt = strExc->getRestartFrequency();

  *maxTime = strExc->getMaxTime();

}

//------------------------------------------------------------------------------

double HeatTransferHandler::update(bool* lastIt, int it, DistVec<double>& T)
{

  int algNum = strExc->getAlgorithmNumber();
  double dt = strExc->getTimeStep();
  if (steady)
    dt = 0.0;

  if (algNum == 6 && it == 0) 
    dt *= 0.5;

  if (it > it0) {
    strExc->sendHeatPower(P);

    if (steady) {
      strExc->negotiateStopping(lastIt);
      if (*lastIt) 
	return 0.0;
    }
  }
  strExc->getTemperature(T);
  if (*lastIt) 
    return 0.0;

  if (algNum == 1) 
    *lastIt = true;

  return dt;

}

//------------------------------------------------------------------------------
double HeatTransferHandler::updateStep1(bool* lastIt, int it, DistVec<double>& T)
{

  int algNum = strExc->getAlgorithmNumber();
  double dt = strExc->getTimeStep();
  if (steady)
    dt = 0.0;

  if (algNum == 6 && it == 0)
    dt *= 0.5;
  if (it > it0) {
    strExc->sendHeatPower(P);

   // Vamshi -> This is needed so as to synchronize the fluid with thermal in Quasistatic Thermal simulation
   /*   if (steady) {
         strExc->negotiateStopping(lastIt);
         if (*lastIt)
          return 0.0;
        }
   */

  }

  return dt;

}

//------------------------------------------------------------------------------
double HeatTransferHandler::updateStep2(bool* lastIt, int it, DistVec<double>& T)
{
  int algNum = strExc->getAlgorithmNumber();
  double dt = strExc->getTimeStep();
  if (steady)
    dt = 0.0;

  if (algNum == 6 && it == 0)
    dt *= 0.5;

  if (it > it0) {
    if (steady) {
      strExc->negotiateStopping(lastIt);
      if (*lastIt)
        return 0.0;
    }
  }

  strExc->getTemperature(T);

  if (*lastIt)
    return 0.0;

  if (algNum == 1)
    *lastIt = true;

  return dt;

}

//------------------------------------------------------------------------------
