#include <cmath>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif

#include <KspConvCriterion.h>

//------------------------------------------------------------------------------

KspConvCriterion::KspConvCriterion(KspData &data)
{

  typeEpsFormula = data.epsFormula;

  eps0 = data.eps;

  nlits = 0;

}

//------------------------------------------------------------------------------

double KspConvCriterion::compute(int nlit, int nlmaxits, double fnorm)
{

  double eps;

  if (nlit == 0 && nlmaxits > 1) nlits = 0;

  if (typeEpsFormula == KspData::CONSTANT)
    eps = eps0;
  else if (typeEpsFormula == KspData::EISENSTADT)
    eps = computeEpsilonEisenstadt(fnorm);

  return eps;

}

//------------------------------------------------------------------------------

double KspConvCriterion::computeEpsilonEisenstadt(double fnorm)
{

  double eps, eps1;

  double gamma  = 0.1;
  double alpha  = 2.0;
  double delta  = 0.1;

  if (nlits == 0)
    eps = eps0;
  else {

    eps = gamma * pow(fnorm/fnormPrev, alpha);

    eps1 = gamma * pow(epsPrev, alpha);

    if (eps1 < delta) eps = min(eps0, eps);
    else eps = min(eps0, max(eps, eps1));

  }

  fnormPrev = fnorm;
  epsPrev = eps;
  ++nlits;

  return eps;

}

//------------------------------------------------------------------------------
