#include <cstdlib>
#include <cmath>

#include <CorotSolver.h>

#include <MatchNode.h>
#include <Domain.h>
#include <Vector3D.h>
#include <DenseMatrixOps.h>

#include <BCApplier.h>

//#define HB_COROTSOLVER_DEBUG

//------------------------------------------------------------------------------

CorotSolver::CorotSolver(DefoMeshMotionData &data, MatchNodeSet **matchNodes, Domain *dom)
  : domain(dom), X0(dom->getNodeDistInfo()), Xtilde(dom->getNodeDistInfo())
{

  numLocSub = domain->getNumLocSub();

  com = domain->getCommunicator();

  domain->getReferenceMeshPosition(X0);

  domain->getNdAeroLists(nInterfNd, interfNd, nInfNd, infNd, nInternalNd, internalNd, matchNodes); //HB

  computeCG(X0, cg0);
#ifdef HB_COROTSOLVER_DEBUG
  com->fprintf(stderr," -> In CorotSolver::CorotSolver, cg0 = %.2e  %.2e  %.2e\n",cg0[0],cg0[1],cg0[2]);
#endif

  cgN[0] = cg0[0];
  cgN[1] = cg0[1];
  cgN[2] = cg0[2];
  
  typedef double (*ddummy)[3];
  gapVec = new ddummy[numLocSub];

  if (matchNodes) 
    locAllocGap = false;
  else
    locAllocGap = true;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    if (matchNodes) 
      gapVec[iSub] = matchNodes[iSub]->getGap(nInterfNd[iSub], interfNd[iSub]);
    else {
      gapVec[iSub] = new double[ nInterfNd[iSub] ][3];
      for (int i=0; i<nInterfNd[iSub]; ++i)
	for (int j=0; j<3; ++j)
	  gapVec[iSub][i][j] = 0.0;
    }
  }

  double zeroRot[3] = {0.0, 0.0, 0.0};
  computeRotMat(zeroRot, R);

  //HB: look if a symmetry plane was specified in the input file
  double nx = data.symmetry.nx;
  double ny = data.symmetry.ny;
  double nz = data.symmetry.nz;
  double nrm= sqrt(nx*nx+ny*ny+nz*nz);
  SymAxis = CorotSolver::NONE;
  if(nrm!=0.0){
    nx /= nrm; ny /= nrm; nz /= nrm; 
    if((fabs(nx)==1.0) & (ny==0.0) & (nz==0.0))
      SymAxis = CorotSolver::AXIS_X;
    else if((nx==0.0) & (fabs(ny)==1.0) & (nz==0.0))
      SymAxis = CorotSolver::AXIS_Y;
    else if((nx==0.0) & (ny==0.0) & (fabs(nz)==1.0))
      SymAxis = CorotSolver::AXIS_Z;
    else {
      com->fprintf(stderr," *** ERROR: corotational solver only supports a canonical plane as a symmetry plane.\n");
      exit(-1);
    }
  } 
  switch(SymAxis) {
    case(CorotSolver::NONE):
      com->fprintf(stderr," ... No symmetry plane is used in the corotational solver.\n");
      com->fprintf(stderr,"     -> the 3 rotations axis (X,Y,Z) & 3 translation axis (X,Y,Z) are used.\n");
      break;
    case(CorotSolver::AXIS_X):
      com->fprintf(stderr," ... Symmetry plane of normal X is used in the corotational solver.\n");
      com->fprintf(stderr,"     -> only rotation around axis X & translations in the Y-Z plane are allowed.\n");
      break;
    case(CorotSolver::AXIS_Y):
      com->fprintf(stderr," ... Symmetry plane of normal Y is used in the corotational solver.\n");
      com->fprintf(stderr,"     -> only rotation around axis Y & translations in the X-Z plane are allowed.\n");
      break;
    case(CorotSolver::AXIS_Z):
      com->fprintf(stderr," ... Symmetry plane of normal Z is used in the corotational solver.\n");
      com->fprintf(stderr,"     -> only rotation around axis Z & translations in the X-Y plane are allowed.\n");
      break;
  }
}

//------------------------------------------------------------------------------

CorotSolver::~CorotSolver()
{

  if (gapVec) {
    if (locAllocGap) {
#pragma omp parallel for
      for (int iSub = 0; iSub < numLocSub; ++iSub)
	delete [] gapVec[iSub];
    }

    delete [] gapVec;
  }

}

//------------------------------------------------------------------------------

void CorotSolver::computeRotGradAndJac(DistSVec<double,3> &X, double RR[3][3], 
				       double cg1[3], double grad[3], double jac[3][3])
{

  int iSub, i, j;

  const DistInfo &distInfo = X.info();

  double (*allgrad)[3] = reinterpret_cast<double (*)[3]>
                         (alloca(sizeof(double) * 3*distInfo.numGlobSub));
  double (*alljac)[3][3] = reinterpret_cast<double (*)[3][3]>
                           (alloca(sizeof(double) * 9*distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) {
    for (i=0; i<3; ++i) {
      allgrad[iSub][i] = 0.0;
      for (j=0; j<3; ++j)
        alljac[iSub][i][j] = 0.0;
    }
  }

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    int globSub = distInfo.locSubToGlobSub[iSub];
    computeSubGradAndJac(X(iSub), RR, cg1, allgrad[globSub], alljac[globSub], iSub);
  }

  com->globalSum(3*distInfo.numGlobSub, reinterpret_cast<double *>(allgrad));
  com->globalSum(9*distInfo.numGlobSub, reinterpret_cast<double *>(alljac));

  for (i=0; i<3; ++i) {
    grad[i] = 0.0;
    for (j=0; j<3; ++j)
      jac[i][j] = 0.0;
  }

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) {
    for (i=0; i<3; ++i) {
      grad[i] += allgrad[iSub][i];
      for (j=0; j<3; ++j)
        jac[i][j] += alljac[iSub][i][j];
    }
  }

}

//------------------------------------------------------------------------------

inline
void invRotLocVec(double mat[3][3], double v[3]) 
{

  double c[3];
  for (int j = 0; j < 3; j++)
    c[j] = mat[0][j]*v[0] + mat[1][j]*v[1] + mat[2][j]*v[2];

  v[0] = c[0];
  v[1] = c[1];
  v[2] = c[2];

}

//------------------------------------------------------------------------------

double CorotSolver::computeSubGradAndJac(SVec<double,3> &subNd, double RR[3][3], 
					 double cg1[3], double grad[3], double jac[3][3],
					 int iSub)
{
  int i, j;
  double d2 = 0.0;

  for (i = 0; i < 3; i++)  {
    grad[i] = 0.0;
    for (j = 0; j < 3; j++)
      jac[i][j] = 0.0;
  }

  Vec3D rotGrad[3];
  SVec<double,3> &x0 = X0(iSub);

  for (i = 0; i < nInterfNd[iSub]; i++) {
    double rd[3];
    // rotate the local vectors using R(n-1)
    rd[0] = x0[interfNd[iSub][i]][0] - cg0[0] - gapVec[iSub][i][0];
    rd[1] = x0[interfNd[iSub][i]][1] - cg0[1] - gapVec[iSub][i][1];
    rd[2] = x0[interfNd[iSub][i]][2] - cg0[2] - gapVec[iSub][i][2];
    rotLocVec(RR, rd);

    Vec3D eVec;
    //compute freq. used values
    eVec[0] = subNd[interfNd[iSub][i]][0] - cg1[0] - rd[0] 
	    - gapVec[iSub][i][0];
    eVec[1] = subNd[interfNd[iSub][i]][1] - cg1[1] - rd[1] 
	    - gapVec[iSub][i][1];
    eVec[2] = subNd[interfNd[iSub][i]][2] - cg1[2] - rd[2]
	    - gapVec[iSub][i][2];

    d2 += eVec*eVec;

    rotGrad[0][0] = -rd[1];
    rotGrad[0][1] = rd[0];
    rotGrad[0][2] = 0;
  
    rotGrad[1][0] = rd[2];
    rotGrad[1][1] = 0;
    rotGrad[1][2] = -rd[0];

    rotGrad[2][0] = 0;
    rotGrad[2][1] = -rd[2];
    rotGrad[2][2] = rd[1];

    grad[0] += -2*(eVec * rotGrad[0]);
    grad[1] += -2*(eVec * rotGrad[1]);
    grad[2] += -2*(eVec * rotGrad[2]);

    jac[0][0] +=  2 * (rd[0]*eVec[0]
                    +  rd[1]*eVec[1]
                    +  rotGrad[0][0]*rotGrad[0][0]
                    +  rotGrad[0][1]*rotGrad[0][1]);

    jac[0][1] += -2 * (rd[2]*eVec[1]
                    - rotGrad[0][0]*rotGrad[1][0]);

    jac[0][2] += -2 * (rd[2]*eVec[0]
                    -  rotGrad[0][1]*rotGrad[2][1]);

    jac[1][1] += 2 * (rd[0]*eVec[0]
                    +  rd[2]*eVec[2]
                    +  rotGrad[1][0]*rotGrad[1][0]
                    +  rotGrad[1][2]*rotGrad[1][2]);

    jac[1][2] += -2 * (rd[1]*eVec[0]
                    -   rotGrad[2][2]*rotGrad[1][2]);

    jac[2][2] += 2 * (rd[2]*eVec[2]
                    +  rd[1]*eVec[1]
                    +  rotGrad[2][1]*rotGrad[2][1]
                    +  rotGrad[2][2]*rotGrad[2][2]);
  }

  // Fill in symmetric terms
  jac[1][0] = jac[0][1];
  jac[2][0] = jac[0][2];
  jac[2][1] = jac[1][2];

  return d2;
}

//------------------------------------------------------------------------------

void CorotSolver::rotLocVec(double mat[3][3], double v[3]) 
{
  double c[3];
  for (int j = 0; j < 3; j++)
    c[j] = mat[j][0]*v[0] +
           mat[j][1]*v[1] +
           mat[j][2]*v[2];
  
  v[0] = c[0];
  v[1] = c[1];
  v[2] = c[2];
}

//------------------------------------------------------------------------------
//computes delta R

void CorotSolver::computeRotMat(double *angle, double mat[3][3])
{
  // trig functions of angles
  double c1 = cos(angle[0]);
  double s1 = sin(angle[0]);
  double c2 = cos(angle[1]);
  double s2 = sin(angle[1]);
  double c3 = cos(angle[2]);
  double s3 = sin(angle[2]);

  /* 
     compute rotation matrix computed as R1.R2.R3
     where R1 is rotation about z
           R2 is rotation about y
           R3 is rotation about x
  */

  mat[0][0] = c1*c2;
  mat[0][1] = c1*s2*s3 - c3*s1;
  mat[0][2] = c1*c3*s2 + s1*s3;

  mat[1][0] = c2*s1;
  mat[1][1] = c1*c3+s1*s2*s3;
  mat[1][2] = c3*s1*s2 - c1*s3;

  mat[2][0] = -s2;
  mat[2][1] = c2*s3;
  mat[2][2] = c2*c3;
}

//------------------------------------------------------------------------------

void CorotSolver::printRotMat(double mat[3][3])
{
  com->fprintf(stderr," Rotation matrix = \n");
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++)
      com->fprintf(stderr," %e  ",mat[i][j]);
    com->fprintf(stderr,"\n");
  }
}

//------------------------------------------------------------------------------

void CorotSolver::computeCG(DistSVec<double,3> &X, double cg[3])
{
  int iSub, i, j, size = 0;
  const DistInfo &distInfo = X.info();

  double (*allcg)[3] = reinterpret_cast<double (*)[3]>
                       (alloca(sizeof(double) * 3*distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub)
    for (j=0; j<3; ++j)
      allcg[iSub][j] = 0.0;

#pragma omp parallel for private(i,j) reduction(+: size)
  for (iSub=0; iSub<numLocSub; ++iSub) {

    double (*locX)[3] = X.subData(iSub);
    bool *flag = X.getMasterFlag(iSub);

    int locsize = 0;
    double loccg[3] = {0.0, 0.0, 0.0};

    for (i=0; i<nInterfNd[iSub]; ++i) {
      if (flag[ interfNd[iSub][i] ]) {
	++locsize;
	for (j=0; j<3; ++j)
	  loccg[j] += locX[ interfNd[iSub][i] ][j];
      }
    }

    size += locsize;

    for (j=0; j<3; ++j)
      allcg[distInfo.locSubToGlobSub[iSub]][j] = loccg[j];

  }

  com->globalSum(1, &size);
  com->globalSum(3*distInfo.numGlobSub, reinterpret_cast<double *>(allcg));

  for (j=0; j<3; ++j)
    cg[j] = 0.0;

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub)
    for (j=0; j<3; ++j)
      cg[j] += allcg[iSub][j];

  double invTotNd = 1.0 / double(size);

  for (j=0; j<3; ++j)
    cg[j] *= invTotNd;
}

//------------------------------------------------------------------------------

void CorotSolver::computeInfNodeRot(double RR[3][3], DistSVec<double,3> &X, 
				    double cg00[3], double cg1[3])
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; iSub++)  {
    SVec<double, 3> &subNd = X(iSub);
    SVec<double, 3> &subX0 = X0(iSub);
    double x[3];
    for (int i = 0; i < nInfNd[iSub]; i++)  {
      x[0] = subX0[infNd[iSub][i]][0] - cg00[0];
      x[1] = subX0[infNd[iSub][i]][1] - cg00[1];
      x[2] = subX0[infNd[iSub][i]][2] - cg00[2];

      rotLocVec(RR, x);

      subNd[infNd[iSub][i]][0] = x[0] + cg1[0];
      subNd[infNd[iSub][i]][1] = x[1] + cg1[1];
      subNd[infNd[iSub][i]][2] = x[2] + cg1[2];

    }
  }
}

//------------------------------------------------------------------------------
// Create Xtilde{b,i} by rotating Xn, and then update dX so that
// Xtilde+dx = Xn+1

void CorotSolver::computeNodeRot(double dRot[3][3], DistSVec<double,3> &X, 
				 DistSVec<double,3> &deltaX, double deltaT[3])
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)  {

    double (*locX)[3] = X.subData(iSub);
    double (*locdX)[3] = deltaX.subData(iSub);

    int i;
    for (i = 0; i < nInterfNd[iSub]; ++i)  {
      locdX[interfNd[iSub][i]][0] += locX[interfNd[iSub][i]][0];
      locdX[interfNd[iSub][i]][1] += locX[interfNd[iSub][i]][1];
      locdX[interfNd[iSub][i]][2] += locX[interfNd[iSub][i]][2];
      
      rotLocVec(dRot, locX[interfNd[iSub][i]]);

      locX[interfNd[iSub][i]][0] += deltaT[0];
      locX[interfNd[iSub][i]][1] += deltaT[1];
      locX[interfNd[iSub][i]][2] += deltaT[2];

      locdX[interfNd[iSub][i]][0] -= locX[interfNd[iSub][i]][0];
      locdX[interfNd[iSub][i]][1] -= locX[interfNd[iSub][i]][1];
      locdX[interfNd[iSub][i]][2] -= locX[interfNd[iSub][i]][2];
    }
    for (i = 0; i < nInternalNd[iSub]; ++i)  {
      rotLocVec(dRot, locX[internalNd[iSub][i]]);

      locX[internalNd[iSub][i]][0] += deltaT[0];
      locX[internalNd[iSub][i]][1] += deltaT[1];
      locX[internalNd[iSub][i]][2] += deltaT[2];
    }
  }
}

//------------------------------------------------------------------------------

void CorotSolver::solveDeltaRot(DistSVec<double,3> &X, double deltaRot[3][3], double cg1[3])
{
  double jac[3][3], grad[3];
  double dRot[3][3];

  // Initialize deltaRot to Identity Matrix
  double zeroRot[3] = {0.0, 0.0, 0.0};
  computeRotMat(zeroRot, deltaRot);

  int maxits = 10;
  double atol = 1.e-12;
  double rtol = 1.e-10;
  double res0, res, target;

  for (int iter = 0; iter < maxits; ++iter) {

    // compute rotation gradients and derivatives
    computeRotGradAndJac(X, R, cg1, grad, jac);

    //HB: zero terms depending on the axis of the plane of symmetry
    switch(SymAxis) {
      case(CorotSolver::AXIS_X): 
        grad[1] = grad[2] = 0.0;
        jac[0][1] = jac[1][0] = jac[0][2] = jac[2][0] = 0.0;
        jac[1][2] = jac[2][1] = 0.0;
        jac[1][1] = jac[2][2] = 1.0;
        break;
      
      case(CorotSolver::AXIS_Y):
        grad[0] = grad[2] = 0.0;
        jac[0][1] = jac[1][0] = jac[0][2] = jac[2][0] = 0.0;
        jac[1][2] = jac[2][1] = 0.0;
        jac[0][0] = jac[2][2] = 1.0;
        break;

      case(CorotSolver::AXIS_Z):
        grad[0] = grad[1] = 0.0;
        jac[0][1] = jac[1][0] = jac[0][2] = jac[2][0] = 0.0;
        jac[1][2] = jac[2][1] = 0.0;
        jac[0][0] = jac[1][1] = 1.0;
        break;
    }

    res = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);

#ifdef HB_COROTSOLVER_DEBUG
  com->printf(1, "*** incremental rotation solver residual at it %d = %.2e\n", iter,res);
#endif

    if (iter == 0) { res0 = res; target = rtol*res0; } 
    
    //if (res == 0.0 || res <= target) break;
    if (res <= atol || res <= target) break; //HB: add absolute tol. to avoid iterations
                                             //if the initial residual is already small

    // rotation results come back in grad
    solveRotMat(jac, grad);

    grad[0] *= -1.0;
    grad[1] *= -1.0;
    grad[2] *= -1.0;

    // update dRot
    computeRotMat(grad, dRot);

    denseMatrixTimesDenseMatrixInPlace(dRot, R);
    denseMatrixTimesDenseMatrixInPlace(dRot, deltaRot);

  }

  if (res>target & res>atol) {
    com->printf(1, "*** Warning: incremental rotation solver reached %d its", maxits);
    com->printf(1, " (initial res = %.2e, final res=%.2e, target=%.2e)\n", res0, res, target);
  }

#ifdef HB_COROTSOLVER_DEBUG
  printRotMat(deltaRot);
#endif
}

//------------------------------------------------------------------------------

void CorotSolver::solveRotMat(double m[3][3], double v[3])  
{
  int i,j,k;
  for (i = 0; i < 2; i++)
    for (j = i+1; j < 3; ++j) {
      double coef = m[j][i]/m[i][i];
      for (k = i+1; k < 3; ++k)
        m[j][k] -= coef*m[i][k];
      v[j] -= coef*v[i];
    }
  
  for (i=2; i >= 0; i--) {
    for (j=2; j > i; j--)
      v[i] -= m[i][j] * v[j];
    v[i] /= m[i][i];

  }
}

//------------------------------------------------------------------------------

void CorotSolver::solve(DistSVec<double,3> &dX, DistSVec<double,3> &X, BCApplier* meshMotionBCs)
{
  Xtilde = X + dX;

  // compute cg(n+1)
  double cg1[3];
  computeCG(Xtilde, cg1);

#define SYMMETRIC_CG
#ifdef SYMMETRIC_CG
   //HB: force translation to be parallel to the symmetry plane 
  switch(SymAxis) {
    case(CorotSolver::AXIS_X):
      cg1[0] = cg0[0];
      break;

    case(CorotSolver::AXIS_Y):
      cg1[1] = cg0[1];
      break;

    case(CorotSolver::AXIS_Z):
      cg1[2] = cg0[2];
      break;
  }
#endif

#ifdef HB_COROTSOLVER_DEBUG
  com->fprintf(stderr," -> In CorotSolver::solve, cg1 = %.2e  %.2e  %.2e\n",cg1[0],cg1[1],cg1[2]);
#endif

  // solve for the incremental rotations via Newton-Rhapson
  double deltaRot[3][3];
  solveDeltaRot(Xtilde, deltaRot, cg1);

  // Update the boundary xtilde to include the gap rotation
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    double (*dx)[3] = dX.subData(iSub);
    for (int i = 0; i < nInterfNd[iSub]; ++i) {
      double nd[3];
      nd[0] = dx[ interfNd[iSub][i] ] [0];
      nd[1] = dx[ interfNd[iSub][i] ] [1];
      nd[2] = dx[ interfNd[iSub][i] ] [2];
      nd[0] -= gapVec[iSub][i][0];
      nd[1] -= gapVec[iSub][i][1];
      nd[2] -= gapVec[iSub][i][2];
      double gv[3];
      gv[0] = gapVec[iSub][i][0];
      gv[1] = gapVec[iSub][i][1];
      gv[2] = gapVec[iSub][i][2];
      rotLocVec(R, gv);

      nd[0] += gv[0];
      nd[1] += gv[1];
      nd[2] += gv[2]; 

      dx[ interfNd[iSub][i] ] [0] = nd[0];
      dx[ interfNd[iSub][i] ] [1] = nd[1];
      dx[ interfNd[iSub][i] ] [2] = nd[2];
    }
  }

  // rotate cgN by incremental rotation matrix, deltaRot
  rotLocVec(deltaRot, cgN);

  int i;
  double deltaT[3];
  for (i = 0; i < 3; i++)
    deltaT[i] = cg1[i] - cgN[i];

//#ifdef HB_COROTSOLVER_DEBUG
//  com->fprintf(stderr," -> In CorotSolver::solve, cg1 - cgN = %.2e  %.2e  %.2e\n",deltaT[0],deltaT[1],deltaT[2]);
//#endif

  // compute Xinf(n+1)
  computeInfNodeRot(R, X, cg0, cg1);

  // compute Xboundary(n+1) and Xinternal(n) guess and update dX
  computeNodeRot(deltaRot, X, dX, deltaT);

  // update CG(n)
  for (i = 0; i < 3; i++)
    cgN[i] = cg1[i];

#ifdef HB_COROTSOLVER_DEBUG
  com->fprintf(stderr," -> In CorotSolver::solve, cgN = %.2e  %.2e  %.2e\n",cgN[0],cgN[1],cgN[2]);
#endif

  //HB: we may also want to apply the rotation to the normal directions used in the
  //projections in BCApplier. Not currently done, because it is assumed that in most 
  //practical cases (to date) there would be only one sliding plane, and this sliding 
  //plane would also be the symmetry plane so that the plane normal wouldn't be affected
  //by the rotation of the fluid mesh around the symmetry plane normal ...
  //if(meshMotionBCs) meshMotionBCs->rotateProjNormal(deltaRot);
}

//------------------------------------------------------------------------------

void
CorotSolver::setup(DistSVec<double,3> &X)
{
  // NOTE: This is the restart procedure that is not based on saving the
  // rotation tensor. The result is that we may hit a singularity sometimes
  // compute cg(n+1)
  double cg1[3];
  computeCG(X, cg1);

  // solve for the incremental rotations via Newton-Rhapson
  double deltaRot[3][3];
  solveDeltaRot(X, deltaRot, cg1);
  
  // update CG(n)
  int i;
  for (i = 0; i < 3; i++)
    cgN[i] = cg1[i];
}
