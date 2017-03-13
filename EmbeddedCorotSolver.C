#include <cstdlib>
#include <cmath>

#include <EmbeddedCorotSolver.h>

#include <MatchNode.h>
#include <Domain.h>
#include <Vector3D.h>
#include <DenseMatrixOps.h>

#include <BCApplier.h>

//------------------------------------------------------------------------------

EmbeddedCorotSolver::EmbeddedCorotSolver(IoData &iodata, MatchNodeSet **mns, Domain *dom, double *Xstruct, int nNodes, int (*structElem)[3])
  : iod(iodata), matchNodes(mns), domain(dom), Xs0(Xstruct), stElem(structElem), X0(dom->getNodeDistInfo())
{

  numStNodes = nNodes;

  numLocSub = domain->getNumLocSub();

  com = domain->getCommunicator();

  domain->getReferenceMeshPosition(X0);
  meshMotionBCs = domain->getMeshMotionBCs();

  computeCG(Xs0, cg0);
  
  cgN[0] = cg0[0];
  cgN[1] = cg0[1];
  cgN[2] = cg0[2];

  double zeroRot[3] = {0.0, 0.0, 0.0};
  computeRotMat(zeroRot, R);

  type = EmbeddedCorotSolver::BASIC;
  if (iod.dmesh.type==DefoMeshMotionData::COROTATIONAL)
    type = EmbeddedCorotSolver::COROTATIONAL;

  //HB: look if a symmetry plane was specified in the input file
  n[0] = iod.dmesh.symmetry.nx;
  n[1] = iod.dmesh.symmetry.ny;
  n[2] = iod.dmesh.symmetry.nz;
  double nrm= sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  SymAxis = EmbeddedCorotSolver::NONE;
  if (type == EmbeddedCorotSolver::BASIC) {
    if(nrm!=0.0) {
      n[0] /= nrm; n[1] /= nrm; n[2] /= nrm; 
      if((fabs(n[0])==1.0) & (n[1]==0.0) & (n[2]==0.0))
        SymAxis = EmbeddedCorotSolver::AXIS_X;
      else if((n[0]==0.0) & (fabs(n[1])==1.0) & (n[2]==0.0))
        SymAxis = EmbeddedCorotSolver::AXIS_Y;
      else if((n[0]==0.0) & (n[1]==0.0) & (fabs(n[2])==1.0))
        SymAxis = EmbeddedCorotSolver::AXIS_Z;
      else {
        com->fprintf(stderr," *** ERROR: embedded corotational solver only supports a canonical plane as a symmetry plane.\n");
        exit(-1);
      }
    } 
    switch(SymAxis) {
      case(EmbeddedCorotSolver::NONE):
        com->fprintf(stderr," ... No symmetry plane is used in the embedded corotational solver.\n");
        com->fprintf(stderr,"     -> the 3 rotations axis (X,Y,Z) & 3 translation axis (X,Y,Z) are used.\n");
        break;
      case(EmbeddedCorotSolver::AXIS_X):
        com->fprintf(stderr," ... Symmetry plane of normal X is used in the embedded corotational solver.\n");
        com->fprintf(stderr,"     -> only rotation around axis X & translations in the Y-Z plane are allowed.\n");
        break;
      case(EmbeddedCorotSolver::AXIS_Y):
        com->fprintf(stderr," ... Symmetry plane of normal Y is used in the embedded corotational solver.\n");
        com->fprintf(stderr,"     -> only rotation around axis Y & translations in the X-Z plane are allowed.\n");
        break;
      case(EmbeddedCorotSolver::AXIS_Z):
        com->fprintf(stderr," ... Symmetry plane of normal Z is used in the embedded corotational solver.\n");
        com->fprintf(stderr,"     -> only rotation around axis Z & translations in the X-Y plane are allowed.\n");
        break;
    }
  }

  if (matchNodes) {
    int (*ptr)[2] = new int[numLocSub][2];

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub)
      ptr[iSub][0] = matchNodes[iSub]->size();

    ptr[0][1] = 0;
    for (int iSub=1; iSub<numLocSub; ++iSub)
      ptr[iSub][1] = ptr[iSub - 1][1] + ptr[iSub - 1][0];

    int numCPUMatchedNodes = ptr[numLocSub - 1][1] + ptr[numLocSub - 1][0];

    // gather the global matched point nodes for this fluid CPU

    int (*cpuMatchNodes)[3] = 0;

    if (numCPUMatchedNodes > 0) {
      cpuMatchNodes = new int[numCPUMatchedNodes][3];

#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub)
        matchNodes[iSub]->exportInfo(iSub, cpuMatchNodes + ptr[iSub][1]);
    }

//get embedded surface matcher
    char *embsurfmatch;
    int sp = strlen(iod.input.prefix) + 1;

    embsurfmatch = new char[sp + strlen(iod.input.embsurfmatch)];
    sprintf(embsurfmatch,"%s%s", iod.input.prefix, iod.input.embsurfmatch);

    FILE *fp = fopen(embsurfmatch, "r");

    if (!fp)  {
      fprintf(stderr, "*** ERROR: Embedded surface match file \'%s\' file does not exist\n",embsurfmatch);
      exit(-1);
    }

    int numNodes;
    char line[MAXLINE];
    char *toto = fgets(line, MAXLINE, fp);
    toto = fgets(line, MAXLINE, fp);
    toto = fgets(line, MAXLINE, fp);
    toto = fgets(line, MAXLINE, fp);
    sscanf(line, "%*s %d", &numNodes);

    int *elempos = new int[numNodes];
    double (*pos)[2] = new double[numNodes][2];
  // read match points
    for (int i = 0; i < numNodes; i++) {
      toto = fgets(line, MAXLINE, fp);
      sscanf(line, "%d %lf %lf", &(elempos[i]), &(pos[i][0]), &(pos[i][1]));
      elempos[i]--;
    }
    fclose(fp);

    double (*xs0)[3] = new double[nNodes][3];
    for (int i=0; i<nNodes; i++)
      for (int j=0; j<3; j++)
        xs0[i][j] = Xs0[3*i + j];

    for (int i=0; i<numCPUMatchedNodes; ++i) {
      matchNodes[cpuMatchNodes[i][2]]->setBufferPosition(cpuMatchNodes[i][0], elempos[cpuMatchNodes[i][1]], pos[cpuMatchNodes[i][1]], stElem, xs0);
    }

    delete[] xs0;
    delete[] elempos;
    delete[] pos;
    delete[] cpuMatchNodes;
    delete[] ptr;
    delete[] embsurfmatch;
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

void EmbeddedCorotSolver::computeRotGradAndJac(double *Xs, double RR[3][3], 
					 double cg1[3], double grad[3], double jac[3][3])
{
  int i, j;

  for (i = 0; i < 3; i++)  {
    grad[i] = 0.0;
    for (j = 0; j < 3; j++)
      jac[i][j] = 0.0;
  }

  Vec3D rotGrad[3];

  for (i = 0; i < numStNodes; i++) {
    double rd[3];
    // rotate the local vectors using R(n-1)
    rd[0] = Xs0[3*i+0] - cg0[0];
    rd[1] = Xs0[3*i+1] - cg0[1];
    rd[2] = Xs0[3*i+2] - cg0[2];
    rotLocVec(RR, rd);

    Vec3D eVec;
    //compute freq. used values
    eVec[0] = Xs[3*i+0] - cg1[0] - rd[0];
    eVec[1] = Xs[3*i+1] - cg1[1] - rd[1];
    eVec[2] = Xs[3*i+2] - cg1[2] - rd[2];

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

}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::rotLocVec(double mat[3][3], double v[3]) 
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

void EmbeddedCorotSolver::computeRotMat(double *angle, double mat[3][3])
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

void EmbeddedCorotSolver::printRotMat(double mat[3][3])
{
  com->fprintf(stderr," Rotation matrix = \n");
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++)
      com->fprintf(stderr," %e  ",mat[i][j]);
    com->fprintf(stderr,"\n");
  }
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::computeCG(double *Xs, double cg[3])
{
  
  for (int j=0; j<3; ++j) cg[j] = 0.0;

  for (int i=0; i<numStNodes; ++i) {
    for (int j=0; j<3; ++j)
      cg[j] += Xs[3*i+j];
  }

  double invTotNd = 1.0 / double(numStNodes);

  for (int j=0; j<3; ++j)
    cg[j] *= invTotNd;
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::computeNodeRot(double RR[3][3], DistSVec<double,3> &X, 
				 double cg00[3], double cg1[3])
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)  {

    double (*x)[3]  = X.subData(iSub);
    double (*x0)[3] = X0.subData(iSub);

    for (int i = 0; i < X.subSize(iSub); ++i)  {
      x[i][0] = x0[i][0] - cg00[0];
      x[i][1] = x0[i][1] - cg00[1];
      x[i][2] = x0[i][2] - cg00[2];

      rotLocVec(RR, x[i]);

      x[i][0] = x[i][0] + cg1[0];
      x[i][1] = x[i][1] + cg1[1];
      x[i][2] = x[i][2] + cg1[2];
    }
  }
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::computeDeltaNodeRot(double RR[3][3], DistSVec<double,3> &X, DistSVec<double,3> &dX, 
				 double cg00[3], double cg1[3])
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)  {

    double (*x)[3]  = X.subData(iSub);
    double (*dx)[3] = dX.subData(iSub);

    for (int i = 0; i < X.subSize(iSub); ++i)  {
      dx[i][0] = x[i][0] - cg00[0];
      dx[i][1] = x[i][1] - cg00[1];
      dx[i][2] = x[i][2] - cg00[2];

      rotLocVec(RR, dx[i]);

      dx[i][0] = dx[i][0] + cg1[0];
      dx[i][1] = dx[i][1] + cg1[1];
      dx[i][2] = dx[i][2] + cg1[2];

      dx[i][0] = dx[i][0] - x[i][0];
      dx[i][1] = dx[i][1] - x[i][1];
      dx[i][2] = dx[i][2] - x[i][2];

      x[i][0] = x[i][0] + dx[i][0];
      x[i][1] = x[i][1] + dx[i][1];
      x[i][2] = x[i][2] + dx[i][2];
    }
  }
}

//------------------------------------------------------------------------------
void EmbeddedCorotSolver::solveDeltaRot(double *Xs, double deltaRot[3][3], double cg1[3])
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
    computeRotGradAndJac(Xs, R, cg1, grad, jac);

    //HB: zero terms depending on the axis of the plane of symmetry
    switch(SymAxis) {
      case(EmbeddedCorotSolver::AXIS_X): // grad[2] is rotation angle about X-axis
        grad[0] = grad[1] = 0.0;
        jac[0][1] = jac[1][0] = jac[0][2] = jac[2][0] = 0.0;
        jac[1][2] = jac[2][1] = 0.0;
        jac[0][0] = jac[1][1] = 1.0;
        break;

      case(EmbeddedCorotSolver::AXIS_Y):
        grad[0] = grad[2] = 0.0;
        jac[0][1] = jac[1][0] = jac[0][2] = jac[2][0] = 0.0;
        jac[1][2] = jac[2][1] = 0.0;
        jac[0][0] = jac[2][2] = 1.0;
        break;

      case(EmbeddedCorotSolver::AXIS_Z): 
        grad[1] = grad[2] = 0.0;
        jac[0][1] = jac[1][0] = jac[0][2] = jac[2][0] = 0.0;
        jac[1][2] = jac[2][1] = 0.0;
        jac[1][1] = jac[2][2] = 1.0;
        break;
      
    }

    res = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);

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
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::solveRotMat(double m[3][3], double v[3])  
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
void EmbeddedCorotSolver::applyProjector(DistSVec<double,3> &Xdot)
{
   if(meshMotionBCs) meshMotionBCs->applyP(Xdot);
}

//------------------------------------------------------------------------------
void EmbeddedCorotSolver::fixNodes(double *Xs, int nNodes, DistSVec<double,3> &X, DistSVec<double,3> &dX)
{
  DistSVec<double,3> dXp(dX);
  if (matchNodes) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      double (*x)[3] = X.subData(iSub);
      double (*dx)[3] = dX.subData(iSub);
      double (*dxp)[3] = dXp.subData(iSub);
      double (*xs)[3] = new double[nNodes][3];
      for (int i=0; i<nNodes; i++)
        for (int j=0; j<3; j++)
          xs[i][j] = Xs[3*i + j];

      matchNodes[iSub]->getDisplacement(xs, stElem, x, dx, dxp, iod.ref.rv.tlength);
      delete[] xs;
    }
  }
  if(meshMotionBCs) {
    double meandX[3];
    computeMeanDXForSlidingPlane(dX,meandX);
    meshMotionBCs->applyD2(dX,meandX);
    meshMotionBCs->applyP(dX);
  }
  dX = dX - dXp;
}

//------------------------------------------------------------------------------
void EmbeddedCorotSolver::computeMeanDXForSlidingPlane(DistSVec<double,3> &dX, double meandX[3])
{
  int iSub, i, j, size = 0;
  const DistInfo &distInfo = dX.info();

  int** DofType = (meshMotionBCs) ? meshMotionBCs->getDofType(): 0;

  if(DofType)
  {
    double (*sumdX)[3] = reinterpret_cast<double (*)[3]>
                         (alloca(sizeof(double) * 3*distInfo.numGlobSub));

    for (iSub=0; iSub<distInfo.numGlobSub; ++iSub)
      for (j=0; j<3; ++j)
        sumdX[iSub][j] = 0.0;

    #pragma omp parallel for private(i,j) reduction(+: size)
    for (iSub=0; iSub<numLocSub; ++iSub) {

      double (*dx)[3] = dX.subData(iSub);
      int (*dofType)[3] = reinterpret_cast<int (*)[3]>(DofType[iSub]);
      bool *flag = dX.getMasterFlag(iSub);

      int locsize = 0;
      double locsumdx[3] = {0.0, 0.0, 0.0};

      for(int i=0;i<dX.subSize(iSub); i++) {
        if (flag[i]) {
	  if (dofType[i][0]==BC_MATCHEDSLIDE ||
	      dofType[i][1]==BC_MATCHEDSLIDE ||
	      dofType[i][2]==BC_MATCHEDSLIDE) {
            ++locsize;
            for (j=0; j<3; ++j)
              if (dofType[i][j]==BC_MATCHEDSLIDE) locsumdx[j] += dx[i][j];
          }
        }
      }

      size += locsize;

      for (j=0; j<3; ++j)
        sumdX[distInfo.locSubToGlobSub[iSub]][j] = locsumdx[j];

    }

    com->globalSum(1, &size);
    com->globalSum(3*distInfo.numGlobSub, reinterpret_cast<double *>(sumdX));

    for (j=0; j<3; ++j)
      meandX[j] = 0.0;

    for (iSub=0; iSub<distInfo.numGlobSub; ++iSub)
      for (j=0; j<3; ++j)
        meandX[j] += sumdX[iSub][j];

    double invTotN;
    if (size != 0) invTotN = 1.0 / double(size);
    else invTotN = 0.;

    for (j=0; j<3; ++j)
      meandX[j] *= invTotN;
  }
}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::solve(double *Xtilde, int nNodes, DistSVec<double,3> &X, DistSVec<double,3> &dX)
{

  if(nNodes!=numStNodes) {
    com->fprintf(stderr,"Number of structure nodes has changed!\n");
    exit(-1);
  }

  // compute cg(n+1)
  double cg1[3];
  computeCG(Xtilde, cg1);

  switch(SymAxis) {
    case(EmbeddedCorotSolver::AXIS_X):
      cg1[0] = cg0[0];
      break;

    case(EmbeddedCorotSolver::AXIS_Y):
      cg1[1] = cg0[1];
      break;

    case(EmbeddedCorotSolver::AXIS_Z):
      cg1[2] = cg0[2];
      break;
  }

  // solve for the incremental rotations via Newton-Rhapson
  double deltaRot[3][3];
  solveDeltaRot(Xtilde, deltaRot, cg1);

// Update node
  if (type==EmbeddedCorotSolver::BASIC)
    computeNodeRot(R, X, cg0, cg1);
  else if (type==EmbeddedCorotSolver::COROTATIONAL) {
    computeDeltaNodeRot(deltaRot, X, dX, cgN, cg1);
    fixNodes(Xtilde,nNodes,X,dX);
  }

  cgN[0] = cg1[0];
  cgN[1] = cg1[1];
  cgN[2] = cg1[2];

}

//------------------------------------------------------------------------------

void EmbeddedCorotSolver::setup(double *Xtilde, int nNodes)
{

  if(nNodes!=numStNodes) {
    com->fprintf(stderr,"Number of structure nodes has changed!\n");
    exit(-1);
  }

  // compute cg(n+1)
  double cg1[3];
  computeCG(Xtilde, cg1);

  switch(SymAxis) {
    case(EmbeddedCorotSolver::AXIS_X):
      cg1[0] = cg0[0];
      break;

    case(EmbeddedCorotSolver::AXIS_Y):
      cg1[1] = cg0[1];
      break;

    case(EmbeddedCorotSolver::AXIS_Z):
      cg1[2] = cg0[2];
      break;
  }

  // solve for the incremental rotations via Newton-Rhapson
  double deltaRot[3][3];
  solveDeltaRot(Xtilde, deltaRot, cg1);

  cgN[0] = cg1[0];
  cgN[1] = cg1[1];
  cgN[2] = cg1[2];
}

//------------------------------------------------------------------------------
