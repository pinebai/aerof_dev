#include <cmath>
#include <cstring>

#include <KspSolver.h>

#include <IoData.h>
#include <KspConvCriterion.h>
#include <complex>

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
KspSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::
KspSolver(KspData &data, MatVecProdOp *mvp, PrecOp *pc, IoOp *io)
{

  this->mvpOp = mvp;
  this->pcOp = pc;
  this->ioOp = io;

  this->maxits = data.maxIts;

  this->absoluteEps = data.absoluteEps;

  this->kspConvCriterion = new KspConvCriterion(data);

  this->kspBinaryOutput = NULL;

  if (data.checkFinalRes == KspData::YES)
    this->checkFinalRes = true;
  else
    this->checkFinalRes = false;

  if (strcmp(data.output, "") == 0)
    this->output = 0;
  else if (strcmp(data.output, "stdout") == 0)
    this->output = stdout;
  else if (strcmp(data.output, "stderr") == 0)
    this->output = stderr;
  else {
    this->output = fopen(data.output, "w");
    if (!this->output) {
      this->ioOp->fprintf(stderr, "*** Error: could not open \'%s\'\n", data.output);
      exit(1);
    }
  }
}

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
void
KspSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::setup(int nlit, int nlmaxits, VecType &b)
{

  this->eps = this->kspConvCriterion->compute(nlit, nlmaxits, b.norm());

}

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
void
KspSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::setKspBinaryOutput(KspBinaryOutput<VecType>* kspBinOut)
{

  this->kspBinaryOutput = kspBinOut;

}

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp>
RichardsonSolver<VecType,MatVecProdOp,PrecOp,IoOp>::
RichardsonSolver(const typename VecType::InfoType &info, KspData &data, 
		 MatVecProdOp *A, PrecOp *P, IoOp *IO) : 
  KspSolver<VecType,MatVecProdOp,PrecOp,IoOp>(data, A, P, IO), dx(info), r(info)
{

}

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp>
int
RichardsonSolver<VecType,MatVecProdOp,PrecOp,IoOp>::solve(VecType &b, VecType &x)
{

  double res, target;

  if (this->output) this->ioOp->fprintf(this->output, "Richardson iterations:\n");

  int iter;
  for (iter=0; iter<this->maxits; ++iter) {
    this->mvpOp->apply(x, r);
    r = b - r;
    res = sqrt(r*r);

    if (iter == 0) 
      target = this->eps * b.norm(); // res;

    if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", iter, res, target);

    if (res == 0.0 || res <= target) break;

    this->pcOp->apply(r, dx);
    x += dx;
  }

  if (this->checkFinalRes) {
    this->mvpOp->apply(x, r);
    r = b - r;
    if (this->output) 
      this->ioOp->fprintf(this->output, "  %d %e (actual residual norm)\n", iter, sqrt(r*r));
  }

  this->ioOp->printf(5, "Richardson solver: its=%d, res=%.2e, target=%.2e\n", iter, res, target);
  if (iter == this->maxits) {
    this->ioOp->printf(1, "*** Warning: Richardson solver reached %d its", this->maxits);
    this->ioOp->printf(1, " (res=%.2e, target=%.2e)\n", res, target);
  }

  return iter;

}
//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp>
int
RichardsonSolver<VecType,MatVecProdOp,PrecOp,IoOp>::solveLS(VecType &b, VecType &x)
{
  double res, target;

  if (this->output) this->ioOp->fprintf(this->output, "Richardson iterations:\n");

  int iter;
  for (iter=0; iter<this->maxits; ++iter) {
    this->mvpOp->apply(x, r);
    r = b - r;
    res = sqrt(r*r);

    if (iter == 0)
      target = this->eps * b.norm(); // res;

    if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", iter, res, target);

    if (res == 0.0 || res <= target) break;

    this->pcOp->apply(r, dx);
    x += dx;
  }

  if (this->checkFinalRes) {
    this->mvpOp->apply(x, r);
    r = b - r;
    if (this->output)
      this->ioOp->fprintf(this->output, "  %d %e (actual residual norm)\n", iter, sqrt(r*r));
  }

  this->ioOp->printf(5, "Richardson solver: its=%d, res=%.2e, target=%.2e\n", iter, res, target);
  if (iter == this->maxits) {
    this->ioOp->printf(1, "*** Warning: Richardson solver reached %d its", this->maxits);
    this->ioOp->printf(1, " (res=%.2e, target=%.2e)\n", res, target);
  }

  return iter;

}

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp>
CgSolver<VecType,MatVecProdOp,PrecOp,IoOp>::
CgSolver(const typename VecType::InfoType &info, KspData &data, 
	 MatVecProdOp *A, PrecOp *P, IoOp *IO) : 
  KspSolver<VecType,MatVecProdOp,PrecOp,IoOp>(data, A, P, IO), 
  r(info), Ap(info), y(info), p(info)
{

}

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp>
int 
CgSolver<VecType,MatVecProdOp,PrecOp,IoOp>::solve(VecType &f, VecType &x)
{
  int typePrec = 1;

  this->mvpOp->apply(x, Ap);
  r = f - Ap;

  double res = sqrt(r*r);
  double target = this->eps * sqrt(f*f);

  if (this->output) {
    this->ioOp->fprintf(this->output, "Cg iterations:\n");
    this->ioOp->fprintf(this->output, "  %d %e %e\n", 0, res, target);
  }

  if (res <= target) return 0;
  double initial = res;

  if(typePrec) this->pcOp->apply(r, p);
  else p = r;
 
  this->mvpOp->apply(p, Ap);

  double pAp = p * Ap;
  double nu = p * r / pAp;

  x += nu * p;
  r -= nu * Ap;

  res = sqrt(r*r);

  if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", 1, res, target);

  int iter;
  for (iter = 1; iter < this->maxits && (res > target); ++iter) {

    if(typePrec) this->pcOp->apply(r, y);
    else y = r;

    double alpha = Ap * y / pAp;
    p = y - alpha * p;

    this->mvpOp->apply(p, Ap);
    pAp = p * Ap;
    nu = p * r / pAp;

    x += nu * p;
    r -= nu * Ap;

    res = sqrt(r*r);

    if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", iter+1, res, target);

  }

  this->ioOp->printf(5, "Cg solver: its=%d, res=%.2e, target=%.2e\n", iter, res, target);
  if (iter == this->maxits) {
    this->ioOp->printf(1, "*** Warning: Cg solver reached %d its", this->maxits);
    this->ioOp->printf(1, " (Residual: initial=%.2e, res=%.2e, target=%.2e)\n", initial, res, target);
  }

  return iter;

}

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp>
int
CgSolver<VecType,MatVecProdOp,PrecOp,IoOp>::solveLS(VecType &f, VecType &x)
{
  this->mvpOp->apply(x, Ap);
  r = f - Ap;

  double res = sqrt(r*r);
  double target = this->eps * sqrt(f*f);

  if (this->output) {
    this->ioOp->fprintf(this->output, "Cg iterations:\n");
    this->ioOp->fprintf(this->output, "  %d %e %e\n", 0, res, target);
  }

  if (res <= target) return 0;

  this->pcOp->apply(r, p);
  this->mvpOp->apply(p, Ap);

  double pAp = p * Ap;
  double nu = p * r / pAp;

  x += nu * p;
  r -= nu * Ap;

  res = sqrt(r*r);

  if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", 1, res, target);

  int iter;
  for (iter = 1; iter < this->maxits && (res > target); ++iter) {
    this->pcOp->apply(r, y);
    double alpha = Ap * y / pAp;
    p = y - alpha * p;
    this->mvpOp->apply(p, Ap);
    pAp = p * Ap;
    nu = p * r / pAp;

    x += nu * p;
    r -= nu * Ap;

    res = sqrt(r*r);

    if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", iter+1, res, target);
  }

  this->ioOp->printf(5, "Cg solver: its=%d, res=%.2e, target=%.2e\n", iter, res, target);
  if (iter == this->maxits) {
    this->ioOp->printf(1, "*** Warning: Cg solver reached %d its", this->maxits);
    this->ioOp->printf(1, " (res=%.2e, target=%.2e)\n", res, target);
  }

  return iter;

}

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
GmresSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::
GmresSolver(const typename VecType::InfoType &info, KspData &data, 
	    MatVecProdOp *A, PrecOp *P, IoOp *IO) : 
  KspSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>(data, A, P, IO),
  numVec(data.numVectors), cs(data.numVectors, 2), H(data.numVectors+1, data.numVectors), 
  g(data.numVectors+1), y(data.numVectors), V(data.numVectors+1, info), w(info), r(info)
{

  double sizeMB = (numVec + 1) * V[0].sizeMB();

  outputConvergenceInfo = true;

  this->ioOp->globalSum(1, &sizeMB);

  this->ioOp->printf(2, "Memory for Gmres(%d) solver: %3.2f MB\n", numVec, sizeMB);

}

//------------------------------------------------------------------------------
/*
@ARTICLE{saad-schultz-86,
  author = "Saad, Y. and Schultz, M. H.",
  title = "{GMRES}: a Generalized Minimal Residual Algorithm for Solving
           Nonsymmetric Linear Problems",
  journal = siamjscistat,
  year = 1986,
  volume = 7,
  number = 3,
  pages = "856--869",
} 
*/
template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
void GmresSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::
disableConvergenceInfo() {

  outputConvergenceInfo = false;
}

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
int 
GmresSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::solve(VecType &b, VecType &x)
{

  int typePrec = 2;
  //int typePrec = 0;
  double beta, l2res, target, res0;

  int iter = 0;
  int exitLoop = 0;

  int numOutputVecs = 0;

  if (!this->pcOp)
    typePrec = 0;

  this->ioOp->printf(10, "preconditioner type is %d", typePrec);
  do {

    this->mvpOp->apply(x, w);
    r = b - w;

    if (typePrec == 1) { this->pcOp->apply(r, w); r = w; }
    
    beta = r.norm();

    if (iter == 0) {
      target = this->eps * (res0 = b.norm()); // beta;
      if (this->output) this->ioOp->fprintf(this->output, "Gmres(%d) iterations:\n", numVec);
      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e (%e)  \n", 0, beta, target, this->eps);
    } 
    else
      if (this->output) this->ioOp->fprintf(this->output, "  --- restart ---\n");

    if (beta == 0.0) return 0;

    V[0] = (1.0/beta) * r;

    g = 0.0;
    g[0] = beta;

    int j;
    for (j=0; j<numVec; ++j) {

      switch (typePrec) {
      case 0: { this->mvpOp->apply(V[j], w); } break;
      case 1: { this->mvpOp->apply(V[j], r); this->pcOp->apply(r, w); } break;
      case 2: { this->pcOp->apply(V[j], r); this->mvpOp->apply(r, w); } break;
      }



      for (int i=0; i<=j; ++i) {
        // For complex vectors, w has to not be conjugated in the definition
        // of the dot product
        H[i][j] = w * V[i];
        w -= H[i][j] * V[i];
      }

      H[j+1][j] = w.norm();
    
      if (H[j+1][j] == 0.0) {

        applyPreviousRotations(j, H, cs);
        applyNewRotation(j, H, cs, g);

        ++iter; exitLoop = 1; break;

      }
      
      V[j+1] = (1./H[j+1][j]) * w;

      applyPreviousRotations(j, H, cs);
      applyNewRotation(j, H, cs, g);

      l2res = sqrt(sqNorm(g[j+1]));

      ++iter;

      if (this->output)
	this->ioOp->fprintf(this->output, "  %d %e %e \n", iter, l2res, target);

      if (l2res <= target || iter >= this->maxits || 
          l2res <= this->absoluteEps) { exitLoop = 1; break; }

    }

    if (j == numVec) --j;

    backwardSolve(j, H, g, y);

    w = 0.0;

    for (int m=0; m<=j; ++m) w += y[m] * V[m];

    if (typePrec == 2) { this->pcOp->apply(w, r); w = r; }

    x += w;

    numOutputVecs = j+1;

  } while (exitLoop == 0);


  if (this->checkFinalRes) {
    this->mvpOp->apply(x, w);
    r = b - w;
    if (typePrec == 1) { 
      this->pcOp->apply(r, w); 
      r = w; 
    }
    if (this->output) 
      this->ioOp->fprintf(this->output, "  %d %e (actual residual norm)\n", iter, r.norm());
    if (r.norm() > target)  {
      if (this->output)  
        this->ioOp->fprintf(this->output, "  %d %e (actual residual norm) vs l2res: %e, target:%e\n", iter, r.norm(), l2res,
target);
      iter = -999;
    }
  
  }

  this->ioOp->printf(5, "Gmres(%d) solver: its=%d, res=%.2e, target=%.2e\n", numVec, iter, l2res, target);
  if (iter == this->maxits && l2res > target && outputConvergenceInfo) {
    this->ioOp->printf(1, "*** Warning: Gmres(%d) solver reached %d its", numVec, this->maxits);
    this->ioOp->printf(1, " (initial=%.2e, res=%.2e, target=%.2e, ratio = %.2e)\n", res0, l2res, target, l2res/target);
  }

  if (this->kspBinaryOutput) {
    if (typePrec == 2) {  //apply preconditioner before outputting
      for (int iVec=0; iVec<numOutputVecs; ++iVec) {
        this->pcOp->apply(V[iVec], r);
        V[iVec] = r;
      }
    }
    this->kspBinaryOutput->writeKrylovVectors(V, y, numOutputVecs);
  }

  return iter;

}




template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
int
GmresSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::solveNew(VecType &b, VecType &x)
{

  int typePrec = 2;
  //int typePrec = 0;
  double beta, l2res, target, res0;

  int iter = 0;
  int exitLoop = 0;

  int numOutputVecs = 0;

  if (!this->pcOp)
    typePrec = 0;

  this->ioOp->printf(10, "preconditioner type is %d", typePrec);
  do {

    this->mvpOp->applyTranspose(x, w);
    r = b - w;

    if (typePrec == 1) { this->pcOp->applyT(r, w); r = w; }

    beta = r.norm();

    if (iter == 0) {
      target = this->eps * (res0 = b.norm()); // beta;
      if (this->output) this->ioOp->fprintf(this->output, "Gmres(%d) iterations:\n", numVec);
      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e (%e)  \n", 0, beta, target, this->eps);
    }
    else
      if (this->output) this->ioOp->fprintf(this->output, "  --- restart ---\n");

    if (beta == 0.0) return 0;

    V[0] = (1.0/beta) * r;

    g = 0.0;
    g[0] = beta;

    int j;
    for (j=0; j<numVec; ++j) {

      switch (typePrec) {
      case 0: { this->mvpOp->applyTranspose(V[j], w); } break;
      case 1: { this->mvpOp->applyTranspose(V[j], r); this->pcOp->applyT(r, w); } break;
      case 2: { this->pcOp->applyT(V[j], r); this->mvpOp->applyTranspose(r, w); } break;
      }



      for (int i=0; i<=j; ++i) {
        // For complex vectors, w has to not be conjugated in the definition
        // of the dot product
        H[i][j] = w * V[i];
        w -= H[i][j] * V[i];
      }

      H[j+1][j] = w.norm();

      if (H[j+1][j] == 0.0) {

        applyPreviousRotations(j, H, cs);
        applyNewRotation(j, H, cs, g);

        ++iter; exitLoop = 1; break;

      }

      V[j+1] = (1./H[j+1][j]) * w;

      applyPreviousRotations(j, H, cs);
      applyNewRotation(j, H, cs, g);

      l2res = sqrt(sqNorm(g[j+1]));

      ++iter;

      if (this->output)
	this->ioOp->fprintf(this->output, "  %d %e %e \n", iter, l2res, target);

      if (l2res <= target || iter >= this->maxits ||
          l2res <= this->absoluteEps) { exitLoop = 1; break; }

    }

    if (j == numVec) --j;

    backwardSolve(j, H, g, y);

    w = 0.0;

    for (int m=0; m<=j; ++m) w += y[m] * V[m];

    if (typePrec == 2) { this->pcOp->applyT(w, r); w = r; }

    x += w;

    numOutputVecs = j+1;

  } while (exitLoop == 0);


  if (this->checkFinalRes) {
    this->mvpOp->applyTranspose(x, w);
    r = b - w;
    if (typePrec == 1) {
      this->pcOp->applyT(r, w);
      r = w;
    }
    if (this->output)
      this->ioOp->fprintf(this->output, "  %d %e (actual residual norm)\n", iter, r.norm());
    if (r.norm() > target)  {
      if (this->output)
        this->ioOp->fprintf(this->output, "  %d %e (actual residual norm) vs l2res: %e, target:%e\n", iter, r.norm(), l2res,
target);
      iter = -999;
    }

  }

//  this->ioOp->printf(5, "Gmres(%d) solver: its=%d, res=%.2e, target=%.2e\n", numVec, iter, l2res, target);
//  if (iter == this->maxits && l2res > target && outputConvergenceInfo) {
//    this->ioOp->printf(1, "*** Warning: Gmres(%d) solver reached %d its", numVec, this->maxits);
//    this->ioOp->printf(1, " (initial=%.2e, res=%.2e, target=%.2e, ratio = %.2e)\n", res0, l2res, target, l2res/target);
//  }
//
//  if (this->kspBinaryOutput) {
//    if (typePrec == 2) {  //apply preconditioner before outputting
//      for (int iVec=0; iVec<numOutputVecs; ++iVec) {
//        this->pcOp->applyTranspose(V[iVec], r);
//        V[iVec] = r;
//      }
//    }
//    this->kspBinaryOutput->writeKrylovVectors(V, y, numOutputVecs);
//  }

  return iter;

}




//------------------------------------------------------------------------------
template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
int
GmresSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::solveLS(VecType &b, VecType &x)
{

  int typePrec = 2;
  double beta, l2res, target;

  int iter = 0;
  int exitLoop = 0;

  do {

    this->mvpOp->apply(x, w);
    r = b - w;

    if (typePrec == 1) { this->pcOp->apply(r, w); r = w; }

    beta = r.norm();

    if (iter == 0) {
      target = this->eps * b.norm(); // beta;
      if (this->output) this->ioOp->fprintf(this->output, "LS : Gmres(%d) iterations:\n", numVec);
      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e (%e)\n", 0, beta, target, this->eps);
    }
    else
      if (this->output) this->ioOp->fprintf(this->output, "  --- restart ---\n");

    if (beta == 0.0) return 0;

    V[0] = (1.0/beta) * r;
    g = 0.0;
    g[0] = beta;

    int j;
    for (j=0; j<numVec; ++j) {
      switch (typePrec) {
      case 0: { this->mvpOp->apply(V[j], w); } break;
      case 1: { this->mvpOp->apply(V[j], r); this->pcOp->apply(r, w); } break;
      case 2: { this->pcOp->apply(V[j], r); this->mvpOp->apply(r, w); } break;
      }

      for (int i=0; i<=j; ++i) {
        // For complex vectors, w has to not be conjugated in the definition
        // of the dot product
        H[i][j] = w * V[i];
        w -= H[i][j] * V[i];
      }

      H[j+1][j] = w.norm();

      if (H[j+1][j] == 0.0) {
        applyPreviousRotations(j, H, cs);
        applyNewRotation(j, H, cs, g);
        ++iter; exitLoop = 1; break;
      }

      V[j+1] = (1./H[j+1][j]) * w;

      applyPreviousRotations(j, H, cs);
      applyNewRotation(j, H, cs, g);
      l2res = sqrt(sqNorm(g[j+1]));
      ++iter;

      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", iter, l2res, target);
      if (l2res <= target || iter >= this->maxits) { exitLoop = 1; break; }
    }

    if (j == numVec) --j;
    backwardSolve(j, H, g, y);
    w = 0.0;

    for (int m=0; m<=j; ++m) w += y[m] * V[m];

    if (typePrec == 2) { this->pcOp->apply(w, r); w = r; }

    x += w;

  } while (exitLoop == 0);

  if (this->checkFinalRes) {
    this->mvpOp->apply(x, w);
    r = b - w;
    if (typePrec == 1) {
      this->pcOp->apply(r, w);
      r = w;
    }
    if (this->output)
      this->ioOp->fprintf(this->output, " LS  %d %e (actual residual norm)\n", iter, r.norm());
    if (r.norm() > target)  {
      if (this->output)
        this->ioOp->fprintf(this->output, " LS  %d %e (actual residual norm) vs l2res: %e, target:%e\n", iter, r.norm(), l2res,
target);
      iter = -999;
    }
  }

  this->ioOp->printf(5, "LS : Gmres(%d) solver: its=%d, res=%.2e, target=%.2e\n", numVec, iter, l2res, target);
  if (iter == this->maxits) {
    this->ioOp->printf(1, "*** Warning: LS Gmres(%d) solver reached %d its", numVec, this->maxits);
    this->ioOp->printf(1, " (res=%.2e, target=%.2e)\n", l2res, target);
  }

  return iter;
}                                                                                               

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
int
GmresSolver<VecType,MatVecProdOp,PrecOp,IoOp,ScalarT>::solveT(VecType &b, VecType &x)
{

  int typePrec = 2;//WTF

  double beta, l2res, target, res0;//TODO res0 is new
  res0 = b.norm();

  int iter = 0;
  int exitLoop = 0;

  int numOutputVecs = 0;//TODO new

  do {

    this->mvpOp->applyTranspose(x, w);

    r = b - w;

    if (typePrec == 1) { this->pcOp->applyT(r, w); r = w; }

    beta = r.norm();

    if (iter == 0) {
      //target = this->eps * b.norm(); // TODO BUGHUNT ori
      target = this->eps * beta;
      //target = this->eps * (res0 = b.norm());//TODO BUGHUNT modified
      if (this->output) this->ioOp->fprintf(this->output, "GmRES iterations:\n");
      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", 0, beta, target);
    }
    else
      if (this->output) this->ioOp->fprintf(this->output, "  --- restart ---\n");

    if (beta == 0.0) return 0;

    V[0] = (1.0/beta) * r;

    g = 0.0;
    g[0] = beta;
    int j;
    for (j=0; j<numVec; ++j) {

      switch (typePrec) {
      case 0: { this->mvpOp->applyTranspose(V[j], w); } break;
      case 1: { this->mvpOp->applyTranspose(V[j], r); this->pcOp->applyT(r, w); } break;
      case 2: {this->pcOp->applyT(V[j], r);  this->mvpOp->applyTranspose(r, w);
      } break;

      }

      for (int i=0; i<=j; ++i) {
        H[i][j] = w * V[i];
        w -= H[i][j] * V[i];
      }

      H[j+1][j] = w.norm();

      if (H[j+1][j] == 0.0) {

        applyPreviousRotations(j, H, cs);
        applyNewRotation(j, H, cs, g);

        ++iter; exitLoop = 1; break;

      }

      V[j+1] = (1./H[j+1][j]) * w;

      applyPreviousRotations(j, H, cs);
      applyNewRotation(j, H, cs, g);

      l2res = sqrt(sqNorm(g[j+1]));

      ++iter;

      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", iter, l2res, target);

      if (l2res <= target || iter >= this->maxits || l2res <= this->absoluteEps) { exitLoop = 1; break; }

    }

    if (j == numVec) --j;

    backwardSolve(j, H, g, y);

    w = 0.0;

    for (int m=0; m<=j; ++m) w += y[m] * V[m];

    if (typePrec == 2) { this->pcOp->applyT(w, r); w = r; }

    x += w;

    numOutputVecs = j+1;

  } while (exitLoop == 0);



  // TODO Yungsoo's version
  ///////////////////////////////////////////////////////////////////////////////////////////////
  if (this->checkFinalRes) {

    this->mvpOp->applyTranspose(x, w);
    r = b - w;

    if (typePrec == 1) { this->pcOp->applyT(r, w); r = w; }

    if (this->output)
      this->ioOp->fprintf(this->output, "  %d %e (actual residual norm)\n", iter, r.norm() );

  }

  return iter;
  ///////////////////////////////////////////////////////////////////////////////////////////////


  //TODO BUGHUNT copied from solve
  ////////////////////////////////////////////////////////////////////////////////////////////////
//  if (this->checkFinalRes) {
//      this->mvpOp->applyTranspose(x, w);
//      r = b - w;
//      if (typePrec == 1) {
//        this->pcOp->applyTranspose(r, w);
//        r = w;
//      }
//      if (this->output)
//        this->ioOp->fprintf(this->output, "  %d %e (actual residual norm)\n", iter, r.norm());
//      if (r.norm() > target)  {
//        if (this->output)
//          this->ioOp->fprintf(this->output, "  %d %e (actual residual norm) vs l2res: %e, target:%e\n", iter, r.norm(), l2res,
//  target);
//        iter = -999;
//      }
//
//    }
//
//    this->ioOp->printf(5, "Gmres(%d) solver: its=%d, res=%.2e, target=%.2e\n", numVec, iter, l2res, target);
//    if (iter == this->maxits && l2res > target && outputConvergenceInfo) {
//      this->ioOp->printf(1, "*** Warning: Gmres(%d) solver reached %d its", numVec, this->maxits);
//      this->ioOp->printf(1, " (initial=%.2e, res=%.2e, target=%.2e, ratio = %.2e)\n", res0, l2res, target, l2res/target);
//    }
//
//    if (this->kspBinaryOutput) {
//      if (typePrec == 2) {  //apply preconditioner before outputting
//        for (int iVec=0; iVec<numOutputVecs; ++iVec) {
//          this->pcOp->applyTranspose(V[iVec], r);
//          V[iVec] = r;
//        }
//      }
//      this->kspBinaryOutput->writeKrylovVectors(V, y, numOutputVecs);
//    }
//
//    return iter;
    ////////////////////////////////////////////////////////////////////////////////////////////////


}
//------------------------------------------------------------------------------
inline double myConj(double x)  { return x; }
inline std::complex<double>  myConj(std::complex<double> x)  { return std::conj(x); }
//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
void
GmresSolver<VecType,MatVecProdOp,PrecOp,IoOp,ScalarT>::
applyPreviousRotations(int j, GenFullM<ScalarT> &H, GenFullM<ScalarT> &cs)
{

  for (int i=0; i<j; ++i) {

    //ScalarT hj   = cs[i][1] * H[i][j] + cs[i][0] * H[i+1][j];
    //ScalarT hjp1 = -cs[i][0] * H[i][j] + cs[i][1] * H[i+1][j];

    ScalarT hj   = cs[i][1] * H[i][j] + myConj(cs[i][0]) * H[i+1][j];
    ScalarT hjp1 = -cs[i][0] * H[i][j] + cs[i][1] * H[i+1][j];

    H[i][j] = hj;
    H[i+1][j] = hjp1;

  }

}
//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
void
GmresSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::
applyNewRotation(int j, GenFullM<ScalarT> &H, GenFullM<ScalarT> &cs, Vec<ScalarT> &g)
{

  if (sqNorm(H[j][j]) == 0.0)  {
    cs[j][1] = 0.0;   // cos
    cs[j][0] = 1.0;   // sin
  }
  else  {
    double coef = 1.0 / sqrt( sqNorm(H[j][j]) + sqNorm(H[j+1][j]) );
    cs[j][1] = sqrt(sqNorm(H[j][j])) * coef;  //cos
    cs[j][0] = cs[j][1] * H[j+1][j]/H[j][j];  //sin
  }

  //H[j][j] = cs[j][1] * H[j][j] + cs[j][0] * H[j+1][j];
  H[j][j] = cs[j][1] * H[j][j] + myConj(cs[j][0]) * H[j+1][j];
  H[j+1][j] = 0.0;

  g[j+1] = -cs[j][0] * g[j];
  g[j] = cs[j][1] * g[j];

}

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
void
GmresSolver<VecType,MatVecProdOp,PrecOp,IoOp,ScalarT>::
backwardSolve(int m, GenFullM<ScalarT> &H, Vec<ScalarT> &g, Vec<ScalarT> &y)
{

  for (int i=m; i>=0; --i) {

    y[i] = g[i];

    for (int k=m; k>=i+1; --k) y[i] -= H[i][k] * y[k];

    y[i] /= H[i][i];

  }

}

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
int
GmresSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::solve(VecSet<VecType> &b, 
            VecSet<VecType> &x)
{

  int typePrec = 2;

  double beta, l2res, target;
  int nVec = b.numVectors();

  Vec<double> eps(nVec);
  int iVec;
  for(iVec = 0; iVec < nVec; ++iVec)
   eps[iVec] = this->kspConvCriterion->compute(b[iVec].norm());

  int iter = 0;
  
  // First compute the residual for each b. b is overwritten by it
  for(iVec = 0; iVec < nVec; ++iVec) {
    this->matVecOp->apply(x[iVec], w);
    b[iVec] -= w;
  }
  do {
    this->exitLoop = 1;
    for(iVec = 0; iVec < nVec; ++iVec) {
      r = b[iVec];

      if (typePrec == 1) { this->precOp->apply(r, w); r = w; }

      beta = sqrt(r*r);

      if (iter == 0) {
        target = eps[iVec] * beta;
        if (this->output) this->ioOp->fprintf(this->output, "Gmres iterations:\n");
        if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", 0, beta, target);
      } 
      else
        if (this->output) this->ioOp->fprintf(this->output, "  --- restart ---\n");

      if (beta == 0.0) return 0;

      V[0] = (1./beta) * r;

      g = 0.0;
      g[0] = beta;

      int j;
      for (j=0; j<numVec; ++j) {

        switch (typePrec) {
        case 0: { this->matVecOp->apply(V[j], w); } break;
        case 1: { this->matVecOp->apply(V[j], r); this->precOp->apply(r, w); } break;
        case 2: { this->precOp->apply(V[j], r); this->matVecOp->apply(r, w); } break;
        }
        this->ioOp->fprintf(stderr, "... applied PrecOp-> w = %e\n", w.norm());

        for (int i=0; i<=j; ++i) {
	  H[i][j] = w * V[i];
	  w -= H[i][j] * V[i];
        }

        H[j+1][j] = sqrt(w*w);
    
        if (H[j+1][j] == 0.0) {

	  applyPreviousRotations(j, H, cs);
	  applyNewRotation(j, H, cs, g);

	  ++iter; this->exitLoop = 1; break;

        }
  
        V[j+1] = (1./H[j+1][j]) * w;

        applyPreviousRotations(j, H, cs);
        applyNewRotation(j, H, cs, g);

        l2res = fabs(g[j+1]);

        ++iter;

        if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e\n", iter, l2res, target);

        if (l2res <= target || iter >= this->maxIter) { break; }

      }

      if (j == numVec) {
	 --j;
	 this->exitLoop = 0;
       }

      backwardSolve(j, H, g, y);

      w = 0.0;

      for (int m=0; m<=j; ++m) w += y[m] * V[m];

      if (typePrec == 2) { this->precOp->apply(w, r); w = r; }

      x[iVec] += w;
      
      //update the right hand side
      this->matVecOp->apply(w,r);
      b[iVec] -= r;     
      // now improve the solution of the other systems
      int kVec;
      for(kVec = 0; kVec < nVec; ++kVec) {
	if(kVec == iVec) continue;
	for(int k = 0; k <= j; ++k)
	  g[k] = V[k]*b[iVec];
	//rotate(j, g, cs);
	backwardSolve(j, H, g, y);
	w = 0.0;

        for (int m=0; m<=j; ++m) w += y[m] * V[m];
        if (typePrec == 2) { this->precOp->apply(w, r); w = r; }

        x[kVec] += w;
        //update the right hand side
        this->matVecOp->apply(w,r);
        b[kVec] -= r;     
	
      }
    }

  } while (this->exitLoop == 0);


  if (this->checkFinalRes) {

    for (int k = 0; k < nVec; ++k) {
      
      this->matVecOp->apply(x[k], w);
      r = b[k] - w;

      if (typePrec == 1) { this->precOp->apply(r, w); r = w; }

      if (this->output) 
        this->ioOp->fprintf(this->output, "  %d %e (actual residual norm)\n", iter, r.norm());
    }

  }

  //if (this->kspBinaryOutput) this->kspBinaryOutput->writeKrylovVectors(V, y);

  return iter;

}

//------------------------------------------------------------------------------
                                                        
template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
GcrSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::
GcrSolver(const typename VecType::InfoType &info, KspData &data,
            MatVecProdOp *A, PrecOp *P, IoOp *IO) :
  KspSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>(data, A, P, IO),
  numVec(data.numVectors), p(data.numVectors+1, info), Ap(data.numVectors+1, info), w(info), r(info), R(info), AR(info), temp(info), x0(info), w0(info)
{
                                                        
  ApAp = new ScalarT[data.numVectors];
                                                        
                                                        
  double sizeMB = (numVec + 1) * p[0].sizeMB();
                                                        
                                                        
                                                        
  this->ioOp->globalSum(1, &sizeMB);
                                                        
                                                        
                                                        
  this->ioOp->printf(2, "Memory for Gcr(%d) solver: %3.2f MB\n", numVec, 2*sizeMB);
                                                        
                                                        
                                                        
}
                                                        
                                                        
//------------------------------------------------------------------------------
template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
int
GcrSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::solve(VecType &b, VecType &x)
{
                                                        
  // No preconditioning (typePrec = 0) or Right Preconditioning (typePrec = 1)
                                                        
  int typePrec = 1;
                                                        
  //Modified Gram-Schmidt (MGS = 1)
  int MGS = 1;
                                                        
  double res, target;
                                                        
  int iter = 0;
  int exitLoop = 0;
  do {
                                                        
                                                        
                                                        
    this->mvpOp->apply(x, w);
    //this->ioOp->fprintf(this->output,"debug solve -0.5 : norm(w) = %e  \n",w.norm());
    r = b - w;
                                                        
    switch (typePrec) {
      case 0: { p[0] = r; } break;
      case 1: { p[0] = r; this->pcOp->apply(r, R); } break;
      }
                                                        
    switch (typePrec) {
      case 0: { this->mvpOp->apply(p[0], Ap[0]);} break;
      case 1: { this->mvpOp->apply(R, Ap[0]); } break;
      }
                                                        
                                 
    int j;
    //this->ioOp->fprintf(this->output,"debug solve 0.25  \n");
    //this->ioOp->fprintf(this->output,"debug solve 0.5 : norm(Ap[0]) = %e  \n",ap0);     
    ApAp[0] = Ap[0] * Ap[0];
    //this->ioOp->fprintf(this->output,"debug solve 1  \n");
    res = r.norm();
                                                        
    if (iter == 0) {
      target = this->eps * b.norm();
      if (this->output) this->ioOp->fprintf(this->output, "Gcr(%d) iterations:\n", numVec);
      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e (%e)\n", 0, res, target, this->eps);
    }
    else
      if (this->output) this->ioOp->fprintf(this->output, "  --- restart ---\n");
                                                        
    if (res == 0.0) return 0;
                                                        
                                                        
                                                        
                                                        
    for (j=0; j<numVec; ++j) {
      alpha = r * Ap[j];
      alpha *= 1./ApAp[j];
      if (typePrec == 1) {
        this->pcOp->apply(p[j], temp);
        x += alpha * temp;
      }
      else
        x += alpha * p[j];
                                                        
      r -= alpha * Ap[j];
                                                        
      switch (typePrec) {
      case 0: { this->mvpOp->apply(r, AR);} break;
      case 1: { this->pcOp->apply(r, R); this->mvpOp->apply(R, AR); } break;
      }
                                                        
      p[j+1] = r;
                                                        
      switch (typePrec) {
      case 0: { this->mvpOp->apply(r, Ap[j+1]); } break;      case 1: { Ap[j+1] = AR; } break;
      }
                                                        
      for (int i=0; i<=j; ++i) {
        if (MGS == 1)
          beta = Ap[j+1] * Ap[i];
        else
          beta = AR * Ap[i];
        beta *= (-1.)/ApAp[i];
        p[j+1] += beta * p[i];
        Ap[j+1] += beta * Ap[i];
      }
      ApAp[j+1] = Ap[j+1] * Ap[j+1];
      res = r.norm();
      ++iter;
      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e \n", iter, res, target);
      if (res <= target || iter >= this->maxits) { exitLoop = 1; break; }
                                                                                                                                                                          
    }
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
  } while (exitLoop == 0);
                                                        
                                                        
                                                        
  if (this->checkFinalRes) {
    this->mvpOp->apply(x, w);
    r = b - w;
    if (this->output)
      this->ioOp->fprintf(this->output, "  %d %e (actual residual norm)\n", iter, r.norm());
    iter = -999;
                                                        
                                                        
                                                        
  }
                                                        
  this->ioOp->printf(5, "Gcr(%d) solver: its=%d, res=%.2e, target=%.2e\n", numVec, iter, res, target);
  if (iter == this->maxits) {
    this->ioOp->printf(1, "*** Warning: Gcr(%d) solver reached %d its", numVec, this->maxits);
    this->ioOp->printf(1, " (res=%.2e, target=%.2e)\n", res, target);
  }
  return iter;
                                                        
                                                        
                                                        
}
                                                        
                                                        
//------------------------------------------------------------------------------
template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
int
GcrSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::solveMRhs(VecType &b, VecType &x)
{
                                                                                                                 
  // No preconditioning (typePrec=0) and Right Preconditioning (typePrec = 1) only
  int typePrec = 1;
                                                        
  // Modified Gram-Schmidt
  int MGS = 1;
  double res, target;
                                                        
  int iter = 0;
  int exitLoop = 0;
  int j;
                                                        
  do {
                                                        
                                                        
    this->mvpOp->apply(x, w);
    r = b - w;
                                                        
    w0 = 0.0;
    if (numCalcVec > 0) {
      y = new ScalarT[numCalcVec+1];
                                                        
      for (j=0; j < numCalcVec+1; j++) {
        y[j] = r * Ap[j];
        y[j] *= (1.0 / ApAp[j]);
      }
      temp = 0.0;
                                                        
      for (j=0; j<numCalcVec+1; j++) {
                                                        
        temp +=  (y[j] * p[j]);
      }
      if (typePrec == 1)
        this->pcOp->apply(temp, x0);
      else if (typePrec == 0)
        x0 = temp;
                                                        
      this->mvpOp->apply(x0, w0);
    }
                                                        
    if (numCalcVec > 0) {
      r = r - w0;
                                                        
    }
                                                        
                                                        
    if (numCalcVec == 0) {
      x0 = 0.0;
      switch (typePrec) {
        case 0: { p[0] = r; } break;
        case 1: { p[0] = r; this->pcOp->apply(r, R); } break;
      }
                                                        
                                                        
      switch (typePrec) {
        case 0: { this->mvpOp->apply(p[0], Ap[0]);} break;
        case 1: { this->mvpOp->apply(R, Ap[0]); } break;      }
                                                        
                                                        
      ApAp[0] = Ap[0] * Ap[0];
    }
    else if (numCalcVec > 0) {
                                                        
      switch (typePrec) {
        case 0: { this->mvpOp->apply(r, AR);} break;
        case 1: { this->pcOp->apply(r, R); this->mvpOp->apply(R, AR); } break;
      }
                                                        
                                                        
      p[numCalcVec+1] = r;
                                                        
                                                        
      switch (typePrec) {
        case 0: { this->mvpOp->apply(r, Ap[numCalcVec+1]); } break;
        case 1: { Ap[numCalcVec+1] = AR; } break;
      }
                                                        
                                                        
      for (int i=0; i<=numCalcVec; ++i) {
                                                        
        if (MGS == 1)
          beta = Ap[numCalcVec+1] * Ap[i];
        else
          beta = AR * Ap[i];
        beta *= (-1.)/ApAp[i];
        p[numCalcVec+1] += beta * p[i];
        Ap[numCalcVec+1] += beta * Ap[i];
      }
      ApAp[numCalcVec+1] = Ap[numCalcVec+1] * Ap[numCalcVec+1];
      numCalcVec++;
    }
                                                        
                                                        
    res = r.norm();
                                                                                                                 
    if (iter == 0) {
      target = this->eps * b.norm(); // beta;
      if (this->output) this->ioOp->fprintf(this->output, "Gcr(%d) iterations:\n", numVec);
      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e (%e)\n", 0, res, target, this->eps);
    }
    else
      if (this->output) this->ioOp->fprintf(this->output, "  --- restart ---\n");
                                                                                                                 
    if (res == 0.0) return 0;
                                                                                                                 
    int j;
    for (j=numCalcVec; j<numVec; ++j) {
                                                                                                                 
      alpha = r * Ap[j];
      alpha *= 1./ApAp[j];
      if (typePrec == 1) {
        this->pcOp->apply(p[j], temp);
        x += alpha * temp;
      }
      else if ( typePrec == 0)
        x += alpha * p[j];
                                                        
                                                        
      r -= alpha * Ap[j];
      switch (typePrec) {
      case 0: { this->mvpOp->apply(r, AR);} break;
      case 1: { this->pcOp->apply(r, R); this->mvpOp->apply(R, AR); } break;
      }
                                                        
                                                        
      p[j+1] = r;
                                                                                                                 
      switch (typePrec) {
      case 0: { this->mvpOp->apply(r, Ap[j+1]); } break;      case 1: { Ap[j+1] = AR; } break;
      }
                                                                                                                 
      for (int i=0; i<=j; ++i) {
                                                        
        if (MGS == 1)
          beta = Ap[j+1] * Ap[i];
        else
          beta = AR * Ap[i];
        beta *= (-1.)/ApAp[i];
        p[j+1] += beta * p[i];
        Ap[j+1] += beta * Ap[i];
      }
      ApAp[j+1] = Ap[j+1] * Ap[j+1];
      res = r.norm();
      ++iter;
      ++numCalcVec;
      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e \n", iter, res, target);
      if (res <= target || iter >= this->maxits) { exitLoop = 1; break; }
                                                        
                                                        
                                                        
                                                        
    }
                                                                                                                 
  } while (exitLoop == 0);
                                                        
                                                        
  x = x + x0;
                                                        
                                                        
  if (this->checkFinalRes) {
    this->mvpOp->apply(x, w);
    r = b - w;
    if (this->output)
      this->ioOp->fprintf(this->output, "  %d %e (actual residual norm)\n", iter, r.norm());
    iter = -999;
                                                        
                                                        
                                                        
                                                        
  }

  return iter;
}
//------------------------------------------------------------------------------
/*template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
int
GcrSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT>::solveMRhs(VecType &b, VecType &x)
{
                                                        
  //this->ioOp->fprintf(this->output,"debug 0  \n"); 
  // No preconditioning (typePrec=0) and Right Preconditioning (typePrec = 1) only
  int typePrec = 1;
                                                        
  // Modified Gram-Schmidt
  int MGS = 1;
  double res, target;
                                                        
  //double bn, rn;                                                                                                                  
  int iter = 0;
  int exitLoop = 0;
  int j;
  //this->ioOp->fprintf(this->output,"debug 1 \n");      
  do {
                                                        
    //bn = b.norm();
    //this->ioOp->fprintf(this->output,"b norm = %e \n",bn);
                                                        
    this->mvpOp->apply(x, w);
    r = b - w;
    //this->ioOp->fprintf(this->output,"debug 2  \n");              
    //rn = r.norm();
    //if (rn != bn)
    //  this->ioOp->fprintf(this->output,"r norm = %e \n",rn);
    //fprintf(stderr,"numCalcVec = %d \n",numCalcVec);
    w0 = 0.0;
    if (numCalcVec > 0) {
      y = new ScalarT[numCalcVec+1];
                                                        
      for (j=0; j < numCalcVec+1; j++) {
        y[j] = r * Ap[j];
        y[j] *= (1.0 / ApAp[j]);
        //imag(y[j]) = - imag(y[j]);
        //fprintf(stderr,"y[%d]  = %e + %e i \n",j,real(y[j]),imag(y[j]));
      }
      temp = 0.0;
                                                        
      for (j=0; j<numCalcVec+1; j++) {
        //fprintf(stderr,"norm p[%d] = %e \n",j,(p[j]).norm());
                                                        
        temp +=  (y[j] * p[j]);
      }
      if (typePrec == 1)
        this->pcOp->apply(temp, x0);
      else if (typePrec == 0)
        x0 = temp;
      //fprintf(stderr,"x0 norm = %e \n",x0.norm());
                                                        
      this->mvpOp->apply(x0, w0);
      //this->ioOp->fprintf(this->output,"w0 norm = %e \n",w0.norm());
    }
                                                        
    if (numCalcVec > 0) {
      r = r - w0;
      //rn = r.norm();
      //this->ioOp->fprintf(this->output,"r norm = %e \n",rn);
                                                        
    }
                                                        
    this->ioOp->fprintf(this->output,"debug 3  \n");                                                    
    if (numCalcVec == 0) {
      x0 = 0.0;
      switch (typePrec) {
        case 0: { p[0] = r; } break;
        case 1: { p[0] = r; this->pcOp->apply(r, R); } break;
      }
                                                        
                                                        
      switch (typePrec) {
        case 0: { this->mvpOp->apply(p[0], Ap[0]);} break;
        case 1: { this->mvpOp->apply(R, Ap[0]); } break;  
      }
      //this->mvpOp->apply(R, Ap[0]);//////////////////
      //this->ioOp->fprintf(this->output,"debug 3.5  \n");
      //this->ioOp->fprintf(this->output,"debug 4, norm(Ap[0] = %e  \n",Ap[0].norm());
   
      ApAp[0] = Ap[0] * Ap[0];
      //ApAp[0] = p[0] * p[0];
      this->ioOp->fprintf(this->output,"debug 4.5  \n");
    }
    else if (numCalcVec > 0) {
                                                        
      switch (typePrec) {
        case 0: { this->mvpOp->apply(r, AR);} break;
        case 1: { this->pcOp->apply(r, R); this->mvpOp->apply(R, AR); } break;
      }
                                                        
      this->ioOp->fprintf(this->output,"debug 5  \n");                                                        
      p[numCalcVec+1] = r;
                                                        
                                                        
      switch (typePrec) {
        case 0: { this->mvpOp->apply(r, Ap[numCalcVec+1]); } break;
        case 1: { Ap[numCalcVec+1] = AR; } break;
      }
                                                        
      this->ioOp->fprintf(this->output,"debug 6  \n");                                                  
      for (int i=0; i<=numCalcVec; ++i) {
                                                        
        if (MGS == 1)
          beta = Ap[numCalcVec+1] * Ap[i];
        else
          beta = AR * Ap[i];
        beta *= (-1.)/ApAp[i];
        p[numCalcVec+1] += beta * p[i];
        Ap[numCalcVec+1] += beta * Ap[i];
        Ap[numCalcVec+1] += beta * Ap[i];
      }
      ApAp[numCalcVec+1] = Ap[numCalcVec+1] * Ap[numCalcVec+1];
      numCalcVec++;
    }
                                                        
    this->ioOp->fprintf(this->output,"debug 7  \n");                                                    
    res = r.norm();
                                                                                                                 
    if (iter == 0) {
      target = this->eps * b.norm(); // beta;
      if (this->output) this->ioOp->fprintf(this->output, "Gcr(%d) iterations:\n", numVec);
      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e (%e)\n", 0, res, target, this->eps);
    }
    else
      if (this->output) this->ioOp->fprintf(this->output, "  --- restart ---\n");
                                                                                                                 
    if (res == 0.0) return 0;
                                                                                                                 
    int j;
    ///this->ioOp->fprintf(this->output,"numCalcVec = %d \n", numCalcVec);
    for (j=numCalcVec; j<numVec; ++j) {
                                                                                                                 
      alpha = r * Ap[j];
      alpha *= 1./ApAp[j];
      if (typePrec == 1) {
        this->pcOp->apply(p[j], temp);
        x += alpha * temp;
      }
      else if ( typePrec == 0)
        x += alpha * p[j];
      //fprintf(stderr,"alpha  = %e + %e i \n",real(alpha),imag(alpha));
                                                        
                                                        
      r -= alpha * Ap[j];
      switch (typePrec) {
      case 0: { this->mvpOp->apply(r, AR);} break;
      case 1: { this->pcOp->apply(r, R); this->mvpOp->apply(R, AR); } break;
      }
                                                        
                                                        
      p[j+1] = r;
                                                                                                                 
      switch (typePrec) {
      case 0: { this->mvpOp->apply(r, Ap[j+1]); } break;      case 1: { Ap[j+1] = AR; } break;
      }
                                                                                                                 
      for (int i=0; i<=j; ++i) {
                                                        
        if (MGS == 1)
          beta = Ap[j+1] * Ap[i];
        else
          beta = AR * Ap[i];
        beta *= (-1.)/ApAp[i];
        p[j+1] += beta * p[i];
        Ap[j+1] += beta * Ap[i];
      }
      ApAp[j+1] = Ap[j+1] * Ap[j+1];
      res = r.norm();
      ++iter;
      ++numCalcVec;
      if (this->output) this->ioOp->fprintf(this->output, "  %d %e %e \n", iter, res, target);
      if (res <= target || iter >= this->maxits) { exitLoop = 1; break; }
                                                        
                                                        
                                                        
                                                        
    }
                                                                                                                 
  } while (exitLoop == 0);
                                                        
                                                        
  x = x + x0;
                                                        
                                                        
  if (this->checkFinalRes) {
    this->mvpOp->apply(x, w);
    r = b - w;
    if (this->output)
      this->ioOp->fprintf(this->output, "  %d %e (actual residual norm)\n", iter, r.norm());
    iter = -999;
                                                        
                                                        
                                                        
                                                        
  }
                                                        
                                                        
  this->ioOp->printf(5, "Gcr(%d) solver: its=%d, res=%.2e, target=%.2e\n", numVec, iter, res, target);
  if (iter == this->maxits) {
    this->ioOp->printf(1, "*** Warning: Gcr(%d) solver reached %d its", numVec, this->maxits);
    this->ioOp->printf(1, " (res=%.2e, target=%.2e)\n", res, target);
  }
  return iter;
                                                        
                                                        
                                                        
                                                        
}
                                                        
                                                
*/                                                        

