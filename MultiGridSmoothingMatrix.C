#include <cstdio>

#include <MultiGridSmoothingMatrix.h>
#include <iostream>
#include <DenseMatrixOps.h>

#include <MultiGridLevel.h>

#include <Edge.h>

template <class Scalar, int dim>
MultiGridSmoothingMatrix<Scalar,dim>::
MultiGridSmoothingMatrix(SmoothingMode m, int isub,
                         int nn, int ne, int nBC,
                         MultiGridLevel<Scalar>* mglvl,
                         MultiGridLevel<Scalar>* mglvlR
                          ) : a(nn + 2*(ne+nBC)),
                                                          iSub(isub), indices(0),
                                                          mySmoothingMode(m),
                                                          mgLevel(mglvl),
                                                          mgLevelRefined(mglvlR) {

  n = nn; numEdges = ne; 

  lineJacobiMatrices = NULL; 
 
  switch (mySmoothingMode) {

   case RAS:

    iluA = mgLevelRefined->template createMaskILU<dim>(iSub,0,PcData::RCM, NULL);
    iluJW = new Vec<int>(nn);
    break;

   case LineJacobi:
   case BlockJacobi:

    indices = new SVec<int,dim>(nn);
    break;
   default:
    // Should never be reached in normal operation
    std::cout << "Error: Unsupported smoother in MultiGridSmoothingMatrix()" << std::endl;
    break;
  }
}
 
template <class Scalar, int dim>
MultiGridSmoothingMatrix<Scalar,dim>::
~MultiGridSmoothingMatrix() {

  switch (mySmoothingMode) {

   case RAS:

    delete iluA;
    delete iluJW;
    break;

   case LineJacobi:
   case BlockJacobi:

    delete indices;
    break;
   default:
    // Should never be reached in normal operation
    std::cout << "Error: Unsupported smoother in MultiGridSmoothingMatrix()" << std::endl;
    break;
  }

  if (lineJacobiMatrices) {

    for (int i = 0; i < numLines; ++i)
      delete lineJacobiMatrices[i];
    delete [] lineJacobiMatrices;
  }  

}
 

template <class Scalar, int dim>
void MultiGridSmoothingMatrix<Scalar,dim>::
getData(GenMat<Scalar,dim>& mat) {

  // A coarse level exists
  if (mySmoothingMode == LineJacobi) {

    if (mgLevel) {
      getDataLineJacobi(mat);
    } else {

    }
  } else if (mySmoothingMode == RAS) {

    iluA->getData(mat);
    iluA->numericILU(iluJW->data()); 
  } else {
    
    getDataBlockJacobi(mat);
  }
}

template <class Scalar, int dim>
void MultiGridSmoothingMatrix<Scalar,dim>::
smooth(SVec<Scalar,dim>& r, SVec<Scalar,dim>& du) {

  switch (mySmoothingMode) {

    case BlockJacobi:
      smoothBlockJacobi(r,du);
      break;
    case LineJacobi:
      smoothLineJacobi(r,du);
      break;
    case RAS:
      smoothRAS(r,du);
      break;
    default:
      break;
  }
  
}

template <class Scalar, int dim>
void MultiGridSmoothingMatrix<Scalar,dim>::
getDataBlockJacobi(GenMat<Scalar,dim>& mat) {

  // Load in the diagonal elements, and compute their LU decompositions
  for (int i = 0; i < n; ++i) {
 
    memcpy(a[i], mat.getElem_ii(i), sizeof(Scalar)*dim*dim);
    DenseMatrixOp<Scalar,dim,20>::ludec(a[i],(*indices)[i],1.0,dim);
    //memcpy(a[i+n], mat.getElem_ii(i), sizeof(Scalar)*dim*dim);
  }
/* 
  for (int i = 0; i < numEdges; ++i) {
 
    memcpy(a[2*i+n],mat.getElem_ij(i),sizeof(Scalar)*dim*dim);
    memcpy(a[2*i+n+1],mat.getElem_ji(i),sizeof(Scalar)*dim*dim);
  }
*/
}

template <class Scalar, int dim>
void MultiGridSmoothingMatrix<Scalar,dim>::
getDataLineJacobi(GenMat<Scalar,dim>& mat) {

  DistInfo& nodeDistInfo = mgLevel->getNodeDistInfo();

  // Load in the diagonal elements, and compute their LU decompositions
  bool* masterFlag = mgLevelRefined->getNodeDistInfo().getMasterFlag(iSub);
  int nl = 0,lineid,edgei,edgej,loci,locj;
  numLines = mgLevel->NumLines(iSub);

  if (!lineJacobiMatrices) {
    lineJacobiMatrices = new BlockTridiagonalMatrix<Scalar,dim>*[mgLevel->NumLines(iSub)]; 
    
    for (int i = 0; i < mgLevel->NumLines(iSub); ++i) {
   
      lineJacobiMatrices[i] = new BlockTridiagonalMatrix<Scalar,dim>(mgLevel->lineLength(iSub,i));
    }
  }

  for (int i = 0; i < n; ++i) {

    // If the master flag of this node is not set, then the jacobian (and
    // block inverse) will be computed by another subdomain.  
    if (!masterFlag[i]) 
      continue;
 
    if (mgLevel->isLine(iSub,i,i,&lineid,&loci,&locj)) {

      Scalar (*al)[dim*dim] = lineJacobiMatrices[lineid]->get(loci,locj);
      memcpy(al, mat.getElem_ii(i), sizeof(Scalar)*dim*dim);
    } else {
    
      memcpy(a[i], mat.getElem_ii(i), sizeof(Scalar)*dim*dim);
      DenseMatrixOp<Scalar,dim,20>::ludec(a[i],(*indices)[i],1.0,dim);
    }
  }

  EdgeSet* es = mgLevelRefined->getEdges()[iSub];
  for (int i = 0; i < numEdges; ++i) {
 
    int edgei = es->getPtr()[i][0], edgej = es->getPtr()[i][1];
    if (mgLevel->isLine(iSub,edgei, edgej,&lineid,&loci,&locj)) {
      
      Scalar (*al)[dim*dim] = lineJacobiMatrices[lineid]->get(loci,locj);
      memcpy(al, mat.getElem_ij(i),sizeof(Scalar)*dim*dim);
      
      al = lineJacobiMatrices[lineid]->get(locj,loci);
      memcpy(al, mat.getElem_ji(i),sizeof(Scalar)*dim*dim);
    }
  }
 
  for (int i = 0; i < numLines; ++i) {
   
    lineJacobiMatrices[i]->computeLU();
  }
}

template <class Scalar, int dim>
void MultiGridSmoothingMatrix<Scalar,dim>::
smoothBlockJacobi(SVec<Scalar,dim>& r, SVec<Scalar,dim>& du) {

  du = r;
  //double tmp[dim],dtmp[dim];
  bool* masterFlag = mgLevelRefined->getNodeDistInfo().getMasterFlag(iSub);
  for (int i = 0; i < n; ++i) {

    //memcpy(tmp,du[i],sizeof(Scalar)*dim);
    if (!masterFlag[i]) {

      memset(du[i],0,sizeof(Scalar)*dim);
    }  
    else
      DenseMatrixOp<Scalar,dim,20>::ludfdbksb(a[i],(*indices)[i], du[i], dim);
    /*for (int j = 0 ; j < dim; ++j) {

      dtmp[j] = 0.0;
      for (int k = 0; k < dim; ++k)
        dtmp[j] += a[i+n][j*dim+k]*du[i][k];
      if (fabs(dtmp[j]-tmp[j]) > 1.0e-6) {
        std::cout << "Error: " << std::endl;
        for (int l = 0; l < dim*dim; ++l)
          std::cout << a[i+n][l] << " ";
      }
    }*/
  }
}

template <class Scalar, int dim>
void MultiGridSmoothingMatrix<Scalar,dim>::
smoothLineJacobi(SVec<Scalar,dim>& r, SVec<Scalar,dim>& du) {

  du = 0.0;
  int lineid,loci,locj; 
  bool* masterFlag = mgLevelRefined->getNodeDistInfo().getMasterFlag(iSub);
  for (int i = 0; i < n; ++i) { 
    if (!masterFlag[i])
      continue;
    if (mgLevel->isLine(iSub,i,i,&lineid,&loci,&locj)) {
      continue;
    } 
    memcpy(du[i],r[i],sizeof(Scalar)*dim);
    DenseMatrixOp<Scalar,dim,20>::ludfdbksb(a[i],(*indices)[i], du[i], dim);
  }

  Scalar b[8][dim],x[8][dim];
  for (int i = 0; i < numLines; ++i) {   
 
    int* line_ptr = mgLevel->getLineData(iSub, i);
    int k = 0;
    while (line_ptr[k] >= 0 && k < 8) {
      for (int j = 0; j < dim; ++j) {
        b[k][j] = r[ line_ptr[k] ][j];
        //x[k][j] = du[ line_ptr[k] ][j];
      }
      ++k;
    }
 
    lineJacobiMatrices[i]->solveLU(b,x);  
    k = 0;
    while (line_ptr[k] >= 0 && k < 8) {
      for (int j = 0; j < dim; ++j) {
        du[line_ptr[k]][j] = x[ k ][j];
      }
      ++k;
    }
  }

}

template <class Scalar, int dim>
void MultiGridSmoothingMatrix<Scalar,dim>::
smoothRAS(SVec<Scalar,dim>& r, SVec<Scalar,dim>& du) {

  SVec<Scalar,dim> tmp(r);
  iluA->lusol(tmp,du);
  bool* masterFlag = mgLevelRefined->getNodeDistInfo().getMasterFlag(iSub);
  for (int i = 0; i < n; ++i) { 
    if (masterFlag[i])
      continue;
    memset(du[i],0,sizeof(Scalar)*dim);
  }
}

#define INST_HELPER(type,dim) \
template class MultiGridSmoothingMatrix<type,dim>;

INST_HELPER(double,1);
INST_HELPER(double,2);
INST_HELPER(double,3);
INST_HELPER(double,5);
INST_HELPER(double,6);
INST_HELPER(double,7);

