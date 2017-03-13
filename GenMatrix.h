#ifndef _GEN_MATRIX_H_
#define _GEN_MATRIX_H_

//------------------------------------------------------------------------------

template<class Scalar, int dim>
class GenMat {

public:

  GenMat() {}
  ~GenMat() {}

  virtual GenMat<Scalar,dim> &operator= (const Scalar) = 0;
  virtual GenMat<Scalar,dim> &operator*= (const Scalar) = 0;

  virtual Scalar (*data())[dim*dim] = 0;

  virtual double norm() { fprintf(stderr, "No Implementation of GenMat::norm\n"); return 0; }

  virtual Scalar *getElem_ii(int) = 0;
  virtual Scalar *getElem_ij(int) = 0;
  virtual Scalar *getElem_ji(int) = 0;

  virtual Scalar *getBcElem_ij(int l) {fprintf(stderr, "No Implementation\n"); return 0;}
  virtual Scalar *getBcElem_ji(int l) {fprintf(stderr, "No Implementation\n"); return 0;}
  virtual void addContrib(int, int *, double *) = 0;

  virtual void enableHHTerms(int) { }

  // -------------------------------------------------------------------------
  // Auxilliary terms (for ghost points)
 
  // Return the edge data corresponding to real node i and ghost node j
  virtual Scalar* getRealNodeElem_ij(int i,int j) { return NULL; }

  // Return the edge data correponding to ghost node i and real node j
  virtual Scalar* getGhostNodeElem_ij(int i,int j) { return NULL; }
   
  virtual Scalar* getGhostGhostElem_ij(int i,int j) { return NULL; }
 

  // Return the edge data corresponding to real node i and ghost node j
  virtual Scalar* queryRealNodeElem_ij(int i,int j) { return NULL; }

  // Return the edge data correponding to ghost node i and real node j
  virtual Scalar* queryGhostNodeElem_ij(int i,int j) { return NULL; }
   
  virtual Scalar* queryGhostGhostElem_ij(int i,int j) { return NULL; }

  struct AuxilliaryIterator {
    int row,col;
    Scalar* pData;
  };

  virtual AuxilliaryIterator* begin_realNodes() { return NULL; }
  virtual AuxilliaryIterator* begin_ghostNodes() { return NULL; } 
  virtual AuxilliaryIterator* begin_ghostGhostNodes() { return NULL; } 
  // Returns false if there is no next, true otherwise

  virtual bool next(AuxilliaryIterator*) { return false; }
  virtual void free(AuxilliaryIterator*) { }

  virtual void clearGhost() { }

  // -------------------------------------------------------------------------
  // Auxilliary terms (for HH Boundary Conditions)
 
  virtual Scalar* getElemUH(int) { return NULL; }
  virtual Scalar* getElemHU(int) { return NULL; }
  virtual Scalar* getElemHH(int) { return NULL; }

  virtual SVec<Scalar,dim*3>* getHU() { return NULL; }
  
};

//------------------------------------------------------------------------------

#endif
