#ifndef _MVP_MATRIX_H_
#define _MVP_MATRIX_H_

#include <Vector.h>
#include <GenMatrix.h>
#include <cstdio>

//------------------------------------------------------------------------------

template<class Scalar, int dim>
class MvpMat : public GenMat<Scalar,dim> {

protected:

  int n;
  int numEdges;

  SVec<Scalar,dim*dim> a;

  SVec<Scalar, dim*3>* uh, *hu;
  Vec<Scalar>* hh;

  typedef std::map< std::pair<int,int>, Scalar* > AuxilliaryRows;
  AuxilliaryRows realAuxilliaryRows;
  AuxilliaryRows ghostAuxilliaryRows;
  AuxilliaryRows ghostGhostAuxilliaryRows;

public:

 MvpMat(int nn, int ne, int nBC = 0) : a(nn + 2*(ne+nBC)),uh(0), hu(0), hh(0) { n = nn; numEdges = ne;}
  ~MvpMat() {

    if (uh)
      delete uh;
    if (hu)
      delete hu;
    if (hh)
      delete hh;
  }

  MvpMat<Scalar,dim> &operator= (const Scalar x) { 
    a = x; 
    if (uh)
      *uh = x;
    if (hh)
      *hh = x;
    if (hu)
      *hu = x;
    return *this; 
  }
  MvpMat<Scalar,dim> &operator*= (const Scalar x)  { 
    a *= x; 
    if (uh)
      *uh *= x;
    if (hh)
      *hh *= x;
    if (hu)
      *hu *= x;
    return *this; 
  }

  void enableHHTerms(int subLen) {

    uh = new SVec<Scalar,dim*3>(subLen);
    hu = new SVec<Scalar,dim*3>(subLen);
    hh = new Vec<Scalar>(subLen);
  }

  double norm() {return a.norm();}

  int numNonZeroBlocks() const { return a.size(); }
  Scalar (*data())[dim*dim] { return a.data(); }

  Scalar *getElem_ii(int i) { return *(a.v + i); }
  Scalar *getElem_ij(int l) { return *(a.v + n + 2 * l); }
  Scalar *getElem_ji(int l) { return *(a.v + n + 2 * l + 1); }

  Scalar *getBcElem_ij(int l) { return *(a.v + n + 2 * numEdges + 2 * l); }
  Scalar *getBcElem_ji(int l) { return *(a.v + n + 2 * numEdges + 2 * l + 1); }
  void addContrib(int nnd, int *nd, double *K) {}

  // -------------------------------------------------------------------------
  // Auxilliary terms (for ghost points)
 
  // Return the edge data corresponding to real node i and ghost node j
  Scalar* getRealNodeElem_ij(int i,int j) {
    return getAuxilliaryRow(realAuxilliaryRows,i,j);
  }

  // Return the edge data correponding to ghost node i and real node j
  Scalar* getGhostNodeElem_ij(int i,int j) {
    return getAuxilliaryRow(ghostAuxilliaryRows,i,j);
  }
  
  // Return the edge data correponding to ghost node i and ghost node j
  Scalar* getGhostGhostElem_ij(int i,int j) {
    return getAuxilliaryRow(ghostGhostAuxilliaryRows,i,j);
  }
 
  // Return the edge data corresponding to real node i and ghost node j
  Scalar* queryRealNodeElem_ij(int i,int j) {
    return queryAuxilliaryRow(realAuxilliaryRows,i,j);
  }

  // Return the edge data correponding to ghost node i and real node j
  Scalar* queryGhostNodeElem_ij(int i,int j) {
    return queryAuxilliaryRow(ghostAuxilliaryRows,i,j);
  }
  
  // Return the edge data correponding to ghost node i and ghost node j
  Scalar* queryGhostGhostElem_ij(int i,int j) {
    return queryAuxilliaryRow(ghostGhostAuxilliaryRows,i,j);
  }

  struct MvpAuxilliaryIterator : public GenMat<Scalar,dim>::AuxilliaryIterator {
    typename AuxilliaryRows::iterator itr;
    AuxilliaryRows* map_ptr;
  };

  typename GenMat<Scalar,dim>::AuxilliaryIterator* begin_realNodes() {
    if (realAuxilliaryRows.size() == 0)
      return NULL;

    MvpAuxilliaryIterator* mvpItr = new MvpAuxilliaryIterator;
    mvpItr->itr = realAuxilliaryRows.begin();
    mvpItr->map_ptr = &realAuxilliaryRows;
    updateIterator(mvpItr);
    return mvpItr;
  }

  typename GenMat<Scalar,dim>::AuxilliaryIterator* begin_ghostNodes() {
    if (ghostAuxilliaryRows.size() == 0)
      return NULL;

    MvpAuxilliaryIterator* mvpItr = new MvpAuxilliaryIterator;
    mvpItr->itr = ghostAuxilliaryRows.begin();
    mvpItr->map_ptr = &ghostAuxilliaryRows;
    updateIterator(mvpItr); 
    return mvpItr;
  }
  
  typename GenMat<Scalar,dim>::AuxilliaryIterator* begin_ghostGhostNodes() {
    if (ghostGhostAuxilliaryRows.size() == 0)
      return NULL;

    MvpAuxilliaryIterator* mvpItr = new MvpAuxilliaryIterator;
    mvpItr->itr = ghostGhostAuxilliaryRows.begin();
    mvpItr->map_ptr = &ghostGhostAuxilliaryRows;
    updateIterator(mvpItr); 
    return mvpItr;
  }
  
  bool next(typename GenMat<Scalar,dim>::AuxilliaryIterator* genItr) { 
    MvpAuxilliaryIterator* mvpItr = static_cast< MvpAuxilliaryIterator*>(genItr);
    ++mvpItr->itr;
    if (mvpItr->itr == mvpItr->map_ptr->end())
      return false;
    else {
      updateIterator(mvpItr);
      return true;
    }
  }

  void free(typename GenMat<Scalar,dim>::AuxilliaryIterator* genItr) { 
    delete genItr;
  }
  
  void clearGhost() { 
    
    clearAuxilliary(realAuxilliaryRows);
    clearAuxilliary(ghostAuxilliaryRows);
    clearAuxilliary(ghostGhostAuxilliaryRows);
  }

  Scalar* getElemUH(int i) {
    if (!uh)
      return NULL;
    
    return *(uh->v + i);
  }
    
  Scalar* getElemHU(int i) {
    if (!hu)
      return NULL;
    
    return *(hu->v + i);
  }

  Scalar* getElemHH(int i) {

    if (!hh)
      return NULL;
    return (hh->v + i);
  }

  SVec<Scalar,dim*3>* getHU() {

    return hu;
  }

 protected:

  Scalar* getAuxilliaryRow(AuxilliaryRows& A, int i, int j) {
    std::pair<int,int> ij(i,j);
    typename AuxilliaryRows::iterator itr = A.find( ij );
    if (itr == A.end()) {
      Scalar* s = new Scalar[dim*dim];
      memset(s,0,sizeof(Scalar)*dim*dim);
      A[ij] = s;
      return s;
    }
    else
      return itr->second;
  }

  Scalar* queryAuxilliaryRow(AuxilliaryRows& A, int i, int j) {
    std::pair<int,int> ij(i,j);
    typename AuxilliaryRows::iterator itr = A.find( ij );
    if (itr == A.end()) {
      return NULL;
    }
    else
      return itr->second;
  }

  void updateIterator(MvpAuxilliaryIterator* mvpItr) {

    mvpItr->row = mvpItr->itr->first.first;
    mvpItr->col = mvpItr->itr->first.second;
    mvpItr->pData = mvpItr->itr->second;
  }

  void clearAuxilliary(AuxilliaryRows& A) {
    
    typename AuxilliaryRows::iterator itr = A.begin();
    for (; itr != A.end(); ++itr) {
      delete [] itr->second;
    }
    
    A.clear();
  }
};

//------------------------------------------------------------------------------

#endif
