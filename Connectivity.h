#ifndef _CONNECTIVITY_H_
#define _CONNECTIVITY_H_

#include <cstdio>

#ifdef OLD_STL
#include <map.h>
#else
#include <map>
using std::map;
#endif

#ifdef _FEM_CODE_
class EqNumberer;
#endif
class ElemSet;
class FaceSet;
class BinFileHandler;

// component data structure
struct compStruct {
  int numComp; // number of components
  int *xcomp;  // pointer to renum for the beginning of each component 
  int *order;  // order of the nodes -> order[new] = old
  int *renum;  // renumbering -> renum[old] = new

  compStruct() { xcomp = 0; order = 0; renum = 0; }
  ~compStruct() 
  { 

    if (xcomp) delete [] xcomp;
    if (order) delete [] order;
    if (renum) delete [] renum;
  }

};


class Connectivity {
  int size;           // size of pointer
  int numtarget;      // size of target, number of Connections
  int *pointer;       // pointer to target
  int *target;        // value of the connectivity
  float *weight;      // weights of pointer (or NULL)

protected:

  Connectivity() { pointer = 0; target = 0; weight = 0; }

public:

  Connectivity(ElemSet *);
  Connectivity(FaceSet *fels); //HB
  Connectivity(int _size, int *_pointer, int *_target);
  Connectivity(int _size, int *_count);
  Connectivity(BinFileHandler &);
  
  ~Connectivity();

  int *operator[](int i);
  const int *operator[](int i) const;
  int csize();
  int numConnect(); // Total number of connections
  int num(int) const;
  int num(int, int*);
  int offset(int i) { return pointer[i]; } // Begining of ith part
  int offset(int i,int j); // returns a unique id for connection i to j
  Connectivity* reverse(float *w = 0); // creates t->s from s->t
  Connectivity* transcon( Connectivity* );
  Connectivity* transconOne( Connectivity*);

  void findPseudoDiam(int *n1, int *n2, int *mask=0);
  int  rootLS(int root, int *xls, int *ls, int &w, int *mask=0);

  // Create a rooted level structure
  int *renumSloan(int *mask, int &firstNum, int *ren = 0);
  int *renumRCM(int *mask, int &firstNum, int *ren = 0);
  int *renumSloan();
#ifdef _FEM_CODE_
  compStruct renumByComponent(int);
#else
  compStruct *renumByComponent(int);
#endif
  void print(FILE * = stderr);
  int findMaxDist(int *);
#ifdef _FEM_CODE_
  int findProfileSize( EqNumberer *eqNumber, int unroll=1);
#endif
  int *ptr() { return pointer; }
  Connectivity *merge(Connectivity *cn);

  void write(BinFileHandler &);
  int numNonZeroP();

  template<class Map>
  void renumberTargets(Map &);
};

//------------------------------------------------------------------------------

inline int
Connectivity::csize() { return size; }

inline int
Connectivity::num(int n) const { return (n >= size) ? 0 : pointer[n+1] - pointer[n];  }

inline int *
Connectivity::operator[](int i) { return target +pointer[i] ; }

inline const int *
Connectivity::operator[](int i) const { return target +pointer[i] ; }

inline int
Connectivity::numConnect() { return numtarget; }

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <Connectivity.C>
#endif

#endif
