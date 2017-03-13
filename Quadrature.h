#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

#include "Dunavant.h"

class Quadrature
{
  public:
    int order;          // order
    int n_point;         // no. of quadrature points
    double (*qloc)[3]; // location of quadrature points
    double *weight;  // weight for each quadrature point

    Quadrature () {
     n_point = 0;
     qloc = 0;
     weight = 0;
    }
    
    Quadrature (int qOrder) {
      order = qOrder;
      n_point = dunavant_order_num ( order );
      weight  = new double[n_point];
      qloc = new double [n_point][3];
    
      double* xytab = new double[2*n_point];
      dunavant_rule ( order, n_point, xytab, weight );
    
      for(int i=0; i<n_point; ++i)
      {
        qloc[i][0] = xytab[0+2*i];
        qloc[i][1] = xytab[1+2*i];
        qloc[i][2] = 1.0 - xytab[0+2*i] - xytab[1+2*i];
      }
      delete [] xytab;
    }
    
    // Destructor
    ~Quadrature() {
      if (qloc) delete [] qloc;
      if (weight) delete [] weight;
      n_point = 0;
    }
};

#endif
