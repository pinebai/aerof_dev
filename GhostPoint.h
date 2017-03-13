#ifndef _GHOST_POINT_H_
#define _GHOST_POINT_H_

#include<VarFcn.h>
#include <Vector.h>
#include <Vector3D.h>
#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;

//class VarFcn;

template<class Scalar> class Vec;

template<int dim>
class GhostPoint{
 protected:
  VarFcn *varFcn;

 public:
  double* Vg;	// Sum of weighted states (rho,u,v,w,T) at the ghost-point. 
  		// After population, it is set to the state at the ghost-point.
  double *V;    // Stores the final primitive states
  double* Ws; 	// Sum of the weights 
  int ng; 	// Number of neighbours in the fluid.
  		// After all GP have been populated, ng=0.
  int ghostTag; // We store here the tag of the surrounding nodes. All the tags of the neighbours 
  		// should be the same. In the case of a complex multiphase flow simulation with Fluid/Structure 
  		// Interaction, this might be no longer true. To be done...

  int ghostTag2;

  double* V_tmp;
	  

//  ~GhostPoint();

//=============================================================================

GhostPoint(VarFcn *vf) : varFcn(vf) 
{

    ng = 0;
	Vg = new double[2*dim];
	V  = new double[2*dim];
	Ws = new double[2*dim];

	for(int i=0; i<2*dim; ++i) 
	{
      Vg[i] = 0.0;
      V[i]  = 0.0;
      Ws[i] = 0.0;
    }
    ghostTag = -2; // Inactive nodes tag
	ghostTag2 = -2;
	
	V_tmp = new double[dim];

  }
//=============================================================================

GhostPoint<dim> & operator=(const GhostPoint<dim> &GP) 
{

    varFcn = GP.varFcn;
    Vg = GP.Vg;
    V  = GP.V;
    Ws = GP.Ws;
    ng = GP.ng;
    ghostTag = GP.ghostTag;
	ghostTag2 = GP.ghostTag2;
	
    return *this;
  }
//=============================================================================

GhostPoint<dim> & operator+=(const GhostPoint<dim> &GP) 
{

    if(ghostTag<0) ghostTag = GP.ghostTag;
    else if(ghostTag != GP.ghostTag) 
      {
        fprintf(stderr,"The two ghost States refer to different Fluids\n");
        fprintf(stderr,"ghostTag: %i, GP.ghostTag: %i",ghostTag,GP.ghostTag);
        exit(-1);
      }
    Vg += GP.Vg;
    Ws += GP.Ws;
    ng += GP.ng;
    return *this;
  }
//=============================================================================

void addNeighbour(double *Vi, double *Wi, int tag) 
{

// We want to satisfy interface condition in least squares manner 
	for(int i=0; i<dim; ++i) 
	{
      Vg[i] += Wi[i]*Vi[i];
      Ws[i] += Wi[i];
    }

    // Tag check
    if(ghostTag < 0) 
	 {
      ghostTag = tag;
    }
    else if(ghostTag != tag) 
	 {
      fprintf(stderr,"We have a ghost node here with two active neighbours having different tags\n");
      fprintf(stderr,"ghostTag: %i, neighbourTag: %i",ghostTag,tag);
      exit(-1);
    }

    ng++;

  }

//=============================================================================

void addNeighbour(bool w1, double *Vi_1, int fId1, 
						bool w2, double *Vi_2, int fId2, 
						double *Wi) 
{

	if(!w1 && !w2) return;

	ghostTag = fId1; // tmp

	if(w1)
	{
		for(int i=0; i<dim; ++i) 
		{
			Vg[i] += Wi[i]*Vi_1[i];
			Ws[i] += Wi[i];
		}

		ghostTag = fId1;
	}

	if(w2) 
	{
		for(int i=0; i<dim; ++i) 
		{
			Vg[i+dim] += Wi[i]*Vi_2[i];
			Ws[i+dim] += Wi[i];
		}

		ghostTag2 = fId2;
	}
   
	ng++;

}

//=============================================================================

  double* getState()
  {
    return Vg;    
  }
//=============================================================================

  double* getPrimitiveState()
  {

	for(int k=0; k<dim; k++) V_tmp[k] = V[k];

	return V_tmp;

}
//=============================================================================

double* getPrimitiveState(int dir)
{
	
	if(dir > 0)
	{
		if(Ws[0] <= 0)
		{
			fprintf(stderr, " *** Error: in getting ghost point: dir = %d, w = %f\n", dir, Ws[0]);
			exit(-1);
		}

		for(int k=0; k<dim; k++) V_tmp[k] = V[k];
	}
	else
	{
		if(Ws[dim] <= 0)
		{
			fprintf(stderr, " *** Error: in getting ghost point: dir = %d, w = %f\n", dir, Ws[dim]);
			exit(-1);
		}

		for(int k=0; k<dim; k++) V_tmp[k] = V[k+dim];
	}

	return V_tmp;

}

//=============================================================================

bool getStatus()
{

	bool status = true;

	if(Ws[0] <= 0 && Ws[dim] <= 0) status = false;

	return status;

  }
//=============================================================================

  void reduce()
  {

	/*
	bool valid = true;

	for(int i=0; i<dim; ++i) 
		if(Ws[i] <= 0) valid = false;

	if(!valid) return;

	for(int i=0; i<dim; ++i)
	{
		Vg[i] /= Ws[i]; 
		Ws[i] = 1.0;
	}
	for(int i=0; i<dim; ++i) Ws[i+dim] = 1.0;

	// populate primitive state vector
	for (int i=0; i<dim; ++i) V[i] = Vg[i];

	varFcn->getV4FromTemperature(V, Vg[4], ghostTag);
	*/

	bool s1 = (Ws[0] > 0)   ? true : false;
	bool s2 = (Ws[dim] > 0) ? true : false;
	//std::cout << " " << std::boolalpha << s1 << " " << s2 << "\n";

	if(Ws[0] > 0)
	{
		for(int i=0; i<dim; ++i) 
		{
      Vg[i] /= Ws[i];
      Ws[i] = 1.0;
    }

		for(int i=0; i<dim; ++i) V_tmp[i] = Vg[i];

		double T = Vg[4];
		varFcn->getV4FromTemperature(V_tmp, T, ghostTag);

		for(int i=0; i<dim; ++i) V[i] = V_tmp[i];
  }

	// ------------------

	if(Ws[dim] > 0)
	{
		for(int i=0; i<dim; ++i) 
		{
			Vg[i+dim] /= Ws[i+dim]; 
			Ws[i+dim] = 1.0;
		}

		for(int i=0; i<dim; ++i) V_tmp[i] = Vg[i+dim];

		double T = Vg[4+dim];
		//varFcn->getV4FromTemperature(V_tmp, T, ghostTag2);
		varFcn->getV4FromTemperature(V_tmp, T, ghostTag);

		for(int i=0; i<dim; ++i) V[i+dim] = V_tmp[i];
	}

}

//=============================================================================

void set(double *Vf, int tag)
{

	for(int i=0; i<dim; ++i) V_tmp[i] = Vf[i];

	double T = Vf[4];
	varFcn->getV4FromTemperature(V_tmp, Vf[4], tag);

	for(int i=0; i<dim; ++i) 
	{
		V[i]     = V_tmp[i];
		V[i+dim] = V_tmp[i];

		Ws[i]     = 1.0;
		Ws[i+dim] = 1.0;
	}

	ghostTag  = tag;
	ghostTag2 = tag;
	
}

//=============================================================================

  ~GhostPoint() {
    delete [] Vg;
    delete [] V;
    delete [] Ws;
	delete [] V_tmp;
  }
//=============================================================================
};

#endif
