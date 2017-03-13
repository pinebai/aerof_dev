#include <Extrapolation.h>
#include <IoData.h>
#include <VarFcn.h>
#include <Elem.h>
#include <Vector.h>
#include <Vector3D.h>

#include <cmath>

//------------------------------------------------------------------------------

template<int dim>
Extrapolation<dim>::Extrapolation(IoData& iod, VarFcn *varFcn) : vf(varFcn)
{

  if (iod.schemes.bc.type == BoundarySchemeData::CONSTANT_EXTRAPOLATION) type = 0;
  else if (iod.schemes.bc.type == BoundarySchemeData::LINEAR_EXTRAPOLATION)  type = 1;
  else{    fprintf(stderr, "*** Error: Wrong extrapolation order\n");
	   exit(1);
  }

  fluid = iod.eqs.fluidModel.fluid;
  
  rho = iod.bc.inlet.density;

  extrapolationdata = 0;

}

//------------------------------------------------------------------------------

template<int dim>
Extrapolation<dim>::~Extrapolation()
{
  
  if (extrapolationdata) delete [] extrapolationdata;

}

//------------------------------------------------------------------------------
inline
double pointToPointDistance(Vec3D a, Vec3D b)
{

  double t = (b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) + (b[2]-a[2])*(b[2]-a[2]);
  return sqrt(t);

}

//------------------------------------------------------------------------------
template<int dim>
void Extrapolation<dim>::removeHydroStaticContribution(int i, int n[3],SVec<double,dim> &Ufar, 
				SVec<double,dim> &V, SVec<double,3> &X, double *Vinter, int node)
{
    /* Treatment for pressure/density if gravity is present
     * Only the hydrodynamic part of the pressure/density 
     * is extrapolated, not the hydrostatic part.
     */
  double p[3];

  if(fluid==FluidModelData::PERFECT_GAS ||
     fluid==FluidModelData::STIFFENED_GAS){
     /* perfect gas - stiffened gas EOS
      * Here, ptotal(V[4]) = phydrodynamic(p) + 
      *                      phydrostatic(rho*x*vec_gravity) */
    p[0] = V[n[0]][4] - vf->hydrostaticPressure(rho,X[n[0]]);
    p[1] = V[n[1]][4] - vf->hydrostaticPressure(rho,X[n[1]]);
    p[2] = V[n[2]][4] - vf->hydrostaticPressure(rho,X[n[2]]);
    Vinter[4] = p[2] + extrapolationdata[i][0].r * ( p[0] - p[2] ) +
        extrapolationdata[i][0].t * ( p[1] - p[2] );
    Vinter[4] += vf->hydrostaticPressure(rho,X[node]);

  }else if(fluid==FluidModelData::LIQUID){

    p[0] = V[n[0]][0]-Ufar[n[0]][0];
    p[1] = V[n[1]][0]-Ufar[n[1]][0];
    p[2] = V[n[2]][0]-Ufar[n[2]][0];
    Vinter[0] = p[2] + extrapolationdata[i][0].r * ( p[0] - p[2] ) +
        extrapolationdata[i][0].t * ( p[1] - p[2] );
    double tep = Vinter[0];
    Vinter[0] += Ufar[node][0];

  }

}
//------------------------------------------------------------------------------

template<int dim>
void Extrapolation<dim>::computeFaceInterpolation(int i, bool &master, int node, ElemSet &elems, 
			SVec<double,dim> &Ufar, SVec<double,dim> &V, double* Vinter1, 
			double* Vinter2, int* LinTet, 
			int* locToGlobNodeMap, SVec<double,3>& X)
{

  LinTet[0] = -1;
  LinTet[1] = -1;
  LinTet[2] = -1;

  int tet0 = extrapolationdata[i][0].tet;
  int tet1 = extrapolationdata[i][1].tet;
  int loc[3];
  int n[3];
  int m[3];

  if(tet0 != -1){

    master = true;
    Elem &elem0 = elems[tet0];
    loc[0] = elem0.faceDef(extrapolationdata[i][0].face,0);
    loc[1] = elem0.faceDef(extrapolationdata[i][0].face,1);
    loc[2] = elem0.faceDef(extrapolationdata[i][0].face,2);
    n[0] = elem0.nodeNum(loc[0]);
    n[1] = elem0.nodeNum(loc[1]);
    n[2] = elem0.nodeNum(loc[2]);

    for(int k=0; k<dim; ++k)
      Vinter1[k] = V[n[2]][k] + extrapolationdata[i][0].r * ( V[n[0]][k] - V[n[2]][k] ) +
	extrapolationdata[i][0].t * ( V[n[1]][k] - V[n[2]][k] );

    if(vf->gravity_value()>0.0)
      removeHydroStaticContribution(i,n,Ufar,V,X,Vinter1,node);

    if (tet1 != -1){ // linear extrapolation
      double V2[dim];
      Elem &elem1 = elems[tet1];
      loc[0] = elem1.faceDef(extrapolationdata[i][1].face,0);
      loc[1] = elem1.faceDef(extrapolationdata[i][1].face,1);
      loc[2] = elem1.faceDef(extrapolationdata[i][1].face,2);
      m[0] = elem1.nodeNum(loc[0]);
      m[1] = elem1.nodeNum(loc[1]);
      m[2] = elem1.nodeNum(loc[2]);

      for(int kk=0; kk<dim; kk++)
	V2[kk] = V[m[2]][kk] + extrapolationdata[i][1].r * ( V[m[0]][kk] - V[m[2]][kk] ) +
                  extrapolationdata[i][1].t * ( V[m[1]][kk] - V[m[2]][kk] );

      if(vf->gravity_value()>0.0)
        removeHydroStaticContribution(i,m,Ufar,V,X,V2,node);

      Vec3D X0 = X[node];
      Vec3D temp0 = X[n[0]];
      Vec3D temp1 = X[n[1]];
      Vec3D temp2 = X[n[2]];
      Vec3D X1 = temp2 + extrapolationdata[i][0].r * (temp0 - temp2) +
		 extrapolationdata[i][0].t * (temp1 - temp2 );
      temp0 = X[m[0]];
      temp1 = X[m[1]];
      temp2 = X[m[2]];
      Vec3D X2 = temp2 + extrapolationdata[i][1].r * (temp0 - temp2) +
		 extrapolationdata[i][1].t * (temp1 - temp2 );

      double dx1 = pointToPointDistance(X1, X0);
      double dx2 = pointToPointDistance(X2, X1);

      if(dx2!=0.0){
        for (int kkk=0; kkk<dim; kkk++)
          Vinter2[kkk] = Vinter1[kkk] + dx1/dx2 * (Vinter1[kkk]-V2[kkk]);
        LinTet[0] = m[0];
        LinTet[1] = m[1];
        LinTet[2] = m[2];

      }
      else
        fprintf(stderr, "***Warning: 1st order extrapolation for node %d, dx2=0\n", locToGlobNodeMap[node]+1);


    }else if(type!=0){
      fprintf(stderr, "*** Warning: 1st order extrapolation for node %d, 2nd tet in different subdomain\n", locToGlobNodeMap[node]+1);
    }

  }else
    master = false;
  

}

//------------------------------------------------------------------------------

