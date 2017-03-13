// DebugTools.h
// Class with various tools for debugging
#ifndef _DEBUG_TOOLS_H_
#define _DEBUG_TOOLS_H_

#define GCCBACKTRACE
#ifdef GCCBACKTRACE
#include <unistd.h>
#include <execinfo.h>
#include <cstdio>
#include <cstdlib>
#endif

#include <DistVector.h>
#include <Vector3D.h>


class DebugTools {

  public:

  template <int N,class Oper>
  static void CheckJacobian(const double jac[N*N], const double u[N],
                     const double f[N],const Oper& op, const char* name, double eps = 1.0e-7) {

    fprintf(stderr,"Checking jacobian for %s (output as (calc,act))\n", name);
    
    double fp[N], up[N];
    double Jp[N*N];
   
    for (int j = 0; j < N; ++j) {

      for (int k = 0; k < N; ++k)
        up[k] = u[k]*(1.0 + ((j!=k)?0.0:eps)) + ((j==k)?1.0e-7:0.0);
 
      op.Compute(up,fp);
      for (int k = 0; k < N; ++k)
        Jp[k*N+j] = (fp[k]-f[k])/(up[j]-u[j]);
    }

    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j)
        fprintf(stderr,"[%e %e] ", jac[i*N+j], Jp[i*N+j]);
      fprintf(stderr,"\n");
    }
  }

  static void PrintBacktrace() {

#ifdef GCCBACKTRACE
     void *array[100];
     size_t size;
     char **strings;
     size_t i;
     
     size = backtrace (array, 100);
     strings = backtrace_symbols (array, size);
     
     fprintf (stderr,"Backtrace from current location:\n"
                     "-----------------------------\n"
                     "Obtained %zd stack frames.\n", size);
     
     for (i = 0; i < size; i++)
        fprintf (stderr,"%s\n", strings[i]);
     
     free (strings);
#endif
   }
  static bool TryWaitForDebug() {

#ifdef AEROF_MPI_DEBUG
    int my_pid;// = getpid();
    MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
    char fn[256];
    sprintf(fn, ".aerofdebug.rnk.%d",my_pid);
    FILE* file = fopen(fn,"r"); 
    std::cout << "PID " << getpid() << " is MPI rank " << my_pid << std::endl;
    if (!file)
      return false; 
    bool debug_process = true;
/*    while (!feof(file)) {
      int p;
      fscanf(file, "%d",&p);
      if (p == my_pid) {
        debug_process = true;
        break;
      }
    }
*/
 
    volatile int wait = (debug_process?1:0);
    if (wait) {
 
      std::cout << "Process ID for rank " << my_pid << " is " << getpid() << std::endl;
    }

    while (wait) {

      sleep(5);
    }

    return debug_process;    
#else
    return false;
#endif
  }

  static void SpitRank() {

//#ifdef AEROF_MPI_DEBUG

    int rnk;
    MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
    std::cout << "Rank is " << rnk << std::endl;
//#endif
  }

  template <class Scalar,int dim>
  static void PrintElement(const char* tag, DistSVec<Scalar,dim>& V,int rank,int iSub, int i) {

    int rnk;
    MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
    if (rnk == rank) {

      std::cout << tag << "(" << iSub << ", " << i << ") = {";
      for (int k = 0; k < dim; ++k) 
        std::cout << V(iSub)[i][k] << ((k < dim-1)?", ":"}");
      std::cout << std::endl;
    }
  }

  template <class Scalar>
  static void PrintFluidId(const char* tag, DistVec<int>& fluidId, DistSVec<Scalar,3>& X, Vec3D pos, int rank) {

    if(pos.norm() < 1.0e-8) { //print every node that has fluidId 2
      int numLocSub = fluidId.numLocSub();
#pragma omp parallel for
      for(int iSub=0; iSub<numLocSub; iSub++) {
        double (*Xsub)[3] = X.subData(iSub);
        int *id = fluidId.subData(iSub);
        for(int i=0; i<fluidId.subSize(iSub); i++) {
          if(id[i]==2) {
            Vec3D Xnode = Vec3D(Xsub[i][0], Xsub[i][1], Xsub[i][2]);
            std::cout << tag << "Rank " << rank << ", Sub " << iSub << ", LocId " << i << ", Pos " << Xnode[0] << " " << Xnode[1] << " " << Xnode[2] << ", FluidId " << id[i] << "." << std::endl;
          }
        }
      }
      return;
    }

    double eps = 1.0e-5;
    int numLocSub = fluidId.numLocSub();
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
      double (*Xsub)[3] = X.subData(iSub);
      int *id = fluidId.subData(iSub);
      for(int i=0; i<fluidId.subSize(iSub); i++) {
        Vec3D Xnode = Vec3D(Xsub[i][0], Xsub[i][1], Xsub[i][2]);
        if(Vec3D(Xnode-pos).norm()<1.0e-5) //found it
          std::cout << tag << "Rank " << rank << ", Sub " << iSub << ", LocId " << i << ", Pos " << Xnode[0] << " " << Xnode[1] << " " << Xnode[2] << ", FluidId " << id[i] << "." << std::endl;
      }
    }
  }

};

/*
//----------------------------------------------------------
 int id_var;
 double diff, err_x, err_J;
 double V_bk_i[dim], V_bk_j[dim], V_p[dim], V_m[dim];
 double flux_p[dim], flux_m[dim], Jac_num[dim][dim];

 for(id_var=0; id_var<dim; ++id_var) V_bk_i[id_var] = Vi[id_var];
 for(id_var=0; id_var<dim; ++id_var) V_bk_j[id_var] = Vj[id_var];

 for(int dd=1; dd<=10; ++dd){

   diff = pow(10.0, -dd);

   for(id_var=0; id_var<dim; ++id_var){
     
     //Plus
     //---------------------------------------
     for(k=0; k<dim; ++k) V_p[k] = V_bk_i[k];
     V_p[id_var] = V_p[id_var] + diff;
     fluxFcn[BC_INTERNAL]->compute(length, 0.0, edgeNorm[l], edgeNormVel[l], V_p, Vj, flux_p, fluidId[i]);

     //Minus
     //---------------------------------------
     for(k=0; k<dim; ++k) V_m[k] = V_bk_i[k];
     V_m[id_var] = V_m[id_var] - diff;
     fluxFcn[BC_INTERNAL]->compute(length, 0.0, edgeNorm[l], edgeNormVel[l], V_m, Vj, flux_m, fluidId[i]);
     
     for(k=0; k<dim; ++k)
       Jac_num[k][id_var] = ( flux_p[k] - flux_m[k] ) / (2.0*diff);

   }

   for(id_var=0; id_var<dim; ++id_var){
     for(k=0; k<dim; ++k)
       fprintf(stderr, "%15.5e", Jac_num[id_var][k]);
     fprintf(stderr, " | ");     
     for(k=0; k<dim; ++k)
       fprintf(stderr, "%15.5e", dfdVi[id_var*dim+k] );
     fprintf(stderr, "\n");    
   }

   err_J = 0.0;
   for(id_var=0; id_var<dim; ++id_var){
     for(k=0; k<dim; ++k){
       err_x = sqrt( pow(Jac_num[id_var][k]-dfdVi[id_var*dim+k], 2) );
       err_J = std::max(err_J, err_x);
     }
   }

   fprintf(stderr, "****** %15.5e %24.16e ******* \n", diff, err_J);

 }

 fprintf(stderr, "\n");

//----------------------------------------------------------
*/


#endif
