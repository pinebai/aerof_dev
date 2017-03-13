#ifndef _TS_RESTART_H_
#define _TS_RESTART_H_

#include<Vector.h>
#include<Vector3D.h>

#include <FluidSelector.h>

class IoData;
class RefVal;
class DistGeoState;
class Domain;

template<int dimLS> class LevelSet;
template<class Scalar, int dim> class DistSVec;
template<class Scalar> class DistVec;
template<int dim> class DistTimeState;
template<int dim> class PostOperator;

//------------------------------------------------------------------------------

class TsRestart {

  RefVal *refVal;

  int index;

  // Boolean to flag whether pressure snapshots are needed
  // for a Kirchhoff integral.
  bool writeKPtraces;

public:

  int iteration;
  double etime;
  double residual;
  double energy[2];

  char *solutions[3];
  char *positions[3];
  char *levelsets[3];
  char *cracking[3];
  char *fluidId[3];
  char *data[3];

  char *structPos;

  int frequency;
  double frequency_dt, prtout;

  bool deleteCharStar;

private:
  bool toWrite(int it, bool lastIt, double t);

public:

  TsRestart(IoData &, RefVal *);
  TsRestart();

  void writeRestartFileNames(const char* fn);

  static void readRestartFileNames(const char* fn,
				   char* sols,
				   char* posit,
				   char* ls,
				   char* crk,
				   char* fid,
				   char* dat,
				   char* spos, Communicator* com);

  void updatePrtout(double t);

  template<int dim, int dimLS>
  void writeToDisk(int, bool, int, double, double, 
		   DistTimeState<dim> &, DistGeoState &, LevelSet<dimLS> *levelSet = 0,
                   class DynamicNodalTransfer* = NULL,
		   class FluidSelector* = NULL);

  /** Function to write the structure positions to disk. Used for the embedded-method only. */
  void writeStructPosToDisk(int, bool, Vec<Vec3D>&);

  /** Function to write the cracking data to disk. */
  void writeCrackingDataToDisk(int, bool, int, double, class DynamicNodalTransfer* = NULL);
 
// Included (MB)
  void rstVar(IoData &);
  ~TsRestart()
    {
      int last=deleteCharStar ? 3 : 1;
      for(int i=0;i<last;++i)
	{
	  delete[] data[i];
	  delete[] solutions[i];
	  delete[] positions[i];
	  delete[] levelsets[i];
	}
      delete[] structPos;
    }


  //! Function to write pressure snapshots for a Kirchhoff integral
  ///
  /// This function writes pressure snapshots for a Kirchhoff integral.
  ///
  /// UH (07/2012)
  /// 
  template<int dim>
  void writeKPtracesToDisk
    (
      IoData &iod,
      bool lastIt, int it, double t,
      DistSVec<double,3> &X,
      DistVec<double> &A,
      DistSVec<double,dim> &U,
      DistTimeState<dim> *timeState,
      Domain *domain,
      PostOperator<dim> *postOp
    );


};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <TsRestart.C>
#endif

#endif
