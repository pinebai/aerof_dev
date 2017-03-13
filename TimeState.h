#ifndef _TIME_STATE_H_
#define _TIME_STATE_H_

#include "LevelSet/LevelSetStructure.h"

class TimeData;
class GeoState;
class TimeLowMachPrec;

template<class Scalar> class Vec;
template<class Scalar, int dim> class SVec;
template<class Scalar, int dim> class GenMat;

//------------------------------------------------------------------------------

template<int dim>
class TimeState {

  TimeData &data;

  Vec<double> &dt;
  Vec<double> &idti;
  Vec<double> &idtv;
  Vec<double> &dtau;
  SVec<double,dim> &Un;
  SVec<double,dim> &Unm1;
  SVec<double,dim> &Unm2;
  SVec<double,dim> &Rn;

  Vec<double> *hhn;
  Vec<double> *hhnm1;
  //Vec<double> &Unm2;
  //Vec<double> &Rn;


public:

  TimeState(TimeData &, Vec<double> &, Vec<double> &, Vec<double> &, Vec<double> &, SVec<double,dim> &, 
	    SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);
  ~TimeState() {}

  void attachHH(Vec<double> *hhn,Vec<double> *hhnm1);

  void add_dAW_dt(bool *, GeoState &, Vec<double> &, 
		  SVec<double,dim> &, SVec<double,dim> &, LevelSetStructure *LSS=0);

  void add_dAW_dt_HH(bool *, GeoState &, Vec<double> &, 
		     Vec<double> &, Vec<double> &);

  void add_GASPrec_dAW_dt(bool *, GeoState &, Vec<double> &, 
		          SVec<double,dim> &, SVec<double,dim> &, 
                          double, double, Vec<double> &, 
                          TimeLowMachPrec &, LevelSetStructure *LSS=0);

  void add_LiquidPrec_dAW_dt(bool *, GeoState &, Vec<double> &, VarFcn *,
		             SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, 
                             TimeLowMachPrec &, LevelSetStructure *LSS=0);

  void add_dAW_dtRestrict(bool *, GeoState &, Vec<double> &, 
			  SVec<double,dim> &, SVec<double,dim> &, const std::vector<int> &sampledLocNodes) ;
  template<int dimLS>
  void add_dAW_dtLS(bool *, GeoState &, Vec<double> &, 
		    SVec<double,dimLS> &, SVec<double,dimLS> &, SVec<double,dimLS> &, 
		    SVec<double,dimLS> &, SVec<double,dimLS> &,bool);

  void add_dAW_dtau(bool *, GeoState &, Vec<double> &, 
		  SVec<double,dim> &, SVec<double,dim> &, LevelSetStructure *LSS=0);

  void add_GASPrec_dAW_dtau(bool *, GeoState &, Vec<double> &, 
		          SVec<double,dim> &, SVec<double,dim> &, 
                          double, double, Vec<double> &, 
                          TimeLowMachPrec &, LevelSetStructure *LSS=0);

  void add_LiquidPrec_dAW_dtau(bool *, GeoState &, Vec<double> &, VarFcn *,
		             SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, 
                             TimeLowMachPrec &, LevelSetStructure *LSS=0);

  template<class Scalar, int neq>
  void addToJacobianNoPrec(bool *, Vec<double> &, GenMat<Scalar,neq> &, SVec<double,dim> &,
                     VarFcn *, int*);

  template<class Scalar,int neq>
    void addToJacobianHH(Vec<double>&,  GenMat<Scalar,neq> &, Vec<double>&);

  template<class Scalar, int neq>
  void addToJacobianLS(bool *, Vec<double> &, GenMat<Scalar,neq> &, SVec<double,dim> &,bool);
  
  template<class Scalar, int neq>
  void addToJacobianNoPrecLocal(int, double, SVec<double,dim> &, GenMat<Scalar,neq> &,int);

  template<class Scalar, int neq>
  void addToJacobianGasPrec(bool *, Vec<double> &, GenMat<Scalar,neq> &, SVec<double,dim> &,
                     VarFcn *, double, double, TimeLowMachPrec &, Vec<double> &, int*);
  template<class Scalar, int neq>
  void addToJacobianGasPrecLocal(int, double, double, double, TimeLowMachPrec &, double,
				 SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToJacobianLiquidPrec(bool *, Vec<double> &, GenMat<Scalar,neq> &, SVec<double,dim> &,
                     VarFcn *, TimeLowMachPrec &, Vec<double> &, int*);
  template<class Scalar, int neq>
  void addToJacobianLiquidPrecLocal(int, double, VarFcn *, TimeLowMachPrec &, double,
				    SVec<double,dim> &, GenMat<Scalar,neq> &);
  
  template<class Scalar, int neq>
  void addToH1(bool *, Vec<double> &, GenMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToH1(bool *, Vec<double> &, GenMat<Scalar,neq> &, Scalar);

  template<class Scalar, int neq> 
  void addToH2(bool *, VarFcn *, Vec<double> &,
              SVec<double,dim> &, GenMat<Scalar,neq> &, Scalar , double);
 
  template<class Scalar, int neq>
  void addToH2(bool *, VarFcn *, Vec<double> &, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToH2(bool *, VarFcn *, Vec<double> &, SVec<double,dim> &,
               GenMat<Scalar,neq> &, Scalar);

  template<class Scalar, int neq>
  void addToH2Minus(bool *, VarFcn *, Vec<double> &, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<class Scalar, int neq> 
  void addToH2NoPrec(bool *, VarFcn *, Vec<double> &,
              SVec<double,dim> &, GenMat<Scalar,neq> &);
 
  template<class Scalar, int neq> 
  void addToH2GasPrec(bool *, VarFcn *, Vec<double> &,
              SVec<double,dim> &, GenMat<Scalar,neq> &, 
              double , double, Vec<double> &, TimeLowMachPrec &);
  template<class Scalar, int neq>
  void addToH2GasPrecLocal(int, double, VarFcn *, double, double, TimeLowMachPrec &, 
                           double, SVec<double,dim> &, GenMat<Scalar,neq> &);
 
  template<class Scalar, int neq> 
  void addToH2LiquidPrec(bool *, VarFcn *, Vec<double> &,
              SVec<double,dim> &, GenMat<Scalar,neq> &, 
              Vec<double> &, TimeLowMachPrec &);
  template<class Scalar, int neq>
  void addToH2LiquidPrecLocal(int, double, VarFcn *, TimeLowMachPrec &,
		              double, SVec<double,dim> &, GenMat<Scalar,neq> &);

  void get_dW_dt(bool *, GeoState &, Vec<double> &, SVec<double,dim> &, SVec<double,dim> &);

  void get_dWBar_dt(bool *, GeoState &, Vec<double> &, SVec<double,dim> &, SVec<double,dim> &,
                    SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);

  double getTimeNorm()  {  return dt.norm(); }

  Vec<double>& getDt() { return dt; }

  private:


  struct TimeFDCoefs {
    double c_np1, c_n, c_nm1, c_nm2;
  };
  enum DescriptorCase {
    DESCRIPTOR, HYBRID, NONDESCRIPTOR
  };
  DescriptorCase descriptorCase;

  void computeTimeFDCoefs(GeoState &, TimeFDCoefs &, Vec<double> &, int );
  void computeTimeFDCoefsSpecialBDF(GeoState &, TimeFDCoefs &, Vec<double> &, int );
  void computeDualTimeFDCoefs(GeoState &, TimeFDCoefs &, Vec<double> &, int );
  
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <TimeState.C>
#endif

#endif
