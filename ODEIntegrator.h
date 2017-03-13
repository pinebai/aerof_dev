/* ODEIntegrator.h
 
   Abstracting ODE integration
*/

#ifndef _ODE_INTEGRATOR_H
#define _ODE_INTEGRATOR_H

class ODEIntegrator {

 public:

  ODEIntegrator(double _t0, double _tfinal, int _n) {

    t0 = _t0;
    tfinal = _tfinal;
    N = _n;
    
    dt = (tfinal-t0)/N;
  }

  ~ODEIntegrator() {

  }

  template <class T,class R, class M>
    void integrateFE(T& obj,R& f,R (T::*df)(double t,const R& r, const M&),const M& m ) {

    bool continueCondition = true;
    int it=0;
    double t = t0;
    while(1){
      f += dt*(obj.*df)(t,f,m);
      t += dt;
      ++it;
      if(it==N) break;
    }

  }

  template <class T,class R, class M>
    void integrateFE(T& obj,R& f,R (T::*df)(double t, const M&),const M& m ) {

    bool continueCondition = true;
    int it=0;
    double t = t0;
    while(1){
      f += dt*(obj.*df)(t,m);
      t += dt;
      ++it;
      if(it==N) break;
    }

  }

  template <class T,class R, class M>
    void integrateMidpoint(T& obj,R& f,R (T::*df)(double t,const R& r, const M&),const M& m ) {

    bool continueCondition = true;
    int it=0;
    double t = t0;
    while(1){
      f += dt*(obj.*df)(t+0.5*dt,f+0.5*dt*obj.*df(t,f,m),m);
      t += dt;
      ++it;
      if(it==N) break;
    }
  }

  template <class T,class R, class M>
    void integrateMidpoint(T& obj,R& f,R (T::*df)(double t, const M&),const M& m ) {

    bool continueCondition = true;
    int it=0;
    double t = t0;
    while(1){
      f += dt*(obj.*df)(t+0.5*dt,m);
      t += dt;
      ++it;
      if(it==N) break;
    }
  }

 private:

  // Number of iterations
  int N;

  // Timestep
  double dt;

  // Start,end times
  double t0,tfinal;
};

#endif // _ODE_INTEGRATOR_H
