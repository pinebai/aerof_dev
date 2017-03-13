/* RKIntegrator.h
  
   Class to abstract Runge-Kutta integration
 */

template <class T>
class RKIntegrator {

 public:

  enum Method { FE, RK2, RK4 };

  RKIntegrator( Method method, int N) : myMethod(method), k1(NULL),k2(NULL),
    k3(NULL), k4(NULL), ytemp(NULL) {
        
    ytemp = new T(N);
    k1 = new T(N);
    if (myMethod == RK2 ||
	myMethod == RK4)
      k2 = new T(N);

    if (myMethod == RK4) {
      k3 = new T(N);
      k4 = new T(N);
    }    
    
  }

  ~RKIntegrator() {
      
    if (k1) delete  k1;
    if (k2) delete  k2;
    if (k3) delete  k3;
    if (k4) delete  k4;
    if (ytemp) delete ytemp;
  }

  template <class Obj>
    void integrate(Obj* obj, void (Obj::*F)(double t, T& Y,T& k),
		   T& y,double t0,double h) {

    if (myMethod == FE) {

      (obj->*F)(t0, y, *k1);
      y += h*(*k1);
      
    } else if (myMethod == RK2) {

      (obj->*F)(t0, y, *k1);
      *ytemp = y + 0.5*h*(*k1);
      (obj->*F)(t0+0.5*h, *ytemp, *k2);
      y += h*(*k2);

    } else if (myMethod == RK4) {

      (obj->*F)(t0, y, *k1);
      *ytemp = y + 0.5*h*(*k1);
      (obj->*F)(t0+0.5*h, *ytemp, *k2);
      *ytemp = y + 0.5*h*(*k2);
      (obj->*F)(t0+0.5*h, *ytemp, *k3);
      *ytemp = y + h*(*k3);
      (obj->*F)(t0+h, *ytemp, *k4);
      y += h*( (1.0/6.0)*(*k1) + (1.0/3.0)*(*k2) + 
	       (1.0/3.0)*(*k3) + (1.0/6.0)*(*k4));      
    }

  }

 private:

  T* k1,*k2,*k3,*k4;
  T* ytemp;
  
  Method myMethod;

};
