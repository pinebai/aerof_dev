#ifndef TAYLOR_H_
#define TAYLOR_H_

template <class Scalar, int nvar>
class Taylor2 {
    Scalar f0;
    Scalar a[nvar];
    Scalar b[nvar][nvar];
  public:
    Taylor2() {};
    Taylor2(Scalar v, int idx) { 
      f0 = v; 
      for(int i = 0; i < nvar; ++i) {
        for(int j = 0; j < nvar; ++j)
          b[i][j] = 0;
        a[i] = 0;
      }
      a[idx] = 1.0;
    }
    Scalar &val() { return f0; }
    Scalar &d(int i) { return a[i]; }
    Scalar &d(int i, int j) { return b[i][j]; }
    Scalar val() const { return f0; }
    Scalar d(int i) const { return a[i]; }
    Scalar d(int i, int j) const { return b[i][j]; }
};

template <class Scalar, int nvar>
Taylor2<Scalar, nvar> operator+(const Taylor2<Scalar, nvar> &a, Scalar b) {
  Taylor2<Scalar, nvar> res = a;
  res.val() += b;
  return res;
}

template <class Scalar, int nvar>
Taylor2<Scalar, nvar> operator-(const Taylor2<Scalar, nvar> &a, Scalar b) {
  Taylor2<Scalar, nvar> res = a;
  res.val() -= b;
  return res;
}

template <class Scalar, int nvar>
Taylor2<Scalar, nvar> operator*(const Taylor2<Scalar, nvar> &a, Scalar b) {
  Taylor2<Scalar, nvar> res;
  res.val() = b*a.val();
  for(int i = 0; i < nvar; ++i)
    res.d(i) = b*a.d(i);
  for(int i = 0; i < nvar; ++i)
    for(int j = 0; j < nvar; ++j)
      res.d(i,j) = b*a.d(i,j);
  return res;
}

template <class Scalar, int nvar>
Taylor2<Scalar, nvar> operator+(Scalar b, const Taylor2<Scalar, nvar> &a) {
  Taylor2<Scalar, nvar> res = a;
  res.val() += b;
  return res;
}

template <class Scalar, int nvar>
Taylor2<Scalar, nvar> operator-(Scalar b, const Taylor2<Scalar, nvar> &a) {
  Taylor2<Scalar, nvar> res; 
  res.val() = b-a.val();
  for(int i = 0; i < nvar; ++i)
    res.d(i) = -a.d(i);
  for(int i = 0; i < nvar; ++i)
    for(int j = 0; j < nvar; ++j)
      res.d(i,j) = -a.d(i,j);
  return res;
  return res;
}

template <class Scalar, int nvar>
Taylor2<Scalar, nvar> operator*(Scalar b, const Taylor2<Scalar, nvar> &a) {
  Taylor2<Scalar, nvar> res;
  res.val() = b*a.val();
  for(int i = 0; i < nvar; ++i)
    res.d(i) = b*a.d(i);
  for(int i = 0; i < nvar; ++i)
    for(int j = 0; j < nvar; ++j)
      res.d(i,j) = b*a.d(i,j);
  return res;
}

template <class Scalar, int nvar>
Taylor2<Scalar, nvar> operator+(const Taylor2<Scalar, nvar> &a, const Taylor2<Scalar, nvar> &b) {
  Taylor2<Scalar, nvar> res;
  res.val() = a.val()+b.val();
  for(int i = 0; i < nvar; ++i)
    res.d(i) = a.d(i)+b.d(i);
  for(int i = 0; i < nvar; ++i)
    for(int j = 0; j < nvar; ++j)
      res.d(i,j) = b.d(i,j)+a.d(i,j);
  return res;
}

template <class Scalar, int nvar>
Taylor2<Scalar, nvar> operator-(const Taylor2<Scalar, nvar> &a, const Taylor2<Scalar, nvar> &b) {
  Taylor2<Scalar, nvar> res;
  res.val() = a.val()-b.val();
  for(int i = 0; i < nvar; ++i)
    res.d(i) = a.d(i)-b.d(i);
  for(int i = 0; i < nvar; ++i)
    for(int j = 0; j < nvar; ++j)
      res.d(i,j) = a.d(i,j)-b.d(i,j);
  return res;
}

template <class Scalar, int nvar>
Taylor2<Scalar, nvar> operator*(const Taylor2<Scalar, nvar> &a, const Taylor2<Scalar, nvar> &b) {
  Taylor2<Scalar, nvar> res;
  res.val() = a.val()*b.val();
  for(int i = 0; i < nvar; ++i)
    res.d(i) = b.val()*a.d(i)+a.val()*b.d(i);
  for(int i = 0; i < nvar; ++i)
    for(int j = 0; j < nvar; ++j)
      res.d(i,j) = b.val()*a.d(i,j)+a.val()*b.d(i,j)+0.5*(a.d(i)*b.d(j)+a.d(j)*b.d(i));
  return res;
}


#endif

