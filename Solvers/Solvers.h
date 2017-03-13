#ifndef _SOLVERS_H_

class IoData;
class GeoSource;
class Domain;

template <int dim>
class 
NavierStokesCoupledSolver {
   public:
     static void solve(IoData &ioData, GeoSource &geoSource, Domain &domain);
};

template <>
class
NavierStokesCoupledSolver<5> {
   public:
     static void solve(IoData &ioData, GeoSource &geoSource, Domain &domain);
};

template <>
class
NavierStokesCoupledSolver<6> {
   public:
     static void solve(IoData &ioData, GeoSource &geoSource, Domain &domain);
};

template <>
class
NavierStokesCoupledSolver<7> {
   public:
     static void solve(IoData &ioData, GeoSource &geoSource, Domain &domain);
};

template <int dim, int neq1, int neq2>
class
NavierStokesSegSolver {
   public:
     static void solve(IoData &ioData, GeoSource &geoSource, Domain &domain);
};

template<int dim, int dimLS>
class
LevelSetSolver {
   public:
     static void solve(IoData &ioData, GeoSource &geoSource, Domain &domain);
};

template<int dim>
class
NavierStokesEmbedded {
   public:
     static void solve(IoData &ioData, GeoSource &geoSource, Domain &domain);
};

template<int dim>
class
NavierStokesEmbeddedCoupledSolver {
   public:
     static void solve(IoData &ioData, GeoSource &geoSource, Domain &domain);
};

template<int dim, int neq1, int neq2>
class
NavierStokesEmbeddedSegSolver {
   public:
     static void solve(IoData &ioData, GeoSource &geoSource, Domain &domain);
};

template<int dim, int dimLS>
class
NavierStokesMultiPhysicsEmbedded {
   public:
     static void solve(IoData &ioData, GeoSource &geoSource, Domain &domain);
};

#define _SOLVERS_H_
#endif
