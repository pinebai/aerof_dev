#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>


#include "Communicator.h"
#include "Domain.h"
#include "DistGeoState.h"
#include "DistVector.h"
#include "IoData.h"
#include "KirchhoffIntegrator.h"
#include "SubDomain.h"
#include "Timer.h"

#ifdef AEROACOUSTIC

#include "fftw3.h"
#include "gsl/gsl_sf.h"

#endif // AEROACOUSTIC

/////////////////
#define _UH_DEBUG_
/////////////////

#ifdef AEROACOUSTIC
#pragma message "Compiling KirchhoffIntegrator.C with Aeroacoustic support"
#endif // AEROACOUSTIC

KirchhoffIntegrator::KirchhoffIntegrator
(
 IoData &iod,
 Domain *domain
 ) :
d_iod(iod), d_domain_p(domain), d_X_p((DistSVec<double, 3>*) 0),
d_globalToNodeID(),
d_timer_p(domain->getTimer()),
d_SurfType(SPHERE), d_R(0.0),
d_cyl_axis(),
d_omega_t()
{
  
#ifndef AEROACOUSTIC 
  domain->getCommunicator()->fprintf(stderr,"*** Error: AERO-F was compiled without Aeroacoustic capability."
                                            "  Rerun cmake with -DAEROACOUSTIC=ON to use this feature\n");
  exit(-1);
#endif
 
  d_X_p = new DistSVec<double, 3>(d_domain_p->getNodeDistInfo());
  
  DistVec<double> A(d_domain_p->getNodeDistInfo());
  
  //--- Get the mesh position
  DistGeoState geoState(d_iod, d_domain_p);
  // restart the geoState (positions of the mesh)
  // At return X contains the latest position of the mesh.
  {
    char temp[1]; temp[0] = '\0';
    geoState.setup1(temp, d_X_p, &A);
  }
  
  if (d_iod.surfKI.d_surfaceType == KirchhoffData::SPHERICAL)
    d_SurfType = SPHERE;
  else if (d_iod.surfKI.d_surfaceType == KirchhoffData::CYLINDRICAL)
    d_SurfType = CYLINDER;
  
}


KirchhoffIntegrator::~KirchhoffIntegrator
(
)
{
  
  if (d_X_p)
    delete d_X_p;
  
}


void KirchhoffIntegrator::Compute
(
)
{
#ifdef AEROACOUSTIC
  char prefix[192];
  sprintf(&prefix[0], "%s%s", d_iod.input.prefix, d_iod.input.strKPtraces);
  
  int *vecLen = new int[d_domain_p->getNumLocSub() + 1];
  vecLen[0] = 0;
  
  struct Data {
    int nid;
    double value;
  };
  
  //
  // Step 1 -- Read the pressure snapshots (snapshots in time)
  //
  
  int numSnapshots = 1;
  double Tf = 0.0;
  bool goodFile = true;
  
  SubDomain **SubDomain_p = d_domain_p->getSubDomain();
  Communicator *MyCom_p = d_domain_p->getCommunicator();

#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    char filename[256];
    long pos = 0;
    sprintf(&filename[0], "%s_sub%d", &prefix[0], SubDomain_p[iSub]->getGlobSubNum());
    ifstream pfile(&filename[0], ios::binary);
    pfile.seekg(pos, ios::beg);
    Data dtmp;
    if (!pfile.read((char*) &dtmp, sizeof(Data)))
    {
      fprintf(stderr, "\n !!! File %s can not be read !!! \n\n", filename);
      goodFile = false;
      continue;
    }
    pos += sizeof(Data);
    if (iSub == 0)
      numSnapshots = dtmp.nid;
    vecLen[iSub+1] = (int) dtmp.value;
    pfile.close();
  }
 
  
#ifdef _UH_DEBUG_
  if (goodFile == false)
  {
    assert(0 > 1);
  }
#endif
  
  
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
    vecLen[iSub+1] += vecLen[iSub];
  
  std::vector<int> nodeID(vecLen[d_domain_p->getNumLocSub()]);
  
  std::vector<double> in(numSnapshots * vecLen[d_domain_p->getNumLocSub()]);
  
  d_omega_t.resize((int) (numSnapshots/2) + 1);
  std::vector< std::complex<double> > out(d_omega_t.size()*vecLen[d_domain_p->getNumLocSub()]);

#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    char filename[256];
    long pos = 0;
    sprintf(&filename[0], "%s_sub%d", &prefix[0], SubDomain_p[iSub]->getGlobSubNum());
    ifstream pfile(&filename[0], ios::binary);
    pfile.seekg(pos, ios::beg);
    Data dtmp;
    if (!pfile.read((char*) &dtmp, sizeof(Data)))
    {
      fprintf(stderr, "\n !!! File %s can not be read !!! \n\n", filename);
    }
    pos += sizeof(Data);
    //---
    if (vecLen[iSub] == vecLen[iSub+1])
    {
      pfile.close();
      continue;
    }
    //---
    double ttt = 0.0;
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    Data *myData = new Data[myLen];
    for (int ii = 0; ii < numSnapshots; ++ii)
    {
      pfile.seekg(pos, ios::beg);
      if (!pfile.read((char*) &ttt, sizeof(double)))
      {
        fprintf(stderr, "\n !!! File %s can not be read !!! \n\n", filename);
      }
      pos += sizeof(double);
      pfile.seekg(pos, ios::beg);
      pfile.read((char*) myData, sizeof(Data)*myLen);
      //---
      if (iSub == 0)
      {
        //
        // UH (08/2012)
        // Here we assume that the time step is constant
        //
        if (ii == 1)
        {
          Tf = ttt * numSnapshots;
        }
        //
        if (ii < d_omega_t.size())
          d_omega_t[ii] = (ii == 0) ? 0.0 : 2.0 * M_PI * ii / Tf;
        //
      }
      //
      //---
      // UH (08/2012)
      // myData[jj].nid is the local index of a node
      // It is local per subdomain (and per processor).
      for (int jj = 0; jj < myLen; ++jj)
      {
        if (ii == 0)
        {
          nodeID[jj + vecLen[iSub]] = myData[jj].nid;
        }
        //---
        in[ii + (jj + vecLen[iSub])*numSnapshots] = myData[jj].value;
      }
      pos += sizeof(Data) * myLen;
    } // for (int ii = 0; ii < numSnapshots; ++ii)

    free(myData);
    pfile.close();

    //
    //--- Step 2: Do the FFT
    //
    
    std::vector<double> tmp11(numSnapshots, 0.0);
    std::vector< std::complex<double> > tmp22(d_omega_t.size());
    fftw_plan planbis = fftw_plan_dft_r2c_1d(numSnapshots, &tmp11[0], (fftw_complex*) &tmp22[0], FFTW_MEASURE);
    
    for (int jj = 0; jj < myLen; ++jj)
    {
      memcpy(&tmp11[0], &in[0] + (jj+vecLen[iSub])*numSnapshots, numSnapshots*sizeof(double));
      fftw_execute(planbis);
      //
      // Scale the result
      // In FFTW, the complex DFT transform is unnormalized. 
      // In other words, applying the real-to-complex (forward) 
      // and then the complex-to-real (backward) transform will multiply 
      // the input by the number of samples (or snapshots).
      //
      for (int ii = 0; ii < d_omega_t.size(); ++ii)
      {
        out[jj + vecLen[iSub] * d_omega_t.size() + ii * myLen] = tmp22[ii]/((double) numSnapshots);
      } // for (int ii = 0; ii < d_omega_t.size(); ++ii)
      
    } // for (int jj = 0; jj < myLen; ++jj)
    
    fftw_destroy_plan(planbis);
    
  } // for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  
  
  //
  // Communicate the frequencies to every processor.
  // Else processors not touching the Kirchhoff surface would not have a value.
  //

  MyCom_p->globalMax(d_omega_t.size(), &d_omega_t[0]);
 
  //
  // At this point, the time-snapshots have been converted
  // to frequency-snapshots (with nodal values for each frequency).
  //
  
  //
  // Create map of node numberings
  // (from global to the subdomain to a local numbering on the surface)
  //
  
  int totNumNodes = 0;
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
    totNumNodes += (*d_X_p)(iSub).size();
  
  if (totNumNodes < vecLen[d_domain_p->getNumLocSub()])
  {
    // The input file has non-matching dimension
    std::cerr << "\n !!! The input file with Kirchhoff data does not have matching dimensions !!!";
    std::cerr << "\n\n";
    return;
  }
  
  d_globalToNodeID.resize(totNumNodes);
  for (int jj = 0; jj < totNumNodes; ++jj)
    d_globalToNodeID[jj] = -1;
  
  totNumNodes = 0;
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    for (int jj = 0; jj < myLen; ++jj)
    {
      d_globalToNodeID[totNumNodes + nodeID[jj + vecLen[iSub]]] = jj;
    }
    totNumNodes += (*d_X_p)(iSub).size();
  }

  //
  // Evaluate and output the quantities of interest
  // (pressure at probe locations and far-field pattern)
  //

  double t20 = -d_timer_p->getTime();
  switch (d_SurfType)
  {
    case SPHERE:
      Sphere(out, vecLen);
      break;
    case CYLINDER:
      Cylinder_TensorGrid(out, vecLen);
      break;
    default:
      std::cerr << "\n !!! The Kirchhoff process is not implemented for this surface !!!\n\n";
      break;
  }

  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ...... Time: ";
    std::cout << scientific << (t20 + d_timer_p->getTime()) << std::endl;
  }
  
#endif // AEROACOUSTIC

}


//----------------------------------------------


void KirchhoffIntegrator::Sphere
(
 std::vector< std::complex<double> > &nodal,
 const int *vecLen
 )
{
  
  //
  // Compute the coefficients of p in series
  //
  
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ... Compute the coefficients for the series of p ...\n";
  }
  
  int Nmax = 0;
  std::vector< std::complex<double> > coeff_p;

  double t20 = - d_timer_p->getTime();
  convertToSHSeries(nodal, vecLen, coeff_p, Nmax);
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ...... Time: ";
    std::cout << scientific << (t20 + d_timer_p->getTime()) << std::endl;
  }

  //
  // Get the coefficient for the normal derivative
  //
    
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ... Compute the coefficients for dp/dn ...\n";
  }
  
  std::vector< std::complex<double> > coeff_dpdn(coeff_p.size());
  double t30 = -d_timer_p->getTime();
  getdpdnSHseries(coeff_p, coeff_dpdn, Nmax);
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ...... Time: ";
    std::cout << scientific << (t30 + d_timer_p->getTime()) << std::endl;
  }
  
  //
  // Compute the Kirchhoff integral for all the frequencies
  // and all the observation points.
  //
  
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ... Compute the pressure value at probes ...\n";
  }
  
  double t40 = -d_timer_p->getTime();
  EvaluateSHSatProbes(coeff_p, Nmax);
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ...... Time: ";
    std::cout << scientific << (t40 + d_timer_p->getTime()) << std::endl;
  }

  //
  // Compute the far-field pattern for all the frequencies
  //
  
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ... Compute the far-field pattern ...\n";
  }
  
  double t50 = -d_timer_p->getTime();
  ffpDataOnSphere(&nodal[0], vecLen, &coeff_dpdn[0], Nmax);
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ...... Time: ";
    std::cout << scientific << (t50 + d_timer_p->getTime()) << std::endl;
  }
  
}


//----------------------------------------------


void KirchhoffIntegrator::convertToSHSeries
(
 std::vector< std::complex<double> > &nodal,
 const int *vecLen,
 std::vector< std::complex<double> > &series,
 int &nmax
 )
{
  
  
  // This function computes the series coefficients for a P1 function
  // defined on the Kirchhoff surface (it should be a sphere).
  // The computation exploits the orthogonality of spherical harmonics.
  
  
  SubDomain **SubDomain_p = d_domain_p->getSubDomain();
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  
  std::vector<double> xi, ww;
  getQuadrature(xi, ww, 12);
  
  nmax = d_iod.surfKI.d_nyquist;
  MyCom_p->globalMax(1, &nmax);

  int numFreq = 0;
  if (vecLen[d_domain_p->getNumLocSub()] > 0)
  {
    numFreq = nodal.size() / vecLen[d_domain_p->getNumLocSub()];
  }
  MyCom_p->globalMax(1, &numFreq);
  
  std::vector<double> yNorm((nmax+1)*(nmax+1));
  
  int sLen = numFreq * (nmax + 1) * (nmax + 1);
  series.resize(sLen);
  for (int ik = 0; ik < sLen; ++ik)
    series[ik] = std::complex<double>(0.0, 0.0);
  
#ifdef _UH_DEBUG_
  double surfaceCheck = 0.0;
#endif
  
  //
  // Get the radius of the surface
  //
  
  int count = 0;
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    
    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
    
    for (int j = 0; j < faces.size(); ++j)
    {
      Face &jFace = faces[j];
      if (jFace.getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double r, t, phi;
      for (int jj = 0; jj < 3; ++jj)
      {
        double *MyCoord = XX[jFace[jj]];
        d_R += sqrt(MyCoord[0]*MyCoord[0] + MyCoord[1]*MyCoord[1] + MyCoord[2]*MyCoord[2]);
        count += 1;
      }
      //
    }
    
  }
  
  if (count > 0)
  {
    d_R = d_R / count;
  }
  MyCom_p->globalMax(1, &d_R);

  int totNumNodes = 0;

#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    if (myLen == 0)
      continue;

    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
    std::complex<double> *val = &nodal[0] + numFreq * vecLen[iSub];
    
    for (int j = 0; j < faces.size(); ++j)
    {
      
      Face &jFace = faces[j];
      if (jFace.getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double x[3], y[3], z[3];
      for (int jj = 0; jj < 3; ++jj)
      {
        double *MyCoord = XX[jFace[jj]];
        x[jj] = MyCoord[0];
        y[jj] = MyCoord[1];
        z[jj] = MyCoord[2];
      }
      //
      std::vector< std::complex<double> > p(3*numFreq);
      for (int ik = 0; ik < numFreq; ++ik)
      {
        std::complex<double> *myVal = val + ik * myLen;
        for (int jj = 0; jj < 3; ++jj)
        {
          p[jj + 3*ik] = myVal[d_globalToNodeID[totNumNodes + jFace[jj]]];
        }
      }
      //
      //
      for (int gp = 0; gp < ww.size(); ++gp)
      {
        
        double xi0 = xi[gp], xi1 = xi[gp+ww.size()], xi2 = xi[gp+2*ww.size()];
        double xp = x[0]*xi0 + x[1]*xi1 + x[2]*xi2;
        double yp = y[0]*xi0 + y[1]*xi1 + y[2]*xi2;
        double zp = z[0]*xi0 + z[1]*xi1 + z[2]*xi2;
        
        double cross = dSigma(&x[0], &y[0], &z[0], xp, yp, zp) * ww[gp];
        
        double myrp = 0.0, mytp = 0.0, myphip = 0.0;
        cart2sph_x(xp, yp, zp, myrp, mytp, myphip);
        
#ifdef _UH_DEBUG_
        surfaceCheck += d_R * d_R * cross;
#endif
        
        for (int nn = 0; nn <= nmax; ++nn)
        {
          std::vector< std::complex<double> > Yn;
          sphericalHarmonic(nn, mytp, myphip, Yn);
          //
          for (int mm = 0; mm < Yn.size(); ++mm)
          {
            yNorm[nn*nn + mm] += real( Yn[mm] * std::conj(Yn[mm]) * cross) ;
            Yn[mm] = std::conj(Yn[mm]);
          }
          //
          std::complex<double> pp;
          for (int ik = 0; ik < numFreq; ++ik)
          {
            std::complex<double> *pp_ik = &p[0] + 3*ik;
            pp = pp_ik[0]*xi0 + pp_ik[1]*xi1 + pp_ik[2]*xi2;
            //
#pragma omp critical
            {
              int shift = ik*(nmax+1)*(nmax+1) + nn*nn;
              for (int mm = 0; mm < Yn.size(); ++mm)
                series[shift + mm] += pp * Yn[mm] * cross;
            }
            
          } // for (int ik = 0; ik < numFreq; ++ik)
          
        } // for (int nn = 0; nn < nmax; ++nn)
        
      } // for (int gp = 0; gp < ww.size(); ++gp)
      
    } // for (int j = 0; j < faces.size(); ++j)
    
    totNumNodes += (*d_X_p)(iSub).size();

  } // for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  
  MyCom_p->globalSum(sLen, &series[0]);
  MyCom_p->globalSum(yNorm.size(), &yNorm[0]);
  
#ifdef _UH_DEBUG_
  MyCom_p->globalSum(1, &surfaceCheck);
  if (MyCom_p->cpuNum() == 0)
  {
    fprintf(stdout, " ... Spherical Surface = %4.3e (Error %4.3e)\n",
            surfaceCheck, std::abs(surfaceCheck - 4*M_PI*d_R*d_R)/(4*M_PI*d_R*d_R));
  }
#endif
  
  //
  //--- Filter the coefficients
  //

  double errNorm = 0.0;
  for (int ik = 0; ik < yNorm.size(); ++ik)
  {
    errNorm = (errNorm > std::abs(yNorm[ik]-1.0)) ? errNorm : std::abs(yNorm[ik]-1.0);
  }
  errNorm *= 4.0;

  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ... Filter level (relative) = " << scientific << errNorm << "\n";
  }

  for (int ik = 0; ik < numFreq; ++ik)
  {
    //--- Maximum
    double maxEntry = 0.0;
    for (int in = 0; in < (nmax+1)*(nmax+1); ++in)
    {
      if (std::abs(series[in + ik*(nmax+1)*(nmax+1)]) > maxEntry)
        maxEntry = std::abs(series[in + ik*(nmax+1)*(nmax+1)]);
    }
    //--- Filter
    for (int in = 0; in < (nmax+1)*(nmax+1); ++in)
    {
      if (std::abs(series[in + ik*(nmax+1)*(nmax+1)]) < errNorm*maxEntry)
        series[in + ik*(nmax+1)*(nmax+1)] = std::complex<double>(0.0, 0.0);
    }
  }
  
#ifdef _UH_DEBUG_
  if (MyCom_p->cpuNum() == 0)
  {

    //
    // Note the file 'trace.txt' is complete only when working with 1 proc
    //

    int totNumNodes = 0;
    std::ofstream MyTrace("trace.txt");
    for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
    {
      SVec<double,3> XX = (*d_X_p)(iSub);
      std::complex<double> *val = &nodal[0] + numFreq * vecLen[iSub];
      int myLen = vecLen[iSub+1] - vecLen[iSub];
      for (int ik = 0; ik < numFreq; ++ik)
      {
        for (int i3 = 0; i3 < XX.size(); ++i3)
        {
          if (d_globalToNodeID[totNumNodes + i3] < 0)
            continue;
          MyTrace << i3 << " " << d_globalToNodeID[totNumNodes + i3] << " ";
          double ro = 0.0, to = 0.0, po = 0.0;
          cart2sph_x(XX[i3][0], XX[i3][1], XX[i3][2], ro, to, po);
          MyTrace << XX[i3][0] << " " << XX[i3][1] << " " << XX[i3][2] << " ";
          MyTrace << ro << " ";
          MyTrace << to << " ";
          MyTrace << po << " ";
          std::complex<double> pp = evaluateSHS(&series[0] + ik*(nmax+1)*(nmax+1), nmax, to, po);
          MyTrace << real(pp) << " " << imag(pp) << " ";
          pp = val[ik * myLen + d_globalToNodeID[totNumNodes + i3]];
          MyTrace << real(pp) << " " << imag(pp) << " ";
          MyTrace << std::endl;
        }
      }
      totNumNodes += (*d_X_p)(iSub).size();
    }
    MyTrace.close();

    std::ofstream MyFile("p_coeff.txt");
    for (int ik = 0; ik < numFreq; ++ik)
    {
      MyFile << ik << " " << d_omega_t[ik] << std::endl;
      MyFile << (nmax+1)*(nmax+1) << " 0 " << std::endl;
      for (int in = 0; in < (nmax+1)*(nmax+1); ++in)
      {
        MyFile << scientific << real(series[in + ik*(nmax+1)*(nmax+1)]) << " ";
        MyFile << scientific << imag(series[in + ik*(nmax+1)*(nmax+1)]) << " ";
        MyFile << std::endl;
      }
    }
    MyFile.close();
  }
  
#endif
  
}


//-------------------------------------------


void KirchhoffIntegrator::Cylinder_TensorGrid
(
 std::vector< std::complex<double> > &nodal,
 const int *vecLen
 )
{
  
  
  int nKappa = 0, nTheta = 0;
  bool isTensorGrid = true;
  
  CylinderGrid(isTensorGrid, nKappa, nTheta);

#ifdef _UH_DEBUG_
  std::cout << " nKappa " << nKappa << " nTheta " << nTheta << std::endl;
#endif
  
  if (isTensorGrid == false)
  {
    std::cerr << "\n !!! The mesh on the cylinder does not appear to be 'orthogonal' !!! \n";
    std::cerr << "\n !!! The outputs are not computed !!! \n";
    return;
  }
  
  
  //-------------------------------
  //
  // Integrate p against conj( exp(i*2*M_PI*s1*(y - y0)/L) * exp(i*s2*theta) )
  //
  //-------------------------------
  
  
  int numFreq = d_omega_t.size();
  std::vector< std::complex<double> > pcoeff(numFreq*(2*nKappa+1)*(2*nTheta+1),
                                             std::complex<double>(0.0, 0.0));
  
  SubDomain **SubDomain_p = d_domain_p->getSubDomain();
  std::vector<double> xi, ww;
  getQuadrature(xi, ww, 6);
  
#ifdef _UH_DEBUG_
  double surfaceCheck = 0.0;
#endif
  
  int totNumNodes = 0;
#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    if (myLen == 0)
      continue;
    
    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
    std::complex<double> *val = &nodal[0] + numFreq * vecLen[iSub];
    
    for (int j = 0; j < faces.size(); ++j)
    {
      
      Face &jFace = faces[j];
      
      if (jFace.getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double x[3], y[3], z[3];
      for (int jj = 0; jj < 3; ++jj)
      {
        double *MyCoord = XX[jFace[jj]];
        x[jj] = MyCoord[0];
        y[jj] = MyCoord[1];
        z[jj] = MyCoord[2];
      }
      //
      std::vector< std::complex<double> > p(3*numFreq);
      for (int ik = 0; ik < numFreq; ++ik)
      {
        std::complex<double> *myVal = val + ik * myLen;
        for (int jj = 0; jj < 3; ++jj)
        {
          p[jj + 3*ik] = myVal[d_globalToNodeID[totNumNodes + jFace[jj]]];
        }
      }
      //
      for (int gp = 0; gp < ww.size(); ++gp)
      {
        
        double xp = 0.0, yp = 0.0, zp = 0.0;
        //
        double xi0 = xi[gp], xi1 = xi[gp+ww.size()], xi2 = xi[gp+2*ww.size()];
        xp = x[0]*xi0 + x[1]*xi1 + x[2]*xi2;
        yp = y[0]*xi0 + y[1]*xi1 + y[2]*xi2;
        zp = z[0]*xi0 + z[1]*xi1 + z[2]*xi2;
        //
        double cross = dSigma(&x[0], &y[0], &z[0], xp, yp, zp);
        cross *= ww[gp];
        //
        double myrp = sqrt(xp*xp + zp*zp);
        double mytp = (zp >= 0.0) ? acos(xp/myrp) : 2.0*M_PI - acos(xp/myrp);
        
#ifdef _UH_DEBUG_
        surfaceCheck += cross * d_R;
#endif
        
        for (int na = -nKappa; na <= nKappa; ++na)
        {
          
          double myarg = 2.0*M_PI*na*(yp - d_cyl_axis[0])/(d_cyl_axis[1] - d_cyl_axis[0]);
          std::complex<double> conj_einy = std::complex<double>(cos(myarg), -sin(myarg));
          
          for (int nt = -nTheta; nt <= nTheta; ++nt)
          {
            
            std::complex<double> conj_eint = std::complex<double>(cos(nt*mytp), -sin(nt*mytp));
            
#pragma omp critical
            for (int ik = 0; ik < numFreq; ++ik)
            {
              std::complex<double> pp;
              pp = p[3*ik]*xi0 + p[1+3*ik]*xi1 + p[2+3*ik]*xi2;
              //
              int shift = ik * (2*nKappa+1) * (2*nTheta+1);
              shift += (nt + nTheta) + (na + nKappa) * (2*nTheta + 1);
              //
              pcoeff[shift] += pp * conj_eint * conj_einy * cross;
            } // for (int ik = 0; ik < numFreq; ++ik)
            
          } // for (int nt = -nTheta; nt <= nTheta; ++nt)
          
        } // for (int na = -nKappa; na <= nKappa; ++na)
        
      } // for (int gp = 0; gp < ww.size(); ++gp)
      
    } // for (int j = 0; j < faces.size(); ++j)
    
    totNumNodes += (*d_X_p)(iSub).size();

  }
  
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  MyCom_p->globalSum(pcoeff.size(), &pcoeff[0]);
  
  for (int ii = 0; ii < pcoeff.size(); ++ii)
    pcoeff[ii] *= 0.5/(M_PI * (d_cyl_axis[1] - d_cyl_axis[0]));
  
#ifdef _UH_DEBUG_
  MyCom_p->globalSum(1, &surfaceCheck);
  std::cout << " surfaceCheck " << scientific << surfaceCheck << std::endl;
#endif
  
  
  //-------------------------------
  //
  // Filter the coefficients to keep only propagating waves
  // Compute the coefficients for the normal derivative
  //
  //-------------------------------
  
  
  std::vector< std::complex<double> > dpdn_coeff(pcoeff.size(), std::complex<double>(0.0, 0.0));
  
  for (int na = -nKappa; na <= nKappa; ++na)
  {
    double kappa_a = 2.0*M_PI*na/(d_cyl_axis[1]-d_cyl_axis[0]);
    for (int ik = 0; ik < numFreq; ++ik)
    {
      double omega_t = d_omega_t[ik];
      if (std::abs(omega_t) > std::abs(kappa_a))
      {
        double wk = sqrt(omega_t*omega_t - kappa_a*kappa_a);
        for (int nt = -nTheta; nt <= nTheta; ++nt)
        {
          int shift = ik * (2*nTheta+1) * (2*nKappa + 1);
          shift += (nt + nTheta) + (na + nKappa) * (2*nTheta + 1);
          if (pcoeff[shift] == std::complex<double>(0.0, 0.0))
            continue;
          dpdn_coeff[shift] = pcoeff[shift] * wk * hankel_prime(nt, wk * d_R) / hankel(nt, wk * d_R);
        }
      } // if (std::abs(omega_t) > std::abs(kappa_a))
      else
      {
        for (int nt = -nTheta; nt <= nTheta; ++nt)
        {
          int shift = ik * (2*nTheta+1) * (2*nKappa + 1);
          shift += (nt + nTheta) + (na + nKappa) * (2*nTheta + 1);
          pcoeff[shift] = std::complex<double>(0.0, 0.0);
          dpdn_coeff[shift] = std::complex<double>(0.0, 0.0);
        }
      } // if (std::abs(omega_t) > std::abs(kappa_a))
    }
  } // for (int na = -nKappa; na <= nKappa; ++na)
  
  
  //-------------------------------
  //
  // Compute values of p at probe locations
  //
  //-------------------------------
  
  
  std::vector<double> observations;
  getObservationPoints(observations);
  
  int numProbes = observations.size()/3;
  
  std::vector< std::complex<double> > pnoise(numFreq * numProbes,
                                             std::complex<double>(0.0, 0.0));
  
  ///--- UH_DBG Do a parallel distribution on the frequencies???
  for (int ik = 0; ik < numFreq; ++ik)
  {
    
    double omega_t = d_omega_t[ik];
    if (omega_t == 0.0)
      continue;
    
    for (int id = 0; id < numProbes; ++id)
    {
      
      double ro = 0.0, to = 0.0;
      ro = sqrt(observations[3*id]*observations[3*id]
                + observations[3*id+2]*observations[3*id+2]);

      if (ro < d_R)
      {
        std::cerr << "\n !!! Probe " << id << " can not be evaluated !!!\n\n";
        continue;
      }
      
      if (observations[3*id+2] >= 0.0)
        to = acos(observations[3*id]/ro);
      else
        to = 2.0*M_PI - acos(observations[3*id]/ro);
      double yo = observations[3*id+1];
      
      for (int na = -nKappa; na <= nKappa; ++na)
      {
        
        double kappa_a = 2.0*M_PI*na/(d_cyl_axis[1]-d_cyl_axis[0]);
        
        if (std::abs(omega_t) <= std::abs(kappa_a))
          continue;
        
        double myarg = kappa_a * (yo - d_cyl_axis[0]);
        double wt = sqrt(omega_t*omega_t - kappa_a*kappa_a);
        
        std::complex<double> einy = std::complex<double>(cos(myarg), sin(myarg));
        
        for (int nt = -nTheta; nt <= nTheta; ++nt)
        {
          
          int shift = ik * (2*nTheta+1) * (2*nKappa + 1);
          shift += (nt + nTheta) + (na + nKappa) * (2*nTheta + 1);
          if (pcoeff[shift] == std::complex<double>(0.0, 0.0))
            continue;
          
          std::complex<double> Hn_ro = hankel(nt, wt * ro);
          std::complex<double> Hn_dR = hankel(nt, wt * d_R);
          
          std::complex<double> eint = std::complex<double>(cos(nt*to), sin(nt*to));
          std::complex<double> ctmp = eint * einy * Hn_ro / Hn_dR;
          pnoise[id + ik*numProbes] += ctmp * pcoeff[shift];
        }
        
      }
      
    } // for (int id = 0; id < numProbes; ++id)
    
  } // for (int ik = 1; ik < numFreq; ++ik)
  
  if (MyCom_p->cpuNum() == 0)
    writePValues(pnoise);
  
  
#ifdef _UH_DEBUG_
  //
  // Note that this trace file is complete only with 1-proc
  //
  if (MyCom_p->cpuNum() == 0)
  {
    int totNumNodes = 0;
    std::ofstream MyTrace("trace.txt");
    for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
    {
      SVec<double,3> XX = (*d_X_p)(iSub);
      std::complex<double> *val = &nodal[0] + numFreq * vecLen[iSub];
      int myLen = vecLen[iSub+1] - vecLen[iSub];
      for (int ik = 0; ik < numFreq; ++ik)
      {
        double omega_t = d_omega_t[ik];
        if (std::abs(omega_t) == 0.0)
          continue;
      
        for (int i3 = 0; i3 < (*d_X_p)(iSub).size(); ++i3)
        {
          if (d_globalToNodeID[totNumNodes + i3] < 0)
            continue;
          MyTrace << i3 << " " << d_globalToNodeID[totNumNodes + i3] << " ";

          double *MyNode = XX[i3];
          double ro = 0.0, to = 0.0;
          ro = sqrt(MyNode[0]*MyNode[0] + MyNode[2]*MyNode[2]);
        
          if (MyNode[2] >= 0.0)
            to = acos(MyNode[0]/ro);
          else
            to = 2.0*M_PI - acos(MyNode[0]/ro);
          double yo = MyNode[1];
        
          MyTrace << MyNode[0] << " " << MyNode[1] << " " << MyNode[2] << " ";
        
          std::complex<double> pp = std::complex<double>(0.0, 0.0);
        
          for (int na = -nKappa; na <= nKappa; ++na)
          {
          
            double kappa_a = 2.0*M_PI*na/(d_cyl_axis[1]-d_cyl_axis[0]);
          
            if (std::abs(omega_t) <= std::abs(kappa_a))
              continue;
          
            double myarg = kappa_a * (yo - d_cyl_axis[0]);
            double wt = sqrt(omega_t*omega_t - kappa_a*kappa_a);
          
            std::complex<double> einy = std::complex<double>(cos(myarg), sin(myarg));
          
            for (int nt = -nTheta; nt <= nTheta; ++nt)
            {
            
              int shift = ik * (2*nTheta+1) * (2*nKappa + 1);
              shift += (nt + nTheta) + (na + nKappa) * (2*nTheta + 1);
              if (pcoeff[shift] == std::complex<double>(0.0, 0.0))
                continue;
            
              std::complex<double> Hn_ro = hankel(nt, wt * ro);
              std::complex<double> Hn_dR = hankel(nt, wt * d_R);
            
              std::complex<double> eint = std::complex<double>(cos(nt*to), sin(nt*to));
              pp += eint * einy * Hn_ro / Hn_dR * pcoeff[shift];
            
            }
          
          }
        
          MyTrace << real(pp) << " " << imag(pp) << " ";

          pp = val[ik * myLen + d_globalToNodeID[totNumNodes + i3]];
          MyTrace << real(pp) << " " << imag(pp) << " ";
          MyTrace << std::endl;
        
        }
      }
      totNumNodes += (*d_X_p)(iSub).size();
    }
    MyTrace.close();
  }

  if (MyCom_p->cpuNum() == 0)
  {
    std::ofstream MyFile("p_coeff.txt");
    for (int ik = 0; ik < numFreq; ++ik)
    {
      double omega_t = d_omega_t[ik];
      MyFile << ik << " " << omega_t << " ";
      MyFile << (2*nTheta+1)*(2*nKappa+1) << " 0 " << std::endl;
      for (int na = -nKappa; na <= nKappa; ++na)
      {
        double kappa_a = 2.0*M_PI*na/(d_cyl_axis[1]-d_cyl_axis[0]);
        double wt = (std::abs(omega_t) > std::abs(kappa_a)) ? 
                    sqrt(omega_t*omega_t - kappa_a*kappa_a) : 0.0;
        for (int nt = -nTheta; nt <= nTheta; ++nt)
        {
          MyFile << scientific << wt << " " << nt << " ";
          int shift = ik * (2*nTheta+1) * (2*nKappa + 1);
          shift += (nt + nTheta) + (na + nKappa) * (2*nTheta + 1);
          MyFile << scientific << real(pcoeff[shift]) << " ";
          MyFile << scientific << imag(pcoeff[shift]) << " ";
          MyFile << std::endl;
        }
      }
    }
    MyFile.close();
  }

  if (MyCom_p->cpuNum() == 0)
  {
    std::ofstream MFile("dpdn_coeff.txt");
    for (int ik = 0; ik < numFreq; ++ik)
    {
      double omega_t = d_omega_t[ik];
      MFile << ik << " " << d_omega_t[ik] << std::endl;
      MFile << (2*nTheta+1)*(2*nKappa+1) << " 0 " << std::endl;
      for (int na = -nKappa; na <= nKappa; ++na)
      {
        double kappa_a = 2.0*M_PI*na/(d_cyl_axis[1]-d_cyl_axis[0]);
        double wt = (std::abs(omega_t) > std::abs(kappa_a)) ? 
                    sqrt(omega_t*omega_t - kappa_a*kappa_a) : 0.0;
        for (int nt = -nTheta; nt <= nTheta; ++nt)
        {
          MFile << scientific << wt << " " << nt << " ";
          int shift = ik * (2*nTheta+1) * (2*nKappa + 1);
          shift += (nt + nTheta) + (na + nKappa) * (2*nTheta + 1);
          MFile << scientific << real(dpdn_coeff[shift]) << " ";
          MFile << scientific << imag(dpdn_coeff[shift]) << " ";
          MFile << std::endl;
        }
      }
    }
    MFile.close();
  }

#endif

  
  //-------------------------------
  //
  // Compute far-field patterns
  //
  //-------------------------------
  
  
  int numInc = d_iod.surfKI.d_angularIncrement;
  int numDir = (numInc/2 + 1) * numInc;
  numFreq = d_omega_t.size();
  
  std::vector<double> directions(3*numDir, 0.0);
  numDir = 0;
  for (int ib = 0; ib <= numInc/2; ++ib)
  {
    double beta = ib * M_PI * 2.0 / numInc - 0.5 * M_PI;
    double cb = cos(beta), sb = sin(beta);
    for (int ia = 0; ia < numInc; ++ia)
    {
      double alpha = ia * 2.0 * M_PI / numInc;
      directions[3*numDir] = cos(alpha) * cb;
      directions[3*numDir+1] = sin(alpha) * cb;
      directions[3*numDir+2] = sb;
      numDir += 1;
    }
  }
  
  std::vector< std::complex<double> > ffp(numDir * numFreq, std::complex<double>(0.0, 0.0));
  
#ifdef _UH_DEBUG_
  surfaceCheck = 0.0;
#endif
  
  totNumNodes = 0;
#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    if (myLen == 0)
      continue;
    
    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
    const std::complex<double> *val = &nodal[0] + numFreq * vecLen[iSub];
    
    for (int j = 0; j < faces.size(); ++j)
    {
      
      Face &jFace = faces[j];
      if (jFace.getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double x[3], y[3], z[3];
      for (int jj = 0; jj < 3; ++jj)
      {
        double *MyCoord = XX[jFace[jj]];
        x[jj] = MyCoord[0];
        y[jj] = MyCoord[1];
        z[jj] = MyCoord[2];
      }
      //
      std::vector< std::complex<double> > p(3*numFreq);
      for (int ik = 0; ik < numFreq; ++ik)
      {
        const std::complex<double> *myVal = val + ik * myLen;
        for (int jj = 0; jj < 3; ++jj)
        {
          p[jj + 3*ik] = myVal[d_globalToNodeID[totNumNodes + jFace[jj]]];
        }
      }
      //
      for (int gp = 0; gp < ww.size(); ++gp)
      {
        
        double xi0 = xi[gp], xi1 = xi[gp + ww.size()], xi2 = xi[gp+2*ww.size()];
        double xp = x[0]*xi0 + x[1]*xi1 + x[2]*xi2;
        double yp = y[0]*xi0 + y[1]*xi1 + y[2]*xi2;
        double zp = z[0]*xi0 + z[1]*xi1 + z[2]*xi2;
        
        double cross = dSigma(&x[0], &y[0], &z[0], xp, yp, zp) * ww[gp];
        
        double myrp = sqrt(xp*xp + zp*zp);
        double mytp = acos(xp / myrp);
        mytp = (zp < 0.0) ? 2.0 * M_PI - mytp : mytp;
        
#ifdef _UH_DEBUG_
        surfaceCheck += d_R * cross;
#endif
        
        //
        std::vector< std::complex<double> > dpdn(numFreq);
        for (int ik = 0; ik < numFreq; ++ik)
        {
          double omega_t = d_omega_t[ik];
          dpdn[ik] = std::complex<double>(0.0, 0.0);
          for (int na = -nKappa; na <= nKappa; ++na)
          {
            
            double kappa_a = 2.0*M_PI*na/(d_cyl_axis[1]-d_cyl_axis[0]);
            if (std::abs(omega_t) <= std::abs(kappa_a))
              continue;
            
            double myarg = kappa_a * (yp - d_cyl_axis[0]);
            std::complex<double> einy = std::complex<double>(cos(myarg), sin(myarg));
            
            for (int nt = -nTheta; nt <= nTheta; ++nt)
            {
              int shift = ik * (2*nTheta+1) * (2*nKappa + 1);
              shift += (nt + nTheta) + (na + nKappa) * (2*nTheta + 1);
              if ((std::real(dpdn_coeff[shift]) == 0.0) && 
                  (std::imag(dpdn_coeff[shift]) == 0.0))
                continue;
              
              std::complex<double> eint = std::complex<double>(cos(nt*mytp), sin(nt*mytp));
              dpdn[ik] += eint * einy * dpdn_coeff[shift];
              
            }
            
          }
        }
        //
        std::vector< std::complex<double> > p_p(numFreq);
        for (int ik = 0; ik < numFreq; ++ik)
        {
          p_p[ik] = p[3*ik] * xi0 + p[1+3*ik] * xi1 + p[2 + 3*ik] * xi2;
        }
        //
        
        for (int id = 0; id < numDir; ++id)
        {
          //
          double DirDotX = xp * directions[3*id];
          DirDotX += zp * directions[3*id+2];
          //
          double DirDotNu = DirDotX / myrp;
          //
          DirDotX += yp * directions[3*id+1];
          //
          for (int ik = 0; ik < numFreq; ++ik)
          {
            
            double omega_t = d_omega_t[ik];
            //
            std::complex<double> e_mi_k_x_d(cos(omega_t*DirDotX), -sin(omega_t*DirDotX));
            //
            std::complex<double> cTmp = e_mi_k_x_d * d_R * cross;
#pragma omp critical
            {
              ffp[id + ik * numDir] += cTmp * (dpdn[ik]
                      + p_p[ik] * std::complex<double>(0.0, omega_t * DirDotNu));
            }
            
          } // for (int ik = 0; ik < numFreq; ++ik)
          
        } // for (int id = 0; id < numDir; ++id)
        
      } // for (int gp = 0; gp < ww.size(); ++gp)
      
    } // for (int j = 0; j < faces.size(); ++j)
    
    totNumNodes += (*d_X_p)(iSub).size();

  } // for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  
  
  MyCom_p->globalSum(ffp.size(), &ffp[0]);
  
  
#ifdef _UH_DEBUG_
  MyCom_p->globalSum(1, &surfaceCheck);
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " Surface " << surfaceCheck << " ";
    std::cout << std::abs(surfaceCheck - d_R * (d_cyl_axis[1] - d_cyl_axis[0]) * 2*M_PI)/(d_R * (d_cyl_axis[1] - d_cyl_axis[0]) * 2*M_PI);
    std::cout << std::endl;
  }
#endif
  
  //
  // Scale the computed far-field pattern
  // The minus sign appears because of the normal orientation.
  //
  for (int ii = 0; ii < ffp.size(); ++ii)
    ffp[ii] = -ffp[ii] / (4.0 * M_PI);
  
  if (MyCom_p->cpuNum() == 0)
    writeFFPValues(ffp);
  
}


//-------------------------------------------


double KirchhoffIntegrator::dSigma
(
 const double *x, const double *y, const double *z,
 double &xp, double &yp, double &zp
 ) const
{
  
  double area = 0.0;
  
  if (d_SurfType == SPHERE)
  {
    
    double myrp = sqrt(xp*xp + yp*yp + zp*zp);
    double dMdr[3], dMds[3], tmpdot;
    tmpdot = xp * (x[1] - x[0]) + yp * (y[1] - y[0]) + zp * (z[1] - z[0]);
    dMdr[0] = d_R/myrp * (x[1] - x[0]);
    dMdr[0] -= d_R/(myrp*myrp*myrp) * xp * tmpdot;
    dMdr[1] = d_R/myrp * (y[1] - y[0]);
    dMdr[1] -= d_R/(myrp*myrp*myrp) * yp * tmpdot;
    dMdr[2] = d_R/myrp * (z[1] - z[0]);
    dMdr[2] -= d_R/(myrp*myrp*myrp) * zp * tmpdot;
    tmpdot = xp * (x[2] - x[0]) + yp * (y[2] - y[0]) + zp * (z[2] - z[0]);
    dMds[0] = d_R/myrp * (x[2] - x[0]);
    dMds[0] -= d_R/(myrp*myrp*myrp) * xp * tmpdot;
    dMds[1] = d_R/myrp * (y[2] - y[0]);
    dMds[1] -= d_R/(myrp*myrp*myrp) * yp * tmpdot;
    dMds[2] = d_R/myrp * (z[2] - z[0]);
    dMds[2] -= d_R/(myrp*myrp*myrp) * zp * tmpdot;
    area += std::pow(dMdr[1]*dMds[2] - dMdr[2]*dMds[1], 2.0);
    area += std::pow(dMdr[2]*dMds[0] - dMdr[0]*dMds[2], 2.0);
    area += std::pow(dMdr[0]*dMds[1] - dMdr[1]*dMds[0], 2.0);
    area = sqrt(area) * 0.5 / (d_R * d_R);
   
    //
    // Place the integration point on the sphere
    //
    double scale = d_R / myrp;
    xp *= scale;
    yp *= scale;
    zp *= scale;
    
  }
  else if (d_SurfType == CYLINDER)
  {
    
    double myrp = sqrt(xp*xp + zp*zp);
    //
    double dMdr[3], dMds[3], tmpdot;
    tmpdot = xp * (x[1] - x[0]) +  zp * (z[1] - z[0]);
    dMdr[0] = d_R/myrp * (x[1] - x[0]);
    dMdr[0] -= d_R/(myrp*myrp*myrp) * xp * tmpdot;
    dMdr[1] = y[1] - y[0];
    dMdr[2] = d_R/myrp * (z[1] - z[0]);
    dMdr[2] -= d_R/(myrp*myrp*myrp) * zp * tmpdot;
    //
    tmpdot = xp * (x[2] - x[0]) + zp * (z[2] - z[0]);
    dMds[0] = d_R/myrp * (x[2] - x[0]);
    dMds[0] -= d_R/(myrp*myrp*myrp) * xp * tmpdot;
    dMds[1] = y[2] - y[0];
    dMds[2] = d_R/myrp * (z[2] - z[0]);
    dMds[2] -= d_R/(myrp*myrp*myrp) * zp * tmpdot;
    //
    area += std::pow(dMdr[1]*dMds[2] - dMdr[2]*dMds[1], 2.0);
    area += std::pow(dMdr[2]*dMds[0] - dMdr[0]*dMds[2], 2.0);
    area += std::pow(dMdr[0]*dMds[1] - dMdr[1]*dMds[0], 2.0);
    area = sqrt(area) * 0.5 / (d_R);

    //
    // Place the integration point on the sphere
    //
    double scale = d_R / myrp;
    xp *= scale;
    zp *= scale;

  }
  
  return area;
  
}


//-------------------------------------------


struct ltdouble
{
  bool operator()(const double a, const double b) const
  {
    //
    // This operator compares number with a 'tolerance'
    // It is used to identify coordinates of vertices on a surface
    // If the nodes do not have 'exactly' the same coordinates,
    // the tolerance is used to combine those values into 1 coordinate.
    //
    double tol = 5.0e-03;
    bool a_lt_b = true;
    double small = (std::abs(a) < std::abs(b)) ? std::abs(a) : std::abs(b);
    if ((a == b) || (std::abs(a - b) < tol * small))
    {
      a_lt_b = false;
    }
    else
    {
      a_lt_b = (a * (1.0 + tol) < b * (1.0 - tol));
    }
    return a_lt_b;
  }
};


//-------------------------------------------


void KirchhoffIntegrator::CylinderGrid
(
 bool &isTensorGrid,
 int &nKappa,
 int &nTheta
 )
{
  
  //-------------------------------
  //
  // Gather the coordinates along the axis and the radial direction
  // Count the number of points
  //
  //-------------------------------
  
  
  std::set<double, ltdouble> axis_coord;
  std::set<double, ltdouble> angle_coord;
  int numPoints = 0;
  d_R = 0.0;
  
  int totNumNodes = 0;
#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    SVec<double,3> XX = (*d_X_p)(iSub);
    for (int i3 = 0; i3 < (*d_X_p)(iSub).size(); ++i3)
    {
      if (d_globalToNodeID[totNumNodes + i3] < 0)
        continue;
      double rr = sqrt(XX[i3][0] * XX[i3][0] + XX[i3][2] * XX[i3][2]);
      if (iSub == 0)
        d_R = (rr > d_R) ? rr : d_R;
      double theta;
      if (XX[i3][2] >= 0.0)
        theta = acos(XX[i3][0]/rr);
      else
      {
        theta = 2.0*M_PI - acos(XX[i3][0]/rr);
        if (std::abs(theta - 2.0 * M_PI) < M_PI * 1.0e-04)
          theta = 0.0;
      }
#pragma omp critical
      {
        axis_coord.insert(XX[i3][1]);
        angle_coord.insert(theta);
        numPoints += 1;
      }
    } // for (int i3 = 0; i3 < (*d_X_p)(iSub).size(); ++i3)
    totNumNodes += (*d_X_p)(iSub).size();
  }
  
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  MyCom_p->globalSum(1, &numPoints);
  MyCom_p->globalMax(1, &d_R);
  
  int numAxis = axis_coord.size();
  MyCom_p->globalMax(1, &numAxis);
  RecInfo rInfo;
  d_cyl_axis[0] = 1.0e308; d_cyl_axis[1] = -1.0e308;
  if (MyCom_p->cpuNum() > 0)
  {
    std::vector<double> tmpD(axis_coord.begin(), axis_coord.end());
    MyCom_p->sendTo(0, 11, &tmpD[0], axis_coord.size());
  }
  else
  {
    std::vector<double> tmpD(numAxis);
    for (int ii = 1; ii < MyCom_p->size(); ++ii)
    {
      rInfo = MyCom_p->recFrom(ii, 11, &tmpD[0], numAxis);
      for (int jj = 0; jj < rInfo.len; ++jj)
        axis_coord.insert(tmpD[jj]);
    }
    std::set<double>::iterator it;
    it = axis_coord.begin();
    d_cyl_axis[0] = *it;
    d_cyl_axis[1] = *it;
    for (it = axis_coord.begin(); it != axis_coord.end(); ++it)
    {
      d_cyl_axis[0] = (d_cyl_axis[0] > *it) ? *it : d_cyl_axis[0];
      d_cyl_axis[1] = (d_cyl_axis[1] < *it) ? *it : d_cyl_axis[1];
    }
  }
  
  numAxis = axis_coord.size();
  MyCom_p->broadcast(1, &numAxis);
  MyCom_p->broadcast(2, &d_cyl_axis[0]);
  
  int numAngle = angle_coord.size();
  MyCom_p->globalMax(1, &numAngle);
  if (MyCom_p->cpuNum() > 0)
  {
    std::vector<double> tmpD(angle_coord.begin(), angle_coord.end());
    MyCom_p->sendTo(0, 11, &tmpD[0], angle_coord.size());
  }
  else
  {
    std::vector<double> tmpD(numAngle);
    for (int ii = 1; ii < MyCom_p->size(); ++ii)
    {
      rInfo = MyCom_p->recFrom(ii, 11, &tmpD[0], numAngle);
      for (int jj = 0; jj < rInfo.len; ++jj)
        angle_coord.insert(tmpD[jj]);
    }
  }
  
  numAngle = angle_coord.size();
  MyCom_p->broadcast(1, &numAngle);
  
  if (numPoints != numAngle * numAxis)
  {
    isTensorGrid = false;
  }
  
  //
  // Define the wavenumbers in the angular direction
  // and in the axial direction
  //
  
  //
  // Angular Direction Wavenumber
  //
  // If numAngle is even: -numAngle/2+1, ..., 0, numAngle/2
  // If numAngle is odd : -(numAngle-1)/2, ..., 0, ..., (numAngle-1)/2
  //
  // For simplicity, we will set numAngle to be an odd number
  //
  
  nTheta = (numAngle < d_iod.surfKI.d_nyquist) ? numAngle : d_iod.surfKI.d_nyquist;
  nTheta = (nTheta % 2 == 1) ? nTheta/2 : (nTheta-2)/2;
  
  //
  // Axial Direction Wavenumber
  //
  // If numAxis is even: -numAxis/2+1, ..., 0, numAxis/2
  // If numAxis is odd : -(numAxis-1)/2, ..., 0, ..., (numAxis-1)/2
  //
  // For simplicity, we will set numAxis to be an odd number
  //
  
  nKappa = (numAxis % 2 == 1) ? numAxis/2 : (numAxis-2)/2;
  
  
}


//-------------------------------------------


void KirchhoffIntegrator::getObservationPoints
(
 std::vector<double> &observations
 ) const
{
  
  int numProbes = 0;
  for (int ii = 0; ii < Probes::MAXNODES; ++ii)
  {
    Probes& myProbes = d_iod.output.transient.probes;
    if (myProbes.myNodes[ii].locationX < -1.0e10)
      break;
    numProbes += 1;
  }
  
  observations.resize(3 * numProbes);
  
  for (int ii = 0; ii < numProbes; ++ii)
  {
    Probes& myProbes = d_iod.output.transient.probes;
    //
    observations[3*ii] = myProbes.myNodes[ii].locationX;
    observations[3*ii+1] = myProbes.myNodes[ii].locationY;
    observations[3*ii+2] = myProbes.myNodes[ii].locationZ;
  }
  
}


//-------------------------------------------


void KirchhoffIntegrator::writePValues
(
 std::vector< std::complex<double> > &pnoise
 ) const
{
  
  char *myName = new char[strlen(d_iod.output.transient.probes.prefix) + strlen(d_iod.output.transient.probes.pressure) + 1];
  sprintf(myName, "%s%s", d_iod.output.transient.probes.prefix, d_iod.output.transient.probes.pressure);
  
  int numProbes = pnoise.size() / d_omega_t.size();
  
  // Output PNOISE in ASCII file
  FILE *out = fopen(myName, "w");
  for (int ik = 0; ik < d_omega_t.size(); ++ik)
  {
    fprintf(out, "%d %e ", ik, d_omega_t[ik]);
    for (int id = 0; id < numProbes; ++id)
    {
      fprintf(out, "%e ", std::real(pnoise[id + ik * numProbes]));
      fprintf(out, "%e ", std::imag(pnoise[id + ik * numProbes]));
      fprintf(out, "%e ", std::abs(pnoise[id + ik * numProbes]));
    }
    fprintf(out, "\n");
  }
  fclose(out);
  delete[] myName;
  
}


//-------------------------------------------


void KirchhoffIntegrator::writeFFPValues
(
 std::vector< std::complex<double> > &ffp
 ) const
{
  
  int numInc = d_iod.surfKI.d_angularIncrement;
  int numDir = (numInc/2 + 1) * numInc;
  
  char *myName = new char[strlen(d_iod.output.transient.probes.prefix)
                          + strlen(d_iod.output.transient.probes.farfieldpattern) + 1];
  sprintf(myName, "%s%s", d_iod.output.transient.probes.prefix,
          d_iod.output.transient.probes.farfieldpattern);
  
  // Output FFP in ASCII file
  FILE *out = fopen(myName, "w");
  for (int ik = 0; ik < d_omega_t.size(); ++ik)
  {
    fprintf(out, "%d %e ", ik, d_omega_t[ik]);
    int iDir = 0;
    for (int ib = 0; ib <= numInc/2; ++ib)
    {
      double beta = ib * M_PI * 2.0 / numInc - 0.5 * M_PI;
      for (int ia = 0; ia < numInc; ++ia)
      {
        double alpha = ia * 2.0 * M_PI / numInc;
        fprintf(out, "%e %e ", alpha, beta);
        fprintf(out, "%e ", std::real(ffp[iDir + ik*numDir]));
        fprintf(out, "%e ", std::imag(ffp[iDir + ik*numDir]));
        fprintf(out, "%e ", std::abs(ffp[iDir + ik*numDir]));
        iDir += 1;
      }
    }
    fprintf(out, "\n");
  }
  fclose(out);
  delete[] myName;
  
}


//-------------------------------------------


void KirchhoffIntegrator::getdpdnSHseries
(
 const std::vector< std::complex<double> > &coeff,
 std::vector< std::complex<double> > &coeff_dudn,
 const int nMax
 ) const
{
  
  int numFrequencies = d_omega_t.size();
  
  int count = 0;
  for (int ik = 0; ik < numFrequencies; ++ik)
  {
    
    double kappa = d_omega_t[ik];
    
    if ((ik == 0) || (kappa == 0.0))
    {
      for (int in = 0; in <= nMax; ++in)
      {
        double scalar = -(in + 1.0)/d_R;
        for (int jn = 0; jn <= 2*in; ++jn)
        {
          coeff_dudn[count] = coeff[count] * scalar;
          count += 1;
        } // for (int jn = 0; jn <= 2*in; ++jn)
      } // for (int in = 0; in <= nMax; ++in)
    }
    else
    {
      
      for (int in = 0; in <= nMax; ++in)
      {
        //
        double maxEntry = 0.0;
        for (int jn = 0; jn <= 2*in; ++jn)
          maxEntry = (maxEntry > std::abs(coeff[count+jn])) ? maxEntry : std::abs(coeff[count+jn]);
        if (maxEntry == 0.0)
        {
          count += 2*in+1;
          continue;
        }
        //
        std::complex<double> scalar = kappa * besselh_prime(in, kappa*d_R) / besselh(in, kappa*d_R);
        for (int jn = 0; jn <= 2*in; ++jn)
        {
          coeff_dudn[count] = coeff[count] * scalar;
          count += 1;
        } // for (int jn = 0; jn <= 2*in; ++jn)
        //
      } // for (int in = 0; in <= nMax; ++in)
      
    }
    
  } // for (int ik = 0; ik < numFrequencies; ++ik)
  
  
  ////////////////////
#ifdef _UH_DEBUG_
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  if (MyCom_p->cpuNum() == 0)
  {
    std::ofstream MyFile("dpdn_coeff.txt");
    for (int ik = 0; ik < numFrequencies; ++ik)
    {
      MyFile << ik << " " << d_omega_t[ik] << std::endl;
      MyFile << (nMax+1)*(nMax+1) << " 0 " << std::endl;
      for (int in = 0; in < (nMax+1)*(nMax+1); ++in)
      {
        MyFile << scientific << real(coeff_dudn[in + ik*(nMax+1)*(nMax+1)]) << " ";
        MyFile << scientific << imag(coeff_dudn[in + ik*(nMax+1)*(nMax+1)]) << " ";
        MyFile << std::endl;
      }
    }
    MyFile.close();
  }
#endif
  ///////////////////
  
  
}


//-------------------------------------------


void KirchhoffIntegrator::integrateOnSphere
(
 const std::complex<double> *pvalues,
 const int *vecLen,
 const std::complex<double> *dpdn_coeff,
 int nmax,
 std::complex<double> *p_coeff
 )
{
  
  // This function computes the Kirchhoff integral.
  
  int numFreq = d_omega_t.size();
  
  std::vector<double> observations;
  getObservationPoints(observations);
  
  int numProbes = observations.size()/3;
  
  SubDomain **SubDomain_p = d_domain_p->getSubDomain();
  
  std::vector<double> xi, ww;
  getQuadrature(xi, ww, 12);
  
  int intLen = numFreq * numProbes;
  std::vector< std::complex<double> > pnoise(intLen, std::complex<double>(0.0, 0.0));
  
  int totNumNodes = 0;
  double t10 = -d_timer_p->getTime();
#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    if (myLen == 0)
      continue;
    
    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
    const std::complex<double> *val = pvalues + numFreq * vecLen[iSub];
    
    for (int j = 0; j < faces.size(); ++j)
    {
      Face &jFace = faces[j];
      if (jFace.getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double x[3], y[3], z[3];
      for (int jj = 0; jj < 3; ++jj)
      {
        double *MyCoord = XX[jFace[jj]];
        x[jj] = MyCoord[0];
        y[jj] = MyCoord[1];
        z[jj] = MyCoord[2];
      }
      //
      std::vector< std::complex<double> > p(3*numFreq);
      for (int ik = 0; ik < numFreq; ++ik)
      {
        const std::complex<double> *myVal = val + ik * myLen;
        for (int jj = 0; jj < 3; ++jj)
        {
          p[jj + 3*ik] = myVal[d_globalToNodeID[totNumNodes+faces[j][jj]]];
        }
      }
      //
      for (int gp = 0; gp < ww.size(); ++gp)
      {
        
        double xi0 = xi[gp], xi1 = xi[gp+ww.size()], xi2 = xi[gp+2*ww.size()];
        double xp = x[0]*xi0 + x[1]*xi1 + x[2]*xi2;
        double yp = y[0]*xi0 + y[1]*xi1 + y[2]*xi2;
        double zp = z[0]*xi0 + z[1]*xi1 + z[2]*xi2;
        
        double cross = dSigma(&x[0], &y[0], &z[0], xp, yp, zp) * ww[gp];
        
        double myrp = 0.0, tp = 0.0, phip = 0.0;
        cart2sph_x(xp, yp, zp, myrp, tp, phip);
        
        //
        std::complex<double> pp;
        std::vector< std::complex<double> > dpdn(numFreq);
        for (int ik = 0; ik < numFreq; ++ik)
        {
          dpdn[ik] = evaluateSHS(dpdn_coeff + ik*(nmax+1)*(nmax+1), nmax, tp, phip);
        }
        //
        for (int id = 0; id < numProbes; ++id)
        {
          double rv[3];
          rv[0] = xp - observations[3*id];
          rv[1] = yp - observations[3*id+1];
          rv[2] = zp - observations[3*id+2];
          //
          double rho = sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]);
          //
          double drho_dn = (rv[0] * xp + rv[1] * yp + rv[2] * zp ) / (myrp * rho);
          //
          for (int ik = 0; ik < numFreq; ++ik)
          {
            
            double kappa = d_omega_t[ik];
            
            std::complex<double> Green = 0.25/(M_PI*rho) * exp(std::complex<double>(0.0, kappa*rho));
            std::complex<double> dGreendNu = Green * drho_dn * std::complex<double>(-1.0/rho, kappa);
            
            //
            pp = p[3*ik]*xi[gp];
            pp += p[1+3*ik]*xi[gp+ww.size()];
            pp += p[2+3*ik]*xi[gp+2*ww.size()];
            //
#pragma omp critical
            {
              pnoise[id + ik * numProbes] += (pp * dGreendNu - dpdn[ik] * Green) * cross * d_R * d_R;
            }
            
          } // for (int ik = 0; ik < numFreq; ++ik)
          
        } // for (int id = 0; id < numProbes; ++id)
        
      } // for (int gp = 0; gp < ww.size(); ++gp)
      
    } // for (int j = 0; j < faces.size(); ++j)
    
    totNumNodes += (*d_X_p)(iSub).size();

  } // for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
    
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  MyCom_p->globalSum(intLen, &pnoise[0]);
  
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " ***** Time: " << scientific << (t10 + d_timer_p->getTime()) << std::endl;
  }

  //
  // Write the values to the disk
  //

  if (MyCom_p->cpuNum() == 0)
  {
    writePValues(pnoise);
  } // if (MyCom_p->cpuNum() == 0)
  
  
}


//-------------------------------------------


void KirchhoffIntegrator::EvaluateSHSatProbes
(
 std::vector< std::complex<double> > &p_coeff,
 int nmax
 )
{
  
  Communicator *MyCom_p = d_domain_p->getCommunicator();

  int numFreq = d_omega_t.size();
  
  std::vector<double> observations;
  getObservationPoints(observations);
  
  int numProbes = observations.size()/3;
  
  std::vector< std::complex<double> > pnoise(numFreq*numProbes,std::complex<double>(0.0,0.0));
   
  //--- UHDBG ... This loop could be parallelized????
  for (int ik = 1; ik < numFreq; ++ik)
  {
    for (int id = 0; id < numProbes; ++id)
    {
      double ro = 0.0, to = 0.0, po = 0.0;
      cart2sph_x(observations[3*id], observations[3*id+1], observations[3*id+2], ro, to, po);
      pnoise[id + ik * numProbes] = evaluateSHS(&p_coeff[0] + ik*(nmax+1)*(nmax+1),
                                                nmax, to, po, ro, d_omega_t[ik]);
    }
  }

  if (MyCom_p->cpuNum() == 0)
  {
    writePValues(pnoise);
  } // if (MyCom_p->cpuNum() == 0)
  
  
}


//-------------------------------------------


void KirchhoffIntegrator::ffpDataOnSphere
(
 const std::complex<double> *pvalues,
 const int *vecLen,
 const std::complex<double> *dpdn_coeff,
 const int nmax
 ) const
{
  
  // This function computes the far-field pattern for
  // (p, dp/dn) specified on a sphere.
  
  SubDomain **SubDomain_p = d_domain_p->getSubDomain();
  
  std::vector<double> xi, ww;
  getQuadrature(xi, ww, 12);
  
  int numInc = d_iod.surfKI.d_angularIncrement;
  int numDir = (numInc/2 + 1) * numInc;
  int numFreq = d_omega_t.size();
  
  std::vector<double> directions(3*numDir, 0.0);
  numDir = 0;
  for (int ib = 0; ib <= numInc/2; ++ib)
  {
    double beta = ib * M_PI * 2.0 / numInc - 0.5 * M_PI;
    double cb = cos(beta), sb = sin(beta);
    for (int ia = 0; ia < numInc; ++ia)
    {
      double alpha = ia * 2.0 * M_PI / numInc;
      directions[3*numDir] = cos(alpha) * cb;
      directions[3*numDir+1] = sin(alpha) * cb;
      directions[3*numDir+2] = sb;
      numDir += 1;
    }
  }
  
  int ffpLen = numDir * numFreq;
  std::vector< std::complex<double> > ffp(ffpLen, std::complex<double>(0.0, 0.0));
  
#ifdef _UH_DEBUG_
  double surfaceCheck = 0.0;
#endif
  
  
  int totNumNodes = 0;
#pragma omp parallel for
  for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  {
    
    int myLen = vecLen[iSub+1] - vecLen[iSub];
    if (myLen == 0)
      continue;
    
    FaceSet &faces = SubDomain_p[iSub]->getFaces();
    SVec<double,3> XX = (*d_X_p)(iSub);
    const std::complex<double> *val = pvalues + numFreq * vecLen[iSub];
    
    for (int j = 0; j < faces.size(); ++j)
    {
      
      Face &jFace = faces[j];
      if (jFace.getCode() != BC_KIRCHHOFF_SURFACE)
        continue;
      //
      double x[3], y[3], z[3];
      for (int jj = 0; jj < 3; ++jj)
      {
        double *MyCoord = XX[jFace[jj]];
        x[jj] = MyCoord[0];
        y[jj] = MyCoord[1];
        z[jj] = MyCoord[2];
      }
      //
      std::vector< std::complex<double> > p(3*numFreq);
      for (int ik = 0; ik < numFreq; ++ik)
      {
        const std::complex<double> *myVal = val + ik * myLen;
        for (int jj = 0; jj < 3; ++jj)
        {
          p[jj + 3*ik] = myVal[d_globalToNodeID[totNumNodes + faces[j][jj]]];
        }
      }
      //
      //
      for (int gp = 0; gp < ww.size(); ++gp)
      {
        
        double xi0 = xi[gp], xi1 = xi[gp+ww.size()], xi2 = xi[gp+2*ww.size()];
        double xp = x[0]*xi0 + x[1]*xi1 + x[2]*xi2;
        double yp = y[0]*xi0 + y[1]*xi1 + y[2]*xi2;
        double zp = z[0]*xi0 + z[1]*xi1 + z[2]*xi2;
        
        double cross = dSigma(&x[0], &y[0], &z[0], xp, yp, zp) * ww[gp];
        
        double myrp = 0.0, tp = 0.0, phip = 0.0;
        cart2sph_x(xp, yp, zp, myrp, tp, phip);
        
#ifdef _UH_DEBUG_
        surfaceCheck += d_R * d_R * cross;
#endif
        
        //
        std::vector< std::complex<double> > p_p(numFreq);
        std::vector< std::complex<double> > dpdn(numFreq);
        for (int ik = 0; ik < numFreq; ++ik)
        {
          dpdn[ik] = evaluateSHS(dpdn_coeff + ik*(nmax+1)*(nmax+1), nmax, tp, phip);
          p_p[ik] = p[3*ik]*xi0 + p[1+3*ik]*xi1 + p[2+3*ik]*xi2;
        }
        //
        
        for (int id = 0; id < numDir; ++id)
        {
          //
          double DirDotX = xp * directions[3*id];
          DirDotX += yp * directions[3*id+1];
          DirDotX += zp * directions[3*id+2];
          //
          double DirDotNu = DirDotX / myrp;
          //
          for (int ik = 0; ik < numFreq; ++ik)
          {
            
            double kappa = d_omega_t[ik];
            //
            std::complex<double> e_mi_k_x_d(cos(kappa*DirDotX), -sin(kappa*DirDotX));
            //
            std::complex<double> cTmp = e_mi_k_x_d * d_R * d_R * cross;
#pragma omp critical
            {
              ffp[id + ik * numDir] += (dpdn[ik]
                     + p_p[ik] * std::complex<double>(0.0, kappa * DirDotNu)) * cTmp;
            }
            
          } // for (int ik = 0; ik < numFreq; ++ik)
          
        } // for (int id = 0; id < numDir; ++id)
        
      } // for (int gp = 0; gp < ww.size(); ++gp)
      
    } // for (int j = 0; j < faces.size(); ++j)
    
    totNumNodes += (*d_X_p)(iSub).size();

  } // for (int iSub = 0; iSub < d_domain_p->getNumLocSub(); ++iSub)
  
  
  Communicator *MyCom_p = d_domain_p->getCommunicator();
  MyCom_p->globalSum(ffpLen, &ffp[0]);
  
  
#ifdef _UH_DEBUG_
  MyCom_p->globalSum(1, &surfaceCheck);
  if (MyCom_p->cpuNum() == 0)
  {
    std::cout << " Surface " << surfaceCheck << " ";
    std::cout << std::abs(surfaceCheck - d_R * d_R * 4*M_PI)/(4*d_R * d_R * M_PI);
    std::cout << std::endl;
  }
#endif
  
  //
  // Scale the computed far-field pattern.
  // The minus sign appears because of the normal orientation.
  //
  double dTmp = - 0.25 / M_PI;
  for (int ii = 0; ii < ffpLen; ++ii)
  {
    ffp[ii] = ffp[ii] * dTmp;
  }
  
  if (MyCom_p->cpuNum() == 0)
  {
    writeFFPValues(ffp);
  }
  
}


//-------------------------------------------


std::complex<double> KirchhoffIntegrator::evaluateSHS
(
 const std::complex<double> *coeff,
 int nmax,
 double theta,
 double phi,
 double r,
 double kappa
 ) const
{
  
  std::complex<double> result(0.0, 0.0);
  
  int count = 0;
  std::vector< std::complex<double> > Yn(2*nmax + 1);
  if (r > 0.0)
  {
    for (int n = 0; n <= nmax; ++n)
    {
      //
      double maxEntry = 0.0;
      for (int in = 0; in <= 2*n; ++in)
      {
        maxEntry = (maxEntry > std::abs(coeff[count+in])) ? maxEntry : std::abs(coeff[count+in]);
      }
      if (maxEntry == 0.0)
      {
        count += 2*n+1;
        continue;
      }
      //
      sphericalHarmonic(n, theta, phi, Yn);
      std::complex<double> hn_r = besselh(n, kappa * r);
      std::complex<double> hn_dR = besselh(n, kappa * d_R);
      //
      for (int m = -n; m <= n; ++m)
      {
        result += coeff[count] * Yn[m + n] * hn_r / hn_dR;
        count += 1;
      }
      //
    }
  }
  else {
    for (int n = 0; n <= nmax; ++n)
    {
      //
      double maxEntry = 0.0;
      for (int in = 0; in <= 2*n; ++in)
      {
        maxEntry = (maxEntry > std::abs(coeff[count+in])) ? maxEntry : std::abs(coeff[count+in]);
      }
      if (maxEntry == 0.0)
      {
        count += 2*n+1;
        continue;
      }
      //
      sphericalHarmonic(n, theta, phi, Yn);
      //
      for (int m = -n; m <= n; ++m)
      {
        result += coeff[count] * Yn[m + n];
        count += 1;
      }
      //
    }
  }
  
  return result;
  
}


//-------------------------------------------


void KirchhoffIntegrator::getQuadrature
(
 std::vector<double> &_xigauss,
 std::vector<double> &w,
 int rule
 ) const
{
  
  if (w.size() != rule)
  {
    w.resize(rule);
  }
  
  if (_xigauss.size() != 3*rule)
  {
    _xigauss.resize(3*rule);
  }
  
  if (rule == 3)
  {
    
    w[0] = 1.0/3.0;
    w[1] = w[0];
    w[2] = w[0];
    
    _xigauss[0] = 4.0/6.0;
    _xigauss[1] = 1.0/6.0;
    _xigauss[2] = 1.0/6.0;
    
    _xigauss[3] = 1.0/6.0;
    _xigauss[4] = 4.0/6.0;
    _xigauss[5] = 1.0/6.0;
    
  }
  else if (rule == 6)
  {
    
    w[0] = 0.109951743655322;
    w[1] = 0.109951743655322;
    w[2] = 0.109951743655322;
    
    w[3] = 0.223381589678011;
    w[4] = 0.223381589678011;
    w[5] = 0.223381589678011;
    
    _xigauss[0] = 0.816847572980459;
    _xigauss[1] = 0.091576213509771;
    _xigauss[2] = 0.091576213509771;
    _xigauss[3] = 0.108103018168070;
    _xigauss[4] = 0.445948490915965;
    _xigauss[5] = 0.445948490915965;
    
    _xigauss[6] = 0.091576213509771;
    _xigauss[7] = 0.816847572980459;
    _xigauss[8] = 0.091576213509771;
    _xigauss[9] = 0.445948490915965;
    _xigauss[10] = 0.108103018168070;
    _xigauss[11] = 0.445948490915965;
    
  }
  else if (rule == 12)
  {
    
    w[0] = 0.050844906370207;
    w[1] = 0.050844906370207;
    w[2] = 0.050844906370207;
    w[3] = 0.116786275726379;
    w[4] = 0.116786275726379;
    w[5] = 0.116786275726379;
    w[6] = 0.082851075618374;
    w[7] = 0.082851075618374;
    w[8] = 0.082851075618374;
    w[9] = 0.082851075618374;
    w[10] = 0.082851075618374;
    w[11] = 0.082851075618374;
    
    _xigauss[ 0] = 0.873821971016996;
    _xigauss[ 1] = 0.063089014491502;
    _xigauss[ 2] = 0.063089014491502;
    _xigauss[ 3] = 0.501426509658179;
    _xigauss[ 4] = 0.249286745170910;
    _xigauss[ 5] = 0.249286745170910;
    _xigauss[ 6] = 0.636502499121399;
    _xigauss[ 7] = 0.636502499121399;
    _xigauss[ 8] = 0.310352451033785;
    _xigauss[ 9] = 0.310352451033785;
    _xigauss[10] = 0.053145049844816;
    _xigauss[11] = 0.053145049844816;
    
    _xigauss[12] = 0.063089014491502;
    _xigauss[13] = 0.873821971016996;
    _xigauss[14] = 0.063089014491502;
    _xigauss[15] = 0.249286745170910;
    _xigauss[16] = 0.501426509658179;
    _xigauss[17] = 0.249286745170910;
    _xigauss[18] = 0.310352451033785;
    _xigauss[19] = 0.053145049844816;
    _xigauss[20] = 0.636502499121399;
    _xigauss[21] = 0.053145049844816;
    _xigauss[22] = 0.636502499121399;
    _xigauss[23] = 0.310352451033785;
    
  }
  
  for (int ig = 2*w.size(); ig < 3*w.size(); ++ig)
  {
    _xigauss[ig] = 1.0 - _xigauss[ig-2*w.size()] - _xigauss[ig-w.size()];
  }
  
}


//-------------------------------------------


std::complex<double> KirchhoffIntegrator::hankel
(
 const int n,
 const double x
 ) const
{
  
#ifdef AEROACOUSTIC
  double Hre, Him;
  
  if (n >= 0)
  {
    Hre = gsl_sf_bessel_Jn(n, x);
    Him = gsl_sf_bessel_Yn(n, x);
  }
  else
  {
    double scale = ((-n) % 2 == 0) ? 1.0 : -1.0;
    Hre = scale * gsl_sf_bessel_Jn(-n, x);
    Him = scale * gsl_sf_bessel_Yn(-n, x);
  }
  
  return std::complex<double>(Hre, Him);
#endif // AEROACOUSTIC

  return std::complex<double>(0.0, 0.0);
  
}


//-------------------------------------------


std::complex<double> KirchhoffIntegrator::hankel_prime
(
 const int n,
 const double x
 ) const
{
  
  std::complex<double> Hnp1 = hankel(n+1, x);
  std::complex<double> Hnm1 = hankel(n-1, x);
  
  return 0.5 * (Hnm1 - Hnp1);
  
}


//-------------------------------------------


std::complex<double> KirchhoffIntegrator::besselh
(
 const int n,
 const double x
 ) const
{
  
  double hre, him;
#ifdef AEROACOUSTIC
  if (n >= 0)
  {
    hre = gsl_sf_bessel_jl(n, x);
    him = gsl_sf_bessel_yl(n, x);
  }
  else
  {
    double scale = ((-n) % 2 == 0) ? 1.0 : -1.0;
    hre = scale * gsl_sf_bessel_jl(-n, x);
    him = scale * gsl_sf_bessel_yl(-n, x);
  }
  
  return std::complex<double>(hre, him);
#endif // AEROACOUSTIC
  return std::complex<double>(0.0,0.0);
}


//-------------------------------------------


std::complex<double> KirchhoffIntegrator::besselh_prime
(
 const int n,
 const double x
 ) const
{
  
  std::complex<double> h_nm1 = besselh(n-1, x);
  std::complex<double> h_np1 = besselh(n+1, x);
  
  std::complex<double> hprime;
  hprime = ((double) n)*h_nm1 - (n+1.0)*h_np1;
  hprime *= 1.0/(2.0*n+1.0);
  
  return hprime;
  
}


//-------------------------------------------


void KirchhoffIntegrator::sphericalHarmonic
(
 const int n,
 const double theta,
 const double phi,
 std::vector< std::complex<double> > &Yn
 ) const
{
  
  //
  // UH (08/30/2012)
  // This routine has been verified against MATLAB.
  //
#ifdef AEROACOUSTIC
  if (Yn.size() < 2*n+1)
  {
    Yn.resize(2*n+1);
  }
  
  double cost = cos(theta);
  
  std::vector<double> Pn(n+1, 0.0);
  for (int im = 0; im <= n; ++im)
  {
    Pn[im] = gsl_sf_legendre_sphPlm(n, im, cost);
  }
  
  for (int jj = 0; jj <= n; ++jj)
  {
    Yn[n + jj] = Pn[jj] * std::complex<double>(cos(jj*phi), sin(jj*phi));
  }
  
  double coeff = -1.0;
  for (int jj = 1; jj <= n; ++jj)
  {
    Yn[n - jj] = coeff * conj(Yn[n + jj]);
    coeff *= -1.0;
  }
  
#endif // AEROACOUSTIC
}


//-------------------------------------------


void KirchhoffIntegrator::cart2sph_x
(
 double xo, double yo, double zo,
 double &ro, double &to, double &po
 ) const
{
  
  // UH (09/2012)
  // This function returns spherical coordinates
  // theta is between 0 and pi
  // phi is between 0 and 2*pi
  
  ro = sqrt(xo*xo + yo*yo + zo*zo);
  to = (ro == 0.0) ? 0.0 : acos(zo/ro);
  
  if ((std::abs(xo) < 1.0e-15 * ro) && (std::abs(yo) < 1.0e-15 * ro))
  {
    po = 0.0;
  }
  else if (std::abs(xo) < 1.0e-15 * ro)
  {
    po = (yo < 0.0) ? M_PI * 1.5 : M_PI * 0.5;
  }
  else if (std::abs(yo) < 1.0e-15 * ro)
  {
    po = (xo < 0.0) ? M_PI : 0.0;
  }
  else
  {
    double angle = xo/(ro * sin(to));
    if (std::abs(angle) >= 1.0)
    {
      angle = (angle > 0.0) ? 0.0 : M_PI;
    }
    else
    {
      angle = acos(angle);
    }
    po = (yo < 0.0) ? 2.0*M_PI - angle : angle;
  }
  
}


//-------------------------------------------


void KirchhoffIntegrator::cart2sph_y
(
 double xo, double yo, double zo,
 double &ro, double &to, double &po
 ) const
{
  
  // UH (09/2012)
  // This function returns spherical coordinates
  // theta is between 0 and pi
  // phi is between -pi/2 and 3*pi/2
  
  ro = sqrt(xo*xo + yo*yo + zo*zo);
  if (ro == 0.0)
  {
    to = 0.0;
  }
  else
  {
    to = acos(zo/ro);
  }
  
  if ((std::abs(xo) < 1.0e-15 * ro) && (std::abs(yo) < 1.0e-15 * ro))
  {
    po = 0.0;
  }
  else if (std::abs(xo) < 1.0e-15 * ro)
  {
    po = (yo < 0.0) ? -M_PI * 0.5 : M_PI * 0.5;
  }
  else if (std::abs(yo) < 1.0e-15 * ro)
  {
    po = (xo < 0.0) ? M_PI : 0.0;
  }
  else
  {
    double angle = yo/(ro * sin(to));
    if (std::abs(angle) >= 1.0)
    {
      angle = (angle > 0.0) ? 0.5*M_PI : 1.5*M_PI;
    }
    else
    {
      angle = asin(angle);
    }
    po = (xo > 0.0) ? angle : M_PI - angle;
  }
  
}


