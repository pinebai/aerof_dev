#include <TimeState.h>

#include <TimeData.h>
#include <GeoState.h>
#include <Vector.h>
#include <GenMatrix.h>
#include <BcDef.h>
#include <LowMachPrec.h>

//------------------------------------------------------------------------------

template<int dim>
TimeState<dim>::TimeState(TimeData &_data, Vec<double> &_dt, Vec<double> &_idti, 
			  Vec<double> &_idtv, Vec<double> &_dtau, SVec<double,dim> &_Un,
			  SVec<double,dim> &_Unm1, SVec<double,dim> &_Unm2, 
			  SVec<double,dim> &_Rn) : 
  data(_data), dt(_dt), idti(_idti), idtv(_idtv), dtau(_dtau), Un(_Un), Unm1(_Unm1), Unm2(_Unm2), Rn(_Rn)
{
  if (data.use_modal || data.descriptor_form == 1) {
    descriptorCase = DESCRIPTOR;
  } else if (data.descriptor_form == 2) {
    descriptorCase = HYBRID;
  } else { 
    descriptorCase = NONDESCRIPTOR;
  }
}

template<int dim>
void TimeState<dim>::attachHH(Vec<double> *hhn_,Vec<double> *hhnm1_) {

  hhn = hhn_;
  hhnm1 = hhnm1_;  
}

//------------------------------------------------------------------------------
// Add Time Derivative Term, d(AW)/dt to the flux F
// If running Non-Modal, adds invA*d(AW)/dt to the flux invA*F

template<int dim>
void TimeState<dim>::add_dAW_dt(bool *nodeFlag, GeoState &geoState, 
				Vec<double> &ctrlVol, SVec<double,dim> &Q, 
				SVec<double,dim> &R, LevelSetStructure *LSS)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

	for(int i=0; i<dt.size(); ++i) 
	{
		if(LSS && !(LSS->isActive(0.0,i))) 
		{
      // Node i lies in the structure: Do nothing.
      continue;
    }

    TimeFDCoefs coefs;
    computeTimeFDCoefs(geoState, coefs, ctrlVol, i);

    double invDt = 1.0 / dt[i];

		for(int k=0; k<dim; ++k) 
		{
			double sum = coefs.c_np1 *    Q[i][k] 
				        + coefs.c_n   *   Un[i][k] 
				        + coefs.c_nm1 * Unm1[i][k] 
				        + coefs.c_nm2 * Unm2[i][k];

      double dAWdt = invDt * sum; 
      if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON)
        R[i][k] = dAWdt + 0.5 * (R[i][k] + Rn[i][k]);
      else
        R[i][k] += dAWdt;
    }
  }
}

template<int dim>
void TimeState<dim>::add_dAW_dt_HH(bool *nodeFlag, GeoState &geoState, 
				   Vec<double> &ctrlVol, Vec<double> &Q, 
				   Vec<double> &R)
{
  for (int i=0; i<R.size(); ++i) {

    TimeFDCoefs coefs;
    // we require that non-descriptor be used
    assert(descriptorCase == NONDESCRIPTOR);
    computeTimeFDCoefs(geoState, coefs, ctrlVol, 0);

    double invDt = 1.0 / dt[0];
    double sum = coefs.c_np1*Q[i]+ coefs.c_n*(*hhn)[i] +
      coefs.c_nm1*(*hhnm1)[i];// + coefs.c_nm2*Unm2[i][k];
    double dAWdt = invDt * sum; 
    //if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON)
    //  R[i][k] = dAWdt + 0.5 * (R[i][k] + Rn[i][k]);
    //else
    R[i] += dAWdt;
  }
}

//------------------------------------------------------------------------------
// Add Time Derivative Term, P^-1 x d(AW)/dt to the flux F, 
// where P is the turkel preconditioning matrix.
// If running Non-Modal, adds invA*d(P^-1 x AW)/dt to the flux invA*F

template<int dim>
void TimeState<dim>::add_GASPrec_dAW_dt(bool *nodeFlag, GeoState &geoState, 
				Vec<double> &ctrlVol, SVec<double,dim> &Q, 
				SVec<double,dim> &R, double gam, double pstiff, Vec<double> &irey, TimeLowMachPrec &tprec, LevelSetStructure *LSS)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  for (int i=0; i<dt.size(); ++i) {
    if (LSS && !(LSS->isActive(0.0,i))) {
      // Node i lies in the structure: Do nothing.
      continue;
    }

    TimeFDCoefs coefs;
    computeTimeFDCoefs(geoState, coefs, ctrlVol, i);

    double invDt = 1.0 / dt[i];

    if (dim<5) {
      for (int k=0; k<dim; ++k) {
        double sum = coefs.c_np1*Q[i][k] + coefs.c_n*Un[i][k] +
                     coefs.c_nm1*Unm1[i][k] + coefs.c_nm2*Unm2[i][k];
        double dAWdt = invDt * sum; 
        if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON)
          R[i][k] = dAWdt + 0.5 * (R[i][k] + Rn[i][k]);
        else
          R[i][k] += dAWdt;
      }
    }
    else {
      double ro = Un[i][0];
      double invRho = 1.0/ro;
      double u  = Un[i][1] * invRho;
      double v  = Un[i][2] * invRho;
      double w  = Un[i][3] * invRho;
      double u2 = u*u;
      double v2 = v*v;
      double w2 = w*w;
      double q2 = u2 + v2 + w2;
      double gam1 = gam - 1.0;
      double p  = gam1 * (Un[i][4] - 0.5 * ro * q2) - gam*pstiff;
      double c2 = gam*(p+pstiff)/ro;
      double locMach = sqrt(q2/c2); //local Preconditioning (ARL)
      double beta = tprec.getBeta(locMach, irey[i]);

      double beta2 =   beta * beta;
      double qhat2 = (q2 * gam1)/2.0;
 
      double nu = qhat2/c2;
      double mu = (1.0/beta2) - 1.0;

      double Pinv[5][5] = { {nu*mu + 1.0,  -u*mu*gam1/c2,      -v*mu*gam1/c2,        -w*mu*gam1/c2,       mu*gam1/c2   },
                            {u*nu*mu,     1.0 - u2*mu*gam1/c2, -u*v*mu*gam1/c2,      -u*w*mu*gam1/c2,     u*mu*gam1/c2 },
                            {v*nu*mu,     -u*v*mu*gam1/c2 ,    1.0 - v2*mu*gam1/c2,  -v*w*mu*gam1/c2,     v*mu*gam1/c2 },
                            {w*nu*mu,     -u*w*mu*gam1/c2 ,    -v*w*mu*gam1/c2,      1.0 - w2*mu*gam1/c2, w*mu*gam1/c2 },
       	                  {0.5*mu*(1.0+nu)*q2,    -u*mu*(1+nu), -v*mu*(1+nu), -w*mu*(1+nu), (1.0/beta2)+mu*nu } };


      double dAWdt[dim],Pinv_dAWdt[dim];

      for (int k=0; k<dim; ++k) {
        double sum = coefs.c_np1*Q[i][k] + coefs.c_n*Un[i][k] +
                     coefs.c_nm1*Unm1[i][k] + coefs.c_nm2*Unm2[i][k];
        dAWdt[k] = invDt * sum; 
      }

      for (int k=0; k<5; ++k) {
        Pinv_dAWdt[k] = 0.0;
        for (int l =0; l<5; ++l) {
          Pinv_dAWdt[k] = Pinv_dAWdt[k] + Pinv[k][l]*dAWdt[l];
        }
      }

      //turbulence preconditioning
      if (dim == 6) {
        double t1 = Un[i][5] * invRho;
        double mup = mu*t1*gam1/c2;
        double Pt[6] = {mu*nu*t1, -mup*u, -mup*v, -mup*w, mup, 1.0};
        for (int k=5; k<dim; ++k) {
          Pinv_dAWdt[k] = 0.0;
          for (int l =0; l<dim; ++l) {
            Pinv_dAWdt[k] = Pinv_dAWdt[k] + Pt[l]*dAWdt[l];
          }
        }
      }
      else if (dim == 7) {
        double t1 = Un[i][5] * invRho;
        double t2 = Un[i][6] * invRho;
        double mup1 = mu*t1*gam1/c2;
        double mup2 = mu*t2*gam1/c2;
        double Pt[2][7] = { {mu*nu*t1, -mup1*u, -mup1*v, -mup1*w, mup1, 1.0, 0.0},
                            {mu*nu*t2, -mup2*u, -mup2*v, -mup2*w, mup2, 0.0, 1.0} };
        for (int k=5; k<dim; ++k) {
          Pinv_dAWdt[k] = 0.0;
          for (int l =0; l<dim; ++l) {
            Pinv_dAWdt[k] = Pinv_dAWdt[k] + Pt[k-5][l]*dAWdt[l];
          }
        }
      }

      for (int k=0; k<dim; ++k) {
        if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON)
          R[i][k] = Pinv_dAWdt[k] + 0.5 * (R[i][k] + Rn[i][k]);
        else
          R[i][k] += Pinv_dAWdt[k];
      }
    }
  }
}

//------------------------------------------------------------------------------

// Add Time Derivative Term, P^-1 x d(AW)/dt to the flux F, 
// where P is the preconditioning matrix.
// If running Non-Modal, adds invA*d(P^-1 x AW)/dt to the flux invA*F

template<int dim>
void TimeState<dim>::add_LiquidPrec_dAW_dt(bool *nodeFlag, GeoState &geoState, 
				Vec<double> &ctrlVol, VarFcn *vf, SVec<double,dim> &Q, 
				SVec<double,dim> &R, Vec<double> &irey, 
                                TimeLowMachPrec &tprec, LevelSetStructure *LSS)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  for (int i=0; i<dt.size(); ++i) {
    if (LSS && !(LSS->isActive(0.0,i))) {
      // Node i lies in the structure: Do nothing.
      continue;
    }

    TimeFDCoefs coefs;
    computeTimeFDCoefs(geoState, coefs, ctrlVol, i);

    double invDt = 1.0 / dt[i];

    double dAWdt[dim];

    for (int k=0; k<dim; ++k) {
      double sum = coefs.c_np1*Q[i][k] + coefs.c_n*Un[i][k] +
                   coefs.c_nm1*Unm1[i][k] + coefs.c_nm2*Unm2[i][k];
      dAWdt[k] = invDt * sum; 
    }
    double V[dim];
    vf->conservativeToPrimitive(Un[i],V); // assumption : no steady two-phase flow, hence no phi
    double e = vf->computeRhoEnergy(V)/V[0];
    double pressure = vf->getPressure(V);
    double c = vf->computeSoundSpeed(V);
    double c2 = c*c;
    double locMach = vf->computeMachNumber(V); //local Preconditioning (ARL)
    double beta = tprec.getBeta(locMach,irey[i]);
    double oobeta2 = 1.0/(beta*beta);
    double oobeta2m1 = oobeta2 - 1.0;

    double Pinv[dim];
    for (int j=0; j<dim; j++)
      Pinv[j] = oobeta2m1*V[j];
    Pinv[0] = oobeta2;
    Pinv[4] = oobeta2m1*(e+pressure/V[0] - c2);
 
    for (int k=0; k<dim; ++k) {
      if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON)
        R[i][k] = dAWdt[k] + Pinv[k]*dAWdt[0] + 0.5 * (R[i][k] + Rn[i][k]);
      else
        R[i][k] += dAWdt[k] + Pinv[k]*dAWdt[0];
    }
  }
}

//------------------------------------------------------------------------------
template<int dim>
void TimeState<dim>::add_dAW_dtRestrict(bool *nodeFlag, GeoState &geoState, 
					Vec<double> &ctrlVol, SVec<double,dim> &Q, 
					SVec<double,dim> &R, const std::vector<int> &sampledLocNodes)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  int i;
  for (int iSampledNode=0; iSampledNode<sampledLocNodes.size(); ++iSampledNode) {
    i = sampledLocNodes[iSampledNode];

    TimeFDCoefs coefs;
    computeTimeFDCoefs(geoState, coefs, ctrlVol, i);


    double invDt = 1.0 / dt[i];
    for (int k=0; k<dim; ++k) {
      double sum = coefs.c_np1*Q[i][k] + coefs.c_n*Un[i][k] +
                   coefs.c_nm1*Unm1[i][k] + coefs.c_nm2*Unm2[i][k];
      double dAWdt = invDt * sum;
      if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON)
        R[i][k] = dAWdt + 0.5 * (R[i][k] + Rn[i][k]);
      else
        R[i][k] += dAWdt;
    }
  }
}
//------------------------------------------------------------------------------
// Add Time Derivative Term, d(AW)/dt to the flux F
// If running Non-Modal, adds invA*d(AW)/dt to the flux invA*F

template<int dim>
template<int dimLS>
void TimeState<dim>::add_dAW_dtLS(bool *nodeFlag, GeoState &geoState, 
					Vec<double> &ctrlVol, SVec<double,dimLS> &Q, 
					SVec<double,dimLS> &Qn, SVec<double,dimLS> &Qnm1,
					SVec<double,dimLS> &Qnm2, SVec<double,dimLS> &R,bool requireSpecialBDF)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();
  Vec<double>& ctrlVol_nm1 = geoState.getCtrlVol_nm1();
  Vec<double>& ctrlVol_nm2 = geoState.getCtrlVol_nm2();

  for (int i=0; i<dt.size(); ++i) {

    TimeFDCoefs coefs;
    if (requireSpecialBDF)
      computeTimeFDCoefsSpecialBDF(geoState, coefs, ctrlVol, i);
    else  
      computeTimeFDCoefs(geoState, coefs, ctrlVol, i);

    double invDt = 1.0 / dt[i];

    for (int idim=0; idim<dimLS; idim++){
      double dAWdt;
      if (!requireSpecialBDF) {
        dAWdt = invDt * (coefs.c_np1*Q[i][idim] + coefs.c_n*Qn[i][idim] 
             + coefs.c_nm1*Qnm1[i][idim] + coefs.c_nm2*Qnm2[i][idim]);
      }
      else {
        dAWdt = invDt * ( coefs.c_np1*Q[i][idim] + coefs.c_n*Qn[i][idim]) 
              + coefs.c_nm1*Qnm1[i][idim]; // here Qnm1 is dQdt^n, so no multiplication by invDt
      }

      if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON)
        R[i][idim] = dAWdt + 0.5 * (R[i][idim] + Rn[i][idim]);
      else
        R[i][idim] += dAWdt;
   }

  }
}

//------------------------------------------------------------------------------
template<int dim>
void TimeState<dim>::computeTimeFDCoefs(GeoState &geoState, TimeFDCoefs &coefs, Vec<double> &ctrlVol, int i) 
{ 

  Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();
  Vec<double>& ctrlVol_nm1 = geoState.getCtrlVol_nm1();
  Vec<double>& ctrlVol_nm2 = geoState.getCtrlVol_nm2();

  coefs.c_np1 = data.alpha_np1;
  coefs.c_n   = data.alpha_n * ctrlVol_n[i];
  coefs.c_nm1 = data.alpha_nm1 * ctrlVol_nm1[i];
  coefs.c_nm2 = data.alpha_nm2 * ctrlVol_nm2[i];

	switch (descriptorCase) 
	{
	   case DESCRIPTOR: 
		{
      coefs.c_np1 *= ctrlVol[i];
			break; 
		}
	   case HYBRID:
		{
      double invsqrtCtrlVol = 1.0 / sqrt(ctrlVol[i]);
      coefs.c_np1 *= ctrlVol[i] * invsqrtCtrlVol;
      coefs.c_n   *= invsqrtCtrlVol;
      coefs.c_nm1 *= invsqrtCtrlVol;
      coefs.c_nm2 *= invsqrtCtrlVol;
			break; 
		}
	   case NONDESCRIPTOR: 
		{
      double invCtrlVol = 1.0 / ctrlVol[i];
      coefs.c_n   *= invCtrlVol;
      coefs.c_nm1 *= invCtrlVol;
      coefs.c_nm2 *= invCtrlVol;
			break; 
  }
}

}
//------------------------------------------------------------------------------
template<int dim>
void TimeState<dim>::computeTimeFDCoefsSpecialBDF(GeoState &geoState, TimeFDCoefs &coefs, Vec<double> &ctrlVol, int i) 
{

  Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();
  Vec<double>& ctrlVol_nm1 = geoState.getCtrlVol_nm1();
  Vec<double>& ctrlVol_nm2 = geoState.getCtrlVol_nm2();


  coefs.c_np1 = 2.0;
  coefs.c_n   = -2.0 * ctrlVol_n[i];
  coefs.c_nm1 = -1.0 * ctrlVol_nm1[i];
  coefs.c_nm2 = 0.0;

	switch (descriptorCase) 
	{
	   case DESCRIPTOR: 
		{
      coefs.c_np1 *= ctrlVol[i];
			break; 
		}
  	   case HYBRID:
		{
      double invsqrtCtrlVol = 1.0 / sqrt(ctrlVol[i]);
      coefs.c_np1 *= ctrlVol[i] * invsqrtCtrlVol;
      coefs.c_n   *= invsqrtCtrlVol;
      coefs.c_nm1 *= invsqrtCtrlVol;
			break; 
		}
      case NONDESCRIPTOR: 
		{
      double invCtrlVol = 1.0 / ctrlVol[i];
      coefs.c_n   *= invCtrlVol;
      coefs.c_nm1 *= invCtrlVol;
			break; 
  }
  }

}
//-----------------------------------------------------------------------------

// Add Dual-Time Derivative Term, d(AW)/dtau
// If running Non-Modal, adds invA*d(AW)/dtau

template<int dim>
void TimeState<dim>::add_dAW_dtau(bool *nodeFlag, GeoState &geoState, 
				Vec<double> &ctrlVol, SVec<double,dim> &Q, 
				SVec<double,dim> &R, LevelSetStructure *LSS)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  for (int i=0; i<dt.size(); ++i) {
    if (LSS && !(LSS->isActive(0.0,i))) {
      // Node i lies in the structure: Do nothing.
      continue;
    }

    TimeFDCoefs coefs;
    computeDualTimeFDCoefs(geoState, coefs, ctrlVol, i);

    double invDtau = 1.0 / dtau[i];
    for (int k=0; k<dim; ++k) {
      double sum = coefs.c_np1*Q[i][k] + coefs.c_n*Un[i][k];
      double dAWdtau = invDtau * sum; 
      R[i][k] += dAWdtau;
    }
  }
}

//------------------------------------------------------------------------------
// Add Dual-Time Derivative Term, P^-1 x d(AW)/dtau, 
// where P is the turkel preconditioning matrix.
// If running Non-Modal, adds invA*d(P^-1 x AW)/dtau

template<int dim>
void TimeState<dim>::add_GASPrec_dAW_dtau(bool *nodeFlag, GeoState &geoState, 
				Vec<double> &ctrlVol, SVec<double,dim> &Q, 
				SVec<double,dim> &R, double gam, double pstiff, Vec<double> &irey, TimeLowMachPrec &tprec, LevelSetStructure *LSS)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  for (int i=0; i<dt.size(); ++i) {
    if (LSS && !(LSS->isActive(0.0,i))) {
      // Node i lies in the structure: Do nothing.
      continue;
    }

    TimeFDCoefs coefs;
    computeDualTimeFDCoefs(geoState, coefs, ctrlVol, i);

    double invDtau = 1.0 / dtau[i];

    if (dim<5) {
      for (int k=0; k<dim; ++k) {
        double sum = coefs.c_np1*Q[i][k] + coefs.c_n*Un[i][k];
        double dAWdtau = invDtau * sum; 
        R[i][k] += dAWdtau;
      }
    }
    else {
      double ro = Un[i][0];
      double invRho = 1.0/ro;
      double u  = Un[i][1] * invRho;
      double v  = Un[i][2] * invRho;
      double w  = Un[i][3] * invRho;
      double u2 = u*u;
      double v2 = v*v;
      double w2 = w*w;
      double q2 = u2 + v2 + w2;
      double gam1 = gam - 1.0;
      double p  = gam1 * (Un[i][4] - 0.5 * ro * q2) - gam*pstiff;
      double c2 = gam*(p+pstiff)/ro;
      double locMach = sqrt(q2/c2); //local Preconditioning (ARL)
      double beta = tprec.getBeta(locMach, irey[i]);

      double beta2 =   beta * beta;
      double qhat2 = (q2 * gam1)/2.0;
 
      double nu = qhat2/c2;
      double mu = (1.0/beta2) - 1.0;

      double Pinv[5][5] = { {nu*mu + 1.0,  -u*mu*gam1/c2,      -v*mu*gam1/c2,        -w*mu*gam1/c2,       mu*gam1/c2   },
                            {u*nu*mu,     1.0 - u2*mu*gam1/c2, -u*v*mu*gam1/c2,      -u*w*mu*gam1/c2,     u*mu*gam1/c2 },
                            {v*nu*mu,     -u*v*mu*gam1/c2 ,    1.0 - v2*mu*gam1/c2,  -v*w*mu*gam1/c2,     v*mu*gam1/c2 },
                            {w*nu*mu,     -u*w*mu*gam1/c2 ,    -v*w*mu*gam1/c2,      1.0 - w2*mu*gam1/c2, w*mu*gam1/c2 },
       	                  {0.5*mu*(1.0+nu)*q2,    -u*mu*(1+nu), -v*mu*(1+nu), -w*mu*(1+nu), (1.0/beta2)+mu*nu } };


      double dAWdtau[dim],Pinv_dAWdtau[dim];

      for (int k=0; k<dim; ++k) {
        double sum = coefs.c_np1*Q[i][k] + coefs.c_n*Un[i][k];
        dAWdtau[k] = invDtau * sum; 
      }

      for (int k=0; k<5; ++k) {
        Pinv_dAWdtau[k] = 0.0;
        for (int l =0; l<5; ++l) {
          Pinv_dAWdtau[k] = Pinv_dAWdtau[k] + Pinv[k][l]*dAWdtau[l];
        }
      }

      //turbulence preconditioning
      if (dim == 6) {
        double t1 = Un[i][5] * invRho;
        double mup = mu*t1*gam1/c2;
        double Pt[6] = {mu*nu*t1, -mup*u, -mup*v, -mup*w, mup, 1.0};
        for (int k=5; k<dim; ++k) {
          Pinv_dAWdtau[k] = 0.0;
          for (int l =0; l<dim; ++l) {
            Pinv_dAWdtau[k] = Pinv_dAWdtau[k] + Pt[l]*dAWdtau[l];
          }
        }
      }
      else if (dim == 7) {
        double t1 = Un[i][5] * invRho;
        double t2 = Un[i][6] * invRho;
        double mup1 = mu*t1*gam1/c2;
        double mup2 = mu*t2*gam1/c2;
        double Pt[2][7] = { {mu*nu*t1, -mup1*u, -mup1*v, -mup1*w, mup1, 1.0, 0.0},
                            {mu*nu*t2, -mup2*u, -mup2*v, -mup2*w, mup2, 0.0, 1.0} };
        for (int k=5; k<dim; ++k) {
          Pinv_dAWdtau[k] = 0.0;
          for (int l =0; l<dim; ++l) {
            Pinv_dAWdtau[k] = Pinv_dAWdtau[k] + Pt[k-5][l]*dAWdtau[l];
          }
        }
      }

      for (int k=0; k<dim; ++k) {
        R[i][k] += Pinv_dAWdtau[k];
      }
    }
  }
}

//------------------------------------------------------------------------------

// Add Time Derivative Term, P^-1 x d(AW)/dtau
// where P is the preconditioning matrix.
// If running Non-Modal, adds invA*d(P^-1 x AW)/dtau

template<int dim>
void TimeState<dim>::add_LiquidPrec_dAW_dtau(bool *nodeFlag, GeoState &geoState, 
				Vec<double> &ctrlVol, VarFcn *vf, SVec<double,dim> &Q, 
				SVec<double,dim> &R, Vec<double> &irey, 
                                TimeLowMachPrec &tprec, LevelSetStructure *LSS)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  for (int i=0; i<dt.size(); ++i) {
    if (LSS && !(LSS->isActive(0.0,i))) {
      // Node i lies in the structure: Do nothing.
      continue;
    }

    TimeFDCoefs coefs;
    computeDualTimeFDCoefs(geoState, coefs, ctrlVol, i);

    double invDtau = 1.0 / dtau[i];
    double dAWdtau[dim];

    for (int k=0; k<dim; ++k) {
      double sum = coefs.c_np1*Q[i][k] + coefs.c_n*Un[i][k];
      dAWdtau[k] = invDtau * sum; 
    }
    double V[dim];
    vf->conservativeToPrimitive(Un[i],V); // assumption : no steady two-phase flow, hence no phi
    double e = vf->computeRhoEnergy(V)/V[0];
    double pressure = vf->getPressure(V);
    double c = vf->computeSoundSpeed(V);
    double c2 = c*c;
    double locMach = vf->computeMachNumber(V); //local Preconditioning (ARL)
    double beta = tprec.getBeta(locMach,irey[i]);
    double oobeta2 = 1.0/(beta*beta);
    double oobeta2m1 = oobeta2 - 1.0;

    double Pinv[dim];
    for (int j=0; j<dim; j++)
      Pinv[j] = oobeta2m1*V[j];
    Pinv[0] = oobeta2;
    Pinv[4] = oobeta2m1*(e+pressure/V[0] - c2);
 
    for (int k=0; k<dim; ++k) {
      R[i][k] += dAWdtau[k] + Pinv[k]*dAWdtau[0];
    }
  }
}

//------------------------------------------------------------------------------
template<int dim>
void TimeState<dim>::computeDualTimeFDCoefs(GeoState &geoState, TimeFDCoefs &coefs, Vec<double> &ctrlVol, int i) { 
  Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();

  coefs.c_np1 = 1.0;
  coefs.c_n   = -1.0 * ctrlVol_n[i];
  switch (descriptorCase) {
    case DESCRIPTOR: {
      coefs.c_np1 *= ctrlVol[i];
      break; }
    case HYBRID:{
      double invsqrtCtrlVol = 1.0 / sqrt(ctrlVol[i]);
      coefs.c_np1 *= ctrlVol[i] * invsqrtCtrlVol;
      coefs.c_n   *= invsqrtCtrlVol;
      break; }
    case NONDESCRIPTOR: {
      double invCtrlVol = 1.0 / ctrlVol[i];
      coefs.c_n   *= invCtrlVol;
      break; }
  }
}
//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianNoPrec(bool *nodeFlag, Vec<double> &ctrlVol, GenMat<Scalar,neq> &A,
                                   SVec<double,dim> &U, VarFcn *vf, int* nodeType)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;
  if (data.use_modal == true && data.use_freq == false) A *= 0.5;

  if(!nodeType){
    for (int i=0; i<dt.size(); ++i) 
      addToJacobianNoPrecLocal(i, ctrlVol[i], U, A,i);

  }else{
    for (int i=0; i<dt.size(); ++i) 
      if(!(nodeType[i]==BC_INLET_MOVING || nodeType[i]==BC_OUTLET_MOVING ||
           nodeType[i]==BC_INLET_FIXED || nodeType[i]==BC_OUTLET_FIXED ||
           nodeType[i]==BC_DIRECTSTATE_INLET_MOVING || nodeType[i]==BC_DIRECTSTATE_OUTLET_MOVING ||
           nodeType[i]==BC_DIRECTSTATE_INLET_FIXED || nodeType[i]==BC_DIRECTSTATE_OUTLET_FIXED ||
           nodeType[i]==BC_MASSFLOW_INLET_MOVING || nodeType[i]==BC_MASSFLOW_OUTLET_MOVING ||
           nodeType[i]==BC_MASSFLOW_INLET_FIXED || nodeType[i]==BC_MASSFLOW_OUTLET_FIXED) )
        addToJacobianNoPrecLocal(i, ctrlVol[i], U, A,i);
  }

}

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianLS(bool* nodeFlag,Vec<double> &ctrlVol, GenMat<Scalar,neq> &A,
                                     SVec<double,dim> &U, bool requireSpecialBDF) {

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  double c_np1;
  double coef;
  for (int i=0; i<dt.size(); ++i)  {
    if (requireSpecialBDF)
      coef = 2.0;
    else
      coef = data.alpha_np1; 

   switch (descriptorCase) {
     case DESCRIPTOR: {   
       c_np1 = coef * ctrlVol[i] / dt[i];
       break; }
     case HYBRID: {
       c_np1 = coef * sqrt(ctrlVol[i]) / dt[i];
       break; }
     case NONDESCRIPTOR: {
       c_np1 = coef / dt[i];
       break; }
    }

    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq; ++k)
      Aii[k + k*neq] += c_np1;
  }
}

template<int dim>
template<class Scalar,int neq>
void TimeState<dim>::addToJacobianHH(Vec<double>& ctrlVol,  GenMat<Scalar,neq> & A,
				     Vec<double>& hhn) {

  double c_np1;
  for (int i = 0; i < hhn.size(); ++i) {
  
    Scalar *Aii = A.getElemHH(i);
    switch (descriptorCase) {
      /*case DESCRIPTOR: {
      c_np1 = data.alpha_np1 * vol / dt[dt_i];
      break; }
    case HYBRID: {
      c_np1 = data.alpha_np1 * sqrt(vol) / dt[dt_i];
      break; }*/
    case NONDESCRIPTOR : { 
      c_np1 = data.alpha_np1 / dt[0];
      break; }
    }
    *Aii += c_np1;
  }

}

//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianNoPrecLocal(int i, double vol, 
					SVec<double,dim> &U, GenMat<Scalar,neq> &A,
                                        int dt_i)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  double c_np1;

  switch (descriptorCase) {
    case DESCRIPTOR: {
      c_np1 = data.alpha_np1 * vol / dt[dt_i];
      c_np1 += data.dtau_switch * vol / dtau[i];
      break; }
    case HYBRID: {
      c_np1 = data.alpha_np1 * sqrt(vol) / dt[dt_i];
      c_np1 += data.dtau_switch * sqrt(vol) / dtau[i];
      break; }
    case NONDESCRIPTOR : { 
      c_np1 = data.alpha_np1 / dt[dt_i];
      c_np1 += data.dtau_switch * 1.0 / dtau[i];
      break; }
  }
  Scalar *Aii = A.getElem_ii(i);
  for (int k=0; k<neq; ++k)
    Aii[k + k*neq] += c_np1;
}
//------------------------------------------------------------------------------
//  This Part Performs the Low Mach Steady State Preconditioning
//  Reference: Preconditioning Methods for Low-Speed Flows
//             By Turkel et. al. (ICASE Publication)

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianGasPrec(bool *nodeFlag, Vec<double> &ctrlVol, GenMat<Scalar,neq> &A,
                                   SVec<double,dim> &U, VarFcn *vf, double gam, 
				   double pstiff, TimeLowMachPrec &tprec,
				   Vec<double> &irey, int* nodeType)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;
  if (data.use_modal == true && data.use_freq == false) A *= 0.5;
                                                                                                                           
  double c_np1;
  if(!nodeType){
    for (int i=0; i<dt.size(); ++i) 
      addToJacobianGasPrecLocal(i,ctrlVol[i],gam,pstiff,tprec,irey[i],U,A);

  }else{
    for (int i=0; i<dt.size(); ++i) 
      if(!(nodeType[i]==BC_INLET_MOVING || nodeType[i]==BC_OUTLET_MOVING ||
           nodeType[i]==BC_INLET_FIXED || nodeType[i]==BC_OUTLET_FIXED ||
           nodeType[i]==BC_DIRECTSTATE_INLET_MOVING || nodeType[i]==BC_DIRECTSTATE_OUTLET_MOVING ||
           nodeType[i]==BC_DIRECTSTATE_INLET_FIXED || nodeType[i]==BC_DIRECTSTATE_OUTLET_FIXED ||
           nodeType[i]==BC_MASSFLOW_INLET_MOVING || nodeType[i]==BC_MASSFLOW_OUTLET_MOVING ||
           nodeType[i]==BC_MASSFLOW_INLET_FIXED || nodeType[i]==BC_MASSFLOW_OUTLET_FIXED) )
        addToJacobianGasPrecLocal(i,ctrlVol[i],gam,pstiff,tprec,irey[i],U,A);

  }
}
//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianGasPrecLocal(int i, double vol, double gam, 
				double pstiff, TimeLowMachPrec &tprec,
				double irey, SVec<double,dim> &U, GenMat<Scalar,neq> &A)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  double c_np1;

  switch (descriptorCase) {
    case DESCRIPTOR: { 
      c_np1 = data.alpha_np1 * vol / dt[i];
      c_np1 += data.dtau_switch * vol / dtau[i];
      break; }
    case HYBRID: {
      c_np1 = data.alpha_np1 * sqrt(vol) / dt[i];
      c_np1 += data.dtau_switch * sqrt(vol) / dtau[i];
      break; }
    case NONDESCRIPTOR: {
      c_np1 = data.alpha_np1 / dt[i];
      c_np1 += data.dtau_switch * 1.0 / dtau[i];
      break; }
  }

  Scalar *Aii = A.getElem_ii(i);

  if(neq<5){		//turbulence model equation in segregated solver
    for (int k=0; k<neq; ++k)
      Aii[k + k*neq] += c_np1;

  }else{	//Navier-Stokes (part of segregated turb model or alone) or fully coupled

    double ro = Un[i][0];
    double invRho = 1.0/ro;
    double u  = Un[i][1] * invRho;
    double v  = Un[i][2] * invRho;
    double w  = Un[i][3] * invRho;
    double u2 = u*u;
    double v2 = v*v;
    double w2 = w*w;
    double q2 = u2 + v2 + w2;
    double gam1 = gam - 1.0;
    double p  = gam1 * (Un[i][4] - 0.5 * ro * q2) - gam*pstiff;
    double c2 = gam*(p+pstiff)/ro;
    double locMach = sqrt(q2/c2); //local Preconditioning (ARL)
    double beta = tprec.getBeta(locMach, irey);

    double beta2 =   beta * beta;

// Preconditioning in unsteady flow
    double bt =  1.0 /(1.0 + data.dtau_switch*dtau[i]*data.alpha_np1/dt[i]);
    beta2 = beta2/(bt - beta2*(bt - 1.0));

    double qhat2 = (q2 * gam1)/2.0;
 
    double nu = qhat2/c2;
    double mu = (1.0/beta2) - 1.0;

    double Pinv[5][5] = { {nu*mu + 1.0,  -u*mu*gam1/c2,      -v*mu*gam1/c2,        -w*mu*gam1/c2,       mu*gam1/c2   },
                          {u*nu*mu,     1.0 - u2*mu*gam1/c2, -u*v*mu*gam1/c2,      -u*w*mu*gam1/c2,     u*mu*gam1/c2 },
                          {v*nu*mu,     -u*v*mu*gam1/c2 ,    1.0 - v2*mu*gam1/c2,  -v*w*mu*gam1/c2,     v*mu*gam1/c2 },
                          {w*nu*mu,     -u*w*mu*gam1/c2 ,    -v*w*mu*gam1/c2,      1.0 - w2*mu*gam1/c2, w*mu*gam1/c2 },
     	                  {0.5*mu*(1.0+nu)*q2,    -u*mu*(1+nu), -v*mu*(1+nu), -w*mu*(1+nu), (1.0/beta2)+mu*nu } };

    for (int l=0; l<5; ++l)
      for (int m=0; m<5; ++m)
        Aii[l*neq+m] += c_np1*Pinv[l][m];


    //turbulence preconditioning
    if(neq==6){
      double t1 = Un[i][5] * invRho;
      double mup = mu*t1*gam1/c2;
      double Pt[6] = {mu*nu*t1, -mup*u, -mup*v, -mup*w, mup, 1.0};
      for (int k=0; k<6; k++)
        Aii[neq*(neq-1)+k] += c_np1*Pt[k];

    }else if(neq==7){
      double t1 = Un[i][5] * invRho;
      double t2 = Un[i][6] * invRho;
      double mup1 = mu*t1*gam1/c2;
      double mup2 = mu*t2*gam1/c2;
      double Pt[2][7] = { {mu*nu*t1, -mup1*u, -mup1*v, -mup1*w, mup1, 1.0, 0.0},
                          {mu*nu*t2, -mup2*u, -mup2*v, -mup2*w, mup2, 0.0, 1.0} };
      for (int k=0; k<7; k++){
        Aii[neq*(neq-2)+k] += c_np1*Pt[0][k];
        Aii[neq*(neq-1)+k] += c_np1*Pt[1][k];
      }
    }
  }
}
//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianLiquidPrec(bool *nodeFlag, Vec<double> &ctrlVol, GenMat<Scalar,neq> &A,
                                   SVec<double,dim> &U, VarFcn *vf, TimeLowMachPrec &tprec,
				   Vec<double> &irey, int* nodeType)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;
  if (data.use_modal == true && data.use_freq == false) A *= 0.5;
                                                                                                                           
  if(!nodeType){
    for (int i=0; i<dt.size(); ++i) 
      addToJacobianLiquidPrecLocal(i,ctrlVol[i],vf,tprec,irey[i],U,A);

  }else{
    for (int i=0; i<dt.size(); ++i) 
      if(!(nodeType[i]==BC_INLET_MOVING || nodeType[i]==BC_OUTLET_MOVING ||
           nodeType[i]==BC_INLET_FIXED || nodeType[i]==BC_OUTLET_FIXED ||
           nodeType[i]==BC_DIRECTSTATE_INLET_MOVING || nodeType[i]==BC_DIRECTSTATE_OUTLET_MOVING ||
           nodeType[i]==BC_DIRECTSTATE_INLET_FIXED || nodeType[i]==BC_DIRECTSTATE_OUTLET_FIXED ||
           nodeType[i]==BC_MASSFLOW_INLET_MOVING || nodeType[i]==BC_MASSFLOW_OUTLET_MOVING ||
           nodeType[i]==BC_MASSFLOW_INLET_FIXED || nodeType[i]==BC_MASSFLOW_OUTLET_FIXED) )
        addToJacobianLiquidPrecLocal(i,ctrlVol[i],vf,tprec,irey[i],U,A);
  }
}
//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToJacobianLiquidPrecLocal(int i, double vol, VarFcn *vf,
				TimeLowMachPrec &tprec, double irey,
				SVec<double,dim> &U, GenMat<Scalar,neq> &A)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

// ARL : turbulence preconditioning never tested ...
  double c_np1;
  switch (descriptorCase) {
    case DESCRIPTOR: {
      c_np1 = data.alpha_np1 * vol / dt[i];
      c_np1 += data.dtau_switch * vol / dtau[i];
      break; }
    case HYBRID: {
      c_np1 = data.alpha_np1 * sqrt(vol) / dt[i];
      c_np1 += data.dtau_switch * sqrt(vol) / dtau[i];
      break; }
    case NONDESCRIPTOR: {
      c_np1 = data.alpha_np1 / dt[i];
      c_np1 += data.dtau_switch * 1.0 / dtau[i];
      break; }
  }
  Scalar *Aii = A.getElem_ii(i);
  if(neq<5){            //turbulence model equation in segregated solver
    for (int k=0; k<neq; ++k)
      Aii[k + k*neq] += c_np1;
  }else{        //Navier-Stokes (part of segregated turb model or alone) or fully coupled
    double V[dim];
    vf->conservativeToPrimitive(Un[i],V); // assumption : no steady two-phase flow, hence no phi
    double e = vf->computeRhoEnergy(V)/V[0];
    double pressure = vf->getPressure(V);
    double c = vf->computeSoundSpeed(V);
    double c2 = c*c;
    double locMach = vf->computeMachNumber(V); //local Preconditioning (ARL)
    double beta = tprec.getBeta(locMach,irey);
    double beta2 =   beta * beta;

// Preconditioning in unsteady flow
    double bt =  1.0 /(1.0 + data.dtau_switch*dtau[i]*data.alpha_np1/dt[i]);
    beta2 = beta2/(bt - beta2*(bt - 1.0));

    double oobeta2 = 1.0/beta2;
    double oobeta2m1 = oobeta2 - 1.0;

    double Pinv[dim];
    for (int j=0; j<dim; j++)
      Pinv[j] = oobeta2m1*V[j];
    Pinv[0] = oobeta2;
    Pinv[4] = oobeta2m1*(e+pressure/V[0] - c2);
    /* The preconditioning matrix is:
     * Pinv[dim][dim] = { { oobeta2,           0.0, 0.0, 0.0, 0.0 , 0.0, 0.0 },
     *                    {(oobeta2-1.0)*V[1], 1.0, 0.0, 0.0, 0.0 , 0.0, 0.0 },
     *                    {(oobeta2-1.0)*V[2], 0.0, 1.0, 0.0, 0.0 , 0.0, 0.0 },
     *                    {(oobeta2-1.0)*V[3], 0.0, 0.0, 1.0, 0.0 , 0.0, 0.0 },
     *                    {(oobeta2-1.0)*(h-c2), 0.0, 0.0, 0.0, 1.0 , 0.0, 0.0 },
     *                    {(oobeta2-1.0)*V[5], 0.0, 0.0, 0.0, 0.0 , 1.0, 0.0 },
     *                    {(oobeta2-1.0)*V[6], 0.0, 0.0, 0.0, 0.0 , 0.0, 1.0 } };
     * Take the first 5-by-5 matrix to get the Euler preconditioner
     *      the first 6-by-6 matrix to get the "SA"  preconditioner
     *      the whole 7-by-7 matrix to get the "k-e" preconditioner
     */
 
    for (int j=1; j<dim; j++)
      Aii[j + j*neq] += c_np1;
    for (int j=0; j<dim; j++)
      Aii[j*neq] += c_np1*Pinv[j];
  }
}  
//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH1(bool *nodeFlag, Vec<double> &ctrlVol, GenMat<Scalar,neq> &A)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  if (data.use_modal == true && data.use_freq == false) A *= 0.5;

  double c_np1;
  for (int i = 0; i < dt.size(); ++i) {

    if (nodeFlag && !nodeFlag[i]) continue;

    switch (descriptorCase) {
      case DESCRIPTOR: {
        if (data.use_freq == true)
          c_np1 = data.alpha_np1 * ctrlVol[i];
        else {
          c_np1 = data.alpha_np1 * ctrlVol[i] / dt[i];
          c_np1 += data.dtau_switch * ctrlVol[i] / dtau[i]; }
        break; }
      case HYBRID: {
        c_np1 = data.alpha_np1 * sqrt(ctrlVol[i]) / dt[i];
        c_np1 += data.dtau_switch * sqrt(ctrlVol[i]) / dtau[i];
        break; }
      case NONDESCRIPTOR: {
        c_np1 = data.alpha_np1 / dt[i];
        c_np1 += data.dtau_switch * 1.0 / dtau[i];
        break; }
    }
    Scalar *Aii = A.getElem_ii(i);

    for (int k=0; k<neq; ++k)
      Aii[k + k*neq] += c_np1;

  }

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH1(bool *nodeFlag, Vec<double> &ctrlVol,
                GenMat<Scalar,neq> &A, Scalar shift)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

//  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  if (data.use_modal == true && data.use_freq == false) A *= 0.5;

  Scalar c_np1;
  for (int i=0; i<dt.size(); ++i) {

    switch (descriptorCase)  {
      case DESCRIPTOR: {
        if (data.use_freq == true)
          c_np1 = shift * ctrlVol[i];
        else {
          c_np1 = data.alpha_np1 * ctrlVol[i] / dt[i];
          c_np1 += data.dtau_switch * ctrlVol[i] / dtau[i]; }
        break; }
      case HYBRID: {
        c_np1 = data.alpha_np1 * sqrt(ctrlVol[i]) / dt[i];
        c_np1 += data.dtau_switch * sqrt(ctrlVol[i]) / dtau[i];
        break; }
      case NONDESCRIPTOR: {
        if (data.use_freq == true)
          c_np1 = shift;
        else {
          c_np1 = data.alpha_np1 / dt[i];
          c_np1 += data.dtau_switch * 1.0 / dtau[i]; }
        break; }
    }
    Scalar *Aii = A.getElem_ii(i);

    for (int k=0; k<neq; ++k) {
      Aii[k + k*neq] += c_np1;
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH2(bool *nodeFlag, VarFcn *varFcn, Vec<double> &ctrlVol,
			     SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  double dfdUi[dim*dim], dfdVi[dim*dim];

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  double coef = data.alpha_np1;
  if (data.use_modal == true && data.use_freq == false) {
    A *= 2.0;
    coef *= 3.0;
  }
  
  double c_np1;
  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && !nodeFlag[i]) continue;

    if (data.use_freq == true)
      c_np1 = data.alpha_np1 * ctrlVol[i];
    else {
      c_np1 = coef * ctrlVol[i] / dt[i];
      c_np1 += data.dtau_switch * ctrlVol[i] / dtau[i]; }

    int k;
    for (k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;
    for (k=0; k<dim; ++k) dfdUi[k + k*dim] = c_np1;

    varFcn->postMultiplyBydUdV(V[i], dfdUi, dfdVi);
  
    Scalar *Aii = A.getElem_ii(i);

    for (k=0; k<neq*neq; ++k) Aii[k] += dfdVi[k];

  }

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH2(bool *nodeFlag, VarFcn *varFcn,
                Vec<double> &ctrlVol, SVec<double,dim> &V,
                GenMat<Scalar,neq> &A, Scalar shift)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  Scalar dfdUi[dim*dim], dfdVi[dim*dim];

  //if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  if (data.use_modal && data.use_freq == false) A *= 0.5;

  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && nodeFlag[i] == 0) continue;
      Scalar c_np1;
      if (data.use_freq == true)
        c_np1 = shift*ctrlVol[i];
      else {
        c_np1 = data.alpha_np1 * ctrlVol[i] / dt[i];
        c_np1 += data.dtau_switch * ctrlVol[i] / dtau[i]; }

    int k;
    for (k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;
    for (k=0; k<dim; ++k) dfdUi[k + k*dim] = c_np1;

    varFcn->postMultiplyBydUdV(V[i], dfdUi, dfdVi);

    Scalar *Aii = A.getElem_ii(i);

    for (k=0; k<neq*neq; ++k) Aii[k] += dfdVi[k];

  }

}

//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH2(bool *nodeFlag, VarFcn *varFcn, Vec<double> &ctrlVol,
                             SVec<double,dim> &V, GenMat<Scalar,neq> &A, Scalar coefVol, double coefA)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  Scalar dfdUi[dim*dim], dfdVi[dim*dim];

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  //double coef = 1;
//data.alpha_np1;
  if (data.use_modal == true && data.use_freq == false) {
    A *= coefA;
  }

  Scalar c_np1;
  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && !nodeFlag[i]) continue;

    if (data.use_freq == true)
      c_np1 = coefVol * ctrlVol[i];
    else {
      c_np1 = coefVol * ctrlVol[i] / dt[i];
      c_np1 += data.dtau_switch * ctrlVol[i] / dtau[i]; }

    int k;
    for (k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;
    for (k=0; k<dim; ++k) dfdUi[k + k*dim] = c_np1;

    varFcn->postMultiplyBydUdV(V[i], dfdUi, dfdVi);

    Scalar *Aii = A.getElem_ii(i);

    for (k=0; k<neq*neq; ++k) Aii[k] += dfdVi[k];

  }

}



//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH2Minus(bool *nodeFlag, VarFcn *varFcn, Vec<double> &ctrlVol,
                                  SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  double dfdUi[dim*dim], dfdVi[dim*dim];

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  if (data.use_modal == true && data.use_freq == false) A *= -0.5;

  double c_np1;
  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && nodeFlag[i] == 0) continue;

    if (data.use_freq)
      c_np1 = -2.0*data.alpha_np1 * ctrlVol[i];
    else {
      c_np1 = data.alpha_np1 * ctrlVol[i] / dt[i];
      c_np1 += data.dtau_switch * ctrlVol[i] / dtau[i]; }

    int k;
    for (k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;
    for (k=0; k<dim; ++k) dfdUi[k + k*dim] = c_np1;

    varFcn->postMultiplyBydUdV(V[i], dfdUi, dfdVi);

    Scalar *Aii = A.getElem_ii(i);

    for (k=0; k<neq*neq; ++k) Aii[k] += dfdVi[k];

  }

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH2NoPrec(bool *nodeFlag, VarFcn *varFcn, Vec<double> &ctrlVol,
			     SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  double dfdUi[dim*dim], dfdVi[dim*dim];

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  double coef = data.alpha_np1;
  if (data.use_modal == true && data.use_freq == false) {
    A *= 2.0;
    coef *= 3.0;
  }
  
  double c_np1;
  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && !nodeFlag[i]) continue;

    if (data.use_freq == true)
      c_np1 = data.alpha_np1 * ctrlVol[i];
    else {
      c_np1 = coef * ctrlVol[i] / dt[i];
      c_np1 += data.dtau_switch * ctrlVol[i] / dtau[i]; }

    int k;
    for (k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;
    for (k=0; k<dim; ++k) dfdUi[k + k*dim] = c_np1;

    varFcn->postMultiplyBydUdV(V[i], dfdUi, dfdVi);
  
    Scalar *Aii = A.getElem_ii(i);

    for (k=0; k<neq*neq; ++k) Aii[k] += dfdVi[k];
 
  }

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH2GasPrec(bool *nodeFlag, VarFcn *varFcn, Vec<double> &ctrlVol,
			     SVec<double,dim> &V, GenMat<Scalar,neq> &A, 
                             double gam, double pstiff, Vec<double> &irey, 
                             TimeLowMachPrec &tprec)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  if (data.use_modal == true && data.use_freq == false) {
    A *= 2.0;
  }
  
  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && !nodeFlag[i]) continue;

    addToH2GasPrecLocal(i,ctrlVol[i],varFcn,gam,pstiff,tprec,irey[i],V,A);

  }
}
//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH2GasPrecLocal(int i, double vol, VarFcn *vf, double gam, 
				double pstiff, TimeLowMachPrec &tprec,
				double irey, SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  double dfdUi[dim*dim], dfdVi[dim*dim];

  for (int k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;

  double coef = data.alpha_np1;
  if (data.use_modal == true && data.use_freq == false) {
    coef *= 3.0;
  }
 
  double c_np1;
  if (data.use_freq == true)
    c_np1 = data.alpha_np1 * vol;
  else {
    c_np1 = coef * vol / dt[i];
    c_np1 += data.dtau_switch * vol / dtau[i]; }

  if(neq<5){		//turbulence model equation in segregated solver
    for (int k=0; k<neq; ++k)
      dfdUi[k + k*neq] += c_np1;
  }else{	//Navier-Stokes (part of segregated turb model or alone) or fully coupled

    double ro = Un[i][0];
    double invRho = 1.0/ro;
    double u  = Un[i][1] * invRho;
    double v  = Un[i][2] * invRho;
    double w  = Un[i][3] * invRho;
    double u2 = u*u;
    double v2 = v*v;
    double w2 = w*w;
    double q2 = u2 + v2 + w2;
    double gam1 = gam - 1.0;
    double p  = gam1 * (Un[i][4] - 0.5 * ro * q2) - gam*pstiff;
    double c2 = gam*(p+pstiff)/ro;
    double locMach = sqrt(q2/c2); //local Preconditioning (ARL)
    double beta = tprec.getBeta(locMach, irey);

    double beta2 =   beta * beta;

// Preconditioning in unsteady flow
    double bt =  1.0 /(1.0 + data.dtau_switch*dtau[i]*data.alpha_np1/dt[i]);
    beta2 = beta2/(bt - beta2*(bt - 1.0));

    double qhat2 = (q2 * gam1)/2.0;
 
    double nu = qhat2/c2;
    double mu = (1.0/beta2) - 1.0;

    double Pinv[5][5] = { {nu*mu + 1.0,  -u*mu*gam1/c2,      -v*mu*gam1/c2,        -w*mu*gam1/c2,       mu*gam1/c2   },
                          {u*nu*mu,     1.0 - u2*mu*gam1/c2, -u*v*mu*gam1/c2,      -u*w*mu*gam1/c2,     u*mu*gam1/c2 },
                          {v*nu*mu,     -u*v*mu*gam1/c2 ,    1.0 - v2*mu*gam1/c2,  -v*w*mu*gam1/c2,     v*mu*gam1/c2 },
                          {w*nu*mu,     -u*w*mu*gam1/c2 ,    -v*w*mu*gam1/c2,      1.0 - w2*mu*gam1/c2, w*mu*gam1/c2 },
     	                  {0.5*mu*(1.0+nu)*q2,    -u*mu*(1+nu), -v*mu*(1+nu), -w*mu*(1+nu), (1.0/beta2)+mu*nu } };

    for (int l=0; l<5; ++l)
      for (int m=0; m<5; ++m)
        dfdUi[l*neq+m] += c_np1*Pinv[l][m];


    //turbulence preconditioning
    if(neq==6){
      double t1 = Un[i][5] * invRho;
      double mup = mu*t1*gam1/c2;
      double Pt[6] = {mu*nu*t1, -mup*u, -mup*v, -mup*w, mup, 1.0};
      for (int k=0; k<6; k++)
        dfdUi[neq*(neq-1)+k] += c_np1*Pt[k];

    }else if(neq==7){
      double t1 = Un[i][5] * invRho;
      double t2 = Un[i][6] * invRho;
      double mup1 = mu*t1*gam1/c2;
      double mup2 = mu*t2*gam1/c2;
      double Pt[2][7] = { {mu*nu*t1, -mup1*u, -mup1*v, -mup1*w, mup1, 1.0, 0.0},
                          {mu*nu*t2, -mup2*u, -mup2*v, -mup2*w, mup2, 0.0, 1.0} };
      for (int k=0; k<7; k++){
        dfdUi[neq*(neq-2)+k] += c_np1*Pt[0][k];
        dfdUi[neq*(neq-1)+k] += c_np1*Pt[1][k];
      }
    }
  }

  vf->postMultiplyBydUdV(V[i], dfdUi, dfdVi);
  
  Scalar *Aii = A.getElem_ii(i);

  for (int k=0; k<neq*neq; ++k) Aii[k] += dfdVi[k];
 
}
//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH2LiquidPrec(bool *nodeFlag, VarFcn *varFcn, Vec<double> &ctrlVol,
			     SVec<double,dim> &V, GenMat<Scalar,neq> &A, 
                             Vec<double> &irey, TimeLowMachPrec &tprec)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  if (data.typeIntegrator == ImplicitData::CRANK_NICOLSON) A *= 0.5;

  if (data.use_modal == true && data.use_freq == false) {
    A *= 2.0;
  }
  
  for (int i=0; i<dt.size(); ++i) {

    if (nodeFlag && !nodeFlag[i]) continue;

    addToH2LiquidPrecLocal(i, ctrlVol[i], varFcn, tprec, irey[i], V, A);

  }
}
//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void TimeState<dim>::addToH2LiquidPrecLocal(int i, double vol, VarFcn *vf,
				TimeLowMachPrec &tprec, double irey,
				SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{
  if (data.typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  double dfdUi[dim*dim], dfdVi[dim*dim];

  for (int k=0; k<dim*dim; ++k) dfdUi[k] = 0.0;

  double coef = data.alpha_np1;
  if (data.use_modal == true && data.use_freq == false) {
    coef *= 3.0;
  }
 
  double c_np1;
  if (data.use_freq == true)
    c_np1 = data.alpha_np1 * vol;
  else {
    c_np1 = coef * vol / dt[i];
    c_np1 += data.dtau_switch * vol / dtau[i]; }

  if(neq<5){            //turbulence model equation in segregated solver
    for (int k=0; k<neq; ++k)
      dfdUi[k + k*neq] += c_np1;
  }else{        //Navier-Stokes (part of segregated turb model or alone) or fully coupled
    double V[dim];
    vf->conservativeToPrimitive(Un[i],V); // assumption : no steady two-phase flow, hence no phi
    double e = vf->computeRhoEnergy(V)/V[0];
    double pressure = vf->getPressure(V);
    double c = vf->computeSoundSpeed(V);
    double c2 = c*c;
    double locMach = vf->computeMachNumber(V); //local Preconditioning (ARL)
    double beta = tprec.getBeta(locMach,irey);
    double beta2 =   beta * beta;

// Preconditioning in unsteady flow
    double bt =  1.0 /(1.0 + data.dtau_switch*dtau[i]*data.alpha_np1/dt[i]);
    beta2 = beta2/(bt - beta2*(bt - 1.0));

    double oobeta2 = 1.0/beta2;
    double oobeta2m1 = oobeta2 - 1.0;

    double Pinv[dim];
    for (int j=0; j<dim; j++)
      Pinv[j] = oobeta2m1*V[j];
    Pinv[0] = oobeta2;
    Pinv[4] = oobeta2m1*(e+pressure/V[0] - c2);
 
    for (int j=1; j<dim; j++)
      dfdUi[j + j*neq] += c_np1;
    for (int j=0; j<dim; j++)
      dfdUi[j*neq] += c_np1*Pinv[j];
  }

  vf->postMultiplyBydUdV(V[i], dfdUi, dfdVi);
  
  Scalar *Aii = A.getElem_ii(i);

  for (int k=0; k<neq*neq; ++k) Aii[k] += dfdVi[k];
 
}  
//------------------------------------------------------------------------------
                                                                                                                          
template<int dim>
void TimeState<dim>::get_dW_dt(bool *nodeFlag, GeoState &geoState,
                               Vec<double> &ctrlVol, SVec<double,dim> &Q,
                               SVec<double,dim> &R)
{
                                                                                                                          
  Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();
  Vec<double>& ctrlVol_nm1 = geoState.getCtrlVol_nm1();
  Vec<double>& ctrlVol_nm2 = geoState.getCtrlVol_nm2();
                                                                                                                          
  double c_np1, c_n, c_nm1, c_nm2;
                                                                                                                          
  for (int i=0; i<dt.size(); ++i) {
    if (!nodeFlag[i]) {
      double invDt = 1.0 / dt[i];
      c_np1 = data.alpha_np1*ctrlVol[i];
      c_n   = data.alpha_n * ctrlVol_n[i];
      c_nm1 = data.alpha_nm1 * ctrlVol_nm1[i];
      c_nm2 = data.alpha_nm2 * ctrlVol_nm2[i];
                                                                                                                          
      for (int k=0; k<dim; ++k) {
        double dWdt = invDt * (c_np1*Q[i][k] + c_n*Un[i][k] +
                               c_nm1*Unm1[i][k] + c_nm2*Unm2[i][k]);
        R[i][k] += dWdt;
      }
    }
  }
                                                                                                                          
}
                                                                                                                          
//------------------------------------------------------------------------------
                                                                                                                          
template<int dim>
void TimeState<dim>::get_dWBar_dt(bool *nodeFlag, GeoState &geoState,
                                  Vec<double> &ctrlVol, SVec<double,dim> &QBar,
                                  SVec<double,dim> &UnBar, SVec<double,dim> &Unm1Bar,
                                  SVec<double,dim> &Unm2Bar, SVec<double,dim> &R)
{
                                                                                                                          
  Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();
  Vec<double>& ctrlVol_nm1 = geoState.getCtrlVol_nm1();
  Vec<double>& ctrlVol_nm2 = geoState.getCtrlVol_nm2();
                                                                                                                          
  double c_np1, c_n, c_nm1, c_nm2;
                                                                                                                          
  for (int i=0; i<dt.size(); ++i) {
    if (!nodeFlag[i]) {
      double invDt = 1.0 / dt[i];
      c_np1 = data.alpha_np1*ctrlVol[i];
      c_n   = data.alpha_n * ctrlVol_n[i];
      c_nm1 = data.alpha_nm1 * ctrlVol_nm1[i];
      c_nm2 = data.alpha_nm2 * ctrlVol_nm2[i];
                                                                                                                          
      for (int k=0; k<dim; ++k) {
        double dWdt = invDt * (c_np1*QBar[i][k] + c_n*UnBar[i][k] +
                               c_nm1*Unm1Bar[i][k] + c_nm2*Unm2Bar[i][k]);
        R[i][k] += dWdt;
      }
    }
  }
                                                                                                                          
}
                                                                                                                          
//------------------------------------------------------------------------------
