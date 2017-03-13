#include <FluxFcnDesc.h>

//------------------------------------------------------------------------------
// VL and VR are the primitive state variables !!

void FluxFcnVanLeerEuler3D::evalFlux(double vfgam, double vfp, double *normal, double normalVel, 
				     double *V, double *f, int sign) 
{

  //gamma 
  double gam = vfgam;
  double invgam = 1.0/gam;
  double gam1 = gam - 1.0;
  double invgam1 = 1.0/gam1;

  // compute norm to remove area factor from computations
  double norm = sqrt( normal[0]*normal[0]
                    + normal[1]*normal[1]
                    + normal[2]*normal[2] );

  double invNorm = 1.0 / norm;
  	            	     
  // normalize the normal 
  double nVec[3];
  nVec[0] = normal[0]*invNorm;
  nVec[1] = normal[1]*invNorm;
  nVec[2] = normal[2]*invNorm;
  
  //normalize grid velocity
  double gridVel = normalVel * invNorm;

  //compute speed of sound
  double a = vf->computeSoundSpeed(V);
  
  // compute fluid velocity
  double fluidVel = V[1]*nVec[0] + V[2]*nVec[1] + V[3]*nVec[2];

  // compute total velocity --> fluid vel. - grid vel.
  double totVel = fluidVel - gridVel;

  //compute normal mach numbers
  double mach = totVel / a;

  // get conservative variables
  double U[5];
  vf->primitiveToConservative(V, U);
    
  double factor = mach*sign;
  if (factor >= 1.0)  {
    f[0] = norm * U[0] * totVel;

    f[1] = norm * ( U[1] * totVel
                  + V[4] * nVec[0] ); 

    f[2] = norm * ( U[2] * totVel
                  + V[4] * nVec[1] ); 

    f[3] = norm * ( U[3] * totVel
         	  + V[4] * nVec[2] ); 

    f[4] = norm * ( U[4] * totVel
 		  + V[4] * fluidVel );
  }  
  else if (factor > -1.0)  {
    // compute frequently used terms
    double h = 0.25 * norm * sign * V[0] * a * (mach+sign)*(mach+sign);  
    double f1 = (2*a*sign - totVel) * invgam;

    f[0] = h;

    f[1] = h * (V[1] + f1 * nVec[0]);

    f[2] = h * (V[2] + f1 * nVec[1]);

    f[3] = h * (V[3] + f1 * nVec[2]);


    f[4] = h * ( (2*a*a + 2*sign * gam1 * totVel * a 
                 - totVel*totVel*gam1)/(gam*gam - 1.0)
               + .5*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])
               + gridVel*f1 ) ;
  }  
  else {
    for (int i = 0; i < 5; i++)
      f[i] = 0;  
  }

}

//-----------------------------------------------------------------------
// VL and VR are the primitive state variables !!

// Included (MB)
void FluxFcnVanLeerEuler3D::evalDerivativeOfFlux(double vfgam, double vfp, double dvfp, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *V, double *dV, double dM, double *dF, int sign)
{

  //gamma 
  double gam = vfgam;
  double invgam = 1.0/gam;
  double gam1 = gam - 1.0;
  double invgam1 = 1.0/gam1;

  // compute norm to remove area factor from computations
  double norm = sqrt( normal[0]*normal[0]
                    + normal[1]*normal[1]
                    + normal[2]*normal[2] );

  double dNorm = 1.0 / (2.0*norm) * (2.0*normal[0]*dNormal[0]
                                   + 2.0*normal[1]*dNormal[1]
                                   + 2.0*normal[2]*dNormal[2]);


  double invNorm = 1.0 / norm;

  double dInvNorm = -1.0 / (2.0*norm*norm*norm) * (2.0*normal[0]*dNormal[0]
                                                 + 2.0*normal[1]*dNormal[1]
                                                 + 2.0*normal[2]*dNormal[2]);

  // normalize the normal
  double nVec[3];
  nVec[0] = normal[0]*invNorm;
  nVec[1] = normal[1]*invNorm;
  nVec[2] = normal[2]*invNorm;

  double dNVec[3];
  dNVec[0] = dNormal[0]*invNorm + normal[0]*dInvNorm;
  dNVec[1] = dNormal[1]*invNorm + normal[1]*dInvNorm;
  dNVec[2] = dNormal[2]*invNorm + normal[2]*dInvNorm;

  //normalize grid velocity
  double gridVel = normalVel * invNorm;
  double dGridVel = dNormalVel * invNorm + normalVel * dInvNorm;

  //compute speed of sound
  double a = vf->computeSoundSpeed(V);
  double da = vf->computeDerivativeOfSoundSpeed(V, dV, dM);

  // compute fluid velocity
  double fluidVel = V[1]*nVec[0] + V[2]*nVec[1] + V[3]*nVec[2];
  double dFluidVel = dV[1]*nVec[0] + V[1]*dNVec[0] + dV[2]*nVec[1] + V[2]*dNVec[1] + dV[3]*nVec[2] + V[3]*dNVec[2];

  // compute total velocity --> fluid vel. - grid vel.
  double totVel = fluidVel - gridVel;
  double dTotVel = dFluidVel - dGridVel;

  //compute normal mach numbers
  double mach = totVel / a;
  double dMach = ( dTotVel * a - totVel * da ) / ( a * a );

  // get conservative variables
  double U[5], dU[5];
  vf->primitiveToConservative(V, U);
  vf->primitiveToConservativeDerivative(V, dV, U, dU);

  double factor = mach*sign;
  if (factor >= 1.0)  {
    dF[0] = dNorm * U[0] * totVel + norm * dU[0] * totVel + norm * U[0] * dTotVel;

    dF[1] = dNorm * ( U[1] * totVel + V[4] * nVec[0] ) +
           norm * ( dU[1] * totVel + U[1] * dTotVel + dV[4] * nVec[0] + V[4] * dNVec[0] );

    dF[2] = dNorm * ( U[2] * totVel + V[4] * nVec[1] ) +
           norm * ( dU[2] * totVel + U[2] * dTotVel + dV[4] * nVec[1] + V[4] * dNVec[1] );

    dF[3] = dNorm * ( U[3] * totVel + V[4] * nVec[2] ) +
           norm * ( dU[3] * totVel + U[3] * dTotVel + dV[4] * nVec[2] + V[4] * dNVec[2] );

    dF[4] = dNorm * ( U[4] * totVel + V[4] * fluidVel ) +
           norm * ( dU[4] * totVel + U[4] * dTotVel + dV[4] * fluidVel + V[4] * dFluidVel );
  }
  else if (factor > -1.0)  {
    // compute frequently used terms
    double h = 0.25 * norm * sign * V[0] * a * (mach+sign)*(mach+sign);

    double dh = 0.25 * dNorm * sign * V[0] * a * (mach+sign)*(mach+sign) +
               0.25 * norm * sign * dV[0] * a * (mach+sign)*(mach+sign) +
               0.25 * norm * sign * V[0] * a * 2.0 * (mach+sign) * dMach;

    double f1 = invgam * (2*a*sign - totVel);

    double df1 = invgam * (2*da*sign - dTotVel);

    dF[0] = dh;

    dF[1] = dh * (V[1] + f1 * nVec[0]) +
           h * (dV[1] + df1 * nVec[0] + f1 * dNVec[0]);

    dF[2] = dh * (V[2] + f1 * nVec[1]) +
           h * (dV[2] + df1 * nVec[1] + f1 * dNVec[1]);

    dF[3] = dh * (V[3] + f1 * nVec[2]) +
           h * (dV[3] + df1 * nVec[2] + f1 * dNVec[2]);

    dF[4] = dh * ( (2.0*a*a + 2.0*sign * gam1 * totVel * a
                    - totVel*totVel*gam1)/(gam*gam - 1.0)
                    + .5*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3])
                    + gridVel*f1 ) +
            h * ( (4.0*a*da + 2.0*sign * gam1 * dTotVel * a + 2.0*sign * gam1 * totVel * da
                    - 2.0*totVel*dTotVel*gam1)/(gam*gam - 1.0)
                    + .5*(2.0*V[1]*dV[1]+2.0*V[2]*dV[2]+2.0*V[3]*dV[3])
                    + dGridVel*f1 + gridVel*df1 );
  }
  else {
    for (int i = 0; i < 5; i++)
      dF[i] = 0;
  }

}

//------------------------------------------------------------------------------

void FluxFcnVanLeerEuler3D::evalJac(double vfgam, double vfp, double *normal, double normalVel,
				    double *V, double *jac, int sign)

{

  //gamma 
  double gam = vfgam;
  double invgam = 1.0/gam;
  double gam1 = gam - 1.0;
  double invgam1 = 1.0/gam1;

  // compute norm to remove area factor from computations
  double norm = sqrt( normal[0]*normal[0]
                    + normal[1]*normal[1]
                    + normal[2]*normal[2] );

  double invNorm = 1.0 / norm;
 
  // normalize the normal
  double nVec[3];
  nVec[0] = normal[0]*invNorm;
  nVec[1] = normal[1]*invNorm;
  nVec[2] = normal[2]*invNorm;
 
  //normalize grid velocity
  double gridVel = normalVel * invNorm;

  //compute speed of sound
  double a = vf->computeSoundSpeed(V);

  // compute fluid velocity
  double fluidVel = V[1]*nVec[0] + V[2]*nVec[1] + V[3]*nVec[2];

  // compute total velocity --> fluid vel. - grid vel.
  double totVel = fluidVel - gridVel;

  //compute normal mach numbers
  double mach = totVel / a; 

  // get conservative variables
  double U[5];
  vf->primitiveToConservative(V, U);

  // compute frequently used terms
  double q2 = V[1]*V[1]+V[2]*V[2]+V[3]*V[3];
  // compute Energy term
  double E = U[4] / U[0];

  // jacs are computed in the following order: df0/du0, df1/du0, 
  // df2/du0 ... df0/du1, df2/du1, ... df4/du4
  double factor = mach*sign;
  if (factor >= 1.0)  {
    double pOverRho = V[4] / V[0];
    jac[0] = -gridVel;
    jac[1] = norm * nVec[0];
    jac[2] = norm * nVec[1];
    jac[3] = norm * nVec[2];
    jac[4] = 0;

    jac[5] = norm * (0.5*gam1*q2*nVec[0] - V[1]*fluidVel);
    jac[6] = norm * ( totVel + V[1]*nVec[0]*(2.0-gam) );
    jac[7] = norm * ( V[1] * nVec[1] - gam1 * V[2] * nVec[0] );
    jac[8] = norm * ( V[1] * nVec[2] - gam1 * V[3] * nVec[0] );
    jac[9] = norm * gam1 * nVec[0];

    jac[10] = norm * (0.5*gam1*q2*nVec[1] - V[2]*fluidVel);
    jac[11] = norm * ( V[2] * nVec[0] - gam1 * V[1] * nVec[1] );
    jac[12] = norm * ( totVel + V[2]*nVec[1]*(2.0-gam) );
    jac[13] = norm * ( V[2] * nVec[2] - gam1 * V[3] * nVec[1] );
    jac[14] = norm * gam1 * nVec[1];

    jac[15] = norm * (0.5*gam1*q2*nVec[2] - V[3]*fluidVel);
    jac[16] = norm * ( V[3] * nVec[0] - gam1 * V[1] * nVec[2] );
    jac[17] = norm * ( V[3] * nVec[1] - gam1 * V[2] * nVec[2] );
    jac[18] = norm * ( totVel + V[3]*nVec[2]*(2.0-gam) );
    jac[19] = norm * gam1 * nVec[2];

    jac[20] = norm*fluidVel * ( 0.5*+gam1*q2 - E - pOverRho );
    jac[21] = norm * ( nVec[0]*(E + pOverRho) - gam1*V[1]*fluidVel );
    jac[22] = norm * ( nVec[1]*(E + pOverRho) - gam1*V[2]*fluidVel );
    jac[23] = norm * ( nVec[2]*(E + pOverRho) - gam1*V[3]*fluidVel );
    jac[24] = norm * (fluidVel*gam - gridVel);    
  }
  else if (factor > -1.0)  {
    // compute frequently used terms
    double k3 = 1.0 - sign * mach;
    double k4 = 0.5*(1.0 + sign * mach);
    double coef = sign*V[0]*a*k4*k4;
    double f2 = invgam*a*sign;
    double invgamPlus1 = 1.0 / (gam + 1.0);
    double invRho = 1.0 / V[0];
    double invA = 1.0/a;
    double ggam1 = gam*gam1;
 
    double g[4];
    g[0] = V[1] + nVec[0]*f2*(1.0 + k3);
    g[1] = V[2] + nVec[1]*f2*(1.0 + k3);
    g[2] = V[3] + nVec[2]*f2*(1.0 + k3);
    g[3] = f2*gridVel*(1.0 + k3) + 0.5*q2 + a*a*invgam1 
	 - a*a*k3*k3*invgamPlus1; 

    double gGrad[4][5];
    //dg0/du
    gGrad[0][0] = invRho * 
		( nVec[0]*(gridVel*invgam + gam1*invA*sign*E) - g[0] );
    gGrad[0][1] = invRho * ( 1.0 - sign*invgam*nVec[0]
		           * (ggam1*V[1]*invA + sign*nVec[0]) );
    gGrad[0][2] = -sign*invgam*invRho*nVec[0]
        	* ( ggam1*invA*V[2] + sign*nVec[1] );
    gGrad[0][3] = -sign*invgam*invRho*nVec[0]
        	* ( ggam1*invA*V[3] + sign*nVec[2] );
    gGrad[0][4] = sign * invA * invRho * gam1 * nVec[0];

    // dg1/du
    gGrad[1][0] = invRho * ( nVec[1]*(gridVel*invgam + gam1*invA*sign*E)
                           - g[1] );
    gGrad[1][1] = -sign*invgam*invRho*nVec[1]
                * ( ggam1*invA*V[1] + sign*nVec[0] );
    gGrad[1][2] = invRho * ( 1.0 - sign*invgam*nVec[1]
                           * (ggam1*V[2]*invA + sign*nVec[1]) );
    gGrad[1][3] = -sign*invgam*invRho*nVec[1]
                * ( ggam1*invA*V[3] + sign*nVec[2] );
    gGrad[1][4] = sign * invA * invRho * gam1 * nVec[1];

    // dg2/du
    gGrad[2][0] = invRho * ( nVec[2]*(gridVel*invgam + gam1*invA*sign*E)
                           - g[2] );
    gGrad[2][1] = -sign*invgam*invRho*nVec[2]
                * ( ggam1*invA*V[1] + sign*nVec[0] );
    gGrad[2][2] = -sign*invgam*invRho*nVec[2]
                * ( ggam1*invA*V[2] + sign*nVec[1] );
    gGrad[2][3] = invRho * ( 1.0 - sign*invgam*nVec[2]
                           * (ggam1*V[3]*invA + sign*nVec[2]) );
    gGrad[2][4] = sign * invA * invRho * gam1 * nVec[2];

    //dg3/du
    gGrad[3][0] = invRho * ( (gam1*q2 - gam*E)
                           + k3*invgamPlus1 * ( 2*sign*mach*a*a 
                                              + ggam1*(q2-E) 
					      + 2*a*sign*gridVel )
    			   + gridVel * ( gridVel*invgam + sign*E*invA*gam1
                                       - f2*(1.0 + k3) ) );  

    gGrad[3][1] = invRho*( k3*invgamPlus1 * (ggam1*V[1]+2*a*sign*nVec[0])
                         - gridVel * sign * invgam
			 * (ggam1*invA*V[1] + sign*nVec[0]) - gam1*V[1] );
    gGrad[3][2] = invRho*( k3*invgamPlus1 * (ggam1*V[2]+2*a*sign*nVec[1])
                         - gridVel * sign * invgam
		  	 * (ggam1*invA*V[2] + sign*nVec[1]) - gam1*V[2] );
    gGrad[3][3] = invRho*( k3*invgamPlus1 * (ggam1*V[3]+2*a*sign*nVec[2])
                         - gridVel*sign*invgam*(ggam1*invA*V[3]+sign*nVec[2])
                         - gam1*V[3] );
    gGrad[3][4] = invRho*(gam - k3*ggam1*invgamPlus1 + gridVel*invA*sign*gam1); 

    //df/du0
    jac[0] = k4*(0.25*ggam1*sign*E*invA*k3 - gridVel);
    jac[1] = k4*(nVec[0] - 0.25*ggam1*sign*invA*k3*V[1]);
    jac[2] = k4*(nVec[1] - 0.25*ggam1*sign*invA*k3*V[2]);
    jac[3] = k4*(nVec[2] - 0.25*ggam1*sign*invA*k3*V[3]);
    jac[4] = 0.25*k4*k3*ggam1*sign*invA; 
    
    //df1/du through df4/du
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 5; j++)
        jac[5*(i+1) + j] = norm * (jac[j]*g[i] + coef*gGrad[i][j]);

    jac[0] *= norm;
    jac[1] *= norm;
    jac[2] *= norm;
    jac[3] *= norm;
    jac[4] *= norm;
     
  }
  else {
    for (int i = 0; i < 25; i++) 
      jac[i] = 0.0;
  }

}

//------------------------------------------------------------------------------
