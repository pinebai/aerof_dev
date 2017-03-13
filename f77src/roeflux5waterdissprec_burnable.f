       SUBROUTINE ROEFLUX5WATERBURN(type,gamma,Cv,Pr,alpha,beta,
     &     enormal,evitno,Ugr,Ug,Udr,Ud,phi,
     &     betaRef, k1, cmach, irey, prec)
c-----------------------------------------------------------------------
c This routine computes the Flux of Roe for a fluid like water, taken at the vectors Ug, Ud
c U refers to the primitive variables.
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
c
c In this routine, the dissipation term is computed slightly differently:
c   instead of considering |A|*(Wc_j-Wc_i), where Wc refers to conservative
c   variables, the dissipation term is computed by |A|Q*(Wp_j-Wp_i)
c   where Wp = (P, u, v, w, T)
c   |A| and Q are computed with the roe variables, and are computed at the 
c   same time (we don t compute |A| and Q, but rather directly |A|Q since
c   the analytical expression is simple.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), evitno, phi(*)
      REAL*8 Ugr(*), Udr(*)
      REAL*8 Hg, Hd, vitno, updir, gamma
      REAL*8 VdotN , rnorm, invnorm
      REAL*8 flur1 , flur2 , flur3 , flur4 , flur5
      REAL*8 diss1, diss2, diss3, diss4, diss5
      REAL*8 dif1 , dif2 , dif3 , dif4 ,dif5
      REAL*8 cr , cr2 , qir, tet1, tet2, tet3
      REAL*8 vp1 , vp4 , vp5
      REAL*8 uar1 , uar2 , uar3 , uar4 , uar5, uar6
      REAL*8 usro , squsrg , squsrd
      REAL*8 vitg2, vitd2
      REAL*8 Pr, alpha, beta, Cv, Pg, Pd
      REAL*8 beta1, Ugbeta, Udbeta
      REAL*8 coeff1, coeff2, coeff3, coeff4, coeff5, eps
      INTEGER type
c for precontioning
      REAL*8 k1, betaRef, cmach, irey
      REAL*8 vpa, vpb
      REAL*8 oob2, b2, betaPrec, locMach
      REAL*8 Proeor, P1, P2
      REAL*8 r,s,rs,rps,smr,r2,s2,rp,sp
      REAL*8 nx, ny, nz
      INTEGER prec

c
c Initialisation
c

      beta1 = beta - 1.0d0
     
      rnorm = DSQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) +
     &             enormal(3)*enormal(3))
      invnorm = 1.0d0 / rnorm

      normal(1) = enormal(1) * invnorm
      normal(2) = enormal(2) * invnorm
      normal(3) = enormal(3) * invnorm
      nx = normal(1)
      ny = normal(2)
      nz = normal(3)

      vitno = evitno * invnorm
      
c
c Computation of the centred terms
c

      VdotN = Ug(2)*normal(1) + Ug(3)*normal(2) + Ug(4)*normal(3)
      vitg2 = Ug(2)*Ug(2) + Ug(3)*Ug(3) + Ug(4)*Ug(4)
      Ugbeta = Ug(1)**beta
      Pg = Pr + alpha * Ugbeta
      Hg = Cv*Ug(5) + Pg/Ug(1) + 0.5d0*vitg2
      phi(1) = Ug(1)*(VdotN - vitno)
      phi(2) = phi(1)*Ug(2) + Pg*normal(1)
      phi(3) = phi(1)*Ug(3) + Pg*normal(2)
      phi(4) = phi(1)*Ug(4) + Pg*normal(3)
      phi(5) = phi(1)*Hg + Pg*vitno

      VdotN = Ud(2)*normal(1) + Ud(3)*normal(2) + Ud(4)*normal(3)
      VdotN = VdotN - vitno
      vitd2 = Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4)
      Udbeta = Ud(1)**beta
      Pd = Pr + alpha * Udbeta
      Hd = Cv*Ud(5) + Pd/Ud(1) +0.5d0*vitd2
      phi(1) = phi(1) + Ud(1)*VdotN
      phi(2) = phi(2) + Ud(1)*Ud(2)*VdotN + Pd*normal(1)
      phi(3) = phi(3) + Ud(1)*Ud(3)*VdotN + Pd*normal(2)
      phi(4) = phi(4) + Ud(1)*Ud(4)*VdotN + Pd*normal(3)
      phi(5) = phi(5) + Ud(1)*VdotN*Hd + Pd*vitno

c
c Computation of the Roe-averaged state
c
      squsrg = DSQRT(Ug(1))
      squsrd = DSQRT(Ud(1))
      usro  = 1.0d0/(squsrg + squsrd)


c With the analytical implementation, oscillations appeared as 
c soon as the relative difference between Ud(1) and Ug(1) was
c less than 1e-7.so a first order expansion is enough to get 1e-14 
c precision


      if (DABS(Ud(1)-Ug(1)).le.0.005d0*Ug(1)) then
       coeff1 = 0.5d0
       coeff2 = (beta-2.0d0)/24.0d0
       coeff3 = -(beta-2.0d0)/48.0d0
       coeff4 = -(2.0d0*beta**3-3.0d0*beta**2-78.0d0*beta+152.0d0)
     &              / 5760.0d0
       coeff5 = (2.0d0*beta**3-3.0d0*beta**2-38.0d0*beta+72.0d0)
     &              /3840.0d0
       eps                       = (Ud(1)-Ug(1))/Ug(1)
       uar1                      =  1.0d0 + coeff1*eps 
     &                             + coeff2*eps**2+ coeff3*eps**3
     &                             + coeff4*eps**4+ coeff5*eps**5
       uar1                      = Ug(1)*uar1
      else
       uar1                      = ( (Udbeta - Ugbeta) / 
     &                            ( beta * (Ud(1) - Ug(1)) ) ) **
     &                            ( 1.0d0 / beta1 )
      endif

      uar2                      = (squsrg*Ug(2) +
     &                                squsrd*Ud(2))*usro

      uar3                      = (squsrg*Ug(3) +
     &                                squsrd*Ud(3))*usro

      uar4                      = (squsrg*Ug(4) +
     &                                squsrd*Ud(4))*usro

      uar5                      = (squsrg*Hg + 
     &                                squsrd*Hd)*usro
     
      uar6                      = Pr +alpha*uar1**beta  



c
c Computation of the dissipation term
c

c  the dif term uses (P,u,v,w,T)
      dif1 = -alpha*(Ugr(1)**beta-Udr(1)**beta)
      dif2 = -Ugr(2) + Udr(2)
      dif3 = -Ugr(3) + Udr(3)
      dif4 = -Ugr(4) + Udr(4)
      dif5 = -Ugr(5) + Udr(5)

c  values needed by the different matrices
      VdotN = normal(1)*uar2 + normal(2)*uar3 + normal(3)*uar4
      cr2 = alpha*beta*( uar1**beta1 )
      cr  = DSQRT(cr2)
      
c  set preconditioner beta
      if(prec.eq.0) then
        betaPrec = 1.0d0
      else
        locMach = sqrt((uar2*uar2+uar3*uar3+uar4*uar4)/cr2)
        betaPrec = MAX(betaRef, k1*locMach)
        betaPrec = MIN((1.0d0+DSQRT(irey))*betaPrec, cmach)
      endif

      b2 = betaPrec**2
      oob2 = 1.0d0/b2

      tet1 = normal(3)*uar3 - normal(2)*uar4
      tet2 = normal(1)*uar4 - normal(3)*uar2
      tet3 = normal(2)*uar2 - normal(1)*uar3

c  eigenvalues - Preconditioned
      vp1 = VdotN
      vpa = 0.5d0*(1.0d0+b2)*VdotN
      vpb = 0.5d0*dsqrt(4.0d0*b2*cr2+(VdotN*(1.0d0-b2))**2)
      vp4 = vpa+vpb
      vp5 = vpa-vpb

c Dynamic Mesh inclusion
      vp1 = vp1 - vitno
      vp4 = vp4 - vitno
      vp5 = vp5 - vitno
      
      
c  useful values
      r = VdotN - vp4
      s = VdotN - vp5
      rs = -b2*cr2
      rps = VdotN*(1.0d0-b2)
      smr = s-r
      r2 = r*r
      s2 = s*s
      rp = uar1*r/(s*smr)
      sp = -uar1*s/(r*smr)
      P1 = uar6/(uar1**2*rs*Cv)
      P2 = P1*uar1*rps
      Proeor = uar6/uar1
      

      flur1 = dabs(vp1)*(
     &          P1*nx        *dif1  
     &         +(P2*nx*nx   )*dif2
     &         +(P2*nx*ny-nz)*dif3
     &         +(P2*nx*nz+ny)*dif4
     &         +nx           *dif5 )

      flur2 = dabs(vp1)*(
     &          P1*ny        *dif1  
     &         +(P2*ny*nx+nz)*dif2
     &         +(P2*ny*ny   )*dif3
     &         +(P2*ny*nz-nx)*dif4
     &         +ny           *dif5 )
     
      flur3 = dabs(vp1)*(
     &          P1*nz        *dif1  
     &         +(P2*nz*nx-ny)*dif2
     &         +(P2*nz*ny+nx)*dif3
     &         +(P2*nz*nz   )*dif4
     &         +nz           *dif5 )
    
      flur4 = dabs(vp4)/cr2*(
     &        -dif1/(r*smr)+sp*(nx*dif2+ny*dif3+nz*dif4) )

      flur5 = dabs(vp5)/cr2*(
     &         dif1/(s*smr)+rp*(nx*dif2+ny*dif3+nz*dif4) )



      diss1 = (r2*flur4+s2*flur5)*oob2

      diss2 = uar1*(nz*flur2-ny*flur3)
     &      + (oob2*uar2*r2-r*cr2*nx)*flur4
     &      + (oob2*uar2*s2-s*cr2*nx)*flur5
   
      diss3 = uar1*(nx*flur3-nz*flur1)
     &      + (oob2*uar3*r2-r*cr2*ny)*flur4
     &      + (oob2*uar3*s2-s*cr2*ny)*flur5

      diss4 = uar1*(ny*flur1-nx*flur2)
     &      + (oob2*uar4*r2-r*cr2*nz)*flur4
     &      + (oob2*uar4*s2-s*cr2*nz)*flur5

      diss5 = uar1*Cv*(flur1*nx + flur2*ny + flur3*nz)
     &      - uar1*(tet1*flur1 + tet2*flur2 + tet3*flur3)
     &      + (oob2*uar5*r2+Proeor*(cr2-oob2*r2)
     &                - VdotN*cr2*r)*flur4
     &      + (oob2*uar5*s2+Proeor*(cr2-oob2*s2)
     &                - VdotN*cr2*s)*flur5

      phi(1) = phi(1) - gamma*diss1
      phi(2) = phi(2) - gamma*diss2
      phi(3) = phi(3) - gamma*diss3 
      phi(4) = phi(4) - gamma*diss4
      phi(5) = phi(5) - gamma*diss5

      phi(1) = phi(1)*0.5d0*rnorm  
      phi(2) = phi(2)*0.5d0*rnorm
      phi(3) = phi(3)*0.5d0*rnorm
      phi(4) = phi(4)*0.5d0*rnorm
      phi(5) = phi(5)*0.5d0*rnorm

      END


