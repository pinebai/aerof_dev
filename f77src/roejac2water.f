      SUBROUTINE ROEFLUX1WATER(gamma,Cv,Pr,alpha,beta,enormal,
     &               evitno,Ug,Ud,phi1)
c-----------------------------------------------------------------------
c This routine is only used for TURBULENCE and UNCOUPLED SOLVERS
c-----------------------------------------------------------------------
c This routine computes the Flux of Roe taken at the vectors Ug, Ud
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
c-----------------------------------------------------------------------
c This routine computes the same thing as ROEFLUX6 except that it only 
c computes the flux in density and with Ug=Ugr , Ud=Udr (these two equalities
c are always true in our code!)
c-----------------------------------------------------------------------


      IMPLICIT NONE
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), evitno, phi1
      REAL*8 Hg, Hd, vitno, updir, gamma
      REAL*8 VdotN , rnorm, invnorm
      REAL*8 invnormal1
      REAL*8 flur1 , flur2 , flur3 , flur4 , flur5
      REAL*8 dif1 , dif2 , dif3 , dif4 ,dif5
      REAL*8 cr , cr2 , qir, tet1, tet2, tet3, u2mh
      REAL*8 vp1 , vp4 , vp5
      REAL*8 ener1 , ener2
      REAL*8 uar1 , uar2 , uar3 , uar4 , uar5
      REAL*8 usro , squsr1 , squsr2
      REAL*8 vitg2, vitd2
      REAL*8 Pr, alpha, beta, Cv, Pg, Pd
      REAL*8 beta1, Ugbeta, Udbeta
      REAL*8 coeff1, coeff2, coeff3, coeff4, coeff5, eps
      INTEGER type



c
c Initialisation
c

      beta1 = beta - 1.0d0

      rnorm = DSQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) +
     &             enormal(3)*enormal(3))
      invnorm = 1.0d0 / rnorm
c
      normal(1) = enormal(1) * invnorm
      normal(2) = enormal(2) * invnorm
      normal(3) = enormal(3) * invnorm

      vitno = evitno * invnorm

c
c Computation of the centred terms
c

      VdotN = Ug(2)*normal(1) + Ug(3)*normal(2) + Ug(4)*normal(3)
      vitg2 = Ug(2)*Ug(2) + Ug(3)*Ug(3) + Ug(4)*Ug(4)
      Ugbeta = Ug(1)**beta
      Pg = Pr + alpha*Ugbeta
      Hg = Cv*Ug(5) + Pg/Ug(1) + 0.5*vitg2
      phi1 = Ug(1)*(VdotN - vitno)
  
c
      VdotN = Ud(2)*normal(1)+Ud(3)*normal(2)+Ud(4)*normal(3)
      VdotN = VdotN - vitno
      vitd2 = Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4)
      Udbeta = Ud(1)**beta
      Pd = Pr + alpha*Udbeta
      Hd = Cv*Ud(5) + Pd/Ud(1) +0.5*vitd2
      phi1 = phi1 + Ud(1)*VdotN
 
c
c Computation of the Roe-averaged state
c
      squsr1                    = DSQRT(Ug(1))
      squsr2                    = DSQRT(Ud(1))
c     
      ener1                     = Ug(1) * ( Cv*Ug(5) + 0.5d0*vitg2 ) 
c
      ener2                     = Ud(1) * ( Cv*Ud(5) + 0.5d0*vitd2 ) 
c
      usro                      = 1.0d0/(squsr1 + squsr2)
c
      
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
c
      uar2                      = (squsr1*Ug(2) +
     &                                squsr2*Ud(2))*usro
c
      uar3                      = (squsr1*Ug(3) +
     &                                squsr2*Ud(3))*usro
c
      uar4                      = (squsr1*Ug(4) +
     &                                squsr2*Ud(4))*usro
c
      uar5                      = (squsr1*Hg + 
     &                                squsr2*Hd)*usro
c
c Computation of the dissipation term
c
      VdotN                     = normal(1)*uar2 + normal(2)*uar3 +
     &                               normal(3)*uar4
c
      qir                       = 0.5d0*(uar2*uar2 + uar3*uar3 +
     &                                    uar4*uar4)
c
      cr2                       = alpha*beta*(uar1**beta1)
      cr                        = DSQRT(cr2)
c
      tet1 = normal(3)*uar3 - normal(2)*uar4
      tet2 = normal(1)*uar4 - normal(3)*uar2
      tet3 = normal(2)*uar2 - normal(1)*uar3
c
      u2mh  = 2.0d0*qir - uar5

 
      dif1                     = - Ug(1) + Ud(1)
      dif2                     = - Ug(1)*Ug(2) + Ud(1)*Ud(2)
      dif3                     = - Ug(1)*Ug(3) + Ud(1)*Ud(3)
      dif4                     = - Ug(1)*Ug(4) + Ud(1)*Ud(4)
      dif5                     = - ener1 + ener2
c
      vp1                       = gamma * DABS(VdotN - vitno)
      vp4                       = gamma * DABS((VdotN + cr) - vitno)
      vp5                       = gamma * DABS((VdotN - cr) - vitno)
      flur4                  = vp4*0.5d0*
     &                       (    (1.0d0 - VdotN/cr) * dif1 +
     &                             normal(1)/cr  * dif2 +
     &                             normal(2)/cr  * dif3 +
     &                             normal(3)/cr  * dif4  )
    
c
      flur5                  = vp5*0.5d0*
     &                       (    (1.0d0 + VdotN/cr) * dif1 -
     &                             normal(1)/cr  * dif2 -
     &                             normal(2)/cr  * dif3 -
     &                             normal(3)/cr  * dif4  )

c
      phi1 =  phi1 - (flur4 + flur5)

c 
      phi1 = phi1*0.5d0*rnorm  
c
      END






      SUBROUTINE ROEJAC2WATER(type,gamma,gam,enormal,evitno
     &                       ,Ug,Ud,jacg,jacd)
c-----------------------------------------------------------------------
c This routine computes the Flux of Roe taken at the vectors Ug, Ud
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 gam, enormal(3), evitno, Ug(*), Ud(*), gamma
      REAL*8 jacg(*), jacd(*), phi1, updir, jacleft, jacright
      INTEGER type
c
      CALL ROEFLUX1(gamma, gam, enormal, evitno, Ug, Ud, phi1)

      updir = 0.5d0 + dsign(0.5d0, phi1)

      jacleft = phi1 * updir / Ug(1)
      jacright = phi1 * (1.0d0 - updir) / Ud(1)

c      if (type.eq.1) then
c         jacg(1) = jacleft
c         jacd(1) = jacright
c      else if (type.eq.2) then
c         jacg(1) = jacleft
c         jacg(2) = 0.0
c         jacg(3) = 0.0
c         jacg(4) = jacg(1)

c         jacd(1) = jacright
c         jacd(2) = 0.0
c         jacd(3) = 0.0
c         jacd(4) = jacd(1)
c      endif

      END
