      SUBROUTINE ROEFLUX1(gamma,gam,pstiff,enormal,evitno,
     &     Ug,Ud,phi1,mach,k1,cmach,irey,prec)
c-----------------------------------------------------------------------
c This routine computes the Flux of Roe taken at the vectors Ug, Ud
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
c gamma is the dissipation coefficient
c gam is the ratio of cp/cv
c The Roe-Turkel Preconditioning is applied for LowMach Simulations
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), evitno
      REAL*8 phi1
      REAL*8 energ, enerd
      REAL*8 H, vitno, updir, gamma 
      REAL*8 VdotN , rnorm, invnorm
      REAL*8 flur1 , flur2 , flur3 , flur4 , flur5
      REAL*8 dif1 , dif2 , dif3 , dif4 ,dif5
      REAL*8 cr , cr2 , qir
      REAL*8 vp1 , vp4 , vp5
      REAL*8 ener1 , ener2
      REAL*8 uar1 , uar2 , uar3 , uar4 , uar5
      REAL*8 usro , squsr1 , squsr2
      REAL*8 tet1 , tet2 , tet3
      REAL*8 gam , gam1, vitg2, vitd2, pstiff
      REAL*8 r,s,t,A,B,H1,beta, mach00,beta2
      REAL*8 maxu2,maxrho,minpres,mach,k1
      REAL*8 locMach, cmach, irey
      INTEGER type, prec

c
c Initialisation
c
      gam1 = gam - 1.d0      
c
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
      phi1 = Ug(1)*(VdotN - vitno)
c
      VdotN = Ud(2)*normal(1)+Ud(3)*normal(2)+Ud(4)*normal(3)
      VdotN = VdotN - vitno
      vitd2 = Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4)
      phi1 = phi1 + Ud(1)*VdotN


c
c Computation of the Roe-averaged state
c

      squsr1                    = DSQRT(Ug(1))
      squsr2                    = DSQRT(Ud(1))
c     
      ener1                     = (Ug(5)+gam*pstiff)/gam1 +
     &                            0.5d0*Ug(1)*vitg2
c
      ener2                     = (Ud(5)+gam*pstiff)/gam1 +
     &                            0.5d0*Ud(1)*vitd2
c
      usro                      = 1.d0/(squsr1 + squsr2)
c
      uar1                      = (squsr1*Ug(1) +
     &                                  squsr2*Ud(1))*usro
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
      uar5                      = ((ener1 + Ug(5))/
     &                                squsr1 +
     &                                (ener2 + Ud(5))/
     &                                squsr2)*usro


c
c Computation of the dissipation term 
c if prec = 1  then the dissipation is preconditioned
c else if prec = 0 then it is not preconditioned
c
c Reference: Implicit Upwind Schemes for Lowmach number Compressible Flows
c            By Cecile Viozat (INRIA Publication)


      VdotN                     = normal(1)*uar2 + normal(2)*uar3 +
     &                               normal(3)*uar4
c
      qir                       = 0.5d0*(uar2*uar2 + uar3*uar3 +
     &                                    uar4*uar4)
c
      tet1                      = normal(3)*uar3 - normal(2)*uar4
      tet2                      = normal(1)*uar4 - normal(3)*uar2
      tet3                      = normal(2)*uar2 - normal(1)*uar3
c
      cr2                       = gam1*(uar5 - qir)
      cr                        = DSQRT(cr2)
      cr2                       = 1.d0/cr2
c
        
      if (prec .eq. 0) then
        beta = 1.d0
      else
c
c     Better way to compute beta
c        local Preconditioning (ARL)
c        beta = MIN(MAX(0.000001d0,beta),1.d0)
      locMach = DSQRT(2.0d0*qir*cr2)
      beta = MAX(k1*locMach, mach)
      beta = (1.0d0+DSQRT(irey))*beta
      beta = MIN(beta, cmach)
c
c     Crude way to compute beta
c
c        beta = mach 
      end if
      
      beta2 = beta * beta 

      energ = (Ug(5)+gam*Pstiff)/gam1 + 
     &     0.5d0*Ug(1)*(Ug(2)*Ug(2) + Ug(3)*Ug(3) + Ug(4)*Ug(4))
      enerd = (Ud(5)+gam*Pstiff)/gam1 + 
     &     0.5d0*Ud(1)*(Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4))


      dif1                     = - Ug(1)+Ud(1)
      dif2                     = - Ug(1)*Ug(2) + Ud(1)*Ud(2)
      dif3                     = - Ug(1)*Ug(3) + Ud(1)*Ud(3)
      dif4                     = - Ug(1)*Ug(4) + Ud(1)*Ud(4)
      dif5                     = - energ + enerd

      vp1 = VdotN
      vp4 = 0.5d0*((1.d0+beta2)*VdotN +
     &      DSQRT(((1.d0-beta2)*VdotN)**2 +
     &      4.d0*beta2*cr**2))
      vp5 = 0.5d0*((1.d0+beta2)*VdotN - 
     &      DSQRT(((1.d0-beta2)*VdotN)**2 +
     &      4.d0*beta2*cr**2))

c Roe-Turkel coefficients

      r = vp4 - vp1*beta2
      s = vp5 - vp1*beta2
      t = 0.5d0*(vp5-vp4) 

c Dynamic mesh inclusion 

      vp1 = (vp1-vitno)
      vp4 = (vp4-vitno)
      vp5 = (vp5-vitno)


c
      flur1                     = DABS(vp1)*
     &              ((normal(1)*(1.d0 - gam1*qir*cr2) - tet1)*dif1 +
     &               (normal(1)*gam1*uar2*cr2)*dif2  +
     &               (normal(3)  + (normal(1)*gam1*uar3*cr2))*dif3   +
     &               (-normal(2) + (normal(1)*gam1*uar4*cr2))*dif4   -
     &               (normal(1)*gam1*cr2)*dif5)
c
      flur2                     = DABS(vp1)*
     &              ((normal(2)*(1.d0 - gam1*qir*cr2) - tet2)*dif1 +
     &               (-normal(3) + (normal(2)*gam1*uar2*cr2))*dif2   +
     &               (normal(2)*gam1*uar3*cr2)*dif3  +
     &               (normal(1)  + (normal(2)*gam1*uar4*cr2))*dif4   -
     &               (normal(2)*gam1*cr2)*dif5)
c
      flur3                     = DABS(vp1)*
     &              ((normal(3)*(1.d0 - gam1*qir*cr2) - tet3)*dif1 +
     &               (normal(2)  + (normal(3)*gam1*uar2*cr2))*dif2   +
     &               (-normal(1) + (normal(3)*gam1*uar3*cr2))*dif3   +
     &               (normal(3)*gam1*uar4*cr2)*dif4  -
     &               (normal(3)*gam1*cr2)*dif5)

c Roe-Turkel
       
       flur4                     = DABS(vp4)* 
     &          ((cr**2*VdotN/t + gam1*qir*s/(t*beta2))*dif1  -
     &           (cr**2*normal(1)/t + gam1*uar2*s/(t*beta2))*dif2 -
     &           (cr**2*normal(2)/t + gam1*uar3*s/(t*beta2))*dif3 -
     &           (cr**2*normal(3)/t + gam1*uar4*s/(t*beta2))*dif4 +
     &            gam1*s/(t*beta2)*dif5)
      

 
       flur5                     = DABS(vp5)*  
     &          (-(cr**2*VdotN/t  + gam1*qir*r/(beta2*t))*dif1 +
     &           (cr**2*normal(1)/t + gam1*uar2*r/(beta2*t))*dif2 +
     &           (cr**2*normal(2)/t + gam1*uar3*r/(beta2*t))*dif3 +
     &           (cr**2*normal(3)/t + gam1*uar4*r/(beta2*t))*dif4 -
     &            gam1*r/(beta2*t)*dif5)


c
c Final phi including the numerical viscosity parameter
c

      phi1 =  phi1 - gamma*(normal(1)*flur1 +
     &     normal(2)*flur2 + normal(3)*flur3 +
     &     0.5d0*(flur4 + flur5)*cr2)


      phi1 = phi1*0.5d0*rnorm


      END


      SUBROUTINE ROEJAC2(type,gamma,gam,pstiff,enormal,evitno,
     &     Ug,Ud,jacg,jacd,mach,k1,cmach,irey,prec)
c-----------------------------------------------------------------------
c This routine computes the Flux of Roe taken at the vectors Ug, Ud
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 gam, pstiff, enormal(3), evitno, Ug(*), Ud(*), gamma
      REAL*8 mach, k1, cmach, irey
      REAL*8 jacg(*), jacd(*), phi1, updir, jacleft, jacright
      INTEGER prec
      INTEGER type
c
c      CALL ROEFLUX1(gamma, gam, enormal, evitno, Ug, Ud, phi1)
      CALL ROEFLUX1(gamma,gam,pstiff,enormal,evitno,
     &     Ug,Ud,phi1,mach,k1,cmach,irey,prec)

      updir = 0.5d0 + dsign(0.5d0, phi1)

      jacleft = phi1 * updir / Ug(1)
      jacright = phi1 * (1.0d0 - updir) / Ud(1)

      if (type.eq.1) then
         jacg(1) = jacleft
         jacd(1) = jacright
      else if (type.eq.2) then
         jacg(1) = jacleft
         jacg(2) = 0.0
         jacg(3) = 0.0
         jacg(4) = jacg(1)

         jacd(1) = jacright
         jacd(2) = 0.0
         jacd(3) = 0.0
         jacd(4) = jacd(1)
      endif

      END
