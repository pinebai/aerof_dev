      SUBROUTINE ROEFLUX5JWL(type,gamma,omega,A1,A2,R1r,R2r,
     &     enormal,evitno,
     &     Ugr,Ug,Udr,Ud,phi,mach,k1,cmach,shockreducer,
     &     irey,length,prec)
c-----------------------------------------------------------------------
c This routine computes the Flux of Roe taken at the vectors Ug, Ud
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
c gamma is the dissipation coefficient
c The Roe-Turkel Preconditioning is applied for LowMach Simulations
c The equation of state is the JWL one for burned gas
c The flux is computed by phi = 0.5(F(Ug)+F(Ud) - P^{-1}|PA|DeltaW) 
c where W are the conservative variable and 
c       P is the low-mach preconditioner
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 omega,A1,A2,R1r,R2r
      REAL*8 ooomega,omegap1,frhog,frhod,frhogr,frhodr
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), evitno, phi(*)
      REAL*8 Ugr(*), Udr(*), energ, enerd
      REAL*8 Hg, Hd, vitno, updir, gamma 
      REAL*8 VdotN , rnorm, invnorm
      REAL*8 flur1 , flur2 , flur3 , flur4 , flur5
      REAL*8 dif1 , dif2 , dif3 , dif4 ,dif5
      REAL*8 cr , cr2 , qir
      REAL*8 qirmfpuar1, omegacr2
      REAL*8 vp1 , vp4 , vp5
      REAL*8 uar2 , uar3 , uar4 , uar5, fpuar1
      REAL*8 usro , squsr1 , squsr2
      REAL*8 tet1 , tet2 , tet3
      REAL*8 vitg2, vitd2
      REAL*8 r,s,t,beta,beta2
      REAL*8 mach,k1,shock,locMach, cmach, irey, length
      REAL*8 temp1, temp2, shockreducer
      INTEGER type, prec

c
c Initialisation
c
      rnorm = DSQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) + 
     &             enormal(3)*enormal(3))
      
      invnorm = 1.0d0 / rnorm
c
      normal(1) = enormal(1) * invnorm
      normal(2) = enormal(2) * invnorm
      normal(3) = enormal(3) * invnorm

      vitno = evitno * invnorm

      ooomega = 1.0/omega
      omegap1 = omega + 1.0

      frhog =A1*(1.0-omega*Ug(1)/R1r)*exp(-R1r/Ug(1))
     &      +A2*(1.0-omega*Ug(1)/R2r)*exp(-R2r/Ug(1))
      frhod =A1*(1.0-omega*Ud(1)/R1r)*exp(-R1r/Ud(1))
     &      +A2*(1.0-omega*Ud(1)/R2r)*exp(-R2r/Ud(1))
      frhogr =A1*(1.0-omega*Ugr(1)/R1r)*exp(-R1r/Ugr(1))
     &       +A2*(1.0-omega*Ugr(1)/R2r)*exp(-R2r/Ugr(1))
      frhodr =A1*(1.0-omega*Udr(1)/R1r)*exp(-R1r/Udr(1))
     &       +A2*(1.0-omega*Udr(1)/R2r)*exp(-R2r/Udr(1))
      
c
c Computation of the centred terms
c
      VdotN = Ug(2)*normal(1) + Ug(3)*normal(2) + Ug(4)*normal(3)
      vitg2 = Ug(2)*Ug(2) + Ug(3)*Ug(3) + Ug(4)*Ug(4)
      Hg = (omegap1*Ug(5)+0.5d0*omega*Ug(1)*vitg2-frhog)/(omega*Ug(1))
      phi(1) = Ug(1)*(VdotN - vitno)
      phi(2) = phi(1)*Ug(2) + Ug(5)*normal(1)
      phi(3) = phi(1)*Ug(3) + Ug(5)*normal(2)
      phi(4) = phi(1)*Ug(4) + Ug(5)*normal(3)
      phi(5) = phi(1)*Hg + Ug(5)*vitno
c
      VdotN = Ud(2)*normal(1)+Ud(3)*normal(2)+Ud(4)*normal(3)
      VdotN = VdotN - vitno
      vitd2 = Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4)
      Hd = (omegap1*Ud(5)+0.5d0*omega*Ud(1)*vitd2-frhod)/(omega*Ud(1))
      phi(1) = phi(1) + Ud(1)*VdotN
      phi(2) = phi(2) + Ud(1)*Ud(2)*VdotN +Ud(5)*normal(1)
      phi(3) = phi(3) + Ud(1)*Ud(3)*VdotN +Ud(5)*normal(2)
      phi(4) = phi(4) + Ud(1)*Ud(4)*VdotN +Ud(5)*normal(3)
      phi(5) = phi(5) + Ud(1)*VdotN*Hd +Ud(5)*vitno
c
c Computation of the Roe-averaged state
c

      squsr1   = DSQRT(Ug(1))
      squsr2   = DSQRT(Ud(1))
c     
      usro     = 1.d0/(squsr1 + squsr2)
c
      uar2     = (squsr1*Ug(2) + squsr2*Ud(2))*usro
c
      uar3     = (squsr1*Ug(3) + squsr2*Ud(3))*usro
c
      uar4     = (squsr1*Ug(4) + squsr2*Ud(4))*usro
c
      uar5     = (squsr1*Hg    + squsr2*Hd   )*usro
c uar1 needs not be computed since we only need f'(uar1) 
c (see computation below)


c
c Computation of the dissipation term 
c if prec = 1  then the dissipation is preconditioned
c else if prec = 0 then it is not preconditioned
c
c Reference: Implicit Upwind Schemes for Lowmach number Compressible Flows
c            By Cecile Viozat (INRIA Publication)


      VdotN    = normal(1)*uar2 + normal(2)*uar3 + normal(3)*uar4
c
      qir      = 0.5d0*(uar2*uar2 + uar3*uar3 + uar4*uar4)
c
      tet1     = normal(3)*uar3 - normal(2)*uar4
      tet2     = normal(1)*uar4 - normal(3)*uar2
      tet3     = normal(2)*uar2 - normal(1)*uar3
c
      temp1 = 0.0
      temp2 = 0.0
      if (ABS(Ug(1) - Ud(1)) .lt. 1.0e-14*ABS(Ud(1))) then
        temp1 = A1*(-omega/R1r+(1.0-omega*Ug(1)/R1r)*R1r/Ug(1)**2)
     &             *exp(-R1r/Ug(1))
     &         + A2*(-omega/R2r+(1.0-omega*Ug(1)/R2r)*R2r/Ug(1)**2)
     &             *exp(-R2r/Ug(1))
        fpuar1 = temp1
      else
        temp2 = 0.0
        temp2 = (frhod-frhog)/(Ud(1)-Ug(1))
        fpuar1 = temp2
      endif
c
      cr2      = omega*(uar5-qir)+fpuar1
c      if (cr2 .le. 0.0) then
c        write(*,*) 'rho is ', Ug(1), Ud(1), Ud(1)-Ug(1)
c        write(*,*) 'done computing frho', frhog, frhod, frhogr, frhodr
c        write(*,*) 'pressure is ', Ug(5), Ud(5)
c        write(*,*) 'computing cr (cr2 =', cr2, ', fpuar1 =', fpuar1,')'
c        write(*,*) 'temp1 = ', temp1, ' and temp2 = ', temp2
c      endif
      cr       = DSQRT(cr2)
      cr2      = 1.d0/cr2
c
      qirmfpuar1 = qir - fpuar1*ooomega
      omegacr2   = omega*cr2
        
      if (prec .eq. 0) then
        beta = 1.d0
      else
c       local Preconditioning (ARL)
        shock = DABS(Ugr(5) - Udr(5))/(Ugr(5)+Udr(5))/length
        locMach = DSQRT(2.0d0*qir*cr2)
        beta = MAX(k1*locMach, mach)
        beta = (1.0d0+DSQRT(irey))*beta+shockreducer*shock
        beta = MIN(beta, cmach)
      end if
      
      beta2 = beta * beta 

      energ = (Ugr(5)-frhogr)*ooomega + 
     &     0.5d0*Ugr(1)*(Ugr(2)*Ugr(2) + Ugr(3)*Ugr(3) + Ugr(4)*Ugr(4))
      enerd = (Udr(5)-frhodr)*ooomega + 
     &     0.5d0*Udr(1)*(Udr(2)*Udr(2) + Udr(3)*Udr(3) + Udr(4)*Udr(4))


      dif1   = - Ugr(1) + Udr(1)
      dif2   = - Ugr(1)*Ugr(2) + Udr(1)*Udr(2)
      dif3   = - Ugr(1)*Ugr(3) + Udr(1)*Udr(3)
      dif4   = - Ugr(1)*Ugr(4) + Udr(1)*Udr(4)
      dif5   = - energ + enerd

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

c Temporary variables temp1 and temp2 to avoid recomputation
      temp1 = uar5 - 2.0d0*qir

      flur1                  = DABS(vp1)*
     &           ((normal(1)*(omegacr2*temp1)-tet1)*dif1 +
     &            (normal(1)*omegacr2*uar2)*dif2  +
     &            (normal(3)  + (normal(1)*omegacr2*uar3))*dif3   +
     &            (-normal(2) + (normal(1)*omegacr2*uar4))*dif4   -
     &            (normal(1)*omegacr2)*dif5)

      flur2                  = DABS(vp1)*
     &           ((normal(2)*(omegacr2*temp1)-tet2)*dif1 +
     &            (-normal(3) + (normal(2)*omegacr2*uar2))*dif2   +
     &            (normal(2)*omegacr2*uar3)*dif3  +
     &            (normal(1)  + (normal(2)*omegacr2*uar4))*dif4   -
     &            (normal(2)*omegacr2)*dif5)

      flur3                  = DABS(vp1)*
     &           ((normal(3)*(omegacr2*temp1)-tet3)*dif1 +
     &            (normal(2)  + (normal(3)*omegacr2*uar2))*dif2   +
     &            (-normal(1) + (normal(3)*omegacr2*uar3))*dif3   +
     &            (normal(3)*omegacr2*uar4)*dif4  -
     &            (normal(3)*omegacr2)*dif5)

      temp1 = omega*s/beta2
       
       flur4                     = DABS(vp4)* 
     &          ((cr**2*VdotN + s/beta2*(omega*qir+fpuar1))*dif1  -
     &           (cr**2*normal(1) + temp1*uar2)*dif2 -
     &           (cr**2*normal(2) + temp1*uar3)*dif3 -
     &           (cr**2*normal(3) + temp1*uar4)*dif4 +
     &            temp1*dif5)/t

      temp1 = omega*r/beta2
 
       flur5                     = DABS(vp5)*  
     &          (-(cr**2*VdotN + r/beta2*(omega*qir+fpuar1))*dif1 +
     &           (cr**2*normal(1) + temp1*uar2)*dif2 +
     &           (cr**2*normal(2) + temp1*uar3)*dif3 +
     &           (cr**2*normal(3) + temp1*uar4)*dif4 -
     &            temp1*dif5)/t

c
c Final phi including the numerical viscosity parameter
c
      temp1 = flur4+flur5
      temp2 = r*flur4+s*flur5

      phi(1) =  phi(1) - gamma*(normal(1)*flur1 +
     &     normal(2)*flur2 + normal(3)*flur3 +
     &     0.5d0*temp1*cr2)
c

      phi(2) = phi(2) - gamma*(  
     &              (uar2*normal(1))*flur1 +
     &              (uar2*normal(2) - normal(3))*flur2 +
     &              (uar2*normal(3) + normal(2))*flur3 +
     &              (normal(1)*temp2 + uar2*temp1)*0.5d0*cr2)

c
      phi(3) = phi(3) - gamma*(
     &              (uar3*normal(1) + normal(3))*flur1 +
     &              (uar3*normal(2))*flur2 +
     &              (uar3*normal(3) - normal(1))*flur3 +
     &              (normal(2)*temp2 + uar3*temp1)*0.5d0*cr2)

c
      phi(4) = phi(4) - gamma*( 
     &              (uar4*normal(1) - normal(2))*flur1 +
     &              (uar4*normal(2) + normal(1))*flur2 +
     &              (uar4*normal(3))*flur3 +
     &              (normal(3)*temp2 + uar4*temp1)*0.5d0*cr2)

c
      phi(5) = phi(5) - gamma*(
     &              (qirmfpuar1*normal(1) + tet1)*flur1 +
     &              (qirmfpuar1*normal(2) + tet2)*flur2 +
     &              (qirmfpuar1*normal(3) + tet3)*flur3 +
     &              (VdotN*temp2 + uar5*temp1)*0.5d0*cr2)

      phi(1) = phi(1)*0.5d0*rnorm
      phi(2) = phi(2)*0.5d0*rnorm
      phi(3) = phi(3)*0.5d0*rnorm
      phi(4) = phi(4)*0.5d0*rnorm
      phi(5) = phi(5)*0.5d0*rnorm

c
c For one and two equation turbulence models
c

      if (type.eq.1) then
         updir = 0.5d0 + dsign(0.5d0, phi(1))
         phi(6) = phi(1) * (updir * Ugr(6) + (1.0d0 - updir) * Udr(6))
      else if (type.eq.2) then
         updir = 0.5d0 + dsign(0.5d0, phi(1))
         phi(6) = phi(1) * (updir * Ugr(6) + (1.0d0 - updir) * Udr(6))
         phi(7) = phi(1) * (updir * Ugr(7) + (1.0d0 - updir) * Udr(7))
      endif
      
      END
