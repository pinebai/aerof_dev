      SUBROUTINE ROEJAC6JWL(type,gamma,omega,A1,A2,R1r,R2r,
     &                    enormal,evitno,Ug,Ud,jac,sw,
     &                    mach,k1,cmach,shockreducer,
     &                    irey,length,prec)
c-----------------------------------------------------------------------
c This routine computes the jacobian of the flux of Roe (defect correction) 
c with respect to the conservative variables taken at the vectors Ug, Ud.
c normal is the normal of the boundary concerned by the flux.
c The JWL EOS for burned gas is considered here.
c The Roe-Turkel Low-Mach preconditioning can be applied.
c-----------------------------------------------------------------------
c    INCLUDE 'Param3D.h'
      implicit none
c-----------------------------------------------------------------------
      REAL*8 omega,A1,A2,R1r,R2r
      REAL*8 ooomega,omegap1,frhog,frhod,fprho
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), H(5,5), jac(*)
      REAL*8 tmw(5,5),dmat(5,5)
      REAL*8 VdotN , rnorm, evitno, vitno, invnorm
      REAL*8 g_VdotN(5)
      REAL*8 flur1 , flur2 , flur3 , flur4 , flur5
      REAL*8 dif1 , dif2 , dif3 , dif4 ,dif5
      REAL*8 cr , cr2 , qir, fpuar1
      REAL*8 qirmfpuar1, omegacr2
      REAL*8 vp1 , vp4 , vp5 , svp1 , svp4 , svp5
      REAL*8 enerd, energ, enth
      REAL*8 uar2 , uar3 , uar4 , uar5
      REAL*8 usro , squsr1 , squsr2
      REAL*8 tet1 , tet2 , tet3
      REAL*8 Hg,Hd,usc2,temp1
      INTEGER k,i,type,dim,j
      REAL*8 gamma
      REAL*8 r,s,t,beta,beta2
      REAL*8 mach,k1,locMach,cmach,irey,shock,shockreducer,length
      INTEGER prec,sw

c
c Initialisation
c
      rnorm = DSQRT(enormal(1)**2 + enormal(2)**2 + enormal(3)**2)
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
      fprho = A1*(-omega/R1r+(1.0-omega*Ug(1)/R1r)*R1r/Ug(1)**2)
     &          *exp(-R1r/Ug(1))
     &      + A2*(-omega/R2r+(1.0-omega*Ug(1)/R2r)*R2r/Ug(1)**2)
     &          *exp(-R2r/Ug(1))
      
c
c Computation of the centred terms
c
      VdotN = Ug(2)*normal(1) + Ug(3)*normal(2) + Ug(4)*normal(3)
      qir = 0.5d0*(Ug(2)**2 + Ug(3)**2 + Ug(4)**2)
      enth = (omegap1*Ug(5)+omega*Ug(1)*qir-frhog)/(omega*Ug(1))
       
c
      H(1,1) = - vitno
      H(1,2) = normal(1)
      H(1,3) = normal(2)
      H(1,4) = normal(3)
      H(1,5) = 0.0d0 
c
      H(2,1) = (omega*qir+fprho)*normal(1) - Ug(2)*VdotN 
      H(2,2) = (1.d0-omega)*normal(1)*Ug(2) + VdotN - vitno 
      H(2,3) = -omega*normal(1)*Ug(3) + Ug(2)*normal(2)
      H(2,4) = -omega*normal(1)*Ug(4) + Ug(2)*normal(3) 
      H(2,5) = omega*normal(1)
c
      H(3,1) = (omega*qir+fprho)*normal(2) - Ug(3)*VdotN 
      H(3,2) = -omega*normal(2)*Ug(2) + Ug(3)*normal(1)
      H(3,3) = (1.d0-omega)*normal(2)*Ug(3) + VdotN - vitno 
      H(3,4) = -omega*normal(2)*Ug(4) + Ug(3)*normal(3)
      H(3,5) = omega*normal(2)
c
      H(4,1) = (omega*qir+fprho)*normal(3) - Ug(4)*VdotN 
      H(4,2) = -omega*normal(3)*Ug(2) + Ug(4)*normal(1)
      H(4,3) = -omega*normal(3)*Ug(3) + Ug(4)*normal(2)
      H(4,4) = (1.d0-omega)*normal(3)*Ug(4) + VdotN - vitno 
      H(4,5) = omega*normal(3) 
c
      H(5,1) = (omega*qir-enth)*VdotN
      H(5,2) = enth*normal(1) - omega*Ug(2)*VdotN
      H(5,3) = enth*normal(2) - omega*Ug(3)*VdotN
      H(5,4) = enth*normal(3) - omega*Ug(4)*VdotN
      H(5,5) = (1.d0+omega)*VdotN - vitno 

c
c Computation of the Roe-averaged state
c
      squsr1 = DSQRT(Ug(1))
      squsr2 = DSQRT(Ud(1))
c
      energ = (Ug(5)-frhog)*ooomega +
     &     0.5d0*Ug(1)*(Ug(2)*Ug(2) + Ug(3)*Ug(3) + Ug(4)*Ug(4))
      enerd = (Ud(5)-frhod)*ooomega +
     &     0.5d0*Ud(1)*(Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4))
c    
      usro   = 1.d0/(squsr1 + squsr2)
c
      uar2   = (squsr1*Ug(2) + squsr2*Ud(2))*usro
c
      uar3   = (squsr1*Ug(3) + squsr2*Ud(3))*usro
c
      uar4   = (squsr1*Ug(4) + squsr2*Ud(4))*usro
c
      uar5   = ((energ + Ug(5))/squsr1 +
     &          (enerd + Ud(5))/squsr2)*usro
c
      squsr1 = 0.5d0*usro/squsr1
c
      enth    = uar5 
c
      if (ABS(Ug(1) - Ud(1)) .lt. 1.0e-14*ABS(Ud(1))) then
        fpuar1 = A1*(-omega/R1r+(1.0-omega*Ug(1)/R1r)*R1r/Ug(1)**2)
     &             *exp(-R1r/Ug(1))
     &         + A2*(-omega/R2r+(1.0-omega*Ug(1)/R2r)*R2r/Ug(1)**2)
     &             *exp(-R2r/Ug(1))
      else
        fpuar1 = (frhod-frhog)/(Ud(1)-Ug(1))
      endif

c
c Computation of the dissipation term
c if prec = 1 then the dissipation is preconditioned
c else if prec = 0 then it is not preconditioned
c
c Reference: Implicit Upwind Schemes for Lowmach number Compressible Flows
c            By Cecile Viozat (INRIA Publication)
c

      VdotN = normal(1)*uar2 + normal(2)*uar3 + normal(3)*uar4
c      
      tet1     = normal(3)*uar3 - normal(2)*uar4
      tet2     = normal(1)*uar4 - normal(3)*uar2
      tet3     = normal(2)*uar2 - normal(1)*uar3
c
      qir   = 0.5d0*(uar2*uar2 + uar3*uar3 + uar4*uar4)
c
c
      cr2   = omega*(uar5-qir)+fpuar1
      cr    = DSQRT(cr2)
      cr2   = 1.d0/cr2

      qirmfpuar1 = qir - fpuar1*ooomega
      omegacr2   = omega*cr2
      usc2 = 0.5d0*cr2

      if (prec .eq. 0) then
        beta = 1.d0
      else
        shock = DABS(Ug(5) - Ud(5))/(Ug(5)+Ud(5))/length
        locMach = DSQRT(2.0d0*qir*cr2)
        beta = MAX(k1*locMach, mach)
        beta = (1.0d0+DSQRT(irey))*beta+shockreducer*shock
        beta = MIN(beta, cmach)
      end if
  
      beta2 = beta * beta


      vp1 = VdotN
      vp4 = 0.5d0*((1.d0+beta2)*VdotN +
     &      DSQRT(((1.d0-beta2)*VdotN)**2 +
     &      4.d0*beta2*cr**2))
      vp5 = 0.5d0*((1.d0+beta2)*VdotN -
     &      DSQRT(((1.d0-beta2)*VdotN)**2 +
     &      4.d0*beta2*cr**2))

c Roe-Turkel

      r = vp4 - vp1*beta2
      s = vp5 - vp1*beta2
      t = 0.5d0*(vp5-vp4)

c
c Dynamic mesh and gamma addition 
c
c Might be required to remove gamma while doing lowmach
c ns or les simulations else it might lead to lot of
c suprious oscillations
c

      vp1 = DABS(vp1-vitno)*gamma
      vp4 = DABS(vp4-vitno)*gamma
      vp5 = DABS(vp5-vitno)*gamma


      temp1 = uar5 - 2.0d0*qir

      tmw(1,1) = normal(1)*omegacr2*temp1-tet1
      tmw(1,2) = normal(1)*omegacr2*uar2
      tmw(1,3) = normal(3)  + normal(1)*omegacr2*uar3
      tmw(1,4) = -normal(2) + normal(1)*omegacr2*uar4
      tmw(1,5) = -normal(1)*omegacr2

      tmw(2,1) = normal(2)*omegacr2*temp1-tet2
      tmw(2,2) = -normal(3) + normal(2)*omegacr2*uar2
      tmw(2,3) = normal(2)*omegacr2*uar3
      tmw(2,4) = normal(1)  + normal(2)*omegacr2*uar4
      tmw(2,5) = -normal(2)*omegacr2

      tmw(3,1) = normal(3)*omegacr2*temp1-tet3
      tmw(3,2) = normal(2)  + normal(3)*omegacr2*uar2
      tmw(3,3) = -normal(1) + normal(3)*omegacr2*uar3
      tmw(3,4) = normal(3)*omegacr2*uar4
      tmw(3,5) = -normal(3)*omegacr2

      temp1 = omega*s/beta2
       
      tmw(4,1) = (cr**2*VdotN + s/beta2*(omega*qir+fpuar1))/t
      tmw(4,2) = -(cr**2*normal(1) + temp1*uar2)/t
      tmw(4,3) = -(cr**2*normal(2) + temp1*uar3)/t
      tmw(4,4) = -(cr**2*normal(3) + temp1*uar4)/t
      tmw(4,5) = temp1/t

      temp1 = omega*r/beta2
 
      tmw(5,1) = -(cr**2*VdotN + r/beta2*(omega*qir+fpuar1))/t
      tmw(5,2) = (cr**2*normal(1) + temp1*uar2)/t
      tmw(5,3) = (cr**2*normal(2) + temp1*uar3)/t
      tmw(5,4) = (cr**2*normal(3) + temp1*uar4)/t
      tmw(5,5) = -temp1/t

c
      DO 70 j=1,5
c
        dmat(1,j)   = vp1*(normal(1)*tmw(1,j) +
     &                     normal(2)*tmw(2,j) +
     &                     normal(3)*tmw(3,j))
     &              + vp4*usc2*tmw(4,j)
     &              + vp5*usc2*tmw(5,j)

        dmat(2,j) = vp1*(normal(1)*uar2*tmw(1,j) +
     &                  (normal(2)*uar2 - normal(3))*tmw(2,j) +
     &                  (normal(3)*uar2 + normal(2))*tmw(3,j))
     &            + vp4*usc2*(uar2 + r*normal(1))*tmw(4,j)
     &            + vp5*usc2*(uar2 + s*normal(1))*tmw(5,j)

        dmat(3,j) = vp1*(normal(2)*uar3*tmw(2,j)   +
     &                  (normal(1)*uar3 + normal(3))*tmw(1,j) +
     &                  (normal(3)*uar3 - normal(1))*tmw(3,j))
     &            + vp4*usc2*(uar3 + r*normal(2))*tmw(4,j)
     &            + vp5*usc2*(uar3 + s*normal(2))*tmw(5,j)

        dmat(4,j) = vp1*(normal(3)*uar4*tmw(3,j) +
     &                  (normal(1)*uar4 - normal(2))*tmw(1,j) +
     &                  (normal(2)*uar4 + normal(1))*tmw(2,j))
     &            + vp4*usc2*(uar4 + r*normal(3))*tmw(4,j)
     &            + vp5*usc2*(uar4 + s*normal(3))*tmw(5,j)

        dmat(5,j) = vp1*((normal(1)*qirmfpuar1 + tet1)*tmw(1,j) +
     &                   (normal(2)*qirmfpuar1 + tet2)*tmw(2,j) +
     &                   (normal(3)*qirmfpuar1 + tet3)*tmw(3,j))
     &            + vp4*usc2*(uar5 + VdotN*r)*tmw(4,j)
     &            + vp5*usc2*(uar5 + VdotN*s)*tmw(5,j)

70       CONTINUE


c
c Final operation
c
      if (sw.eq.1) then
      do i=1,5
c
         H(1,i) = H(1,i)+dmat(1,i) 
c
         H(2,i) = H(2,i)+dmat(2,i)
c
         H(3,i) = H(3,i)+dmat(3,i)

         H(4,i) = H(4,i)+dmat(4,i)
c
         H(5,i) = H(5,i)+dmat(5,i)
c
      end do
      end if

      if (sw.eq.2) then
      do i=1,5
c
         H(1,i) = H(1,i)-dmat(1,i)
c
         H(2,i) = H(2,i)-dmat(2,i)
c
         H(3,i) = H(3,i)-dmat(3,i)
c
         H(4,i) = H(4,i)-dmat(4,i)

         H(5,i) = H(5,i)-dmat(5,i)
c
      end do
      end if

c
c "Normalization" of the flux jacobian
c
      dim = 5 + type
      do i = 1,5
         do k = 1,5
            jac(dim*(i-1) + k) = 0.5d0*H(i,k)*rnorm
         enddo
         do k = 6,dim
            jac(dim*(i-1) + k) = 0.d0
         enddo
      enddo
      do i = 6,dim
         do k = 1,dim
            jac(dim*(i-1) + k) = 0.d0
         enddo
      enddo
c
      end
