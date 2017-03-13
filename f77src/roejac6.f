      SUBROUTINE ROEJAC6(type,gamma,gam,pstiff,enormal,evitno,Ug,Ud,
     &                    jac,sw,mach,k1,cmach,shockreducer,
     &                    irey,length,prec)
c-----------------------------------------------------------------------
c This routine computes the jacobian of the flux of Roe (defect correction) 
c with respect to the conservative variables taken at the vectors Ug, Ud.
c normal is the normal of the boundary concerned by the flux.
c-----------------------------------------------------------------------
c    INCLUDE 'Param3D.h'
      implicit none
c-----------------------------------------------------------------------
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), H(5,5), jac(*)
      REAL*8 tmw(5,5),dmat(5,5)
      REAL*8 VdotN , rnorm, evitno, vitno, invnorm
      REAL*8 g_VdotN(5)
      REAL*8 flur1 , flur2 , flur3 , flur4 , flur5
      REAL*8 dif1 , dif2 , dif3 , dif4 ,dif5
      REAL*8 cr , cr2 , qir
      REAL*8 vp1 , vp4 , vp5 , svp1 , svp4 , svp5
      REAL*8 ener1 , ener2
      REAL*8 uar2 , uar3 , uar4 , uar5
      REAL*8 usro , squsr1 , squsr2
      REAL*8 tet1 , tet2 , tet3
      REAL*8 enth,gc,usc2
      INTEGER k,i,type,dim,j
      REAL*8 gam, gam1, gamma, pstiff
      REAL*8 rRT,sRT,tRT,beta, betaRT2
      REAL*8 mach,k1,locMach,cmach,irey,shock,shockreducer,length
      INTEGER prec,sw

c
c Initialisation
c
      gam1 = gam - 1.d0

      rnorm = DSQRT(enormal(1)**2 + enormal(2)**2 + enormal(3)**2)
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
      qir = 0.5d0*(Ug(2)**2 + Ug(3)**2 + Ug(4)**2)
      enth = gam*(Ug(5)+pstiff)+gam1*Ug(1)*qir
      enth = enth/(gam1*Ug(1))
       
c
      H(1,1) = - vitno
      H(1,2) = normal(1)
      H(1,3) = normal(2)
      H(1,4) = normal(3)
      H(1,5) = 0.0d0 
c
      H(2,1) = gam1*normal(1)*qir - Ug(2)*VdotN 
      H(2,2) = (1.d0-gam1)*normal(1)*Ug(2) + VdotN - vitno 
      H(2,3) = -gam1*normal(1)*Ug(3) + Ug(2)*normal(2)
      H(2,4) = -gam1*normal(1)*Ug(4) + Ug(2)*normal(3) 
      H(2,5) = gam1*normal(1)
c
      H(3,1) = gam1*normal(2)*qir-Ug(3)*VdotN
      H(3,2) = -gam1*normal(2)*Ug(2) + Ug(3)*normal(1)
      H(3,3) = (1.d0-gam1)*normal(2)*Ug(3) + VdotN - vitno 
      H(3,4) = -gam1*normal(2)*Ug(4) + Ug(3)*normal(3)
      H(3,5) = gam1*normal(2)
c
      H(4,1) = gam1*normal(3)*qir-Ug(4)*VdotN
      H(4,2) = -gam1*normal(3)*Ug(2) + Ug(4)*normal(1)
      H(4,3) = -gam1*normal(3)*Ug(3) + Ug(4)*normal(2)
      H(4,4) = (1.d0-gam1)*normal(3)*Ug(4) + VdotN - vitno 
      H(4,5) = gam1*normal(3) 
c
      H(5,1) = (gam1*qir-enth)*VdotN
      H(5,2) = enth*normal(1) - gam1*Ug(2)*VdotN
      H(5,3) = enth*normal(2) - gam1*Ug(3)*VdotN
      H(5,4) = enth*normal(3) - gam1*Ug(4)*VdotN
      H(5,5) = (1.d0+gam1)*VdotN - vitno 

c
c Computation of the Roe-averaged state
c
      squsr1 = DSQRT(Ug(1))
      squsr2 = DSQRT(Ud(1))
c
      ener1  = (Ug(5)+gam*pstiff)/gam1 + 0.5*Ug(1)*(Ug(2)**2 +
     &         Ug(3)**2 + Ug(4)**2)
c
      ener2  = (Ud(5)+gam*pstiff)/gam1 + 0.5*Ud(1)*(Ud(2)**2 +
     &         Ud(3)**2 + Ud(4)**2)
c    
      usro   = 1.d0/(squsr1 + squsr2)
c
      uar2   = (squsr1*Ug(2) + squsr2*Ud(2))*usro
c
      uar3   = (squsr1*Ug(3) + squsr2*Ud(3))*usro
c
      uar4   = (squsr1*Ug(4) + squsr2*Ud(4))*usro
c
      uar5   = ((ener1 + Ug(5))/squsr1 +
     &          (ener2 + Ud(5))/squsr2)*usro
c
      squsr1 = 0.5d0*usro/squsr1

c
      enth    = uar5 
c

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
      qir   = 0.5d0*(uar2*uar2 + uar3*uar3 + uar4*uar4)
c
c
      cr2   = gam1*(uar5 - qir)
      cr    = DSQRT(cr2)
      cr2   = 1.d0/cr2

      gc = gam1/(cr*cr)
      usc2 = 0.5d0/(cr*cr) 


      if (prec .eq. 0) then
        beta = 1.d0
      else
        shock = DABS(Ug(5) - Ud(5))/(Ug(5)+Ud(5))/length
        locMach = DSQRT(2.0d0*qir*cr2)
        beta = MAX(k1*locMach, mach)
        beta = (1.0d0+DSQRT(irey))*beta+shockreducer*shock
        beta = MIN(beta, cmach)
      end if
  
      betaRT2 = beta * beta


      vp1 = VdotN
      vp4 = 0.5d0*((1.d0+betaRT2)*VdotN +
     &      DSQRT(((1.d0-betaRT2)*VdotN)**2 +
     &      4.d0*betaRT2*cr**2))
      vp5 = 0.5d0*((1.d0+betaRT2)*VdotN -
     &      DSQRT(((1.d0-betaRT2)*VdotN)**2 +
     &      4.d0*betaRT2*cr**2))

c Roe-Turkel

      rRT = vp4 - vp1*betaRT2
      sRT = vp5 - vp1*betaRT2
      tRT = 0.5d0*(vp5-vp4)

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


c
         tmw(1,1)   = normal(1) + normal(2)*uar4 - uar3*normal(3) 
     &   - gc*normal(1)*qir
         tmw(1,2)   = gc*normal(1)*uar2
         tmw(1,3)   = gc*normal(1)*uar3 + normal(3)
         tmw(1,4)   = gc*normal(1)*uar4 - normal(2)
         tmw(1,5)   =-gc*normal(1)
c
         tmw(2,1)   = normal(2) + normal(3)*uar2 - uar4*normal(1) 
     &   - gc*normal(2)*qir
         tmw(2,2)   = gc*normal(2)*uar2 - normal(3)
         tmw(2,3)   = gc*normal(2)*uar3
         tmw(2,4)   = gc*normal(2)*uar4 + normal(1)
         tmw(2,5)   =-gc*normal(2)
c
         tmw(3,1)   = normal(3) + normal(1)*uar3 - uar2*normal(2)
     &   - gc*normal(3)*qir
         tmw(3,2)   = gc*normal(3)*uar2 + normal(2)
         tmw(3,3)   = gc*normal(3)*uar3 - normal(1)
         tmw(3,4)   = gc*normal(3)*uar4
         tmw(3,5)   =-gc*normal(3)
c     
c Roe-Turkel
c
         tmw(4,1) = gam1*qir*sRT/(betaRT2*tRT) + VdotN*cr**2/tRT
         tmw(4,2) =-gam1*uar2*sRT/(betaRT2*tRT) - cr**2*normal(1)/tRT
         tmw(4,3) =-gam1*uar3*sRT/(betaRT2*tRT) - cr**2*normal(2)/tRT
         tmw(4,4) =-gam1*uar4*sRT/(betaRT2*tRT) - cr**2*normal(3)/tRT
         tmw(4,5) = gam1*sRT/(betaRT2*tRT)
c
c
c Roe-Turkel
c
         tmw(5,1) = -gam1*qir*rRT/(betaRT2*tRT) - VdotN*cr**2/tRT
         tmw(5,2) = gam1*uar2*rRT/(betaRT2*tRT) + cr**2*normal(1)/tRT
         tmw(5,3) = gam1*uar3*rRT/(betaRT2*tRT) + cr**2*normal(2)/tRT
         tmw(5,4) = gam1*uar4*rRT/(betaRT2*tRT) + cr**2*normal(3)/tRT
         tmw(5,5) = -gam1*rRT/(betaRT2*tRT)
c

c
         DO 70 j=1,5
c
            dmat(1,j)              = vp1*(normal(1)*tmw(1,j) +
     &                                    normal(2)*tmw(2,j) +
     &                                    normal(3)*tmw(3,j))
     &                             + vp4*usc2*tmw(4,j)
     &                             + vp5*usc2*tmw(5,j)
c
c Roe-Turkel
c
            dmat(2,j) = vp1*(normal(1)*uar2*tmw(1,j) +
     &                 (normal(2)*uar2 - normal(3))*tmw(2,j) +
     &                 (normal(3)*uar2 + normal(2))*tmw(3,j))
     &                 + vp4*(usc2*(uar2 + rRT*normal(1)))*tmw(4,j)
     &                 + vp5*(usc2*(uar2 + sRT*normal(1)))*tmw(5,j)
c
c Roe-Turkel
c
            dmat(3,j) = vp1*(normal(2)*uar3*tmw(2,j)   +
     &                 (normal(1)*uar3 + normal(3))*tmw(1,j) +
     &                 (normal(3)*uar3 - normal(1))*tmw(3,j))
     &                 + vp4*(usc2*(uar3 + rRT*normal(2)))*tmw(4,j)
     &                 + vp5*(usc2*(uar3 + sRT*normal(2)))*tmw(5,j)

c
c Roe-turkel
c
            dmat(4,j) = vp1*(normal(3)*uar4*tmw(3,j) +
     &                 (normal(1)*uar4 - normal(2))*tmw(1,j) +
     &                 (normal(2)*uar4 + normal(1))*tmw(2,j))
     &                 + vp4*(usc2*(uar4 + rRT*normal(3)))*tmw(4,j)
     &                 + vp5*(usc2*(uar4 + sRT*normal(3)))*tmw(5,j)
c
c Roe-Turkel
c
            dmat(5,j) = vp1*((normal(1)*qir +
     &                 (normal(3)*uar3 - normal(2)*uar4))*tmw(1,j) +
     &                 (normal(2)*qir +
     &                 (normal(1)*uar4 - normal(3)*uar2))*tmw(2,j) +
     &                 (normal(3)*qir +
     &                 (normal(2)*uar2 - normal(1)*uar3))*tmw(3,j))
     &                 + vp4*(usc2*(enth + VdotN*rRT))*tmw(4,j)
     &                 + vp5*(usc2*(enth + VdotN*sRT))*tmw(5,j)
c
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
