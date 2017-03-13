      SUBROUTINE ROEFLUX7(type,gam,pstiff,enormal,evitno,
     &     Ug,Ud,phi)
c-----------------------------------------------------------------------
c This routine computes the Flux of Roe taken at the vectors Ug, Ud
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
c gam is the ratio of cp/cv
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), evitno, phi(*)
      REAL*8 H, vitno, updir
      REAL*8 VdotN, rnorm, invnorm, VdotNt
      REAL*8 flur4, flur5
      REAL*8 dif1, dif2, dif3, dif4, dif5
      REAL*8 cr, cr2, cr2byt, invtbeta2, qir
      REAL*8 vp1, vp4, vp5
      REAL*8 ener1, ener2
      REAL*8 uar2, uar3, uar4, uar5
      REAL*8 usro, squsr1, squsr2
      REAL*8 gam, gam1, invgam1, vitg2, vitd2, pstiff
      REAL*8 r, s, t, beta
      INTEGER type
c
c Initialisation
c
      gam1 = gam - 1.d0      
      invgam1 = 1.d0/gam1
c
      rnorm = DSQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) + 
     &             enormal(3)*enormal(3))
      
      invnorm = 1.0d0 / rnorm
c
      normal(1) = enormal(1) * invnorm
      normal(2) = enormal(2) * invnorm
      normal(3) = enormal(3) * invnorm

      vitno = evitno * invnorm

      vitg2 = Ug(2)*Ug(2) + Ug(3)*Ug(3) + Ug(4)*Ug(4)
      vitd2 = Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4)

c
c Computation of the Roe-averaged state
c

      squsr1                    = DSQRT(Ug(1))
      squsr2                    = DSQRT(Ud(1))
c     
      ener1                     = (Ug(5)+gam*pstiff)*invgam1 +
     &                            0.5d0*Ug(1)*vitg2
c
      ener2                     = (Ud(5)+gam*pstiff)*invgam1 +
     &                            0.5d0*Ud(1)*vitd2
c
      usro                      = 1.d0/(squsr1 + squsr2)
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

      VdotNt                     = normal(1)*uar2 + normal(2)*uar3 +
     &                               normal(3)*uar4

c     
c
c Computation of the centred terms
c
      if (VdotNt- vitno .gt. 0) then
         VdotN = Ug(2)*normal(1) + Ug(3)*normal(2) + Ug(4)*normal(3)
         H = gam*(Ug(5)+pstiff) + 0.5d0*gam1*Ug(1)*vitg2
         H = H/(gam1*Ug(1))
         phi(1) = Ug(1)*(VdotN - vitno)
         phi(2) = phi(1)*Ug(2) + Ug(5)*normal(1)
         phi(3) = phi(1)*Ug(3) + Ug(5)*normal(2)
         phi(4) = phi(1)*Ug(4) + Ug(5)*normal(3)
         phi(5) = phi(1)*H + Ug(5)*vitno
      else
c
         VdotN = Ud(2)*normal(1)+Ud(3)*normal(2)+Ud(4)*normal(3)
         H = gam*(Ud(5)+pstiff)+0.5d0*gam1*Ud(1)*vitd2
         H = H/(gam1*Ud(1))
         phi(1) = Ud(1)*(VdotN - vitno)
         phi(2) = phi(1)*Ud(2) + Ud(5)*normal(1)
         phi(3) = phi(1)*Ud(3) + Ud(5)*normal(2)
         phi(4) = phi(1)*Ud(4) + Ud(5)*normal(3)
         phi(5) = phi(1)*H + Ud(5)*vitno
      end if

c
c Computation of the dissipation term 
c
c Reference: Implicit Upwind Schemes for Low Mach Number Compressible Flows
c            By Cecile Viozat (INRIA Publication)


      VdotN                     = VdotNt
c
      qir                       = 0.5d0*(uar2*uar2 + uar3*uar3 +
     &                                    uar4*uar4)
c
      cr2                       = gam1*(uar5 - qir)
      cr                        = DSQRT(cr2)
      cr2                       = 1.d0/cr2
c
        
      beta = 1.d0

      vp1 = VdotN
      vp4 = VdotN + cr
      vp5 = VdotN - cr

      r = cr
      s = -cr
      t = -cr 
      cr2byt = -cr
      invtbeta2 = -1.d0/cr

      dif1                     = - Ug(1) + Ud(1)
      dif2                     = - Ug(1)*Ug(2) + Ud(1)*Ud(2)
      dif3                     = - Ug(1)*Ug(3) + Ud(1)*Ud(3)
      dif4                     = - Ug(1)*Ug(4) + Ud(1)*Ud(4)
      dif5                     = - ener1 + ener2

c Dynamic mesh inclusion 

      vp1 = (vp1-vitno)
      vp4 = (vp4-vitno)
      vp5 = (vp5-vitno)


      if (VdotNt- vitno .le. 0) then
         flur4                     = max(vp4,0.)* 
     &        ( cr2byt*( VdotN*dif1 - normal(1)*dif2 - normal(2)*dif3
     &                              - normal(3)*dif4 ) + 
     &          gam1*s*invtbeta2*( qir*dif1 - uar2*dif2 - uar3*dif3
     &                                     - uar4*dif4 + dif5 ) )
         flur5 = 0.d0
      else
         flur4 = 0.d0
         flur5                     = -min(vp5,0.)*  
     &        ( -cr2byt*( VdotN*dif1 - normal(1)*dif2 - normal(2)*dif3
     &                              - normal(3)*dif4 ) + 
     &          (-gam1*r*invtbeta2*( qir*dif1 - uar2*dif2 - uar3*dif3
     &                                     - uar4*dif4 + dif5 ) ) )
      end if

c
c Final phi
c

      phi(1) =  phi(1) - 0.5d0*(flur4 + flur5)*cr2

c

      phi(2) = phi(2) - (normal(1)*(r*flur4 + s*flur5) +
     &               uar2*(flur4 + flur5))*0.5d0*cr2


c
      phi(3) = phi(3) - (normal(2)*(r*flur4 + s*flur5) +
     &               uar3*(flur4 + flur5))*0.5d0*cr2


c
      phi(4) = phi(4) - (normal(3)*(r*flur4 + s*flur5) +
     &                    uar4*(flur4 + flur5))*0.5d0*cr2

c
      phi(5) = phi(5) - (VdotN*(r*flur4 + s*flur5) +
     &                    uar5*(flur4 + flur5))*0.5d0*cr2

      phi(1) = phi(1)*rnorm
      phi(2) = phi(2)*rnorm
      phi(3) = phi(3)*rnorm
      phi(4) = phi(4)*rnorm
      phi(5) = phi(5)*rnorm

c
c For one and two equation turbulence models
c

      if (type.eq.1) then
         updir = 0.5d0 + dsign(0.5d0, phi(1))
         phi(6) = phi(1) * (updir * Ug(6) + (1.0d0 - updir) * Ud(6))
      else if (type.eq.2) then
         updir = 0.5d0 + dsign(0.5d0, phi(1))
         phi(6) = phi(1) * (updir * Ug(6) + (1.0d0 - updir) * Ud(6))
         phi(7) = phi(1) * (updir * Ug(7) + (1.0d0 - updir) * Ud(7))
      endif
      
      END
