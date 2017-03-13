      SUBROUTINE ROEFLUX6(type,gamma,gam,pstiff,enormal,evitno,
     &     Ugr,Ug,Udr,Ud,phi)
c-----------------------------------------------------------------------
c This routine computes the Flux of Roe taken at the vectors Ug, Ud
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), evitno, phi(*)
      REAL*8 Ugr(*), Udr(*), energ, enerd
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
      INTEGER type
c
c Initialisation
c
      gam1 = gam - 1.d0      
c
      rnorm = SQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) + 
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
      H = gam*(Ug(5)+pstiff) + 0.5*gam1*Ug(1)*vitg2
      H = H/(gam1*Ug(1))
      phi(1) = Ug(1)*(VdotN - vitno)
      phi(2) = phi(1)*Ug(2) + Ug(5)*normal(1)
      phi(3) = phi(1)*Ug(3) + Ug(5)*normal(2)
      phi(4) = phi(1)*Ug(4) + Ug(5)*normal(3)
      phi(5) = phi(1)*H + Ug(5)*vitno
c
      VdotN = Ud(2)*normal(1)+Ud(3)*normal(2)+Ud(4)*normal(3)
      VdotN = VdotN - vitno
      vitd2 = Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4)
      H = gam*(Ud(5)+pstiff)+0.5*gam1*Ud(1)*vitd2
      H = H/(gam1*Ud(1))
      phi(1) = phi(1) + Ud(1)*VdotN
      phi(2) = phi(2) + Ud(1)*Ud(2)*VdotN +Ud(5)*normal(1)
      phi(3) = phi(3) + Ud(1)*Ud(3)*VdotN +Ud(5)*normal(2)
      phi(4) = phi(4) + Ud(1)*Ud(4)*VdotN +Ud(5)*normal(3)
      phi(5) = phi(5) + Ud(1)*VdotN*H +Ud(5)*vitno
c
c Computation of the Roe-averaged state
c
      squsr1                    = SQRT(Ug(1))
      squsr2                    = SQRT(Ud(1))
c     
      ener1                     = (Ug(5)+gam*pstiff)/gam1 +
     &                            0.5*Ug(1)*vitg2
c
      ener2                     = (Ud(5)+gam*pstiff)/gam1 +
     &                            0.5*Ud(1)*vitd2
c
      usro                      = 1.0/(squsr1 + squsr2)
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
c
      VdotN                     = normal(1)*uar2 + normal(2)*uar3 +
     &                               normal(3)*uar4
c
      qir                       = 0.5*(uar2*uar2 + uar3*uar3 +
     &                                    uar4*uar4)
c
      tet1                      = normal(3)*uar3 - normal(2)*uar4
      tet2                      = normal(1)*uar4 - normal(3)*uar2
      tet3                      = normal(2)*uar2 - normal(1)*uar3
c
      cr2                       = gam1*(uar5 - qir)
      cr                        = SQRT(cr2)
      cr2                       = 1.0/cr2
c
      energ = (Ugr(5)+gam*pstiff)/gam1 + 
     &     0.5*Ugr(1)*(Ugr(2)*Ugr(2) + Ugr(3)*Ugr(3) + Ugr(4)*Ugr(4))
      enerd = (Udr(5)+gam*pstiff)/gam1 + 
     &     0.5*Udr(1)*(Udr(2)*Udr(2) + Udr(3)*Udr(3) + Udr(4)*Udr(4))
      dif1                     = - Ugr(1) + Udr(1)
      dif2                     = - Ugr(1)*Ugr(2) + Udr(1)*Udr(2)
      dif3                     = - Ugr(1)*Ugr(3) + Udr(1)*Udr(3)
      dif4                     = - Ugr(1)*Ugr(4) + Udr(1)*Udr(4)
      dif5                     = - energ + enerd
c
      vp1                       = gamma * ABS(VdotN - vitno)
      vp4                       = gamma * ABS((VdotN + cr) - vitno)
      vp5                       = gamma * ABS((VdotN - cr) - vitno)
c
      flur1                     = vp1*
     &              ((normal(1)*(1.0 - gam1*qir*cr2) - tet1)*dif1 +
     &               (normal(1)*gam1*uar2*cr2)*dif2  +
     &               (normal(3)  + (normal(1)*gam1*uar3*cr2))*dif3   +
     &               (-normal(2) + (normal(1)*gam1*uar4*cr2))*dif4   -
     &               (normal(1)*gam1*cr2)*dif5)
c
      flur2                     = vp1*
     &              ((normal(2)*(1.0 - gam1*qir*cr2) - tet2)*dif1 +
     &               (-normal(3) + (normal(2)*gam1*uar2*cr2))*dif2   +
     &               (normal(2)*gam1*uar3*cr2)*dif3  +
     &               (normal(1)  + (normal(2)*gam1*uar4*cr2))*dif4   -
     &               (normal(2)*gam1*cr2)*dif5)
c
      flur3                     = vp1*
     &              ((normal(3)*(1.0 - gam1*qir*cr2) - tet3)*dif1 +
     &               (normal(2)  + (normal(3)*gam1*uar2*cr2))*dif2   +
     &               (-normal(1) + (normal(3)*gam1*uar3*cr2))*dif3   +
     &               (normal(3)*gam1*uar4*cr2)*dif4  -
     &               (normal(3)*gam1*cr2)*dif5)
c
      flur4                     = vp4*
     &              ((-cr*VdotN   + gam1*qir)*dif1  +
     &               ( cr*normal(1) - gam1*uar2)*dif2 +
     &               ( cr*normal(2) - gam1*uar3)*dif3 +
     &               ( cr*normal(3) - gam1*uar4)*dif4 +
     &               gam1*dif5)
c
      flur5                     = vp5*
     &              (( cr*VdotN  + gam1* qir)*dif1 +
     &               (-cr*normal(1) - gam1*uar2)*dif2 +
     &               (-cr*normal(2) - gam1*uar3)*dif3 +
     &               (-cr*normal(3) - gam1*uar4)*dif4 +
     &               gam1*dif5)
c
      phi(1) =  phi(1) - (normal(1)*flur1 +
     &     normal(2)*flur2 + normal(3)*flur3 +
     &     0.5*(flur4 + flur5)*cr2)
c
      phi(2) = phi(2) - (  
     &              (uar2*normal(1))*flur1 +
     &              (uar2*normal(2) - normal(3))*flur2 +
     &              (uar2*normal(3) + normal(2))*flur3 +
     &              0.5*normal(1)*(flur4 - flur5)/cr +
     &              0.5*uar2*(flur4 + flur5)*cr2)
c
      phi(3) = phi(3) - (
     &              (uar3*normal(1) + normal(3))*flur1 +
     &              (uar3*normal(2))*flur2 +
     &              (uar3*normal(3) - normal(1))*flur3 +
     &              0.5*normal(2)*(flur4 - flur5)/cr +
     &              0.5*uar3*(flur4 + flur5)*cr2)
c
      phi(4) = phi(4) - ( 
     &              (uar4*normal(1) - normal(2))*flur1 +
     &              (uar4*normal(2) + normal(1))*flur2 +
     &              (uar4*normal(3))*flur3 +
     &              0.5*normal(3)*(flur4 - flur5)/cr +
     &              0.5*uar4*(flur4 + flur5)*cr2)
c
      phi(5) = phi(5) - (
     &              (qir*normal(1) + tet1)*flur1 +
     &              (qir*normal(2) + tet2)*flur2 +
     &              (qir*normal(3) + tet3)*flur3 +
     &              0.5*VdotN*(flur4 - flur5)/cr +
     &              0.5*uar5*(flur4 + flur5)*cr2)
c
      phi(1) = phi(1)*0.5*rnorm
      phi(2) = phi(2)*0.5*rnorm
      phi(3) = phi(3)*0.5*rnorm
      phi(4) = phi(4)*0.5*rnorm
      phi(5) = phi(5)*0.5*rnorm

      if (type.eq.1) then
         updir = 0.5d0 + dsign(0.5d0, phi(1))
         phi(6) = phi(1) * (updir * Ugr(6) + (1.0d0 - updir) * Udr(6))
      else if (type.eq.2) then
         updir = 0.5d0 + dsign(0.5d0, phi(1))
         phi(6) = phi(1) * (updir * Ugr(6) + (1.0d0 - updir) * Udr(6))
         phi(7) = phi(1) * (updir * Ugr(7) + (1.0d0 - updir) * Udr(7))
      endif
c
      END
