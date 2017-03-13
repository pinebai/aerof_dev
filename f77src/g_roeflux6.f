      SUBROUTINE GXROEFLUX6(type,gamma,gam,pstiff,dpstiff,
     &                      enormal,denormal,evitno,devitno,
     &                      Ugr,dUgr,Ug,dUg,Udr,dUdr,Ud,dUd,phi,dphi)
c-----------------------------------------------------------------------
c This routine computes the Roe Flux derivative taken at the vectors 
c Ug, dUg, Ud and dUd. normal is the normal of the boundary concerned 
c by the flux. phi stores the resulting flux.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), evitno, phi(*)
      REAL*8 dUg(*), dUd(*), dnormal(3), denormal(3), devitno, dphi(*)
      REAL*8 Ugr(*), Udr(*), energ, enerd
      REAL*8 dUgr(*), dUdr(*), denerg, denerd
      REAL*8 H, vitno, updir, gamma 
      REAL*8 dH, dvitno 
      REAL*8 VdotN, rnorm, invnorm
      REAL*8 dVdotN, drnorm, dinvnorm
      REAL*8 flur1, flur2, flur3, flur4, flur5
      REAL*8 dflur1, dflur2, dflur3, dflur4, dflur5
      REAL*8 dif1, dif2, dif3, dif4, dif5
      REAL*8 ddif1, ddif2, ddif3, ddif4, ddif5
      REAL*8 cr, cr2, qir
      REAL*8 dcr, dcr2, dqir
      REAL*8 vp1, vp4, vp5
      REAL*8 dvp1, dvp4, dvp5
      REAL*8 ener1, ener2
      REAL*8 dener1, dener2
      REAL*8 uar1, uar2, uar3, uar4, uar5
      REAL*8 duar1, duar2, duar3, duar4, duar5
      REAL*8 usro, squsr1, squsr2
      REAL*8 dusro, dsqusr1, dsqusr2
      REAL*8 tet1, tet2, tet3
      REAL*8 dtet1, dtet2, dtet3
      REAL*8 gam, gam1, vitg2, vitd2, pstiff
      REAL*8 dvitg2, dvitd2, dpstiff
      INTEGER type
c
c Initialisation
c
      gam1 = gam - 1.0d0      
c
      rnorm = SQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) + 
     &        enormal(3)*enormal(3))
      drnorm = 1.0d0 / SQRT(enormal(1)*enormal(1) + 
     &         enormal(2)*enormal(2) + enormal(3)*enormal(3)) * 
     &         (enormal(1)*denormal(1) + enormal(2)*denormal(2) + 
     &         enormal(3)*denormal(3))
      invnorm = 1.0d0 / rnorm
      dinvnorm = -1.0d0 / (rnorm*rnorm) * drnorm
c
      normal(1) = enormal(1) * invnorm
      dnormal(1) = denormal(1) * invnorm + enormal(1) * dinvnorm
      normal(2) = enormal(2) * invnorm
      dnormal(2) = denormal(2) * invnorm + enormal(2) * dinvnorm
      normal(3) = enormal(3) * invnorm
      dnormal(3) = denormal(3) * invnorm + enormal(3) * dinvnorm

      vitno = evitno * invnorm
      dvitno = devitno * invnorm + evitno * dinvnorm
c
c Computation of the centred terms
c
      VdotN = Ug(2)*normal(1) + Ug(3)*normal(2) + Ug(4)*normal(3)
      dVdotN = dUg(2)*normal(1) + Ug(2)*dnormal(1) + 
     &         dUg(3)*normal(2) + Ug(3)*dnormal(2) + 
     &         dUg(4)*normal(3) + Ug(4)*dnormal(3)
      vitg2 = Ug(2)*Ug(2) + Ug(3)*Ug(3) + Ug(4)*Ug(4)
      dvitg2 = 2.0d0*(Ug(2)*dUg(2) + Ug(3)*dUg(3) + Ug(4)*dUg(4))
      H = gam*(Ug(5)+pstiff) + 0.5*gam1*Ug(1)*vitg2
      dH = gam*(dUg(5)+dpstiff) + 0.5*gam1*dUg(1)*vitg2 + 
     &     0.5*gam1*Ug(1)*dvitg2
      dH = (dH*(gam1*Ug(1))-H*(gam1*dUg(1)))/((gam1*Ug(1))*(gam1*Ug(1)))
      H = H/(gam1*Ug(1))
      phi(1) = Ug(1)*(VdotN - vitno)
      dphi(1) = dUg(1)*(VdotN - vitno) + Ug(1)*(dVdotN - dvitno)
      phi(2) = phi(1)*Ug(2) + Ug(5)*normal(1)
      dphi(2) = dphi(1)*Ug(2) + phi(1)*dUg(2) + dUg(5)*normal(1) + 
     &          Ug(5)*dnormal(1)
      phi(3) = phi(1)*Ug(3) + Ug(5)*normal(2)
      dphi(3) = dphi(1)*Ug(3) + phi(1)*dUg(3) + dUg(5)*normal(2) + 
     &          Ug(5)*dnormal(2)
      phi(4) = phi(1)*Ug(4) + Ug(5)*normal(3)
      dphi(4) = dphi(1)*Ug(4) + phi(1)*dUg(4) + dUg(5)*normal(3) + 
     &          Ug(5)*dnormal(3)
      phi(5) = phi(1)*H + Ug(5)*vitno
      dphi(5) = dphi(1)*H + phi(1)*dH + dUg(5)*vitno + Ug(5)*dvitno
c
      VdotN = Ud(2)*normal(1)+Ud(3)*normal(2)+Ud(4)*normal(3)
      dVdotN = dUd(2)*normal(1)+Ud(2)*dnormal(1) +
     &         dUd(3)*normal(2)+Ud(3)*dnormal(2) +
     &         dUd(4)*normal(3)+Ud(4)*dnormal(3)
      VdotN = VdotN - vitno
      dVdotN = dVdotN - dvitno
      vitd2 = Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4)
      dvitd2 = 2.0d0*(Ud(2)*dUd(2) + Ud(3)*dUd(3) + Ud(4)*dUd(4))
      H = gam*(Ud(5)+pstiff)+0.5*gam1*Ud(1)*vitd2
      dH = gam*(dUd(5)+dpstiff) + 0.5*gam1*dUd(1)*vitd2 + 
     &     0.5*gam1*Ud(1)*dvitd2
      dH = (dH*(gam1*Ud(1))-H*(gam1*dUd(1)))/((gam1*Ud(1))*(gam1*Ud(1)))
      H = H/(gam1*Ud(1))
      phi(1) = phi(1) + Ud(1)*VdotN
      dphi(1) = dphi(1) + dUd(1)*VdotN + Ud(1)*dVdotN
      phi(2) = phi(2) + Ud(1)*Ud(2)*VdotN +Ud(5)*normal(1)
      dphi(2) = dphi(2) + dUd(1)*Ud(2)*VdotN + Ud(1)*dUd(2)*VdotN + 
     &          Ud(1)*Ud(2)*dVdotN + dUd(5)*normal(1) +Ud(5)*dnormal(1)
      phi(3) = phi(3) + Ud(1)*Ud(3)*VdotN +Ud(5)*normal(2)
      dphi(3) = dphi(3) + dUd(1)*Ud(3)*VdotN + Ud(1)*dUd(3)*VdotN + 
     &          Ud(1)*Ud(3)*dVdotN + dUd(5)*normal(2) + Ud(5)*dnormal(2)
      phi(4) = phi(4) + Ud(1)*Ud(4)*VdotN +Ud(5)*normal(3)
      dphi(4) = dphi(4) + dUd(1)*Ud(4)*VdotN + Ud(1)*dUd(4)*VdotN + 
     &          Ud(1)*Ud(4)*dVdotN + dUd(5)*normal(3) + Ud(5)*dnormal(3)
      phi(5) = phi(5) + Ud(1)*VdotN*H +Ud(5)*vitno
      dphi(5) = dphi(5) + dUd(1)*VdotN*H + Ud(1)*dVdotN*H + 
     &          Ud(1)*VdotN*dH + dUd(5)*vitno + Ud(5)*dvitno
c
c Computation of the Roe-averaged state
c
      squsr1 = SQRT(Ug(1))
      dsqusr1 = 1.0d0/(2.0d0*SQRT(Ug(1))) * dUg(1)
      squsr2 = SQRT(Ud(1))
      dsqusr2 = 1.0d0/(2.0d0*SQRT(Ud(1))) * dUd(1)
c     
      ener1 = (Ug(5)+gam*pstiff)/gam1 + 0.5d0*Ug(1)*vitg2
      dener1 = (dUg(5)+gam*dpstiff)/gam1 + 0.5d0*dUg(1)*vitg2 + 
     &         0.5d0*Ug(1)*dvitg2
c
      ener2 = (Ud(5)+gam*pstiff)/gam1 + 0.5d0*Ud(1)*vitd2
      dener2 = (dUd(5)+gam*dpstiff)/gam1 + 0.5d0*dUd(1)*vitd2 + 
     &         0.5d0*Ud(1)*dvitd2
c
      usro  = 1.0d0/(squsr1 + squsr2)
      dusro = -1.0d0/((squsr1 + squsr2)*(squsr1 + squsr2)) * 
     &        (dsqusr1 + dsqusr2)
c
      uar1 = (squsr1*Ug(1) + squsr2*Ud(1))*usro
      duar1 = (dsqusr1*Ug(1) + squsr1*dUg(1) + dsqusr2*Ud(1) + 
     &        squsr2*dUd(1))*usro + (squsr1*Ug(1) + 
     &        squsr2*Ud(1))*dusro
c
      uar2 = (squsr1*Ug(2) + squsr2*Ud(2))*usro
      duar2 = (dsqusr1*Ug(2) + squsr1*dUg(2) + dsqusr2*Ud(2) +
     &        squsr2*dUd(2))*usro + (squsr1*Ug(2) +
     &        squsr2*Ud(2))*dusro
c
      uar3 = (squsr1*Ug(3) + squsr2*Ud(3))*usro
      duar3 = (dsqusr1*Ug(3) + squsr1*dUg(3) + dsqusr2*Ud(3) +
     &        squsr2*dUd(3))*usro + (squsr1*Ug(3) +
     &        squsr2*Ud(3))*dusro
c
      uar4 = (squsr1*Ug(4) + squsr2*Ud(4))*usro
      duar4 = (dsqusr1*Ug(4) + squsr1*dUg(4) + dsqusr2*Ud(4) + 
     &        squsr2*dUd(4))*usro + (squsr1*Ug(4) +
     &        squsr2*Ud(4))*dusro
c
      uar5 = ((ener1 + Ug(5))/squsr1 + (ener2 + Ud(5))/squsr2)*usro
      duar5 = (((dener1 + dUg(5))*squsr1 - (ener1 + Ug(5))*dsqusr1)/
     &        (squsr1*squsr1) + ((dener2 + dUd(5))*squsr2 -
     &        (ener2 + Ud(5))*dsqusr2)/(squsr2*squsr2))*usro + 
     &        ((ener1 + Ug(5))/squsr1 + (ener2 + Ud(5))/squsr2)*dusro
c
c Computation of the dissipation term
c
      VdotN = normal(1)*uar2 + normal(2)*uar3 + normal(3)*uar4
      dVdotN = dnormal(1)*uar2 + normal(1)*duar2 + dnormal(2)*uar3 + 
     &         normal(2)*duar3 + dnormal(3)*uar4 + normal(3)*duar4
c
      qir = 0.5d0*(uar2*uar2 + uar3*uar3 + uar4*uar4)
      dqir = uar2*duar2 + uar3*duar3 + uar4*duar4
c
      tet1 = normal(3)*uar3 - normal(2)*uar4
      dtet1 = dnormal(3)*uar3 + normal(3)*duar3 - 
     &        dnormal(2)*uar4 - normal(2)*duar4
      tet2 = normal(1)*uar4 - normal(3)*uar2
      dtet2 = dnormal(1)*uar4 + normal(1)*duar4 - 
     &        dnormal(3)*uar2 - normal(3)*duar2
      tet3 = normal(2)*uar2 - normal(1)*uar3
      dtet3 = dnormal(2)*uar2 + normal(2)*duar2 - 
     &        dnormal(1)*uar3 - normal(1)*duar3
c
      cr2 = gam1*(uar5 - qir)
      dcr2 = gam1*(duar5 - dqir)
      cr = SQRT(cr2)
      dcr = 1.0d0/(2.0d0*SQRT(cr2))*dcr2
      dcr2 = -1.0d0/(cr2*cr2)*dcr2
      cr2 = 1.0/cr2
c
      energ = (Ugr(5)+gam*pstiff)/gam1 + 0.5d0*Ugr(1)*(Ugr(2)*Ugr(2) + 
     &        Ugr(3)*Ugr(3) + Ugr(4)*Ugr(4))
      denerg = (dUgr(5)+gam*dpstiff)/gam1 + 0.5d0*dUgr(1)* 
     &         (Ugr(2)*Ugr(2) + Ugr(3)*Ugr(3) + Ugr(4)*Ugr(4)) +
     &         Ugr(1)*(Ugr(2)*dUgr(2) + Ugr(3)*dUgr(3) + Ugr(4)*dUgr(4))
      enerd = (Udr(5)+gam*pstiff)/gam1 + 0.5d0*Udr(1)*(Udr(2)*Udr(2) + 
     &        Udr(3)*Udr(3) + Udr(4)*Udr(4))
      denerd = (dUdr(5)+gam*dpstiff)/gam1 + 0.5d0*dUdr(1)*
     &         (Udr(2)*Udr(2) + Udr(3)*Udr(3) + Udr(4)*Udr(4)) +
     &         Udr(1)*(Udr(2)*dUdr(2) + Udr(3)*dUdr(3) + Udr(4)*dUdr(4))
      dif1 = - Ugr(1) + Udr(1)
      ddif1 = - dUgr(1) + dUdr(1)
      dif2 = - Ugr(1)*Ugr(2) + Udr(1)*Udr(2)
      ddif2 = - dUgr(1)*Ugr(2) - Ugr(1)*dUgr(2) +  
     &        dUdr(1)*Udr(2) + Udr(1)*dUdr(2)
      dif3 = - Ugr(1)*Ugr(3) + Udr(1)*Udr(3)
      ddif3 = - dUgr(1)*Ugr(3) - Ugr(1)*dUgr(3) + 
     &        dUdr(1)*Udr(3) + Udr(1)*dUdr(3)
      dif4 = - Ugr(1)*Ugr(4) + Udr(1)*Udr(4)
      ddif4 = - dUgr(1)*Ugr(4) - Ugr(1)*dUgr(4) + 
     &        dUdr(1)*Udr(4) + Udr(1)*dUdr(4)
      dif5 = - energ + enerd
      ddif5 = - denerg + denerd
c
      vp1 = gamma * ABS(VdotN - vitno)
      if (ABS(VdotN - vitno) .eq. 0.0d0) then
        dvp1 = 0.0d0
      else
	dvp1 = gamma * (VdotN - vitno)/ABS(VdotN - vitno)*
     &        (dVdotN - dvitno)
      end if
      vp4 = gamma * ABS((VdotN + cr) - vitno)
      if (ABS((VdotN + cr) - vitno) .eq. 0.0d0) then
        dvp4 = 0.0d0
      else
	dvp4 = gamma * ((VdotN + cr) - vitno)/ABS((VdotN + cr) - vitno)*
     &        ((dVdotN + dcr) - dvitno)
      end if
      vp5 = gamma * ABS((VdotN - cr) - vitno)
      if (ABS((VdotN - cr) - vitno) .eq. 0.0d0) then
        dvp5 = 0.0d0
      else
	dvp5 = gamma * ((VdotN - cr) - vitno)/ABS((VdotN - cr) - vitno)*
     &        ((dVdotN - dcr) - dvitno)
      end if
c
      flur1 = vp1*
     &        ((normal(1)*(1.0d0 - gam1*qir*cr2) - tet1)*dif1 +
     &        (normal(1)*gam1*uar2*cr2)*dif2 +
     &        (normal(3) + (normal(1)*gam1*uar3*cr2))*dif3 +
     &        (-normal(2) + (normal(1)*gam1*uar4*cr2))*dif4 -
     &        (normal(1)*gam1*cr2)*dif5)
      dflur1 = dvp1*((normal(1)*(1.0d0 - gam1*qir*cr2) - tet1)*dif1 +
     &         (normal(1)*gam1*uar2*cr2)*dif2 +
     &         (normal(3) + (normal(1)*gam1*uar3*cr2))*dif3 +
     &         (-normal(2) + (normal(1)*gam1*uar4*cr2))*dif4 -
     &         (normal(1)*gam1*cr2)*dif5) + vp1*
     &         ((dnormal(1)*(1.0d0 - gam1*qir*cr2) - normal(1)*
     &         (gam1*dqir*cr2+gam1*qir*dcr2) - dtet1)*dif1 +
     &         (normal(1)*(1.0d0 - gam1*qir*cr2) - tet1)*ddif1 +
     &         (dnormal(1)*gam1*uar2*cr2 + normal(1)*gam1*duar2*cr2 + 
     &         normal(1)*gam1*uar2*dcr2)*dif2 +
     &         (normal(1)*gam1*uar2*cr2)*ddif2 +
     &         (dnormal(3) + (dnormal(1)*gam1*uar3*cr2 + normal(1)*
     &         gam1*duar3*cr2 + normal(1)*gam1*uar3*dcr2))*dif3 +
     &         (normal(3) + (normal(1)*gam1*uar3*cr2))*ddif3 +
     &         (-dnormal(2) + (dnormal(1)*gam1*uar4*cr2 + normal(1)*
     &         gam1*duar4*cr2 + normal(1)*gam1*uar4*dcr2))*dif4 +
     &         (-normal(2) + (normal(1)*gam1*uar4*cr2))*ddif4 -
     &         (dnormal(1)*gam1*cr2 + normal(1)*gam1*dcr2)*dif5 - 
     &         (normal(1)*gam1*cr2)*ddif5)
c
      flur2 = vp1*
     &        ((normal(2)*(1.0d0 - gam1*qir*cr2) - tet2)*dif1 +
     &        (-normal(3) + (normal(2)*gam1*uar2*cr2))*dif2 +
     &        (normal(2)*gam1*uar3*cr2)*dif3 +
     &        (normal(1)  + (normal(2)*gam1*uar4*cr2))*dif4 -
     &        (normal(2)*gam1*cr2)*dif5)
      dflur2 = dvp1*
     &         ((normal(2)*(1.0d0 - gam1*qir*cr2) - tet2)*dif1 +
     &         (-normal(3) + (normal(2)*gam1*uar2*cr2))*dif2 +
     &         (normal(2)*gam1*uar3*cr2)*dif3 +
     &         (normal(1)  + (normal(2)*gam1*uar4*cr2))*dif4 -
     &         (normal(2)*gam1*cr2)*dif5) + vp1*
     &         ((dnormal(2)*(1.0d0 - gam1*qir*cr2) - normal(2)*
     &         (gam1*dqir*cr2 + gam1*qir*dcr2) - dtet2)*dif1 + 
     &         (normal(2)*(1.0d0 - gam1*qir*cr2) - tet2)*ddif1 +
     &         (-dnormal(3) + (dnormal(2)*gam1*uar2*cr2 + normal(2)*
     &         gam1*duar2*cr2 + normal(2)*gam1*uar2*dcr2))*dif2 +
     &         (-normal(3) + (normal(2)*gam1*uar2*cr2))*ddif2 +
     &         (dnormal(2)*gam1*uar3*cr2 + normal(2)*gam1*duar3*cr2 + 
     &         normal(2)*gam1*uar3*dcr2)*dif3 +
     &         (normal(2)*gam1*uar3*cr2)*ddif3 +
     &         (dnormal(1) + (dnormal(2)*gam1*uar4*cr2 + normal(2)*
     &         gam1*duar4*cr2 + normal(2)*gam1*uar4*dcr2))*dif4 +
     &         (normal(1) + (normal(2)*gam1*uar4*cr2))*ddif4 -
     &         (dnormal(2)*gam1*cr2 + normal(2)*gam1*dcr2)*dif5 - 
     &         (normal(2)*gam1*cr2)*ddif5)
c
      flur3 = vp1*
     &        ((normal(3)*(1.0d0 - gam1*qir*cr2) - tet3)*dif1 +
     &        (normal(2)  + (normal(3)*gam1*uar2*cr2))*dif2 +
     &        (-normal(1) + (normal(3)*gam1*uar3*cr2))*dif3 +
     &        (normal(3)*gam1*uar4*cr2)*dif4 -
     &        (normal(3)*gam1*cr2)*dif5)
      dflur3 = dvp1*
     &         ((normal(3)*(1.0d0 - gam1*qir*cr2) - tet3)*dif1 +
     &         (normal(2)  + (normal(3)*gam1*uar2*cr2))*dif2 +
     &         (-normal(1) + (normal(3)*gam1*uar3*cr2))*dif3 +
     &         (normal(3)*gam1*uar4*cr2)*dif4 -
     &         (normal(3)*gam1*cr2)*dif5) + vp1*
     &         ((dnormal(3)*(1.0d0 - gam1*qir*cr2) - normal(3)*
     &         (gam1*dqir*cr2 + gam1*qir*dcr2) - dtet3)*dif1 +
     &         (normal(3)*(1.0d0 - gam1*qir*cr2) - tet3)*ddif1 +
     &         (dnormal(2) + (dnormal(3)*gam1*uar2*cr2 + normal(3)*
     &         gam1*duar2*cr2 + normal(3)*gam1*uar2*dcr2))*dif2 +
     &         (normal(2) + (normal(3)*gam1*uar2*cr2))*ddif2 +
     &         (-dnormal(1) + (dnormal(3)*gam1*uar3*cr2 + normal(3)*
     &         gam1*duar3*cr2 + normal(3)*gam1*uar3*dcr2))*dif3 +
     &         (-normal(1) + (normal(3)*gam1*uar3*cr2))*ddif3 +
     &         (dnormal(3)*gam1*uar4*cr2 + normal(3)*gam1*duar4*cr2 + 
     &         normal(3)*gam1*uar4*dcr2)*dif4 +
     &         (normal(3)*gam1*uar4*cr2)*ddif4 -
     &         (dnormal(3)*gam1*cr2 + normal(3)*gam1*dcr2)*dif5 - 
     &         (normal(3)*gam1*cr2)*ddif5)
c
      flur4 = vp4*
     &        ((-cr*VdotN   + gam1*qir)*dif1 +
     &        ( cr*normal(1) - gam1*uar2)*dif2 +
     &        ( cr*normal(2) - gam1*uar3)*dif3 +
     &        ( cr*normal(3) - gam1*uar4)*dif4 +
     &        gam1*dif5)
      dflur4 = dvp4*
     &         ((-cr*VdotN   + gam1*qir)*dif1 +
     &         ( cr*normal(1) - gam1*uar2)*dif2 +
     &         ( cr*normal(2) - gam1*uar3)*dif3 +
     &         ( cr*normal(3) - gam1*uar4)*dif4 +
     &         gam1*dif5) + vp4*
     &         ((-dcr*VdotN - cr*dVdotN + gam1*dqir)*dif1 +
     &         (-cr*VdotN   + gam1*qir)*ddif1 +
     &         ( dcr*normal(1) + cr*dnormal(1) - gam1*duar2)*dif2 +
     &         ( cr*normal(1) - gam1*uar2)*ddif2 +
     &         ( dcr*normal(2) + cr*dnormal(2) - gam1*duar3)*dif3 +
     &         ( cr*normal(2) - gam1*uar3)*ddif3 +
     &         ( dcr*normal(3) + cr*dnormal(3) - gam1*duar4)*dif4 +
     &         ( cr*normal(3) - gam1*uar4)*ddif4 +
     &         gam1*ddif5)
c
      flur5 = vp5*
     &        (( cr*VdotN  + gam1* qir)*dif1 +
     &        (-cr*normal(1) - gam1*uar2)*dif2 +
     &        (-cr*normal(2) - gam1*uar3)*dif3 +
     &        (-cr*normal(3) - gam1*uar4)*dif4 +
     &        gam1*dif5)
      dflur5 = dvp5*
     &         (( cr*VdotN  + gam1* qir)*dif1 +
     &         (-cr*normal(1) - gam1*uar2)*dif2 +
     &         (-cr*normal(2) - gam1*uar3)*dif3 +
     &         (-cr*normal(3) - gam1*uar4)*dif4 +
     &         gam1*dif5) + vp5*
     &         (( dcr*VdotN + cr*dVdotN + gam1* dqir)*dif1 +
     &         ( cr*VdotN  + gam1* qir)*ddif1 +
     &         (-dcr*normal(1) - cr*dnormal(1) - gam1*duar2)*dif2 +
     &         (-cr*normal(1) - gam1*uar2)*ddif2 +
     &         (-dcr*normal(2) - cr*dnormal(2) - gam1*duar3)*dif3 +
     &         (-cr*normal(2) - gam1*uar3)*ddif3 +
     &         (-dcr*normal(3) - cr*dnormal(3) - gam1*duar4)*dif4 +
     &         (-cr*normal(3) - gam1*uar4)*ddif4 +
     &         gam1*ddif5)
c
      phi(1) = phi(1) - (normal(1)*flur1 +
     &         normal(2)*flur2 + normal(3)*flur3 +
     &         0.5d0*(flur4 + flur5)*cr2)
      dphi(1) = dphi(1) - (dnormal(1)*flur1 + normal(1)*dflur1 +
     &          dnormal(2)*flur2 + normal(2)*dflur2 + dnormal(3)*flur3 +
     &          normal(3)*dflur3 + 0.5d0*(dflur4 + dflur5)*cr2 + 
     &          0.5d0*(flur4 + flur5)*dcr2)
c
      phi(2) = phi(2) - (  
     &         (uar2*normal(1))*flur1 +
     &         (uar2*normal(2) - normal(3))*flur2 +
     &         (uar2*normal(3) + normal(2))*flur3 +
     &         0.5d0*normal(1)*(flur4 - flur5)/cr +
     &         0.5d0*uar2*(flur4 + flur5)*cr2)
      dphi(2) = dphi(2) - (  
     &          (duar2*normal(1) + uar2*dnormal(1))*flur1 +
     &          (uar2*normal(1))*dflur1 +
     &          (duar2*normal(2) + uar2*dnormal(2) - dnormal(3))*flur2 +
     &          (uar2*normal(2) - normal(3))*dflur2 +
     &          (duar2*normal(3) + uar2*dnormal(3) + dnormal(2))*flur3 +
     &          (uar2*normal(3) + normal(2))*dflur3 +
     &          0.5d0*dnormal(1)*(flur4 - flur5)/cr + 
     &          0.5d0*normal(1)*(dflur4 - dflur5)/cr -
     &          0.5d0*normal(1)*(flur4 - flur5)*dcr/(cr*cr) +
     &          0.5d0*duar2*(flur4 + flur5)*cr2 + 0.5d0*uar2*
     &          (dflur4 + dflur5)*cr2 + 0.5d0*uar2*(flur4 + flur5)*dcr2)
c
      phi(3) = phi(3) - (
     &         (uar3*normal(1) + normal(3))*flur1 +
     &         (uar3*normal(2))*flur2 +
     &         (uar3*normal(3) - normal(1))*flur3 +
     &         0.5d0*normal(2)*(flur4 - flur5)/cr +
     &         0.5d0*uar3*(flur4 + flur5)*cr2)
      dphi(3) = dphi(3) - (
     &          (duar3*normal(1) + uar3*dnormal(1) + dnormal(3))*flur1 +
     &          (uar3*normal(1) + normal(3))*dflur1 +
     &          (duar3*normal(2) + uar3*dnormal(2))*flur2 +
     &          (uar3*normal(2))*dflur2 +
     &          (duar3*normal(3) + uar3*dnormal(3) - dnormal(1))*flur3 +
     &          (uar3*normal(3) - normal(1))*dflur3 +
     &          0.5d0*dnormal(2)*(flur4 - flur5)/cr + 
     &          0.5d0*normal(2)*(dflur4 - dflur5)/cr - 
     &          0.5d0*normal(2)*(flur4 - flur5)*dcr/(cr*cr) +
     &          0.5d0*duar3*(flur4 + flur5)*cr2 + 0.5d0*uar3*
     &          (dflur4 + dflur5)*cr2 + 0.5d0*uar3*(flur4 + flur5)*dcr2)
c
      phi(4) = phi(4) - ( 
     &         (uar4*normal(1) - normal(2))*flur1 +
     &         (uar4*normal(2) + normal(1))*flur2 +
     &         (uar4*normal(3))*flur3 +
     &         0.5d0*normal(3)*(flur4 - flur5)/cr +
     &         0.5d0*uar4*(flur4 + flur5)*cr2)
      dphi(4) = dphi(4) - ( 
     &          (duar4*normal(1) + uar4*dnormal(1) - dnormal(2))*flur1 +
     &          (uar4*normal(1) - normal(2))*dflur1 +
     &          (duar4*normal(2) + uar4*dnormal(2) + dnormal(1))*flur2 +
     &          (uar4*normal(2) + normal(1))*dflur2 +
     &          (duar4*normal(3) + uar4*dnormal(3))*flur3 +
     &          (uar4*normal(3))*dflur3 +
     &          0.5d0*dnormal(3)*(flur4 - flur5)/cr +
     &          0.5d0*normal(3)*(dflur4 - dflur5)/cr -
     &          0.5d0*normal(3)*(flur4 - flur5)*dcr/(cr*cr) +
     &          0.5d0*duar4*(flur4 + flur5)*cr2 + 0.5d0*uar4*
     &          (dflur4 + dflur5)*cr2 + 0.5d0*uar4*(flur4 + flur5)*dcr2)
c
      phi(5) = phi(5) - (
     &         (qir*normal(1) + tet1)*flur1 +
     &         (qir*normal(2) + tet2)*flur2 +
     &         (qir*normal(3) + tet3)*flur3 +
     &         0.5d0*VdotN*(flur4 - flur5)/cr +
     &         0.5d0*uar5*(flur4 + flur5)*cr2)
      dphi(5) = dphi(5) - (
     &          (dqir*normal(1) + qir*dnormal(1) + dtet1)*flur1 +
     &          (qir*normal(1) + tet1)*dflur1 +
     &          (dqir*normal(2) + qir*dnormal(2) + dtet2)*flur2 +
     &          (qir*normal(2) + tet2)*dflur2 +
     &          (dqir*normal(3) + qir*dnormal(3) + dtet3)*flur3 +
     &          (qir*normal(3) + tet3)*dflur3 +
     &          0.5d0*dVdotN*(flur4 - flur5)/cr + 
     &          0.5d0*VdotN*(dflur4 - dflur5)/cr -
     &          0.5d0*VdotN*(flur4 - flur5)*dcr/(cr*cr) +
     &          0.5d0*duar5*(flur4 + flur5)*cr2 + 0.5d0*uar5*
     &          (dflur4 + dflur5)*cr2 + 0.5d0*uar5*(flur4 + flur5)*dcr2)
c
      dphi(1) = dphi(1)*0.5*rnorm + phi(1)*0.5*drnorm
      phi(1) = phi(1)*0.5*rnorm
      dphi(2) = dphi(2)*0.5*rnorm + phi(2)*0.5*drnorm
      phi(2) = phi(2)*0.5*rnorm
      dphi(3) = dphi(3)*0.5*rnorm + phi(3)*0.5*drnorm
      phi(3) = phi(3)*0.5*rnorm
      dphi(4) = dphi(4)*0.5*rnorm + phi(4)*0.5*drnorm
      phi(4) = phi(4)*0.5*rnorm
      dphi(5) = dphi(5)*0.5*rnorm + phi(5)*0.5*drnorm
      phi(5) = phi(5)*0.5*rnorm

      if (type.eq.1) then
        updir = 0.5d0 + dsign(0.5d0, phi(1))
        phi(6) = phi(1) * (updir * Ugr(6) + (1.0d0 - updir) * Udr(6))
        dphi(6) = dphi(1) * (updir * Ugr(6) + (1.0d0 - updir) * Udr(6))+
     &  	  phi(1) * (updir * dUgr(6) + (1.0d0 - updir) * dUdr(6))
      else if (type.eq.2) then
        updir = 0.5d0 + dsign(0.5d0, phi(1))
        phi(6) = phi(1) * (updir * Ugr(6) + (1.0d0 - updir) * Udr(6))
        dphi(6) = dphi(1) * (updir * Ugr(6) + (1.0d0 - updir) * Udr(6))+
     &            phi(1) * (updir * dUgr(6) + (1.0d0 - updir) * dUdr(6))
        phi(7) = phi(1) * (updir * Ugr(7) + (1.0d0 - updir) * Udr(7))
        dphi(7) = dphi(1) * (updir * Ugr(7) + (1.0d0 - updir) * Udr(7))+
     &            phi(1) * (updir * dUgr(7) + (1.0d0 - updir) * dUdr(7))
      endif
c
      END
