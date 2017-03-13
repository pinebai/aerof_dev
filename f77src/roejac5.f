      SUBROUTINE ROEJAC5(type,gamma,gam,pstiff,enormal,evitno,Ug,Ud,jac)
c-----------------------------------------------------------------------
c This routine computes the derivative of the flux of Roe 
c with respect to the first variable taken at the vectors Ug, Ud.
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
C We use a parabolic regularization for the flux.
c-----------------------------------------------------------------------
c    INCLUDE 'Param3D.h'
      implicit none
c-----------------------------------------------------------------------
      REAL*8 Ug(*), Ud(*), normal(3), enormal(3), H(5,5), jac(*)
      REAL*8 enth, pstiff
      REAL*8 VdotN , rnorm, evitno, vitno, invnorm
      REAL*8 g_VdotN(5)
      REAL*8 flur1 , flur2 , flur3 , flur4 , flur5
      REAL*8 g_flur1(5),g_flur2(5),g_flur3(5),g_flur4(5),g_flur5(5)
      REAL*8 dif1 , dif2 , dif3 , dif4 ,dif5
      REAL*8 g_dif1(5),g_dif2(5),g_dif3(5),g_dif4(5),g_dif5(5)
      REAL*8 cr , cr2 , qir
      REAL*8 g_cr(5),g_cr2(5),cdiff,g_qir(5)
      REAL*8 vp1 , vp4 , vp5 , svp1 , svp4 , svp5
      REAL*8 deltar , switch , deltarinv , g_deltar
      REAL*8 g_vp1(5),g_vp4(5),g_vp5(5) 
      REAL*8 ener1 , ener2
      REAL*8 uar2 , uar3 , uar4 , uar5
      REAL*8 g_uar2(5),g_uar3(5),g_uar4(5),g_uar5(5)
      REAL*8 usro , squsr1 , squsr2
      REAL*8 tet1 , tet2 , tet3
      REAL*8 g_tet1(5),g_tet2(5),g_tet3(5)
      INTEGER k,i,type,dim
      INTEGER ient
      REAL*8 gam, gam1, gamma, epsiim
c
c Initialisation
c
      ient = 0
      gam1 = gam - 1.d0
      epsiim = 0.1d0

      rnorm = DSQRT(enormal(1)**2 + enormal(2)**2 + enormal(3)**2)
      invnorm = 1.0d0 / rnorm
c
      normal(1) = enormal(1) * invnorm
      normal(2) = enormal(2) * invnorm
      normal(3) = enormal(3) * invnorm

      vitno = evitno * invnorm
c
c Computation of the centred term
c
      VdotN = Ug(2)*normal(1) + Ug(3)*normal(2) + Ug(4)*normal(3)
      qir = 0.5*(Ug(2)**2 + Ug(3)**2 + Ug(4)**2)
      enth = gam*(Ug(5)+pstiff)+gam1*Ug(1)*qir
      enth = enth/(gam1*Ug(1))
c
      H(1,1) = VdotN - vitno
      H(1,2) = Ug(1)*normal(1)
      H(1,3) = Ug(1)*normal(2)
      H(1,4) = Ug(1)*normal(3)
      H(1,5) = 0.0 
c
      H(2,1) = Ug(2)*(VdotN-vitno)
      H(2,2) = H(1,2)*Ug(2)+Ug(1)*(VdotN-vitno)
      H(2,3) = H(1,3)*Ug(2)
      H(2,4) = H(1,4)*Ug(2)
      H(2,5) = normal(1)
c
      H(3,1) = Ug(3)*(VdotN-vitno)
      H(3,2) = H(1,2)*Ug(3)
      H(3,3) = H(1,3)*Ug(3)+Ug(1)*(VdotN-vitno)
      H(3,4) = H(1,4)*Ug(3)
      H(3,5) = normal(2)
c
      H(4,1) = Ug(4)*(VdotN-vitno)
      H(4,2) = H(1,2)*Ug(4)
      H(4,3) = H(1,3)*Ug(4)
      H(4,4) = H(1,4)*Ug(4)+Ug(1)*(VdotN-vitno)
      H(4,5) = normal(3)
c
      H(5,1) = qir*(VdotN-vitno)
      H(5,2) = Ug(1)*Ug(2)*(VdotN-vitno)+H(1,2)*enth
      H(5,3) = Ug(1)*Ug(3)*(VdotN-vitno)+H(1,3)*enth
      H(5,4) = Ug(1)*Ug(4)*(VdotN-vitno)+H(1,4)*enth
      H(5,5) = gam*(VdotN-vitno)/gam1+vitno
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
      usro   = 1.0/(squsr1 + squsr2)
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
      squsr1 = 0.5*usro/squsr1
c
      do k=1,5
c
         g_uar2(k) = 0.0
         g_uar3(k) = 0.0
         g_uar4(k) = 0.0
         g_uar5(k) = 0.0
c
      end do
c
      g_uar2(1) = squsr1*(Ug(2)-uar2)
      g_uar2(2) = 2*Ug(1)*squsr1
c
      g_uar3(1) = squsr1*(Ug(3)-uar3)
      g_uar3(3) = g_uar2(2)
c
      g_uar4(1) = squsr1*(Ug(4)-uar4)
      g_uar4(4) = g_uar2(2)
c
      g_uar5(1) = squsr1*(ener1 - uar5*Ug(1) -
     &            (gam+1)*Ug(5)/gam1-2.0*gam*pstiff/gam1)/Ug(1)
      g_uar5(2) = g_uar2(2)*Ug(2)
      g_uar5(3) = g_uar2(2)*Ug(3)
      g_uar5(4) = g_uar2(2)*Ug(4)
      g_uar5(5) = gam*2*squsr1/gam1
c
c Computation of the dissipation term
c
      VdotN = normal(1)*uar2 + normal(2)*uar3 + normal(3)*uar4
c
      qir   = 0.5*(uar2*uar2 + uar3*uar3 + uar4*uar4)
c
      tet1  = normal(3)*uar3 - normal(2)*uar4
      tet2  = normal(1)*uar4 - normal(3)*uar2
      tet3  = normal(2)*uar2 - normal(1)*uar3
c
      cr2   = gam1*(uar5 - qir)
      cr    = DSQRT(cr2)
      cr2   = 1.0/cr2
c
      dif1  = - Ug(1) + Ud(1)
      dif2  = - Ug(1)*Ug(2) + Ud(1)*Ud(2)
      dif3  = - Ug(1)*Ug(3) + Ud(1)*Ud(3)
      dif4  = - Ug(1)*Ug(4) + Ud(1)*Ud(4)
      dif5  = - ener1 + ener2
c
      vp1   = gamma*(VdotN - vitno)
      vp4   = gamma*(VdotN - vitno + cr)
      vp5   = gamma*(VdotN - vitno - cr)
c
      svp1  = DSIGN(1.0d0,vp1)
      svp4  = DSIGN(1.0d0,vp4)
      svp5  = DSIGN(1.0d0,vp5)
c
      do i=1,5
c
         g_VdotN(i) = normal(1)*g_uar2(i) + normal(2)*g_uar3(i) +
     &                normal(3)*g_uar4(i)
c
         g_qir(i) = uar2*g_uar2(i)+uar3*g_uar3(i)+uar4*g_uar4(i)
c
         g_tet1(i) = normal(3)*g_uar3(i)-normal(2)*g_uar4(i)
         g_tet2(i) = normal(1)*g_uar4(i)-normal(3)*g_uar2(i)
         g_tet3(i) = normal(2)*g_uar2(i)-normal(1)*g_uar3(i)
c
         cdiff = gam1*(g_uar5(i)-g_qir(i))
         g_cr(i) = 0.5*cdiff/cr
         g_cr2(i) = -cdiff*cr2*cr2
c
         g_dif1(i) = 0.0
         g_dif2(i) = 0.0
         g_dif3(i) = 0.0
         g_dif4(i) = 0.0
         g_dif5(i) = 0.0
c
         g_vp1(i) = gamma*g_VdotN(i)
         g_vp4(i) = gamma*(g_VdotN(i)+g_cr(i))
         g_vp5(i) = gamma*(g_VdotN(i)-g_cr(i))
c
         if (ient .eq. 1) then
c
            deltar = epsiim*vp4
            deltarinv = 1.0/deltar
            g_deltar = epsiim*g_vp4(i)
c
            switch = 0.5 + DSIGN(0.5d0, svp1*vp1 - deltar)
            g_vp1(i) = 0.5*(g_deltar-g_deltar*(vp1**2)*(deltarinv**2)+
     &              2*g_vp1(i)*vp1*deltarinv)+
     &              switch*(svp1*vp1-deltar)*(
     &               0.5*g_deltar*(deltarinv**2)*(svp1*vp1-deltar)
     &               -deltarinv*(svp1*g_vp1(i)-g_deltar))
c
            switch = 0.5 + DSIGN(0.5d0, svp4*vp4 - deltar)
            g_vp4(i) = 0.5*(g_deltar-g_deltar*(vp4**2)*(deltarinv**2)+
     &              2*g_vp4(i)*vp4*deltarinv)+
     &              switch*(svp4*vp4-deltar)*(
     &               0.5*g_deltar*(deltarinv**2)*(svp4*vp4-deltar)
     &               -deltarinv*(svp4*g_vp4(i)-g_deltar))
c
            switch = 0.5 + DSIGN(0.5d0, svp5*vp5 - deltar)
            g_vp5(i) = 0.5*(g_deltar-g_deltar*(vp5**2)*(deltarinv**2)+
     &              2*g_vp5(i)*vp5*deltarinv)+
     &              switch*(svp5*vp5-deltar)*(
     &               0.5*g_deltar*(deltarinv**2)*(svp5*vp5-deltar)
     &               -deltarinv*(svp5*g_vp5(i)-g_deltar))
c
         else
c
            g_vp1(i) = svp1*g_vp1(i)
            g_vp4(i) = svp4*g_vp4(i)
            g_vp5(i) = svp5*g_vp5(i)
c
         endif
c
      end do
c
      if (ient .eq. 1) then 
c
         deltar = epsiim*vp4
c
         switch = 0.5 + DSIGN(0.5d0, svp1*vp1 - deltar)
         vp1    = switch*svp1*vp1 + 0.5*(1.0 - switch)*
     &               (deltar + vp1*vp1/deltar)
c
         switch = 0.5 + DSIGN(0.5d0, svp4*vp4 - deltar)
         vp4    = switch*svp4*vp4 + 0.5*(1.0 - switch)*
     &               (deltar + vp4*vp4/deltar)
c
         switch = 0.5 + DSIGN( 0.5d0, svp5*vp5 - deltar)
         vp5    = switch*svp5*vp5 + 0.5*(1.0 - switch)*
     &               (deltar + vp5*vp5/deltar)
c
      else
c
         vp1 = vp1*svp1
         vp4 = vp4*svp4
         vp5 = vp5*svp5
c
      endif
c
      g_dif1(1) = - 1.0
c
      g_dif2(1) = - Ug(2)
      g_dif2(2) = - Ug(1)
c
      g_dif3(1) = - Ug(3)
      g_dif3(3) = - Ug(1)
c
      g_dif4(1) = - Ug(4)
      g_dif4(4) = - Ug(1)
c
      g_dif5(1) = - 0.5*(Ug(2)**2+Ug(3)**2+Ug(4)**2)
      g_dif5(2) = - Ug(1)*Ug(2)
      g_dif5(3) = - Ug(1)*Ug(3)
      g_dif5(4) = - Ug(1)*Ug(4)
      g_dif5(5) = - 1.0/gam1
c
      flur1     = vp1*
     &              ((normal(1)*(1.0 - gam1*qir*cr2) - tet1)*dif1 +
     &               (normal(1)*gam1*uar2*cr2)*dif2  +
     &               (normal(3)  + (normal(1)*gam1*uar3*cr2))*dif3   +
     &               (-normal(2) + (normal(1)*gam1*uar4*cr2))*dif4   -
     &               (normal(1)*gam1*cr2)*dif5)
c
      do i=1,5
c
         g_flur1(i) = g_vp1(i)*
     &              ((normal(1)*(1.0 - gam1*qir*cr2) - tet1)*dif1 +
     &               (normal(1)*gam1*uar2*cr2)*dif2  +
     &               (normal(3)  + (normal(1)*gam1*uar3*cr2))*dif3   +
     &               (-normal(2) + (normal(1)*gam1*uar4*cr2))*dif4   -
     &               (normal(1)*gam1*cr2)*dif5)
c
         g_flur1(i) = g_flur1(i) + vp1*(
     &               (normal(1)*(1.0 - gam1*qir*cr2) - tet1)*g_dif1(i)+
     &               dif1*(normal(1)*(-gam1*g_qir(i)*cr2-gam1*qir*
     &               g_cr2(i))-g_tet1(i)))
c
         g_flur1(i) = g_flur1(i) + vp1*(
     &              g_dif2(i)*(normal(1)*gam1*uar2*cr2)+
     &              dif2*(normal(1)*gam1*(g_uar2(i)*cr2+
     &              uar2*g_cr2(i))))
c
         g_flur1(i) = g_flur1(i) + vp1*(
     &              g_dif3(i)*(normal(3)+(normal(1)*gam1*uar3*cr2))+
     &              dif3*(normal(1)*gam1*(g_uar3(i)*cr2+
     &              uar3*g_cr2(i))))
c
         g_flur1(i) = g_flur1(i) + vp1*(
     &              g_dif4(i)*(-normal(2)+(normal(1)*gam1*uar4*cr2))+
     &              dif4*(normal(1)*gam1*(g_uar4(i)*cr2+
     &              uar4*g_cr2(i))))
c
         g_flur1(i) = g_flur1(i) - vp1*normal(1)*gam1*(
     &              g_dif5(i)*cr2+dif5*g_cr2(i))
c
      end do
c
      flur2    = vp1*
     &              ((normal(2)*(1.0 - gam1*qir*cr2) - tet2)*dif1 +
     &               (-normal(3) + (normal(2)*gam1*uar2*cr2))*dif2   +
     &               (normal(2)*gam1*uar3*cr2)*dif3  +
     &               (normal(1)  + (normal(2)*gam1*uar4*cr2))*dif4   -
     &               (normal(2)*gam1*cr2)*dif5)
c
      do i=1,5
c
         g_flur2(i) = g_vp1(i)*
     &              ((normal(2)*(1.0 - gam1*qir*cr2) - tet2)*dif1 +
     &               (-normal(3) + (normal(2)*gam1*uar2*cr2))*dif2   +
     &               (normal(2)*gam1*uar3*cr2)*dif3  +
     &               (normal(1)  + (normal(2)*gam1*uar4*cr2))*dif4   -
     &               (normal(2)*gam1*cr2)*dif5)
c
         g_flur2(i) = g_flur2(i)+vp1*(
     &              (normal(2)*(1.0 - gam1*qir*cr2) - tet2)*g_dif1(i)+
     &              dif1*(normal(2)*(-gam1*g_qir(i)*cr2-gam1*qir*
     &               g_cr2(i))-g_tet2(i)))
c
         g_flur2(i) = g_flur2(i)+vp1*(
     &              (-normal(3) + (normal(2)*gam1*uar2*cr2))*g_dif2(i)+
     &              dif2*(normal(2)*gam1*(g_uar2(i)*cr2+
     &              uar2*g_cr2(i))))
c
         g_flur2(i) = g_flur2(i)+vp1*(
     &              (normal(2)*gam1*uar3*cr2)*g_dif3(i)+
     &              dif3*(normal(2)*gam1*(g_uar3(i)*cr2+
     &              uar3*g_cr2(i))))
c
         g_flur2(i) = g_flur2(i)+vp1*(
     &              (normal(1)  + (normal(2)*gam1*uar4*cr2))*g_dif4(i)+
     &              dif4*(normal(2)*gam1*(g_uar4(i)*cr2+
     &              uar4*g_cr2(i))))
c
         g_flur2(i) = g_flur2(i)-vp1*normal(2)*gam1*(
     &              cr2*g_dif5(i)+dif5*g_cr2(i))
c
      end do
c
      flur3    = vp1*
     &              ((normal(3)*(1.0 - gam1*qir*cr2) - tet3)*dif1 +
     &               (normal(2)  + (normal(3)*gam1*uar2*cr2))*dif2   +
     &               (-normal(1) + (normal(3)*gam1*uar3*cr2))*dif3   +
     &               (normal(3)*gam1*uar4*cr2)*dif4  -
     &               (normal(3)*gam1*cr2)*dif5)
c
      do i=1,5
c
         g_flur3(i)=g_vp1(i)*
     &              ((normal(3)*(1.0 - gam1*qir*cr2) - tet3)*dif1 +
     &               (normal(2)  + (normal(3)*gam1*uar2*cr2))*dif2   +
     &               (-normal(1) + (normal(3)*gam1*uar3*cr2))*dif3   +
     &               (normal(3)*gam1*uar4*cr2)*dif4  -
     &               (normal(3)*gam1*cr2)*dif5)
c
         g_flur3(i) = g_flur3(i)+vp1*(
     &              g_dif1(i)*(normal(3)*(1.0 - gam1*qir*cr2) - tet3)+
     &              dif1*(normal(3)*gam1*(-g_qir(i)*cr2-qir*g_cr2(i))-
     &              g_tet3(i)))
c
         g_flur3(i)=g_flur3(i)+vp1*(
     &              g_dif2(i)*(normal(2)  + (normal(3)*gam1*uar2*cr2))+
     &              dif2*(normal(3)*gam1*(g_uar2(i)*cr2+
     &              uar2*g_cr2(i))))
c
         g_flur3(i)=g_flur3(i)+vp1*(
     &              g_dif3(i)*(-normal(1) + (normal(3)*gam1*uar3*cr2))+
     &              dif3*(normal(3)*gam1*(g_uar3(i)*cr2+
     &              uar3*g_cr2(i))))
c
         g_flur3(i)=g_flur3(i)+vp1*(
     &              g_dif4(i)*(normal(3)*gam1*uar4*cr2)+
     &              dif4*(normal(3)*gam1*(g_uar4(i)*cr2+
     &              uar4*g_cr2(i))))
c
         g_flur3(i)=g_flur3(i)-vp1*normal(3)*gam1*(
     &             g_dif5(i)*cr2+dif5*g_cr2(i))
c 
      end do
c
      flur4                     = vp4*
     &              ((-cr*VdotN   + gam1*qir)*dif1  +
     &               ( cr*normal(1) - gam1*uar2)*dif2 +
     &               ( cr*normal(2) - gam1*uar3)*dif3 +
     &               ( cr*normal(3) - gam1*uar4)*dif4 +
     &               gam1*dif5)
c
      do i=1,5
c
         g_flur4(i)=g_vp4(i)*
     &              ((-cr*VdotN   + gam1*qir)*dif1  +
     &               ( cr*normal(1) - gam1*uar2)*dif2 +
     &               ( cr*normal(2) - gam1*uar3)*dif3 +
     &               ( cr*normal(3) - gam1*uar4)*dif4 +
     &               gam1*dif5)
c
         g_flur4(i)=g_flur4(i)+vp4*(
     &              g_dif1(i)*(-cr*VdotN   + gam1*qir)+
     &              dif1*(gam1*g_qir(i)-cr*g_VdotN(i)-
     &              g_cr(i)*VdotN))
c
         g_flur4(i)=g_flur4(i)+vp4*(
     &              g_dif2(i)*( cr*normal(1) - gam1*uar2)+
     &              dif2*(normal(1)*g_cr(i)-gam1*g_uar2(i)))
c
         g_flur4(i)=g_flur4(i)+vp4*(
     &              g_dif3(i)*( cr*normal(2) - gam1*uar3)+
     &              dif3*(normal(2)*g_cr(i)-gam1*g_uar3(i)))
c
         g_flur4(i)=g_flur4(i)+vp4*(
     &              g_dif4(i)*( cr*normal(3) - gam1*uar4)+
     &              dif4*(normal(3)*g_cr(i)-gam1*g_uar4(i)))
c
         g_flur4(i)=g_flur4(i)+vp4*gam1*g_dif5(i)
c 
      end do
c
      flur5                     = vp5*
     &              (( cr*VdotN  + gam1* qir)*dif1 +
     &               (-cr*normal(1) - gam1*uar2)*dif2 +
     &               (-cr*normal(2) - gam1*uar3)*dif3 +
     &               (-cr*normal(3) - gam1*uar4)*dif4 +
     &               gam1*dif5)
c
      do i=1,5
c
         g_flur5(i)=g_vp5(i)*
     &              ((cr*VdotN   + gam1*qir)*dif1  +
     &               (-cr*normal(1) - gam1*uar2)*dif2 +
     &               (-cr*normal(2) - gam1*uar3)*dif3 +
     &               (-cr*normal(3) - gam1*uar4)*dif4 +
     &               gam1*dif5)
c
         g_flur5(i)=g_flur5(i)+vp5*(
     &              g_dif1(i)*(cr*VdotN + gam1*qir)+
     &              dif1*(gam1*g_qir(i)+cr*g_VdotN(i)+
     &              g_cr(i)*VdotN))
c
         g_flur5(i)=g_flur5(i)+vp5*(
     &              g_dif2(i)*(-cr*normal(1) - gam1*uar2)+
     &              dif2*(-normal(1)*g_cr(i)-gam1*g_uar2(i)))
c
         g_flur5(i)=g_flur5(i)+vp5*(
     &              g_dif3(i)*(-cr*normal(2) - gam1*uar3)+
     &              dif3*(-normal(2)*g_cr(i)-gam1*g_uar3(i)))
c
         g_flur5(i)=g_flur5(i)+vp5*(
     &              g_dif4(i)*(-cr*normal(3) - gam1*uar4)+
     &              dif4*(-normal(3)*g_cr(i)-gam1*g_uar4(i)))
c
         g_flur5(i)=g_flur5(i)+vp5*gam1*g_dif5(i)
c
      end do
c
c Final operation
c
      do i=1,5
c
         H(1,i) = H(1,i) - (normal(1)*g_flur1(i)+
     &       normal(2)*g_flur2(i) + normal(3)*g_flur3(i) +
     &       0.5*(g_flur4(i)+g_flur5(i))*cr2+0.5*(flur4+flur5)*g_cr2(i))
c
         H(2,i) = H(2,i)-(g_uar2(i)*normal(1)*flur1+
     &        uar2*normal(1)*g_flur1(i))
         H(2,i) = H(2,i)-(g_flur2(i)*
     &        (uar2*normal(2) - normal(3))+flur2*normal(2)*
     &        g_uar2(i))
         H(2,i) = H(2,i)-(g_flur3(i)*(uar2*normal(3) + normal(2))
     &        +flur3*normal(3)*g_uar2(i))
         H(2,i) = H(2,i)-0.5*normal(1)*((g_flur4(i)-g_flur5(i))/cr-
     &          (flur4-flur5)*cr2*g_cr(i))
         H(2,i) = H(2,i)-0.5*(g_uar2(i)*(flur4+flur5)*cr2+
     &          uar2*(g_flur4(i)+g_flur5(i))*cr2+uar2*g_cr2(i)*
     &          (flur4+flur5))
c
         H(3,i) = H(3,i)-(g_uar3(i)*normal(1)*flur1+
     &          g_flur1(i)*(uar3*normal(1) + normal(3)))
         H(3,i) = H(3,i)-normal(2)*(g_uar3(i)*flur2+
     &          uar3*g_flur2(i))
         H(3,i) = H(3,i)-(g_uar3(i)*normal(3)*flur3+
     &          g_flur3(i)*(uar3*normal(3) - normal(1)))
         H(3,i) = H(3,i)-0.5*normal(2)*((g_flur4(i)-g_flur5(i))/cr-
     &          (flur4-flur5)*cr2*g_cr(i))
         H(3,i) = H(3,i)-0.5*(g_uar3(i)*(flur4+flur5)*cr2+
     &          uar3*(g_flur4(i)+g_flur5(i))*cr2+uar3*g_cr2(i)*
     &          (flur4+flur5))
c
         H(4,i) = H(4,i)-(g_uar4(i)*normal(1)*flur1+
     &          g_flur1(i)*(uar4*normal(1) - normal(2)))
         H(4,i) = H(4,i)-(g_uar4(i)*normal(2)*flur2+
     &          g_flur2(i)*(uar4*normal(2) + normal(1)))
         H(4,i) = H(4,i)-normal(3)*(g_uar4(i)*flur3+
     &          uar4*g_flur3(i))
         H(4,i) = H(4,i)-0.5*normal(3)*((g_flur4(i)-g_flur5(i))/cr-
     &          (flur4-flur5)*cr2*g_cr(i))
         H(4,i) = H(4,i)-0.5*(g_uar4(i)*(flur4+flur5)*cr2+
     &          uar4*(g_flur4(i)+g_flur5(i))*cr2+uar4*g_cr2(i)*
     &          (flur4+flur5))
c
         H(5,i) = H(5,i)-(g_flur1(i)*(qir*normal(1) + tet1)+
     &          flur1*(g_tet1(i)+g_qir(i)*normal(1)))
         H(5,i) = H(5,i)-(g_flur2(i)*(qir*normal(2) + tet2)+
     &          flur2*(g_tet2(i)+g_qir(i)*normal(2)))
         H(5,i) = H(5,i)-(g_flur3(i)*(qir*normal(3) + tet3)+
     &          flur3*(g_tet3(i)+g_qir(i)*normal(3)))
         H(5,i) = H(5,i)-0.5*(g_VdotN(i)*(flur4 - flur5)/cr+
     &          VdotN*(g_flur4(i)-g_flur5(i))/cr - VdotN*
     &          (flur4 - flur5)*cr2*g_cr(i))
         H(5,i) = H(5,i)-0.5*(g_uar5(i)*(flur4 + flur5)*cr2+
     &          uar5*(g_flur4(i)+g_flur5(i))*cr2+uar5*
     &          (flur4 + flur5)*g_cr2(i))
c
      end do
c
c "Normalization" of the flux
c
      dim = 5 + type
      do i = 1,5
         do k = 1,5
            jac(dim*(i-1) + k) = H(i,k)*0.5*rnorm
         enddo
         do k = 6,dim
            jac(dim*(i-1) + k) = 0.0
         enddo
      enddo
      do i = 6,dim
         do k = 1,dim
            jac(dim*(i-1) + k) = 0.0
         enddo
      enddo
c
c$$$      normal(1) = rnorm * normal(1)
c$$$      normal(2) = rnorm * normal(2)
c$$$      normal(3) = rnorm * normal(3)
c
      end
