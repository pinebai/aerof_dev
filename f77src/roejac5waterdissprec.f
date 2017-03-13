       SUBROUTINE ROEJAC5WATERDISSPREC(type,gamma,Cp,Pr,alpha,beta,
     &     enormal,evitno,Ug,Ud,jac, betaRef, k1, cmach, irey,prec)
c-----------------------------------------------------------------------
c This routine computes the jacobian of the Flux of Roe for a fluid like water,
c taken at the vectors Ug, Ud
c U refers to the primitive variables.
c normal is the normal of the boundary concerned by the flux.
c jac stores the resulting jacobian.
c
c EOS is
c   Pressure = Pr + alpha*Density^beta
c   h = cp*T
c   where h is the internal enthalpy
c         Pr, alpha, beta and cp are constants
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
      REAL*8 Ug(*), Ud(*), enormal(3), evitno, jac(*)
      REAL*8 gamma
      REAL*8 Cp, Pr, alpha, beta
      INTEGER type
      
c indices and dimension
      INTEGER i, k, dim

c for general computation
      REAL*8 nx, ny, nz, rnorm, invnorm, vitno
      REAL*8 beta1, oobeta, oobeta1
      REAL*8 VdotN, vitg2, vitd2, Ugbeta, Pg, Pd
      REAL*8 enth, c2, c
      REAL*8 Hg, Hd
      REAL*8 dif1, dif2, dif3, dif4, dif5
      REAL*8 tet1, tet2, tet3
      REAL*8 flur1, flur2, flur3, flur4, flur5
      REAL*8 vp1, vp4, vp5, signvp1,signvp4, signvp5
      REAL*8 H(5,5)

c for Roe Variables
      REAL*8 uar1, uar2, uar3, uar4, uar5
      REAL*8 squsrg, squsrd, usro, difro, difroB, pre
      REAL*8 coeff1, coeff2, coeff3, coeff4
      REAL*8 coeff5, coeff6, coeff7, eps
      REAL*8 cr2, cr
      

c for differentiated variables    
      REAL*8 g_VdotN(5), g_cr(5)
      REAL*8 g_uar1(5), g_uar2(5), g_uar3(5), g_uar4(5), g_uar5(5)
      REAL*8 g_dif1(5), g_dif2(5), g_dif3(5), g_dif4(5), g_dif5(5)
      REAL*8 g_tet1(5), g_tet2(5), g_tet3(5)
      REAL*8 g_vp1(5), g_vp4(5), g_vp5(5)
      REAL*8 g_flur1(5),g_flur2(5),g_flur3(5),g_flur4(5),g_flur5(5)
      


c for precontioning
      REAL*8 k1, betaRef, cmach, irey
      REAL*8 vpa, vpb
      REAL*8 oob2, b2, betaPrec, locMach
      REAL*8 P1, P2
      REAL*8 r,s,rs,rps,smr,r2,s2,rp,sp
      INTEGER prec
      
c for differentiated preconditioned variables
      REAL*8 g_vpa(5), g_vpb(5)
      REAL*8 g_P1(5), g_P2(5)
      REAL*8 g_r(5),g_s(5), g_r2(5),g_s2(5)
      REAL*8 g_rs(5),g_rps(5),g_smr(5)
      REAL*8 g_rp(5), g_sp(5)
      

c
c Initialisation
c

      beta1 = beta - 1.0d0
      oobeta = 1.0d0/beta
      oobeta1 = 1.0d0/beta1
     
      rnorm = DSQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) +
     &             enormal(3)*enormal(3))
      invnorm = 1.0d0 / rnorm

      nx = enormal(1) * invnorm
      ny = enormal(2) * invnorm
      nz = enormal(3) * invnorm

      vitno = evitno * invnorm
      

c
c Computation of the centred terms
c

      VdotN = Ug(2)*nx + Ug(3)*ny + Ug(4)*nz
      vitg2 = Ug(2)*Ug(2) + Ug(3)*Ug(3) + Ug(4)*Ug(4)
      Ugbeta = Ug(1)**beta
      Pg = Pr + alpha * Ugbeta
      enth = Cp*Ug(5) + 0.5d0*vitg2
      c2 = alpha*beta*(Ug(1)**beta1)
      c = DSQRT(c2)

c Jacobian w.r.t conservative variables
      H(1,1) =  vitno
      H(1,2) = nx
      H(1,3) = ny
      H(1,4) = nz
      H(1,5) = 0.0d0

      H(2,1) = - (VdotN-vitno)*Ug(2) + c2*nx
      H(2,2) = Ug(2)*nx + VdotN
      H(2,3) = Ug(2)*ny
      H(2,4) = Ug(2)*nz
      H(2,5) = 0.0d0

      H(3,1) = - (VdotN-vitno)*Ug(3) + c2*ny
      H(3,2) = Ug(3)*nx
      H(3,3) = Ug(3)*ny + VdotN
      H(3,4) = Ug(3)*nz
      H(3,5) = 0.0d0

      H(4,1) = - (VdotN-vitno)*Ug(4) + c2*nz
      H(4,2) = Ug(4)*nx
      H(4,3) = Ug(4)*ny
      H(4,4) = Ug(4)*nz + VdotN
      H(4,5) = 0.0d0

      H(5,1) = (c2-enth)*VdotN
      H(5,2) = enth*nx
      H(5,3) = enth*ny
      H(5,4) = enth*nz
      H(5,5) = VdotN   

c
c Computation of the Roe-averaged state
c
      squsrg = DSQRT(Ug(1))
      squsrd = DSQRT(Ud(1))
      usro  = 1.0d0/(squsrg + squsrd)
      Hg = enth
      vitd2 = Ud(2)*Ud(2) + Ud(3)*Ud(3) + Ud(4)*Ud(4)
      Pd = Pr + alpha * Ud(1)**beta
      Hd = Cp*Ud(5) + 0.5d0*vitd2


c With the analytical implementation, oscillations appeared as 
c soon as the relative difference between Ud(1) and Ug(1) was
c less than 1e-7.so a first order expansion is enough to get 1e-14 
c precision

      difro = Ud(1)-Ug(1)
      difroB = Ud(1)**beta - Ug(1)**beta
      pre = 0.5d0*usro/squsrg

      if (DABS(Ud(1)-Ug(1)).le.0.005d0*Ug(1)) then
       coeff1 = 0.5d0
       coeff2 = (beta-2.0d0)/24.0d0
       coeff3 = -(beta-2.0d0)/48.0d0
       coeff4 = -(2.0d0*beta**3-3.0d0*beta**2-78.0d0*beta+152.0d0)
     &              / 5760.0d0
       coeff5 = (2.0d0*beta**3-3.0d0*beta**2-38.0d0*beta+72.0d0)
     &              /3840.0d0
       eps = (Ud(1)-Ug(1))/Ug(1)
       uar1 =  1.0d0 + coeff1*eps 
     &        + coeff2*eps**2+ coeff3*eps**3
     &        + coeff4*eps**4+ coeff5*eps**5
       uar1 = Ug(1)*uar1
      else
       uar1  = ( difroB / 
     &        ( beta * (difro) ) ) **
     &        ( 1.0d0 / beta1 )
      endif


      uar2 = (squsrg*Ug(2) + squsrd*Ud(2))*usro
      uar3 = (squsrg*Ug(3) + squsrd*Ud(3))*usro
      uar4 = (squsrg*Ug(4) + squsrd*Ud(4))*usro
      uar5 = (squsrg*Hg + squsrd*Hd)*usro

c differentiation of roe variables
      do k=1,5
         g_uar1(k) = 0.0d0
         g_uar2(k) = 0.0d0
         g_uar3(k) = 0.0d0
         g_uar4(k) = 0.0d0
         g_uar5(k) = 0.0d0
      end do
      
      if (DABS(Ud(1)-Ug(1)).le.0.005d0*Ug(1)) then
       coeff1 = 0.5d0
       coeff2 = (beta-2.0d0)/24.0d0
       coeff3 = -(beta-2.0d0)/48.0d0
       coeff4 = -(2.0d0*beta**3-3.0d0*beta**2-78.0d0*beta+152.0d0)
     &              / 5760.0d0
       coeff5 = (2.0d0*beta**3-3.0d0*beta**2-38.0d0*beta+72.0d0)
     &              /3840.0d0
       coeff6 = (16.0d0*beta**5-26.0d0*beta**4-1755.0d0*beta**3
     &            +2620.0d0*beta**2+22444.0d0*beta-41424.0d0)
     &          / 483840.0d0
       coeff7 = -(16.0d0*beta**5-26.0d0*beta**4-747.0d0*beta**3
     &             +1108.0d0*beta**2+7324.0d0*beta-13200.0d0)
     &          /165888.0d0
       eps                       = (Ud(1)-Ug(1))/Ug(1)
       g_uar1(1)                =  coeff1 + 2.0d0*coeff2*eps
     &                         + 3.0d0*coeff3*eps**2+4.0d0*coeff4*eps**3
     &                         + 5.0d0*coeff5*eps**4+coeff6*eps**5
     &                         + coeff7*eps**6
       g_uar1(1) = uar1/Ug(1)-Ud(1)/Ug(1)*g_uar1(1)
      else
       g_uar1(1) = -oobeta*oobeta1*uar1**(2.0d0-beta)
     &           *(beta*(Ug(1)**beta1)*difro-difroB)/difro**2
      endif


      g_uar2(1) = -pre*(Ug(2)+uar2)
c      g_uar2(1) = -squsrd*Ug(2)*pre
      g_uar2(2) = 2.0d0*pre

      g_uar3(1) = -pre*(Ug(3)+uar3)
c      g_uar3(1) = -squsrd*Ug(3)*pre
      g_uar3(3) = g_uar2(2)

      g_uar4(1) = -pre*(Ug(4)+uar4)
c      g_uar4(1) = -squsrd*Ug(4)*pre
      g_uar4(4) = g_uar2(2)

      g_uar5(1) = pre*(2.0d0*c2-enth-uar5)
      g_uar5(5) = 2.0d0*pre
      
c
c Computation of the dissipation term
c

c  the dif term uses (P,u,v,w,T)
      dif1 = -alpha*(Ug(1)**beta-Ud(1)**beta)
      dif2 = -Ug(2) + Ud(2)
      dif3 = -Ug(3) + Ud(3)
      dif4 = -Ug(4) + Ud(4)
      dif5 = -Ug(5) + Ud(5)


c  values needed by the different matrices
      VdotN = nx*uar2 + ny*uar3 + nz*uar4
      cr2 = alpha*beta*uar1**beta1
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
      

      tet1 = nz*uar3 - ny*uar4
      tet2 = nx*uar4 - nz*uar2
      tet3 = ny*uar2 - nx*uar3
    

c  eigenvalues - Preconditioned
      vp1 = VdotN
      vpa = 0.5d0*(1.0d0+b2)*VdotN
      vpb = 0.5d0*dsqrt(4.0d0*b2*cr2+(VdotN*(1.0d0-b2))**2)
      vp4 = vpa+vpb
      vp5 = vpa-vpb
      
      signvp1 = DSIGN(1.0d0, vp1)
      signvp4 = DSIGN(1.0d0, vp4)
      signvp5 = DSIGN(1.0d0, vp5)
      
      
c  useful values
      r = VdotN - vp4
      s = VdotN - vp5
      rs = -b2*cr2
      rps = VdotN*(1.0d0-b2)
      smr = vp4 - vp5
      r2 = r*r
      s2 = s*s
      rp = uar1*r/(s*smr)
      sp = -uar1*s/(r*smr)
      P1 = cr2/(uar1*rs*Cp)
      P2 = P1*uar1*rps
      
c differentiation of all the previous values
      do i = 1,5
        g_dif1(i) = 0.0d0
        g_dif2(i) = 0.0d0
        g_dif3(i) = 0.0d0
        g_dif4(i) = 0.0d0
        g_dif5(i) = 0.0d0
      
        g_VdotN(i) = nx*g_uar2(i)+ny*g_uar3(i)+nz*g_uar4(i)
        g_cr(i) = 0.5d0*alpha*beta*beta1*(uar1**(beta-2.0d0))
     &              *g_uar1(i)/cr
        g_tet1(i) = nz*g_uar3(i) - ny*g_uar4(i)
        g_tet2(i) = nx*g_uar4(i) - nz*g_uar2(i)
        g_tet3(i) = ny*g_uar2(i) - nx*g_uar3(i)

        g_vp1(i) = g_VdotN(i)
        g_vpa(i) = 0.5d0*(1.0d0+b2)*g_VdotN(i)
        g_vpb(i) = 0.25d0*(4.0d0*b2*cr*g_cr(i)
     &            + (1.0d0-b2)**2.0d0*VdotN*g_VdotN(i))/vpb
        g_vp4(i) = (g_vpa(i)+g_vpb(i))
        g_vp5(i) = (g_vpa(i)-g_vpb(i))

        g_r(i) = g_VdotN(i) - g_vp4(i)
        g_r2(i) = 2.0d0*r*g_r(i)
        g_s(i) = g_VdotN(i) - g_vp5(i)
        g_s2(i) = 2.0d0*s*g_s(i)
        g_rs(i) = -2.0d0*b2*cr*g_cr(i)
        g_rps(i) = g_VdotN(i)*(1.0d0-b2)
        g_smr(i) = g_s(i)-g_r(i)
        g_rp(i) = ((g_uar1(i)*r+uar1*g_r(i))*s*smr 
     &            - uar1*r*(s*g_smr(i)+g_s(i)*smr))
     &            / (s2*smr*smr)
        g_sp(i) = -((g_uar1(i)*s+uar1*g_s(i))*r*smr 
     &            - uar1*s*(r*g_smr(i)+g_r(i)*smr))
     &            / (r2*smr*smr)
c        g_P1(i) = (uar1*(rs*g_uar6(i)-uar6*g_rs(i))
c     &                     -2.0d0*uar6*rs*g_uar1(i))
c     &                    /(Cp*uar1**3*rs**2)
        g_P1(i) = (2.0d0*cr*g_cr(i)*uar1*rs
     &           -cr2*(g_uar1(i)*rs+uar1*g_rs(i))) 
     &           /(Cp*uar1**2*rs**2) 
        g_P2(i) = g_P1(i)*uar1*rps+P1*(g_uar1(i)*rps+uar1*g_rps(i))

      end do
      
      vp1 = signvp1*vp1
      vp4 = signvp4*vp4
      vp5 = signvp5*vp5
      do i = 1,5
        g_vp1(i) = signvp1*g_vp1(i)
        g_vp4(i) = signvp4*g_vp4(i)
        g_vp5(i) = signvp5*g_vp5(i)
      end do
      
      g_dif1(1) = -c2
      g_dif2(1) = Ug(2)/Ug(1)
      g_dif2(2) = -1.0d0/Ug(1)
      g_dif3(1) = Ug(3)/Ug(1)
      g_dif3(3) = -1.0d0/Ug(1)
      g_dif4(1) = Ug(4)/Ug(1)
      g_dif4(4) = -1.0d0/Ug(1)
      g_dif5(1) = (enth-c2-vitg2)/(Ug(1)*Cp)
      g_dif5(2) = Ug(2)/(Ug(1)*Cp)
      g_dif5(3) = Ug(3)/(Ug(1)*Cp)
      g_dif5(4) = Ug(4)/(Ug(1)*Cp)
      g_dif5(5) = -1.0d0/(Ug(1)*Cp)


      flur1 = vp1*(
     &          P1*nx        *dif1  
     &         +(P2*nx*nx   )*dif2
     &         +(P2*nx*ny-nz)*dif3
     &         +(P2*nx*nz+ny)*dif4
     &         +nx           *dif5 )
     
      do i=1,5
        g_flur1(i) = g_vp1(i)*(
     &          P1*nx        *dif1  
     &         +(P2*nx*nx   )*dif2
     &         +(P2*nx*ny-nz)*dif3
     &         +(P2*nx*nz+ny)*dif4
     &         +nx           *dif5 )
     
        g_flur1(i) = g_flur1(i)+vp1*(
     &          g_P1(i)*nx   *dif1  
     &         + P1*nx       *g_dif1(i)
     &         +(g_P2(i)*nx*nx)*dif2
     &         + P2*nx*nx    *g_dif2(i)
     &         +(g_P2(i)*nx*ny)*dif3
     &         + (P2*nx*ny-nz)*g_dif3(i)
     &         +(g_P2(i)*nx*nz)*dif4
     &         + (P2*nx*nz+ny)*g_dif4(i)
     &         +nx           *g_dif5(i) )
     
      end do

      flur2 = vp1*(
     &          P1*ny        *dif1  
     &         +(P2*ny*nx+nz)*dif2
     &         +(P2*ny*ny   )*dif3
     &         +(P2*ny*nz-nx)*dif4
     &         +ny           *dif5 )
     
      do i=1,5
        g_flur2(i) = g_vp1(i)*(
     &          P1*ny        *dif1  
     &         +(P2*ny*nx+nz)*dif2
     &         +(P2*ny*ny   )*dif3
     &         +(P2*ny*nz-nx)*dif4
     &         +ny           *dif5 )
     
        g_flur2(i) = g_flur2(i) + vp1*(
     &          g_P1(i)*ny   *dif1  
     &         +P1*ny        *g_dif1(i)
     &         +g_P2(i)*ny*nx*dif2
     &         +(P2*ny*nx+nz)*g_dif2(i)
     &         +(g_P2(i)*ny*ny   )*dif3
     &         +(P2*ny*ny   )*g_dif3(i)
     &         +g_P2(i)*ny*nz*dif4
     &         +(P2*ny*nz-nx)*g_dif4(i)
     &         +ny           *g_dif5(i) )
      end do
     
      flur3 = vp1*(
     &          P1*nz        *dif1  
     &         +(P2*nz*nx-ny)*dif2
     &         +(P2*nz*ny+nx)*dif3
     &         +(P2*nz*nz   )*dif4
     &         +nz           *dif5 )
     
      do i=1,5
        g_flur3(i) = g_vp1(i)*(
     &          P1*nz        *dif1  
     &         +(P2*nz*nx-ny)*dif2
     &         +(P2*nz*ny+nx)*dif3
     &         +(P2*nz*nz   )*dif4
     &         +nz           *dif5 )
     
        g_flur3(i) = g_flur3(i)+vp1*(
     &          g_P1(i)*nz   *dif1 
     &         + P1*nz       *g_dif1(i)
     &         +g_P2(i)*nz*nx*dif2
     &         +(P2*nz*nx-ny)*g_dif2(i)
     &         +g_P2(i)*nz*ny*dif3
     &         +(P2*nz*ny+nx)*g_dif3(i)
     &         +g_P2(i)*nz*nz*dif4
     &         +(P2*nz*nz   )*g_dif4(i)
     &         +nz           *g_dif5(i) )
      end do
    
      flur4 = vp4/cr2*(
     &        -dif1/(r*smr)+sp*(nx*dif2+ny*dif3+nz*dif4) )
     
      do i = 1,5
      g_flur4(i) = (g_vp4(i)*cr-2.0d0*vp4*g_cr(i))
     &              /(cr*cr2)*(
     &    -dif1/(r*smr)+sp*(nx*dif2+ny*dif3+nz*dif4) )
     
      g_flur4(i) = g_flur4(i) + vp4/cr2*(
     &       -(g_dif1(i)*r*smr-dif1*(g_r(i)*smr+r*g_smr(i)))/(r2*smr**2)
     &          +g_sp(i)*(nx*dif2+ny*dif3+nz*dif4)
     &          +sp*(nx*g_dif2(i)+ny*g_dif3(i)+nz*g_dif4(i)) )
      end do

      flur5 = vp5/cr2*(
     &         dif1/(s*smr)+rp*(nx*dif2+ny*dif3+nz*dif4) )

      do i = 1,5
      g_flur5(i) = (g_vp5(i)*cr-2.0d0*vp5*g_cr(i))
     &              /(cr*cr2)*(
     &    dif1/(s*smr)+rp*(nx*dif2+ny*dif3+nz*dif4) )
     
      g_flur5(i) = g_flur5(i) + vp5/cr2*(
     &      (g_dif1(i)*s*smr-dif1*(g_s(i)*smr+s*g_smr(i)))/(s2*smr**2)
     &          +g_rp(i)*(nx*dif2+ny*dif3+nz*dif4)
     &          +rp*(nx*g_dif2(i)+ny*g_dif3(i)+nz*g_dif4(i)) )
      end do


c Final operation

      do i=1,5
      H(1,i) = H(1,i) - gamma*oob2*(
     &         g_r2(i)*flur4+g_s2(i)*flur5
     &          + r2*g_flur4(i) + s2*g_flur5(i) )
     
      H(2,i) = H(2,i) - gamma*(
     &          uar1 *(nz*g_flur2(i)- ny*g_flur3(i)) 
     &      + (oob2*uar2*r2-r*cr2*nx)*g_flur4(i)
     &      + (oob2*uar2*s2-s*cr2*nx)*g_flur5(i) )
    
      H(2,i) = H(2,i) - gamma*(
     &      g_uar1(i)*(nz*flur2-ny*flur3)
     &      + (oob2*(g_uar2(i)*r2+uar2*g_r2(i))
     &          -(g_r(i)*cr2+r*2.0d0*cr*g_cr(i))*nx)*flur4
     &      + (oob2*(g_uar2(i)*s2+uar2*g_s2(i))
     &          -(g_s(i)*cr2+s*2.0d0*cr*g_cr(i))*nx)*flur5 )
    
      H(3,i) = H(3,i) - gamma*(
     &          uar1*(nx*g_flur3(i)-nz*g_flur1(i))
     &      + (oob2*uar3*r2-r*cr2*ny)*g_flur4(i)
     &      + (oob2*uar3*s2-s*cr2*ny)*g_flur5(i) )
    
      H(3,i) = H(3,i) - gamma*(
     &      g_uar1(i)*(nx*flur3-nz*flur1)
     &      + (oob2*(g_uar3(i)*r2+uar3*g_r2(i))
     &          -(g_r(i)*cr2+r*2.0d0*cr*g_cr(i))*ny)*flur4
     &      + (oob2*(g_uar3(i)*s2+uar3*g_s2(i))
     &          -(g_s(i)*cr2+s*2.0d0*cr*g_cr(i))*ny)*flur5 )

      H(4,i) = H(4,i) - gamma*(
     &          uar1*(ny*g_flur1(i)-nx*g_flur2(i))
     &      + (oob2*uar4*r2-r*cr2*nz)*g_flur4(i)
     &      + (oob2*uar4*s2-s*cr2*nz)*g_flur5(i) )
     
      H(4,i) = H(4,i) - gamma*(
     &      g_uar1(i)*(ny*flur1-nx*flur2)
     &      + (oob2*(g_uar4(i)*r2+uar4*g_r2(i))
     &          -(g_r(i)*cr2+r*2.0d0*cr*g_cr(i))*nz)*flur4
     &      + (oob2*(g_uar4(i)*s2+uar4*g_s2(i))
     &          -(g_s(i)*cr2+s*2.0d0*cr*g_cr(i))*nz)*flur5 )
     
      H(5,i) = H(5,i) - gamma*(
     &     uar1*Cp*(g_flur1(i)*nx + g_flur2(i)*ny + g_flur3(i)*nz)
     &     -uar1*(tet1*g_flur1(i)+tet2*g_flur2(i)+tet3*g_flur3(i))
     &     + (oob2*(uar5-cr2)*r2+cr2*(cr2-VdotN*r)
     &               )*g_flur4(i)
     &     + (oob2*(uar5-cr2)*s2+cr2*(cr2-VdotN*s)
     &               )*g_flur5(i) )
     
      H(5,i) = H(5,i) - gamma*(
     &         g_uar1(i)*Cp*(flur1*nx + flur2*ny + flur3*nz)
     &      - g_uar1(i)*(tet1*flur1 + tet2*flur2 + tet3*flur3)
     &      - uar1*(g_tet1(i)*flur1 + g_tet2(i)*flur2 + g_tet3(i)*flur3)
     &      + (oob2*((g_uar5(i)-2.0d0*cr*g_cr(i))*r2
     &                +(uar5-cr2)*g_r2(i))
     &         +2.0d0*cr*g_cr(i)*(cr2-VdotN*r)
     &         + cr2*(2.0d0*cr*g_cr(i)-VdotN*g_r(i)-g_VdotN(i)*r))
     &                               *flur4
     &      + (oob2*((g_uar5(i)-2.0d0*cr*g_cr(i))*s2
     &                +(uar5-cr2)*g_s2(i))
     &         +2.0d0*cr*g_cr(i)*(cr2-VdotN*s)
     &         + cr2*(2.0d0*cr*g_cr(i)-VdotN*g_s(i)-g_VdotN(i)*s))
     &                               *flur5 )
     
      end do
    
    



c
c     NORMALIZATION of the FLUX
c
      dim = 5 + type
      do i = 1,5
         do k = 1,5
            jac(dim*(i-1) + k) = H(i,k)*0.5d0*rnorm
         enddo
         do k = 6,dim
            jac(dim*(i-1) + k) = 0.0d0
         enddo
      enddo
      do i = 6,dim
         do k = 1,dim
            jac(dim*(i-1) + k) = 0.0d0
         enddo
      enddo


      END 
