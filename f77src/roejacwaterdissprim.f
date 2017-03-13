      SUBROUTINE ROEJACWATERDISSPRIM(type,gamma,Cv,Pref,alpha,beta,
     &                         enormal,evitno,Ug,Ud,jac)
c
c--------------------------------------------------------------------
c  This routine computes the derivative of the Roe flux w.r.t
c  the 1st variable taken at the vectors Ug, Ud.
c  normal is the normal of the boundary concerned by the flux.
c  jac stores the resulting jacobian.
c--------------------------------------------------------------------


      implicit none

      REAL*8 Ug(*),Ud(*),enormal(3),evitno,jac(*)
      REAL*8 gamma
      REAL*8 Cv,Pref,alpha,beta
      INTEGER type

      REAL*8 normal(3),rnorm,invnorm,vitno
      REAL*8 VdotN,qir,P,e,enth,c2,c
      REAL*8 H(5,5)
      REAL*8 pre,difro,difroB,beta1,oobeta,oobeta1,ener1,ener2
      REAL*8 squsr1,squsr2,usro,h1,h2,uar1,uar2,uar3,uar4,uar5,uar6
      REAL*8 cr, cr2, u2mh
      REAL*8 g_uar1(5),g_uar2(5),g_uar3(5),
     &       g_uar4(5),g_uar5(5),g_uar6(5),
     &       g_uar6ouar1cr2(5)
      REAL*8 g_tet1(5), g_tet2(5), g_tet3(5)
      REAL*8 dif1,dif2,dif3,dif4,dif5
      REAL*8 vp1,vp4,vp5,signvp1,signvp4,signvp5
      REAL*8 g_VdotN(5),g_qir(5),g_cr(5), g_u2mh(5)
      REAL*8 g_dif1(5),g_dif2(5),g_dif3(5),g_dif4(5),g_dif5(5)
      REAL*8 g_vp1(5),g_vp4(5),g_vp5(5)
      REAL*8 flur1,flur2,flur3,flur4,flur5
      REAL*8 g_flur1(5),g_flur2(5),g_flur3(5),g_flur4(5),g_flur5(5)
      INTEGER i,k,dim
      REAL*8 tet1, tet2, tet3
      REAL*8 coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, eps
      REAL*8 coeff7



c
c     INITIALISATION
c

      beta1 = beta -1.0d0
      oobeta = 1.0d0/beta
      oobeta1 = 1.0d0/beta1

      rnorm = DSQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) +
     &              enormal(3)*enormal(3))
      invnorm = 1.0d0 / rnorm

      normal(1) = enormal(1) * invnorm
      normal(2) = enormal(2) * invnorm
      normal(3) = enormal(3) * invnorm

      vitno = evitno * invnorm

c
c     COMPUTATION OF CENTERED TERM
c

      VdotN = Ug(2)*normal(1) + Ug(3)*normal(2) + Ug(4)*normal(3)
      qir = 0.5d0*(Ug(2)*Ug(2) + Ug(3)*Ug(3) + Ug(4)*Ug(4))
      P = Pref + alpha*(Ug(1)**beta)
      e = Cv*Ug(5) + qir
      enth = e + P/Ug(1)
      c2 = alpha*beta*(Ug(1)**beta1)
      c = DSQRT(c2)

      H(1,1) = VdotN - vitno
      H(1,2) = Ug(1)*normal(1)
      H(1,3) = Ug(1)*normal(2)
      H(1,4) = Ug(1)*normal(3)
      H(1,5) = 0.0d0

      H(2,1) = H(1,1)*Ug(2) + c2*normal(1)
      H(2,2) = H(1,2)*Ug(2) + Ug(1)*H(1,1)
      H(2,3) = H(1,3)*Ug(2)
      H(2,4) = H(1,4)*Ug(2)
      H(2,5) = 0.0d0

      H(3,1) = H(1,1)*Ug(3) + c2*normal(2)
      H(3,2) = H(1,2)*Ug(3)
      H(3,3) = H(1,3)*Ug(3) + Ug(1)*H(1,1)
      H(3,4) = H(1,4)*Ug(3)
      H(3,5) = 0.0d0

      H(4,1) = H(1,1)*Ug(4) + c2*normal(3)
      H(4,2) = H(1,2)*Ug(4)
      H(4,3) = H(1,3)*Ug(4)
      H(4,4) = H(1,4)*Ug(4) + Ug(1)*H(1,1)
      H(4,5) = 0.0d0

      H(5,1) = (e+c2)*H(1,1)
      H(5,2) = H(1,2)*enth + Ug(1)*Ug(2)*H(1,1)
      H(5,3) = H(1,3)*enth + Ug(1)*Ug(3)*H(1,1)
      H(5,4) = H(1,4)*enth + Ug(1)*Ug(4)*H(1,1)
      H(5,5) = Ug(1)*Cv*H(1,1)


c
c     COMPUTATION OF THE ROE-AVERAGED STATE
c

      squsr1 = DSQRT(Ug(1))
      squsr2 = DSQRT(Ud(1))
      usro = 1.0d0/(squsr1+squsr2)

      ener1 = e
      ener2 = Cv*Ud(5) + 0.5d0*(Ud(2)*Ud(2)+Ud(3)*Ud(3)+Ud(4)*Ud(4))

      h1 = enth
      h2 = ener2 + (Pref+alpha*(Ud(1)**beta))/Ud(1)

      difro = Ud(1)-Ug(1)
      difroB = Ud(1)**beta - Ug(1)**beta

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
       uar1                      = ( difroB /
     &                            ( beta * (difro) ) ) **
     &                            ( 1.0d0 / beta1 )
      endif

      uar2 = (squsr1*Ug(2) + squsr2*Ud(2))*usro
      uar3 = (squsr1*Ug(3) + squsr2*Ud(3))*usro
      uar4 = (squsr1*Ug(4) + squsr2*Ud(4))*usro
      uar5 = (squsr1*h1 + squsr2*h2)*usro
      uar6 = Pref +alpha*uar1**beta  

      pre = 0.5d0*usro/squsr1

      do k=1,5
         g_uar1(k) = 0.0d0
         g_uar2(k) = 0.0d0
         g_uar3(k) = 0.0d0
         g_uar4(k) = 0.0d0
         g_uar5(k) = 0.0d0
	 g_uar6(k) = 0.0d0
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
     &          /165888
       eps                       = (Ud(1)-Ug(1))/Ug(1)
       g_uar1(1)                =  coeff1 + 2.0d0*coeff2*eps
     &                         + 3.0d0*coeff3*eps**2+4.0d0*coeff4*eps**3
     &                         + 5.0d0*coeff5*eps**4+coeff6*eps**5
     &                         + coeff7*eps**6
       g_uar1(1) = uar1/Ug(1)-Ud(1)/Ug(1)*g_uar1(1) 
      else
        g_uar1(1) = oobeta*oobeta1*uar1**(2.0d0-beta)
     &           *(beta*(Ug(1)**beta1)*difro-difroB)/difro**2
      endif

      g_uar2(1) = pre*(Ug(2)-uar2)
      g_uar2(2) = 2.0d0*Ug(1)*pre

      g_uar3(1) = pre*(Ug(3)-uar3)
      g_uar3(3) = g_uar2(2)

      g_uar4(1) = pre*(Ug(4)-uar4)
      g_uar4(4) = g_uar2(2)

      g_uar5(1) = pre*(h1-uar5) + 2.0d0*pre*(c2- P/Ug(1))
      g_uar5(2) = g_uar2(2)*Ug(2)
      g_uar5(3) = g_uar2(2)*Ug(3)
      g_uar5(4) = g_uar2(2)*Ug(4)
      g_uar5(5) = g_uar2(2)*Cv
      
      g_uar6(1) = alpha*beta*uar1**beta1*g_uar1(1)

c
c     COMPUTATION OF THE DISSIPATION TERM
c

      VdotN = normal(1)*uar2 + normal(2)*uar3 + normal(3)*uar4
      qir   = 0.5d0*(uar2*uar2 + uar3*uar3 + uar4*uar4)
      u2mh  = 2.0d0*qir - uar5

      tet1  = normal(3)*uar3 - normal(2)*uar4
      tet2  = normal(1)*uar4 - normal(3)*uar2
      tet3  = normal(2)*uar2 - normal(1)*uar3


      cr2 = alpha*beta*(uar1**beta1)
      cr  = DSQRT(cr2)

c  the dif term uses (P,u,v,w,T)
      dif1 = -alpha*(Ug(1)**beta-Ud(1)**beta)
      dif2 = -Ug(2) + Ud(2)
      dif3 = -Ug(3) + Ud(3)
      dif4 = -Ug(4) + Ud(4)
      dif5 = -Ug(5) + Ud(5)
      
      vp1 = gamma*(VdotN-vitno)
      vp4 = gamma*(VdotN-vitno+cr)
      vp5 = gamma*(VdotN-vitno-cr)

      signvp1 = DSIGN(1.0d0, vp1)
      signvp4 = DSIGN(1.0d0, vp4)
      signvp5 = DSIGN(1.0d0, vp5)

      do i=1,5
         g_VdotN(i) = normal(1)*g_uar2(i) + normal(2)*g_uar3(i) +
     &                normal(3)*g_uar4(i)

         g_qir(i) = uar2*g_uar2(i)+uar3*g_uar3(i)+uar4*g_uar4(i)

         g_u2mh(i) = 2.0d0*g_qir(i) - g_uar5(i)

         g_tet1(i) = normal(3)*g_uar3(i)-normal(2)*g_uar4(i)
         g_tet2(i) = normal(1)*g_uar4(i)-normal(3)*g_uar2(i)
         g_tet3(i) = normal(2)*g_uar2(i)-normal(1)*g_uar3(i)

         g_cr(i) = 0.5d0*alpha*beta*beta1*(uar1**(beta-2.0d0))
     &              *g_uar1(i)/cr
     
         g_uar6ouar1cr2(i) = (g_uar6(i)*uar1*cr2-
     &                  uar6*(g_uar1(i)*cr2+2.0d0*cr*g_cr(i)*uar1))
     &                  /(uar1**2 * cr2**2)

c     initialization of g_dif#
         g_dif1(i) = 0.0d0
         g_dif2(i) = 0.0d0
         g_dif3(i) = 0.0d0
         g_dif4(i) = 0.0d0
         g_dif5(i) = 0.0d0

         g_vp1(i) = signvp1*gamma*g_VdotN(i)
         g_vp4(i) = signvp4*gamma*(g_VdotN(i)+g_cr(i))
         g_vp5(i) = signvp5*gamma*(g_VdotN(i)-g_cr(i))

      end do

      vp1 = signvp1*vp1
      vp4 = signvp4*vp4
      vp5 = signvp5*vp5

      g_dif1(1) = -alpha*beta*Ug(1)**beta1
      g_dif2(2) = -1.0d0
      g_dif3(3) = -1.0d0
      g_dif4(4) = -1.0d0
      g_dif5(5) = -1.0d0


      flur1 = vp1*(
     &         -uar6*normal(1)/(uar1*cr2) * dif1
     &         -uar1*normal(3) * dif3
     &         +uar1*normal(2) * dif4
     &         +uar1*Cv*normal(1) * dif5 )


      do i=1,5

         g_flur1(i) = g_vp1(i)*(
     &         -uar6*normal(1)/(uar1*cr2) * dif1
     &         -uar1*normal(3) * dif3
     &         +uar1*normal(2) * dif4
     &         +uar1*Cv*normal(1) * dif5 )

         g_flur1(i) = g_flur1(i) + vp1*(
     &         -uar6*normal(1)/(uar1*cr2) * g_dif1(i)
     &         -uar1*normal(3) * g_dif3(i)
     &         +uar1*normal(2) * g_dif4(i)
     &         +uar1*Cv*normal(1) * g_dif5(i) )
     
         g_flur1(i) = g_flur1(i) + vp1*(
     &         -g_uar6ouar1cr2(i)*normal(1) * dif1
     &         -g_uar1(i)*normal(3) * dif3
     &         +g_uar1(i)*normal(2) * dif4
     &         +g_uar1(i)*Cv*normal(1) * dif5 )

      end do


      flur2 = vp1*(
     &         -uar6*normal(2)/(uar1*cr2) * dif1
     &         +uar1*normal(3) * dif2
     &         -uar1*normal(1) * dif4
     &         +uar1*Cv*normal(2) * dif5 )

      

      do i=1,5

         g_flur2(i) = g_vp1(i)*(
     &         -uar6*normal(2)/(uar1*cr2) * dif1
     &         +uar1*normal(3) * dif2
     &         -uar1*normal(1) * dif4
     &         +uar1*Cv*normal(2) * dif5 )
     
         g_flur2(i) = g_flur2(i) + vp1*(
     &         -uar6*normal(2)/(uar1*cr2) * g_dif1(i)
     &         +uar1*normal(3) * g_dif2(i)
     &         -uar1*normal(1) * g_dif4(i)
     &         +uar1*Cv*normal(2) * g_dif5(i) )
     
         g_flur2(i) = g_flur2(i) + vp1*(
     &         -g_uar6ouar1cr2(i)*normal(2) * dif1
     &         +g_uar1(i)*normal(3) * dif2
     &         -g_uar1(i)*normal(1) * dif4
     &         +g_uar1(i)*Cv*normal(2) * dif5 )

      end do

      


      flur3 = vp1*(
     &         -uar6*normal(3)/(uar1*cr2) * dif1
     &         -uar1*normal(2) * dif2
     &         +uar1*normal(1) * dif3
     &         +uar1*Cv*normal(3) * dif5 )

      do i=1,5

         g_flur3(i) = g_vp1(i)*(
     &         -uar6*normal(3)/(uar1*cr2) * dif1
     &         -uar1*normal(2) * dif2
     &         +uar1*normal(1) * dif3
     &         +uar1*Cv*normal(3) * dif5 )
     
         g_flur3(i) = g_flur3(i) + vp1*(
     &         -uar6*normal(3)/(uar1*cr2) * g_dif1(i)
     &         -uar1*normal(2) * g_dif2(i)
     &         +uar1*normal(1) * g_dif3(i)
     &         +uar1*Cv*normal(3) * g_dif5(i) )
     
         g_flur3(i) = g_flur3(i) + vp1*(
     &         -g_uar6ouar1cr2(i)*normal(3) * dif1
     &         -g_uar1(i)*normal(2) * dif2
     &         +g_uar1(i)*normal(1) * dif3
     &         +g_uar1(i)*Cv*normal(3) * dif5 )
     



      end do


       flur4 = 0.5d0*vp4*(
     &          1.0d0/cr2 * dif1
     &         +uar1*(normal(1)*dif2+normal(2)*dif3+normal(3)*dif4)/cr )


      do i=1,5

         g_flur4(i) = 0.5d0*g_vp4(i)*(
     &          1.0d0/cr2 * dif1
     &         +uar1*(normal(1)*dif2+normal(2)*dif3+normal(3)*dif4)/cr )
     
         g_flur4(i) = g_flur4(i) + 0.5d0*vp4*(
     &          1.0d0/cr2 * g_dif1(i)
     &   +uar1/cr*
     &   (normal(1)*g_dif2(i)+normal(2)*g_dif3(i)+normal(3)*g_dif4(i)) )
     
         g_flur4(i) = g_flur4(i) + 0.5d0*vp4*(
     &          -2.0d0*g_cr(i)/(cr*cr2) * dif1
     &   +(g_uar1(i)*cr-uar1*g_cr(i))/cr2
     &       *(normal(1)*dif2+normal(2)*dif3+normal(3)*dif4) )
     
      end do


      flur5 = 0.5d0*vp5*(
     &          1.0d0/cr2 * dif1
     &         -uar1*(normal(1)*dif2+normal(2)*dif3+normal(3)*dif4)/cr )

      do i=1,5

         g_flur5(i) = 0.5d0*g_vp5(i)*(
     &          1.0d0/cr2 * dif1
     &         -uar1*(normal(1)*dif2+normal(2)*dif3+normal(3)*dif4)/cr )
     
         g_flur5(i) = g_flur5(i) + 0.5d0*vp5*(
     &          1.0d0/cr2 * g_dif1(i)
     &   -uar1/cr*
     &   (normal(1)*g_dif2(i)+normal(2)*g_dif3(i)+normal(3)*g_dif4(i)) )
     
         g_flur5(i) = g_flur5(i) + 0.5d0*vp5*(
     &          -2.0d0*g_cr(i)/(cr*cr2) * dif1
     &   -(g_uar1(i)*cr-uar1*g_cr(i))/cr2
     &       *(normal(1)*dif2+normal(2)*dif3+normal(3)*dif4) )
      end do
 

c
c     FINAL OPERATION
c

      do i=1,5

         H(1,i) = H(1,i) - (g_flur4(i)+g_flur5(i))

         H(2,i) = H(2,i) - (-normal(2)*g_flur3(i)
     &                      +normal(3)*g_flur2(i)
     &                      +uar2 * (g_flur4(i)+g_flur5(i))
     &                      +cr*normal(1) * (g_flur4(i)-g_flur5(i)))

         H(2,i) = H(2,i) - ( g_uar2(i) * (flur4+flur5)
     &                      +g_cr(i)*normal(1) * (flur4-flur5))

         H(3,i) = H(3,i) - ( normal(1)*g_flur3(i)
     &                      -normal(3)*g_flur1(i)
     &                      +uar3 * (g_flur4(i)+g_flur5(i))
     &                      +cr*normal(2) * (g_flur4(i)-g_flur5(i)))

         H(3,i) = H(3,i) - ( g_uar3(i) * (flur4+flur5)
     &                      +g_cr(i)*normal(2) * (flur4-flur5))

         H(4,i) = H(4,i) - (-normal(1)*g_flur2(i)
     &                      +normal(2)*g_flur1(i)
     &                      +uar4 * (g_flur4(i)+g_flur5(i))
     &                      +cr*normal(3) * (g_flur4(i)-g_flur5(i)))

         H(4,i) = H(4,i) - ( g_uar4(i) * (flur4+flur5)
     &                      +g_cr(i)*normal(3) * (flur4-flur5))

         H(5,i) = H(5,i) - ( (normal(1)-tet1)*g_flur1(i)
     &                      +(normal(2)-tet2)*g_flur2(i)
     &                      +(normal(3)-tet3)*g_flur3(i)
     &                      + uar5*(g_flur4(i)+g_flur5(i))
     &                      + cr*VdotN*(g_flur4(i)-g_flur5(i)))

        H(5,i) = H(5,i) - ( (-g_tet1(i))*flur1
     &                      +(-g_tet2(i))*flur2
     &                      +(-g_tet3(i))*flur3
     &                      + g_uar5(i)*(flur4+flur5)
     &                   +(g_cr(i)*VdotN+cr*g_VdotN(i))*(flur4-flur5))


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



      end
     
