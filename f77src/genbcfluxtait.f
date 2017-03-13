      SUBROUTINE GENBCFLUXTAIT(type,cp,ps,a,b,enormal,evitno,
     &     V,Uinf,phi)
c-----------------------------------------------------------------------
c This routine computes the Flux at the boundary using
c left and right eigenvectors. See Ghidaglia
c V and Vinf are values of primitive variables at node and infinity
c Vinf is obtained from the values of the conservative variables at
c   infinity Uinf, given as input
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
c gamma is the dissipation coefficient
c EOS is
c   Pressure = ps + a*Density^b
c   h = cp*T
c   where h is the internal enthalpy
c         ps, a, b and cp are constants
c 
c No Low Mach preconditioner is applied so far
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER type, i
      REAL*8 V(*), Uinf(*), normal(3), enormal(3), evitno, phi(*)
      REAL*8 vitno, updir
      REAL*8 Vinf(10)
      REAL*8 VdotN , VdotNinf, rnorm, invnorm
      REAL*8 flur1 , flur2 , flur3 , flur4 , flur5
      REAL*8 r1(5) , r2(5) , r3(5) , r4(5) , r5(5)
      REAL*8 temp(5), vp(5)
      REAL*8 enerinf , ener, enthinf, u2mhinf
      REAL*8 presinf , pres
      REAL*8 cp,oocp,ps,a,b,bm1
      REAL*8 tet1 , tet2 , tet3
      REAL*8 vit2, vitinf2
      REAL*8 flux(5), fluxinf(5)
      REAL*8 c, c2, q

c
c Initialisation
c
      bm1 = b-1.0d0
      oocp = 1.0d0/cp

      rnorm = DSQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) + 
     &             enormal(3)*enormal(3))
      invnorm = 1.0d0 / rnorm
      normal(1) = enormal(1) * invnorm
      normal(2) = enormal(2) * invnorm
      normal(3) = enormal(3) * invnorm
      vitno = evitno * invnorm

      presinf  = ps+a*Uinf(1)**b

      Vinf(1) = Uinf(1)
      Vinf(2) = Uinf(2)/Uinf(1)
      Vinf(3) = Uinf(3)/Uinf(1)
      Vinf(4) = Uinf(4)/Uinf(1)
      Vinf(5) = oocp*(Uinf(5)/Uinf(1)-0.5d0*
     &      (Vinf(2)**2+Vinf(3)**2+Vinf(4)**2) + presinf/Uinf(1))

c A few useful values for the computation
      VdotNinf = Vinf(2)*normal(1)+Vinf(3)*normal(2)+Vinf(4)*normal(3)
      VdotN    = V(2)*normal(1)+V(3)*normal(2)+V(4)*normal(3)
      vitinf2  = Vinf(2)**2+Vinf(3)**2+Vinf(4)**2
      vit2     = V(2)**2+V(3)**2+V(4)**2
      enerinf  = Uinf(5)
      pres     = ps+a*V(1)**b
      ener     = 0.5d0*V(1)*vit2 + V(1)*cp*V(5) - pres
      enthinf  = (enerinf+presinf)/Vinf(1)
      u2mhinf  = vitinf2-enthinf

c Computation of the physical fluxes
      flux(1) = V(1)*(VdotN-vitno)
      flux(2) = flux(1)*V(2)+pres*normal(1)
      flux(3) = flux(1)*V(3)+pres*normal(2)
      flux(4) = flux(1)*V(4)+pres*normal(3)
      flux(5) = (ener+pres)*(VdotN-vitno) + pres*vitno

      fluxinf(1) = Vinf(1)*(VdotNinf-vitno)
      fluxinf(2) = fluxinf(1)*Vinf(2)+presinf*normal(1)
      fluxinf(3) = fluxinf(1)*Vinf(3)+presinf*normal(2)
      fluxinf(4) = fluxinf(1)*Vinf(4)+presinf*normal(3)
      fluxinf(5) = (enerinf+presinf)*(VdotNinf-vitno)
     &           + presinf*vitno

c Again, a few useful values
      tet1 = normal(3)*Vinf(3) - normal(2)*Vinf(4)
      tet2 = normal(1)*Vinf(4) - normal(3)*Vinf(2)
      tet3 = normal(2)*Vinf(2) - normal(1)*Vinf(3)

      c2  = a*b*Vinf(1)**bm1
      c   = DSQRT(c2)

c Computation of righteigenvectors using values at infinity
c formulae found in Arthur's notes (Research/latex/notes/)
      r1(1) = 0.0d0
      r1(2) = 0.0d0
      r1(3) = -normal(3)
      r1(4) =  normal(2)
      r1(5) =  normal(1) - tet1

      r2(1) = 0.0d0
      r2(2) =  normal(3)
      r2(3) = 0.0d0
      r2(4) = -normal(1)
      r2(5) =  normal(2) - tet2

      r3(1) = 0.0d0
      r3(2) = -normal(2)
      r3(3) = normal(1)
      r3(4) = 0.0d0
      r3(5) =  normal(3) - tet3

      r4(1) = 1.0d0
      r4(2) = Vinf(2)+c*normal(1)
      r4(3) = Vinf(3)+c*normal(2)
      r4(4) = Vinf(4)+c*normal(3)
      r4(5) = enthinf+c*VdotNinf

      r5(1) = 1.0d0
      r5(2) = Vinf(2)-c*normal(1)
      r5(3) = Vinf(3)-c*normal(2)
      r5(4) = Vinf(4)-c*normal(3)
      r5(5) = enthinf-c*VdotNinf

c Computation of eigenvalues using value at infinity
      vp(1) = VdotNinf-vitno
      vp(2) = vp(1)
      vp(3) = vp(1)
      vp(4) = VdotNinf+c-vitno
      vp(5) = VdotNinf-c-vitno

c scalar values of lefteigenvector*physical flux
c where physical flux depends on sign of the eigenvalue
c associated with the lefteigenvector
c check that if temp = r_i then i get flur_j*temp = delta_ij

c first three eigenvalues: u.n
c exiting domain
      if ( vp(1) .ge. 0.0d0 ) then
        do i=1,5
          temp(i) = flux(i)
        enddo
c incoming domain
      else
        do i=1,5
          temp(i) = fluxinf(i)
        enddo
      endif


      flur1  = 
     &        (u2mhinf*normal(1)+tet1)*temp(1) -
     &         Vinf(2)*normal(1)*temp(2) -
     &        (Vinf(3)*normal(1)+normal(3))*temp(3) -
     &        (Vinf(4)*normal(1)-normal(2))*temp(4) +
     &        normal(1)*temp(5) 

      flur2  =
     &        (u2mhinf*normal(2)+tet2)*temp(1) -
     &        (Vinf(2)*normal(2)-normal(3))*temp(2) -
     &         Vinf(3)*normal(2)*temp(3) -
     &        (Vinf(4)*normal(2)+normal(1))*temp(4) +
     &        normal(2)*temp(5) 
     
      flur3  =
     &        (u2mhinf*normal(3)+tet3)*temp(1) -
     &        (Vinf(2)*normal(3)+normal(2))*temp(2) -
     &        (Vinf(3)*normal(3)-normal(1))*temp(3) -
     &         Vinf(4)*normal(3)*temp(4) +
     &        normal(3)*temp(5) 
       
c fourth eigenvalue: u.n + c
c exiting domain
      if ( vp(4) .ge. 0.0d0 ) then
        do i=1,5
          temp(i) = flux(i)
        enddo
c incoming domain
      else
        do i=1,5
          temp(i) = fluxinf(i)
        enddo
      endif
      

      flur4   =
     &         0.5d0*(1.0d0-VdotNinf/c)*temp(1) +
     &         0.5d0*normal(1)/c * temp(2) +
     &         0.5d0*normal(2)/c * temp(3) +
     &         0.5d0*normal(3)/c * temp(4) 

c fifth eigenvalue
c exiting domain
      if ( vp(5) .ge. 0.0d0 ) then
        do i=1,5
          temp(i) = flux(i)
        enddo
c incoming domain
      else
        do i=1,5
          temp(i) = fluxinf(i)
        enddo
      endif

      flur5   =
     &         0.5d0*(1.0d0+VdotNinf/c)*temp(1)- 
     &         0.5d0*normal(1)/c * temp(2) -
     &         0.5d0*normal(2)/c * temp(3) -
     &         0.5d0*normal(3)/c * temp(4) 


c
c Final phi 
c

      DO i=1,5
        phi(i) =
     &      flur1*r1(i)+flur2*r2(i)+flur3*r3(i)+flur4*r4(i)+flur5*r5(i)
      ENDDO

      phi(1) = phi(1)*rnorm
      phi(2) = phi(2)*rnorm
      phi(3) = phi(3)*rnorm
      phi(4) = phi(4)*rnorm
      phi(5) = phi(5)*rnorm

c
c For one and two equation turbulence models
c

      if (type.eq.1) then
         Vinf(6) = Uinf(6) / Vinf(1)
         updir = 0.5d0 + dsign(0.5d0, phi(1))
         phi(6) = phi(1) * (updir * V(6) + (1.0d0 - updir) * Vinf(6))
      else if (type.eq.2) then
         Vinf(6) = Uinf(6) / Vinf(1)
         Vinf(7) = Uinf(7) / Vinf(1)
         updir = 0.5d0 + dsign(0.5d0, phi(1))
         phi(6) = phi(1) * (updir * V(6) + (1.0d0 - updir) * Vinf(6))
         phi(7) = phi(1) * (updir * V(7) + (1.0d0 - updir) * Vinf(7))
      endif
      
      END
