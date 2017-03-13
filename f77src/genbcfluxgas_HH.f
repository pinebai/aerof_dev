      SUBROUTINE GENBCFLUXGAS_HH(type,gam,pstiff,enormal,evitno,
     &     U,Uinf,phi,xface,snew,sold,dt)
c-----------------------------------------------------------------------
c This routine computes the Flux at the boundary using
c left and right eigenvectors. See Ghidaglia & Hagstorm,Hariharan
c U and Vinf are values of primitive variables at node and infinity
c normal is the normal of the boundary concerned by the flux.
c phi stores the resulting flux.
c gamma is the dissipation coefficient
c gam is the ratio of cp/cv
c No Low Mach preconditioner is applied so far

c Last Modified 07/18/2012 by karthik <dkarthik@stanford.edu>
c Notes:
c      1. Eigensystem (based on interior)
c      2. Assumes outflow BC (Hagstorm, Hariharan)
c      3. Have to check if grid velocity will play a role
c      4. Assumes isentropicity
c The following are new inputs:
c xface : coordinates of face center (required to compute radius)
c snew  : New value of incoming Riemann invariant 
c sold  : Old value of incoming Riemann invariant
c dt    : Physical time step size.
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER type, i
      REAL*8 U(*), Uinf(*), normal(3), enormal(3), evitno, phi(*)
      REAL*8 vitno, updir
      REAL*8 Vinf(10)
      REAL*8 VdotN , VdotNinf, rnorm, invnorm
      REAL*8 flur1 , flur2 , flur3 , flur4 , flur5
      REAL*8 r1(5) , r2(5) , r3(5) , r4(5) , r5(5)
      REAL*8 temp(5), vp(5)
      REAL*8 enerinf , ener
      REAL*8 tet1 , tet2 , tet3
      REAL*8 gam , gam1, pstiff
      REAL*8 vit2, vitinf2
      REAL*8 flux(5), fluxinf(5)
      REAL*8 c, c2, q
      REAL*8 Rplus, Entropy, radius, cinf, Rinf, vbar, cbar
      REAL*8 xface(3), snew, sold, dt, rnx, rny, rnz
      REAL*8 kt,ku,dSdt1,dSdt2,dSdt3

c
c Initialisation
c
      gam1 = gam - 1.d0      

      rnorm = DSQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) + 
     &             enormal(3)*enormal(3))
      invnorm = 1.0d0 / rnorm
      normal(1) = enormal(1) * invnorm
      normal(2) = enormal(2) * invnorm
      normal(3) = enormal(3) * invnorm
      vitno = evitno * invnorm

      Vinf(1) = Uinf(1)
      Vinf(2) = Uinf(2)/Uinf(1)
      Vinf(3) = Uinf(3)/Uinf(1)
      Vinf(4) = Uinf(4)/Uinf(1)
      Vinf(5) = gam1*(Uinf(5)-
     &              0.5d0*(Uinf(2)**2+Uinf(3)**2+Uinf(4)**2)/Uinf(1))
     &              -gam*pstiff


c A few useful values for the computation
      VdotNinf = Vinf(2)*normal(1)+Vinf(3)*normal(2)+Vinf(4)*normal(3)
      VdotN    = U(2)*normal(1)+U(3)*normal(2)+U(4)*normal(3)
      vitinf2  = Vinf(2)**2+Vinf(3)**2+Vinf(4)**2
      vit2     = U(2)**2+U(3)**2+U(4)**2
      enerinf  = 0.5d0*Vinf(1)*vitinf2 + (Vinf(5)+gam*pstiff)/gam1
      ener     = 0.5d0*U(1)*vit2 + (U(5)+gam*pstiff)/gam1

c Again, a few useful values
      tet1 = normal(3)*U(3) - normal(2)*U(4)
      tet2 = normal(1)*U(4) - normal(3)*U(2)
      tet3 = normal(2)*U(2) - normal(1)*U(3)

      c2  = gam*(U(5)+pstiff)/U(1)
      c   = DSQRT(c2)
      c2  = 1.d0/c2
      q   = 0.5d0*vit2

!================ Start Hagstorm & Hariharan BC =====================

c Compute Riemann Invariants and use them to infer truncated boundary conditions
      radius = DSQRT(xface(1)*xface(1)+xface(2)*xface(2)
     & +xface(3)*xface(3)) ! Assuming center is at 0,0,0 

      ! Choose the direction for 1D analysis. It makes more sense to use 
      ! radius vector rather than face normal and this fact is 
      ! confirmed by results
      rnx = xface(1)/radius
      rny = xface(2)/radius
      rnz = xface(3)/radius

      !rnx = normal(1)
      !rny = normal(2)
      !rnz = normal(3)
 

      !if(U(2)*rnx+U(3)*rny+U(4)*rnz-c.ge.0) then
      !Snew    =    U(2)*rnx   +U(3)*rny+   U(4)*rnz - 2*c/gam1

      !else

      cinf     = sqrt(gam*(Vinf(5)+pstiff)/Vinf(1))
      Rplus    =    U(2)*rnx   +U(3)*rny+   U(4)*rnz + 2*c/gam1
      Rinf     = Vinf(2)*rnx+Vinf(3)*rny+Vinf(4)*rnz + 2*cinf/gam1
      Entropy  = (U(5)+pstiff)/U(1)**gam

      !Tuning is not required, so these lines can be ignored
      !kt       = 1.0
      !ku       = 1.0
!      dSdt1    = cinf*(Rplus-Rinf)/(2.*radius)
!      dSdt2    = (U(2)*rnx   +U(3)*rny+   U(4)*rnz)*cinf/radius
      dSdt3    = 2*(c-cinf)*cinf/(radius*gam1)

!COMMENTED OUT BY KW
!     if(abs(sold).le.1.e-7) sold=-2.*cinf/gam1

c      Snew     = Sold + dt*dSdt3
c      vbar = 0.5*(Rplus+Snew)
c      cbar = 0.25*gam1*(Rplus-Snew)
      vbar = 0.5*(Rplus+Sold)
      cbar = 0.25*gam1*(Rplus-Sold)

      !if(xface(2).gt.20.and.abs(c-cinf).gt.1.e-6) then
!       if(abs(xface(1)-1.20335056676986).le.1.e-6.and.
!     <	  abs(xface(2)-20.3196242693497).le.1.e-6.and.
!     <    abs(xface(3)-23.4989073713820).le.1.e-6) then
!       print*,xface(1),xface(2),xface(3)
!       print*,U(2)*rnx+U(3)*rny+U(4)*rnz,rnx*normal(1)
!     < +rny*normal(2)+rnz*normal(3)
!       print*,dsdt1,dsdt2,dsdt3
!       print*,c,cinf,cbar
!       print*,Rplus,Snew,Sold
!       print*,Vinf(1),(cbar*cbar/(Entropy*gam))**(1./gam1) 
!       print*,dt,dSdt3
!       print*,'=========================='
!      endif

!     Trick Ghidaglia into thinking projected solution is the Infinity condition 
!     In practice, this may not be required and phi=F(Vinf) may be sufficient

      Vinf(1) = (cbar*cbar/(Entropy*gam))**(1./gam1)
      Vinf(5) = Entropy*Vinf(1)**gam - pstiff
      Vinf(2) = vbar*rnx
      Vinf(3) = vbar*rny
      Vinf(4) = vbar*rnz

      VdotNinf = Vinf(2)*normal(1)+Vinf(3)*normal(2)+Vinf(4)*normal(3)
      vitinf2  = Vinf(2)**2+Vinf(3)**2+Vinf(4)**2
      enerinf  = 0.5d0*Vinf(1)*vitinf2 + (Vinf(5)+gam*pstiff)/gam1

      !endif
!================ End Hagstorm & Hariharan BC =====================

c Computation of the physical fluxes
      flux(1) = U(1)*(VdotN-vitno)
      flux(2) = flux(1)*U(2)+U(5)*normal(1)
      flux(3) = flux(1)*U(3)+U(5)*normal(2)
      flux(4) = flux(1)*U(4)+U(5)*normal(3)
      flux(5) = (ener+U(5))*(VdotN-vitno) + U(5)*vitno

      fluxinf(1) = Vinf(1)*(VdotNinf-vitno)
      fluxinf(2) = fluxinf(1)*Vinf(2)+Vinf(5)*normal(1)
      fluxinf(3) = fluxinf(1)*Vinf(3)+Vinf(5)*normal(2)
      fluxinf(4) = fluxinf(1)*Vinf(4)+Vinf(5)*normal(3)
      fluxinf(5) = (enerinf+Vinf(5))*(VdotNinf-vitno)
     &           + Vinf(5)*vitno

c Computation of righteigenvectors using values at Interior
c formulae found in Hirsh and Arthurs notes (Research/latex/notes/
c preconditioning_stiffenedgas.pdf and RoeTurkel_SG_Flux.pdf)
      r1(1) = normal(1)
      r1(2) = U(2)*normal(1)
      r1(3) = U(3)*normal(1)+normal(3)
      r1(4) = U(4)*normal(1)-normal(2)
      r1(5) = q*normal(1) + tet1

      r2(1) = normal(2)
      r2(2) = U(2)*normal(2)-normal(3)
      r2(3) = U(3)*normal(2)
      r2(4) = U(4)*normal(2)+normal(1)
      r2(5) = q*normal(2) + tet2

      r3(1) = normal(3)
      r3(2) = U(2)*normal(3)+normal(2)
      r3(3) = U(3)*normal(3)-normal(1)
      r3(4) = U(4)*normal(3)
      r3(5) = q*normal(3) + tet3

      r4(1) = 0.5d0*c2
      r4(2) = r4(1)*(U(2)+c*normal(1))
      r4(3) = r4(1)*(U(3)+c*normal(2))
      r4(4) = r4(1)*(U(4)+c*normal(3))
      r4(5) = r4(1)*(ener/U(1)+U(5)/U(1)+c*VdotN)

      r5(1) = 0.5d0*c2
      r5(2) = r4(1)*(U(2)-c*normal(1))
      r5(3) = r4(1)*(U(3)-c*normal(2))
      r5(4) = r4(1)*(U(4)-c*normal(3))
      r5(5) = r4(1)*(ener/U(1)+U(5)/U(1)-c*VdotN)

c Computation of eigenvalues using value at Interior
      vp(1) = VdotN-vitno
      vp(2) = vp(1)
      vp(3) = vp(1)
      vp(4) = VdotN+c-vitno
      vp(5) = VdotN-c-vitno

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
     &        ((normal(1)*(1.d0 - gam1*q*c2) - tet1)*temp(1) +
     &        (normal(1)*gam1*U(2)*c2)*temp(2)  +
     &        (normal(3)  + (normal(1)*gam1*U(3)*c2))*temp(3)   +
     &        (-normal(2) + (normal(1)*gam1*U(4)*c2))*temp(4)   -
     &        (normal(1)*gam1*c2)*temp(5))

      flur2  =
     &        ((normal(2)*(1.d0 - gam1*q*c2) - tet2)*temp(1) +
     &        (-normal(3) + (normal(2)*gam1*U(2)*c2))*temp(2)   +
     &        (normal(2)*gam1*U(3)*c2)*temp(3)  +
     &        (normal(1)  + (normal(2)*gam1*U(4)*c2))*temp(4)   -
     &        (normal(2)*gam1*c2)*temp(5))
     
      flur3  =
     &        ((normal(3)*(1.d0 - gam1*q*c2) - tet3)*temp(1) +
     &        (normal(2)  + (normal(3)*gam1*U(2)*c2))*temp(2)   +
     &        (-normal(1) + (normal(3)*gam1*U(3)*c2))*temp(3)   +
     &        (normal(3)*gam1*U(4)*c2)*temp(4)  -
     &        (normal(3)*gam1*c2)*temp(5))
       
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
     &         ((-c*VdotN   + gam1*q)*temp(1)  +
     &         ( c*normal(1) - gam1*U(2))*temp(2) +
     &         ( c*normal(2) - gam1*U(3))*temp(3) +
     &         ( c*normal(3) - gam1*U(4))*temp(4) +
     &         gam1*temp(5))

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
     &         (( c*VdotN  + gam1*q)*temp(1) +
     &         (-c*normal(1) - gam1*U(2))*temp(2) +
     &         (-c*normal(2) - gam1*U(3))*temp(3) +
     &         (-c*normal(3) - gam1*U(4))*temp(4) +
     &         gam1*temp(5))


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

!     Reset Vinf to original values

      Vinf(1) = Uinf(1)
      Vinf(2) = Uinf(2)/Uinf(1)
      Vinf(3) = Uinf(3)/Uinf(1)
      Vinf(4) = Uinf(4)/Uinf(1)
      Vinf(5) = gam1*(Uinf(5)-
     &              0.5d0*(Uinf(2)**2+Uinf(3)**2+Uinf(4)**2)/Uinf(1))
     &              -gam*pstiff

c
c For one and two equation turbulence models
c

      if (type.eq.1) then
         Vinf(6) = Uinf(6) / Uinf(1)
         updir = 0.5d0 + dsign(0.5d0, phi(1))
         phi(6) = phi(1) * (updir * U(6) + (1.0d0 - updir) * Vinf(6))
      else if (type.eq.2) then
         Vinf(6) = Uinf(6) / Uinf(1)
         Vinf(7) = Uinf(7) / Uinf(1)
         updir = 0.5d0 + dsign(0.5d0, phi(1))
         phi(6) = phi(1) * (updir * U(6) + (1.0d0 - updir) * Vinf(6))
         phi(7) = phi(1) * (updir * U(7) + (1.0d0 - updir) * Vinf(7))
      endif
      
      END
