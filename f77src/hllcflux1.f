      SUBROUTINE HLLCFLUX1(gamma,gam,pstiff,enormal,evitno,
     &     Ug,Ud,phi1,mach,k1,cmach,irey,prec)


c---------------------------------------------------------------------   
c Computes the convective fluxes using the approximate Riemann
c solver of HLLC with optional low Mach preconditionning
c
c Reference: Extension of Harten-Lax-van Leer Scheme
c            for Flows at All Speeds
c            By H. Luo, J. D. Baum, and R. Lohner
c            (AIAA Journal - Vol. 43, No. 6, June 2005)
c---------------------------------------------------------------------   

      IMPLICIT NONE
      REAL*8 Ug(*), Ud(*), normal(0:2), enormal(3), evitno, phi1
      REAL*8 VdotN , rnorm, invnorm, updir, vitno
      REAL*8 energ, enerd
      REAL*8 flu(0:5),solLft(0:6),solRgt(0:6)
      REAL*8 vnLft,vnRgt,cLft,cRgt,rhoLft,rhoRgt,HLft,HRgt
      REAL*8 vpLft,vpRgt,vpLftRoe,vpRgtRoe
      REAL*8 rhoInv,roeMoy(0:5),qRoe,vnRoe,cRoe,cRoe2,pStar,vnStar
      REAL*8 SLft,SRgt,SStar,uStar(0:5)
      REAL*8 gam , gam1, vitg2, vitd2, pstiff
      REAL*8 locMach, cmach, irey
      REAL*8 gamma, mach, k1, shock, beta, beta2
      INTEGER prec

c

      gam1      = gam - 1.d0

      rnorm     = DSQRT(enormal(1)*enormal(1) + enormal(2)*enormal(2) +
     &             enormal(3)*enormal(3))

      invnorm   = 1.0d0 / rnorm
c
      normal(0) = enormal(1) * invnorm
      normal(1) = enormal(2) * invnorm
      normal(2) = enormal(3) * invnorm
c
      vitno = evitno * invnorm
c
      solLft(0) = Ug(1)
      solLft(1) = Ug(2)
      solLft(2) = Ug(3)
      solLft(3) = Ug(4)
      solLft(4) = Ug(5)
      solLft(5) = (Ug(5)+gam*pstiff)/gam1+0.5d0*Ug(1)*(Ug(2)**2+
     &     Ug(3)**2+Ug(4)**2)
c
      solRgt(0) = Ud(1)
      solRgt(1) = Ud(2)
      solRgt(2) = Ud(3)
      solRgt(3) = Ud(4)
      solRgt(4) = Ud(5)
      solRgt(5) = (Ud(5)+gam*pstiff)/gam1+0.5d0*Ud(1)*(Ud(2)**2+
     &     Ud(3)**2+Ud(4)**2)

c
c     Celerity and normal component of the velocity
c
      cLft  = sqrt(abs(gam*(solLft(4)+pstiff)/solLft(0)))
      cRgt  = sqrt(abs(gam*(solRgt(4)+pstiff)/solRgt(0)))
      vnLft = normal(0)*solLft(1)+normal(1)*solLft(2)
     &    +normal(2)*solLft(3)
      vnRgt = normal(0)*solRgt(1)+normal(1)*solRgt(2)
     &    +normal(2)*solRgt(3)

c
c     Roe averaging
c
      rhoLft = sqrt(solLft(0))
      rhoRgt = sqrt(solRgt(0))
      HLft   = (solLft(5)+solLft(4))/solLft(0) ! enthalpy, H
      HRgt   = (solRgt(5)+solRgt(4))/solRgt(0) ! enthalpy, H

c     Roe's average to get SLft and SRgt
      rhoInv    = 1./(rhoLft+rhoRgt)
      rhoLft    = rhoLft*rhoInv
      rhoRgt    = rhoRgt*rhoInv
      roeMoy(1) = (rhoLft*solLft(1)+rhoRgt*solRgt(1)) ! u tilde
      roeMoy(2) = (rhoLft*solLft(2)+rhoRgt*solRgt(2)) ! v tilde
      roeMoy(3) = (rhoLft*solLft(3)+rhoRgt*solRgt(3)) ! w tilde
      roeMoy(4) = (rhoLft*HLft+rhoRgt*HRgt)           ! H tilde
      qRoe      = 0.5d0*(roeMoy(1)*roeMoy(1)
     & +roeMoy(2)*roeMoy(2)+roeMoy(3)*roeMoy(3))      ! q=.5 |u|^2=Ec/rho
c
c     H-q = (1/rho) (E + p) - Ec/rho  ==> p = rho*gam1/gam (H-q)
c
      cRoe      = sqrt(abs(gam1*(roeMoy(4)-qRoe)))   ! c tilde
      vnRoe     = normal(0)*roeMoy(1)+normal(1)*roeMoy(2)
     &           +normal(2)*roeMoy(3)

c  Computing the preconditionning coefficient

      cRoe2 = cRoe*cRoe

      if (prec .eq. 0) then
         beta = 1.d0
      else
c     local Preconditioning (ARL)
         locMach = DSQRT(2.0d0*qRoe/cRoe2)
         beta = MIN((1.0d0+DSQRT(irey))*MAX(k1*locMach, mach),cmach)
      end if
c
      beta2 = beta*beta

c
c  Wave speed SLft, Srgt and SStar if necessary, with preconditionning if necessary
c

c Computing the preconditionned eigenvalues
      vpLft = 0.5d0*((1.d0+beta2)*vnLft -
     &      DSQRT(((1.d0-beta2)*vnLft)**2 +
     &      4.d0*beta2*cLft**2))
      vpLftRoe =  0.5d0*((1.d0+beta2)*vnRoe -
     &      DSQRT(((1.d0-beta2)*vnRoe)**2 +
     &      4.d0*beta2*cRoe**2))
c
      vpRgt = 0.5d0*((1.d0+beta2)*vnRgt +
     &      DSQRT(((1.d0-beta2)*vnRgt)**2 +
     &      4.d0*beta2*cRgt**2))
      vpRgtRoe = 0.5d0*((1.d0+beta2)*vnRoe +
     &      DSQRT(((1.d0-beta2)*vnRoe)**2 +
     &      4.d0*beta2*cRoe**2))

c Computing the signal velocities
      SLft = min(vpLft,vpLftRoe)
      SRgt = max(vpRgt,vpRgtRoe)

c  Flft : left flux by splitting
      if ( SLft-vitno .ge. 0 ) then

         flu(0) = solLft(0)*(vnLft-vitno)                           
         flu(1) = flu(0)*solLft(1)+solLft(4)*normal(0)
         flu(2) = flu(0)*solLft(2)+solLft(4)*normal(1)
         flu(3) = flu(0)*solLft(3)+solLft(4)*normal(2)
         flu(4) = solLft(5)*(vnLft-vitno)+solLft(4)*vnLft

         phi1   = 2.d0*flu(0)

         go to 1000

      endif
c
c  Frgt : right flux by splitting
c
      if ( SRgt-vitno .le. 0 ) then

         flu(0) = solRgt(0)*(vnRgt-vitno)                           
         flu(1) = flu(0)*solRgt(1)+solRgt(4)*normal(0)
         flu(2) = flu(0)*solRgt(2)+solRgt(4)*normal(1)
         flu(3) = flu(0)*solRgt(3)+solRgt(4)*normal(2)
         flu(4) = solRgt(5)*(vnRgt-vitno)+solRgt(4)*vnRgt

         phi1   = 2.d0*flu(0)

         go to 1000
c
      endif
c
c  Calcul de SStar, contact wave speed
c
      SStar =
     &        solRgt(4)-solLft(4) +
     &        vnLft*solLft(0)*(SLft-vnLft) -
     &        vnRgt*solRgt(0)*(SRgt-vnRgt)
c
      SStar = SStar/(solLft(0)*(SLft-vnLft) - solRgt(0)*(SRgt-vnRgt))
c
      if ( SStar-vitno .ge. 0 ) then
c     ========================
c
        pStar    = solLft(0)*(SLft-vnLft)*(SStar-vnLft)+solLft(4)
        uStar(0) = solLft(0)*(SLft-vnLft)/(SLft-SStar)
        uStar(1) = uStar(0)*(solLft(1)+(SStar-vnLft)*normal(0))
        uStar(2) = uStar(0)*(solLft(2)+(SStar-vnLft)*normal(1))
        uStar(3) = uStar(0)*(solLft(3)+(SStar-vnLft)*normal(2))
        uStar(4) = uStar(0)*(solLft(5)/solLft(0)+
     &   (SStar-vnLft)*(SStar+solLft(4)/(solLft(0)*(SLft-vnLft))))
        vnStar   =
     &   (uStar(1)*normal(0)+uStar(2)*normal(1)
     &      +uStar(3)*normal(2))/uStar(0)
c
        flu(0) = uStar(0)*(vnStar-vitno)
        flu(1) = uStar(1)*(vnStar-vitno)+pStar*normal(0)
        flu(2) = uStar(2)*(vnStar-vitno)+pStar*normal(1)
        flu(3) = uStar(3)*(vnStar-vitno)+pStar*normal(2)
        flu(4) = uStar(4)*(vnStar-vitno)+pStar*vnStar

        phi1   = 2.d0*flu(0)

        go to 1000
c
c     ******
c
      else     ! if ( SStar-vitno < 0 )
c     =====
c
c  Fhllc = Frgt + Srgt*(UStar-Urgt) = F(Ustar)
c
        pStar    = solRgt(0)*(SRgt-vnRgt)*(SStar-vnRgt)+solRgt(4)
        uStar(0) = solRgt(0)*(SRgt-vnRgt)/(SRgt-SStar)
        uStar(1) = uStar(0)*(solRgt(1)+(SStar-vnRgt)*normal(0))
        uStar(2) = uStar(0)*(solRgt(2)+(SStar-vnRgt)*normal(1))
        uStar(3) = uStar(0)*(solRgt(3)+(SStar-vnRgt)*normal(2))
        uStar(4) = uStar(0)*(solRgt(5)/solRgt(0)+
     &   (SStar-vnRgt)*(SStar+solRgt(4)/(solRgt(0)*(SRgt-vnRgt))))
        vnStar   =
     &   (uStar(1)*normal(0)+uStar(2)*normal(1)+uStar(3)*normal(2))
     &     /uStar(0) 

        flu(0) = uStar(0)*(vnStar-vitno)
        flu(1) = uStar(1)*(vnStar-vitno)+pStar*normal(0)
        flu(2) = uStar(2)*(vnStar-vitno)+pStar*normal(1)
        flu(3) = uStar(3)*(vnStar-vitno)+pStar*normal(2)
        flu(4) = uStar(4)*(vnStar-vitno)+pStar*vnStar

        phi1   = 2.d0*flu(0)

        go to 1000
c
      endif    ! if ( SStar <> 0 )
c     *****
c
 1000 continue

      phi1 = phi1*0.5d0*rnorm


      END
