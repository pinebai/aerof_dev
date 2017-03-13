      SUBROUTINE HLLCJAC(type,gamma,gam,pstiff,enormal,evitno,
     &     Ug,Ud,jacL,jacR,mach,k1,cmach,irey,prec)

c---------------------------------------------------------------------   
c This routine computes the jacobian of the hllc flux with respect to
c the conservative variables taken at the vectors Ug, Ud.
c normal is the normal of the boundary concerned by the flux.
c
c Reference: Average-State Jacobians and Implicit Methods for Compressible
c            Viscous and Turbulent Flows
c            By P. Batten, M. A. Leschizner, and U. C. Goldberg
c            (Journal of Computational Physics - Vol. 137, pp. 38-78, 1997)
c---------------------------------------------------------------------   

      IMPLICIT NONE
      REAL*8 Ug(*), Ud(*), normal(0:2), enormal(3), evitno
      REAL*8 jacL(*), jacR(*), HL(5,5), HR(5,5)
      REAL*8 rnorm, invnorm, updir, vitno
      REAL*8 energ, enerd
      REAL*8 solLft(0:5), solRgt(0:5), dSdULft(5), dSdURgt(5)
      REAL*8 dpdULft(5), dpdURgt(5), drdULft(5), drdURgt(5)
      REAL*8 drudULft(5), drudURgt(5), drvdULft(5), drvdURgt(5)
      REAL*8 drwdULft(5), drwdURgt(5), dedULft(5), dedURgt(5)
      REAL*8 rtilde, omLft, omRgt
      REAL*8 vnLft,vnRgt,cLft,cRgt,rhoLft,rhoRgt,HLft,HRgt
      REAL*8 qLft, qRgt, phiLft, phiRgt
      REAL*8 vpLft,vpRgt,vpLftRoe,vpRgtRoe
      REAL*8 rhoInv,roeMoy(0:4),qRoe,vnRoe,cRoe,cRoe2,pStar,vnStar
      REAL*8 SLft,SRgt,SStar,uStar(0:4)
      REAL*8 gam , gam1, vitg2, vitd2, pstiff
      REAL*8 locMach, cmach, irey
      REAL*8 gamma, mach, k1,  beta, beta2
      INTEGER k, i, dim
      INTEGER type, prec

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
c     Celerity, normal component of the velocity and modulus
c
      cLft  = sqrt(abs(gam*(solLft(4)+pstiff)/solLft(0)))
      cRgt  = sqrt(abs(gam*(solRgt(4)+pstiff)/solRgt(0)))
      vnLft = normal(0)*solLft(1)+normal(1)*solLft(2)
     &    +normal(2)*solLft(3)
      vnRgt = normal(0)*solRgt(1)+normal(1)*solRgt(2)
     &    +normal(2)*solRgt(3)
      qLft  = 0.5d0*(solLft(1)**2+solLft(2)**2+solLft(3)**2)
      qRgt  = 0.5d0*(solRgt(1)**2+solRgt(2)**2+solRgt(3)**2)

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

c----------------------------------------------------------------------
c Supersonic case #1 (flow to the right), computing jacobians HL and HR
c----------------------------------------------------------------------
      if ( SLft-vitno .ge. 0 ) then

         CALL JACSUPERSONIC(normal, vitno, gam, Ug, vnLft,
     &        qLft, HLft, HL)

         do i = 1,5
            do k = 1,5
               HR(i,k) = 0.d0
            enddo
         enddo

         go to 1000

      endif

c---------------------------------------------------------------------
c Supersonic case #2 (flow to the left), computing jacobians HL and HR
c---------------------------------------------------------------------
      if ( SRgt-vitno .le. 0 ) then

         CALL JACSUPERSONIC(normal, vitno, gam, Ud, vnRgt,
     &        qRgt, HRgt, HR)

         do i = 1,5
            do k = 1,5
               HL(i,k) = 0.d0
            enddo
         enddo

         go to 1000
c
      endif

c---------------------------------------------------------------------
c Subsonic cases - Computing SStar (contact wave speed)
c---------------------------------------------------------------------
      SStar = solRgt(4)-solLft(4) +
     &        vnLft*solLft(0)*(SLft-vnLft) -
     &        vnRgt*solRgt(0)*(SRgt-vnRgt)
c
      SStar = SStar/(solLft(0)*(SLft-vnLft) - solRgt(0)*(SRgt-vnRgt))
c
      omLft = 1.d0/(SLft-SStar)
      omRgt = 1.d0/(SRgt-SStar)
c
      rtilde = Ud(1)*(SRgt-vnRgt)-Ug(1)*(SLft-vnLft)
      rtilde = 1.d0/rtilde
c
      phiLft = 2.d0*qLft
      phiRgt = 2.d0*qRgt
c
c Computing derivatives of the contact wave speed with respect to Ug and Ud
c
      CALL SMDERIVATIVE(normal, gam, Ug, vnLft, phiLft, rtilde,
     &     SLft, SStar, dSdULft)
      CALL SMDERIVATIVE(normal, gam, Ud, vnRgt, phiRgt, -rtilde,
     &     SRgt, SStar, dSdURgt)
c
c Computing derivatives of the star state pressure
c
      do i = 1,5
         dpdULft(i)  = Ud(1)*(SRgt - vnRgt)*dSdULft(i)
         dpdURgt(i)  = Ug(1)*(SLft - vnLft)*dSdURgt(i)
      enddo


      if ( SStar-vitno .ge. 0 ) then
c     ========================
c
        pStar    = solLft(0)*(SLft-vnLft)*(SStar-vnLft)+solLft(4)
        uStar(0) = solLft(0)*(SLft-vnLft)/(SLft-SStar)
        uStar(1) = uStar(0)*(solLft(1)+(SStar-vnLft)*normal(0))
        uStar(2) = uStar(0)*(solLft(2)+(SStar-vnLft)*normal(1))
        uStar(3) = uStar(0)*(solLft(3)+(SStar-vnLft)*normal(2))
        uStar(4) = uStar(0)*(solLft(5)/solLft(0)+
     &       (SStar-vnLft)*(SStar+solLft(4)/(solLft(0)*(SLft-vnLft))))
        vnStar   = (uStar(1)*normal(0)+uStar(2)*normal(1)
     &       +uStar(3)*normal(2))/uStar(0)
c
c Derivatives of the density (normal and cross)
c
        CALL RDERIVATIVE(normal, uStar(0), omLft, SLft,
     &       dSdULft, drdULft)
        CALL RCROSS(uStar(0), omLft, dSdURgt, drdURgt)
c
c Derivatives of the left star momentum
c
        CALL MOMDERIVATIVE(normal, gam, Ug, uStar, vnLft, phiLft,
     &       omLft, SLft, dpdULft, dSdULft, drudULft, 1)
        CALL MOMDERIVATIVE(normal, gam, Ug, uStar, vnLft, phiLft,
     &       omLft, SLft, dpdULft, dSdULft, drvdULft, 2)
        CALL MOMDERIVATIVE(normal, gam, Ug, uStar, vnLft, phiLft,
     &       omLft, SLft, dpdULft, dSdULft, drwdULft, 3)
c
c Cross derivatives of the left star momentum
c
        CALL MOMCROSS(normal, uStar, omLft, dpdURgt, dSdURgt,
     &       drudURgt, 1)
        CALL MOMCROSS(normal, uStar, omLft, dpdURgt, dSdURgt,
     &       drvdURgt, 2)
        CALL MOMCROSS(normal, uStar, omLft, dpdURgt, dSdURgt,
     &       drwdURgt, 3)
c
c Derivatives of the energy (normal and cross)
c
        CALL EDERIVATIVE(normal, gam, Ug, vnLft, phiLft, omLft, 
     &       solLft(5), solLft(4), uStar(4), pStar, SLft, SStar,
     &       dpdULft, dSdULft, dedULft)
        CALL ECROSS(omLft, uStar(4), pStar, SStar, dpdURgt,
     &       dSdURgt, dedURgt)
c
c Final computation of the Jacobian matrices HL and HR
c
       do k = 1,5
          HL(1,k) = (SStar-vitno)*drdULft(k) + uStar(0)*dSdULft(k)
          HR(1,k) = (SStar-vitno)*drdURgt(k) + uStar(0)*dSdURgt(k)
       enddo
c
       do k = 1,5
          HL(2,k) = (SStar-vitno)*drudULft(k) + uStar(1)*dSdULft(k) +
     &         normal(0)*dpdULft(k)
          HR(2,k) = (SStar-vitno)*drudURgt(k) + uStar(1)*dSdURgt(k) +
     &         normal(0)*dpdURgt(k)
       enddo
c
       do k = 1,5
          HL(3,k) = (SStar-vitno)*drvdULft(k) + uStar(2)*dSdULft(k) +
     &         normal(1)*dpdULft(k)
          HR(3,k) = (SStar-vitno)*drvdURgt(k) + uStar(2)*dSdURgt(k) +
     &         normal(1)*dpdURgt(k)
       enddo
c
       do k = 1,5
          HL(4,k) = (SStar-vitno)*drwdULft(k) + uStar(3)*dSdULft(k) +
     &         normal(2)*dpdULft(k)
          HR(4,k) = (SStar-vitno)*drwdURgt(k) + uStar(3)*dSdURgt(k) +
     &         normal(2)*dpdURgt(k)
       enddo
c
       do k = 1,5
          HL(5,k) = (SStar-vitno)*dedULft(k) + SStar*dpdULft(k) +
     &         (uStar(4) + pStar)*dSdULft(k)
          HR(5,k) = (SStar-vitno)*dedURgt(k) + SStar*dpdURgt(k) +
     &         (uStar(4) + pStar)*dSdURgt(k)
       enddo

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
     &       (SStar-vnRgt)*(SStar+solRgt(4)/(solRgt(0)*(SRgt-vnRgt))))
        vnStar   =  (uStar(1)*normal(0)+uStar(2)*normal(1)+uStar(3)*
     &       normal(2))/uStar(0) 
c
c Derivatives of the density (normal and cross)
c
        CALL RDERIVATIVE(normal, uStar(0), omRgt, SRgt,
     &       dSdURgt, drdURgt)
        CALL RCROSS(uStar(0), omRgt, dSdULft, drdULft)
c
c Derivatives of the left star momentum
c
        CALL MOMDERIVATIVE(normal, gam, Ud, uStar, vnRgt, phiRgt,
     &       omRgt, SRgt, dpdURgt, dSdURgt, drudURgt, 1)
        CALL MOMDERIVATIVE(normal, gam, Ud, uStar, vnRgt, phiRgt,
     &       omRgt, SRgt, dpdURgt, dSdURgt, drvdURgt, 2)
        CALL MOMDERIVATIVE(normal, gam, Ud, uStar, vnRgt, phiRgt,
     &       omRgt, SRgt, dpdURgt, dSdURgt, drwdURgt, 3)
c
c Cross derivatives of the left star momentum
c
        CALL MOMCROSS(normal, uStar, omRgt, dpdULft, dSdULft,
     &       drudULft, 1)
        CALL MOMCROSS(normal, uStar, omRgt, dpdULft, dSdULft,
     &       drvdULft, 2)
        CALL MOMCROSS(normal, uStar, omRgt, dpdULft, dSdULft,
     &       drwdULft, 3)
c
c Derivatives of the energy (normal and cross)
c
        CALL EDERIVATIVE(normal, gam, Ud, vnRgt, phiRgt, omRgt, 
     &       solRgt(5), solRgt(4), uStar(4), pStar, SRgt, SStar,
     &       dpdURgt, dSdURgt, dedURgt)
        CALL ECROSS(omRgt, uStar(4), pStar, SStar, dpdULft,
     &       dSdULft, dedULft)
c
c Final computation of the Jacobian matrices HL and HR
c
       do k = 1,5
          HL(1,k) = (SStar-vitno)*drdULft(k) + uStar(0)*dSdULft(k)
          HR(1,k) = (SStar-vitno)*drdURgt(k) + uStar(0)*dSdURgt(k)
       enddo
c
       do k = 1,5
          HL(2,k) = (SStar-vitno)*drudULft(k) + uStar(1)*dSdULft(k) +
     &         normal(0)*dpdULft(k)
          HR(2,k) = (SStar-vitno)*drudURgt(k) + uStar(1)*dSdURgt(k) +
     &         normal(0)*dpdURgt(k)
       enddo
c
       do k = 1,5
          HL(3,k) = (SStar-vitno)*drvdULft(k) + uStar(2)*dSdULft(k) +
     &         normal(1)*dpdULft(k)
          HR(3,k) = (SStar-vitno)*drvdURgt(k) + uStar(2)*dSdURgt(k) +
     &         normal(1)*dpdURgt(k)
       enddo
c
       do k = 1,5
          HL(4,k) = (SStar-vitno)*drwdULft(k) + uStar(3)*dSdULft(k) +
     &         normal(2)*dpdULft(k)
          HR(4,k) = (SStar-vitno)*drwdURgt(k) + uStar(3)*dSdURgt(k) +
     &         normal(2)*dpdURgt(k)
       enddo
c
       do k = 1,5
          HL(5,k) = (SStar-vitno)*dedULft(k) + SStar*dpdULft(k) +
     &         (uStar(4) + pStar)*dSdULft(k)
          HR(5,k) = (SStar-vitno)*dedURgt(k) + SStar*dpdURgt(k) +
     &         (uStar(4) + pStar)*dSdURgt(k)
       enddo


        go to 1000
c
      endif    ! if ( SStar <> 0 )
c     *****
c
 1000 continue

c---------------------------------------------------------------------
c Updating the input jacobian matrices from matrices HL and HR
c---------------------------------------------------------------------
      dim = 5 + type
      do i = 1,5
         do k = 1,5
            jacL(dim*(i-1) + k) = rnorm*HL(i,k)
            jacR(dim*(i-1) + k) = rnorm*HR(i,k)
         enddo
         do k = 6,dim
            jacL(dim*(i-1) + k) = 0.d0
            jacR(dim*(i-1) + k) = 0.d0
         enddo
      enddo
      do i = 6,dim
         do k = 1,dim
            jacL(dim*(i-1) + k) = 0.d0
            jacR(dim*(i-1) + k) = 0.d0
         enddo
      enddo

      END

c=====================================================================
c=====================================================================
c Following: subroutines used to compute derivatives of conservative
c variables
c=====================================================================
c=====================================================================


      SUBROUTINE JACSUPERSONIC(normal,vitno,gam,U,vn,q,H,jac)

c---------------------------------------------------------------------   
c This routine computes the trivial jacobian in case of a supersonic
c flow.
c---------------------------------------------------------------------

      IMPLICIT NONE
      REAL*8 U(*), jac(5,5), normal(0:2), vitno
      REAL*8 gam, gam1, vn, q, H

c
      gam1     = gam - 1.d0
c
      jac(1,1) = - vitno
      jac(1,2) = normal(0)
      jac(1,3) = normal(1)
      jac(1,4) = normal(2)
      jac(1,5) = 0.0d0 
c     
      jac(2,1) = gam1*normal(0)*q - U(2)*vn
      jac(2,2) = (1.d0 - gam1)*normal(0)*U(2) + vn - vitno
      jac(2,3) = -gam1*normal(0)*U(3) + U(2)*normal(1)
      jac(2,4) = -gam1*normal(0)*U(4) + U(2)*normal(2) 
      jac(2,5) = gam1*normal(0)
c     
      jac(3,1) = gam1*normal(1)*q - U(3)*vn
      jac(3,2) = -gam1*normal(1)*U(2) + U(3)*normal(0)
      jac(3,3) = (1.d0 - gam1)*normal(1)*U(3) + vn - vitno 
      jac(3,4) = -gam1*normal(1)*U(4) + U(3)*normal(2)
      jac(3,5) = gam1*normal(1)
c     
      jac(4,1) = gam1*normal(2)*q - U(4)*vn
      jac(4,2) = -gam1*normal(2)*U(2) + U(4)*normal(0)
      jac(4,3) = -gam1*normal(2)*U(3) + U(4)*normal(1)
      jac(4,4) = (1.d0 - gam1)*normal(2)*U(4) + vn - vitno 
      jac(4,5) = gam1*normal(2)
c     
      jac(5,1) = (gam1*q - H)*vn
      jac(5,2) = H*normal(0) - gam1*U(2)*vn
      jac(5,3) = H*normal(1) - gam1*U(3)*vn
      jac(5,4) = H*normal(2) - gam1*U(4)*vn
      jac(5,5) = (1.d0 + gam1)*vn - vitno 

      END



      SUBROUTINE SMDERIVATIVE(normal,gam,U,vn,phi,rtilde,S,SStar,der)

c---------------------------------------------------------------------   
c This routine computes the derivative of the contact speed SStar
c with respect to U
c---------------------------------------------------------------------

      IMPLICIT NONE
      REAL*8 U(*), der(5), normal(0:2)
      REAL*8 gam, gam1, vn, phi, rtilde, S, SStar

c
      gam1 = gam - 1.d0
c
      der(1) = rtilde*(-vn**2 + 
     &     phi*gam1/2.d0 + SStar*S)
      der(2) = rtilde*(normal(0)*(2.d0*vn - 
     &     S - SStar) - gam1*U(2))
      der(3) = rtilde*(normal(1)*(2.d0*vn - 
     &     S - SStar) - gam1*U(3))
      der(4) = rtilde*(normal(2)*(2.d0*vn - 
     &     S - SStar) - gam1*U(4))
      der(5) = rtilde*gam1

      END



      SUBROUTINE RDERIVATIVE(normal,rhoStar,om,S,SStarder,der)

c---------------------------------------------------------------------   
c This routine computes the derivative of the density with respect
c to U of the opposite side
c---------------------------------------------------------------------

      IMPLICIT NONE
      REAL*8 normal(0:2), rhoStar
      REAL*8 SStarder(5), der(5)
      REAL*8 om, S

c
      der(1) = om*(S + rhoStar*SStarder(1))
      der(2) = om*(-normal(0) + rhoStar*SStarder(2))
      der(3) = om*(-normal(1) + rhoStar*SStarder(3))
      der(4) = om*(-normal(2) + rhoStar*SStarder(4))
      der(5) = om*rhoStar*SStarder(5)

      END



      SUBROUTINE RCROSS(rhoStar,om,SStarder,der)

c---------------------------------------------------------------------   
c This routine computes the derivative of the density with respect
c to U of the opposite side
c---------------------------------------------------------------------

      IMPLICIT NONE
      REAL*8 rhoStar
      REAL*8 SStarder(5), der(5)
      REAL*8 om
      INTEGER i

c
      do i = 1,5
         der(i) = om*rhoStar*SStarder(i)
      enddo

      END



      SUBROUTINE MOMDERIVATIVE(normal,gam,U,uStar,vn,phi,om,S,pder,
     &     SStarder,der,sw)

c---------------------------------------------------------------------   
c This routine computes the derivative of the component of the momentum
c definied by sw with respect to U
c---------------------------------------------------------------------

      IMPLICIT NONE
      REAL*8 U(*), normal(0:2), uStar(0:4)
      REAL*8 pder(5), SStarder(5), der(5)
      REAL*8 gam, gam1, vn, phi,om, S
      INTEGER sw

c
      gam1 = gam - 1.d0
c
      der(1) = om*(vn*U(sw+1) - normal(sw-1)*
     &     (phi*gam1/2.d0 - pder(1)) + uStar(sw)*SStarder(1))
c
      if ( sw .eq. 1 ) then
         der(2) = om*(S - vn + normal(sw-1)*((gam1 - 
     &        1.d0)*U(sw+1) + pder(2)) + uStar(sw)*SStarder(2))
      else
         der(2) = om*(-U(sw+1)*normal(0) + normal(sw-1)*(gam1
     &        *U(2) + pder(2)) + uStar(sw)*SStarder(2))
      endif
c
      if ( sw .eq. 2 ) then
         der(3) = om*(S - vn + normal(sw-1)*((gam1 - 
     &        1.d0)*U(sw+1) + pder(3)) + uStar(sw)*SStarder(3))
      else
         der(3) = om*(-U(sw+1)*normal(1) + normal(sw-1)*(gam1
     &        *U(3) + pder(3)) + uStar(sw)*SStarder(3))
      endif
c
      if ( sw .eq. 3 ) then
         der(4) = om*(S - vn + normal(sw-1)*((gam1 - 
     &        1.d0)*U(sw+1) + pder(4)) + uStar(sw)*SStarder(4))
      else
         der(4) = om*(-U(sw+1)*normal(2) + normal(sw-1)*(gam1
     &        *U(4) + pder(4)) + uStar(sw)*SStarder(4))
      endif
c
      der(5) = om*(normal(sw-1)*(-gam1
     &     + pder(5)) + uStar(sw)*SStarder(5))
      
      END



      SUBROUTINE MOMCROSS(normal,uStar,om,pder,SStarder,der,sw)

c---------------------------------------------------------------------   
c This routine computes the derivative of the component of the momentum
c definied by sw with respect to U of the opposite side
c---------------------------------------------------------------------

      IMPLICIT NONE
      REAL*8 normal(0:2), uStar(0:4)
      REAL*8 pder(5), SStarder(5), der(5)
      REAL*8 om
      INTEGER sw,i

c
      do i = 1,5
         der(i) = om*(normal(sw-1)*pder(i) + uStar(sw)*SStarder(i))
      enddo

      END



      SUBROUTINE EDERIVATIVE(normal,gam,U,vn,phi,om,ener,p,eStar,
     $     pStar,S,SStar,pder,SStarder,der)

c---------------------------------------------------------------------   
c This routine computes the derivative of the contact speed SStar
c with respect to U
c---------------------------------------------------------------------

      IMPLICIT NONE
      REAL*8 U(*), normal(0:2)
      REAL*8 pder(5), SStarder(5), der(5)
      REAL*8 gam, gam1, vn, phi, om
      REAL*8 S, SStar, eStar, pStar, ener, p

c
      gam1 = gam - 1.d0
c
       der(1) = om*((ener + p)*vn/(U(1)) -
     &       vn*phi*gam1/2.d0 + SStar*pder(1) + (pStar +
     &       eStar)*SStarder(1))
       der(2) = om*(-normal(0)*(ener + p)/(U(1)) +
     &       gam1*U(2)*vn + SStar*pder(2) + (pStar +
     &       eStar)*SStarder(2))
       der(3) = om*(-normal(1)*(ener + p)/(U(1)) +
     &       gam1*U(3)*vn + SStar*pder(3) + (pStar +
     &       eStar)*SStarder(3))
       der(4) = om*(-normal(2)*(ener + p)/(U(1)) +
     &       gam1*U(4)*vn + SStar*pder(4) + (pStar +
     &       eStar)*SStarder(4))
       der(1) = om*(S - vn*gam + SStar*pder(5) +
     &      (pStar + eStar)*SStarder(5))

      END



      SUBROUTINE ECROSS(om,eStar,pStar,SStar,pder,SStarder,der)

c---------------------------------------------------------------------   
c This routine computes the derivative of the contact speed SStar
c with respect to U
c---------------------------------------------------------------------

      IMPLICIT NONE
      REAL*8 pder(5), SStarder(5), der(5)
      REAL*8 om
      REAL*8 SStar, eStar, pStar
      INTEGER i

c
      do i = 1,5
         der(i) = om*(SStar*pder(i) + (pStar + eStar)*SStarder(i))
      enddo

      END
