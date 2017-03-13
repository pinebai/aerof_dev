       SUBROUTINE ERIEMANNGG(DL,UL,PL,DR,UR,PR,
     &                      PM,UM,RIL,RIR,
     &                      gamL,PREFL,gamR,PREFR,errcod,
     &                      pcutl,pcutr,rcutl,rcutr,
     &                      TOLPRE,NRITER)
C     PREFL  p_c (L)
*     
*-----------------------------------------------------------------*
*                                                                 *
C     Exact Riemann Solver for the Time-Dependent                 *
C     One Dimensional Euler Equations                             *
C     ref is Int. Journal For Num. Meth. In Fluids 28 (1998)      *
C     pp 395-418 by Ivings, Causon and Toro			  *
*                                                                 *
*-----------------------------------------------------------------*
*     
      IMPLICIT NONE
     
C     Declaration of variables:
     
      INTEGER i, NRITER,errcod

      REAL*8 DL,UL,PL,DR,UR,PR
      REAL*8 PM,UM,RIL,RIR
      REAL*8 gamR, PREFR, gamL, PREFL

      REAL*8 G1R, G2R, G3R, G4R, G5R
      REAL*8 G1L, G2L, G3L, G4L, G5L
      REAL*8 CL, CR, PLL, PRR
      REAL*8 FL, FR, FLD, FRD, UDIFF
      REAL*8 POLD, PRATIO, CHANGE
      REAL*8 POLDL,POLDR
      REAL*8 AK, BK, QRT
      REAL*8 TOLPRE
      REAL*8 pcutl,rcutl,pcutr,rcutr
      REAL*8 pcut

      errcod = 0

C     TOLPRE is precision tolerance for convergence of newton-raphson
C     iteration and NRITER is max # of newton-raphson iterations
c     TOLPRE = 1.0E-03
c     NRITER = 100

      PRR = PR + PREFR
      PLL = PL + PREFL

      PCUT = MAX(pcutl,pcutr)
*     
C     Constants for SG EOS R
*     
      G1R = (gamR - 1.0)/(2.0*gamR)
      G2R = 2.0/(gamR - 1.0)
      G3R = 2.0/(gamR + 1.0)
      G4R = G3R/G2R
      G5R = G1R - 1.0
*     
C     Constants for SG EOS L
*     
      G1L = (gamL - 1.0)/(2.0*gamL)
      G2L = 2.0/(gamL - 1.0)
      G3L = 2.0/(gamL + 1.0)
      G4L = G3L/G2L
      G5L = G1L - 1.0

*     
C     Compute sound speeds
*     
      CL = DSQRT(gamL*PLL/DL)
      CR = DSQRT(gamR*PRR/DR)
      IF(gamR*PRR/DR.lt.0.0 .or. gamL*PLL/DL.lt.0.0)THEN
         WRITE(*,*) 'Negative sound speeds'
         WRITE(*,*) 'CR, CL = ', gamR*PRR/DR, gamL*PLL/DL
      ENDIF
      

      UDIFF = UR - UL

C     SOLVER FOR EACH CASE


C     INITIALIZE NEWTON ITERATION
      
      PM = MIN(PR,PL)
      IF(PM .lt. 0.0) THEN
        PM = TOLPRE
      ENDIF
      POLD = PM
         
*     
*     
      DO 10 I = 1, NRITER
*     
C     FL RELATIONS
*     
         POLDL = POLD + PREFL
         IF(POLD.LE.PL)THEN
*     
C     Rarefaction wave
*     
            PRATIO = POLDL/PLL
            FL     = G2L*CL*(PRATIO**G1L - 1.0)
            FLD    = PRATIO**G5L/(DL*CL)
         ELSE
*     
C     Shock wave
*     
            AK   = G3L/DL
            BK   = G4L*PLL
            QRT  = SQRT(AK/(BK + POLDL))
            FL   = (POLDL - PLL)*QRT
            FLD  = (1.0 - 0.5*(POLDL - PLL)/(BK + POLDL))*QRT
         ENDIF


*
C     FR RELATIONS
*
         POLDR = POLD + PREFR
         IF(POLD.LE.PR)THEN
*     
C     Rarefaction wave
*     
            PRATIO = POLDR/PRR
            FR     = G2R*CR*(PRATIO**G1R - 1.0)
            FRD    = PRATIO**G5R/(DR*CR)
         ELSE
*     
C     Shock wave
*     
            AK   = G3R/DR
            BK   = G4R*PRR
            QRT  = SQRT(AK/(BK + POLDR))
            FR   = (POLDR - PRR)*QRT
            FRD  = (1.0 - 0.5*(POLDR - PRR)/(BK + POLDR))*QRT
         ENDIF

         PM     = POLD - (FL + FR + UDIFF)/(FLD + FRD)
         CHANGE = 2.0*ABS((PM- POLD)/(PM+ POLD))

         IF (PM.LT.PCUT) PM = PCUT

         IF(CHANGE.LE.TOLPRE)GOTO 20
         IF(PM.LT.0.0) PM = TOLPRE*0.001;
         POLD  = PM



 10   CONTINUE
*     
      WRITE(*,*) ' *** Warning: '
      WRITE(*,*) 'Newton-Raphson reached max num. iterations ', NRITER
      WRITE(*,*) 'without converging to the desired tolerance', TOLPRE
      WRITE(*,*) PM, PL, PR
      WRITE(*,*) 'INPUT RIEMANN (DL,UL,PL,DR,UR,PR)'
      WRITE(*,*) DL,UL,PL,DR,UR,PR
      WRITE(*,*) 'PARAMETERS (GAML,PREFL,GAMR,PREFR)'
      WRITE(*,*) GAML,PREFL,GAMR,PREFR
      WRITE(*,*) 'OUTPUT (PM,UM,RIL,RIR)'
      WRITE(*,*) PM,UM,RIL,RIR
      WRITE(*,*) ' *** '
      errcod = 1
*     
 20   CONTINUE

      UM = 0.5*(UL + UR + FR - FL)

      IF (PM .LE. PL) THEN
         RIL  = DL*((PM+PREFL)/PLL)**(1.0/gamL)
C        by isentropic relationship
      ELSE
         RIL  = DL* ((PM+PREFL)/PLL+G4L)/((PM+PREFL)*G4L/PLL + 1.0)
C        by rankine hugoniot relationship (24) of reference.
      ENDIF

      IF (PM .LE. PR) THEN
         RIR  = DR*((PM+PREFR)/PRR)**(1.0/gamR)
      ELSE
         RIR  = DR* ((PM+PREFR)/PRR+G4R)/((PM+PREFR)*G4R/PRR + 1.0)
      ENDIF

      IF (RIL.LT.rcutl) RIL = rcutl
      IF (RIR.LT.rcutr) RIR = rcutr

      END
