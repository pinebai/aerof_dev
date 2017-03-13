      SUBROUTINE ERIEMANNGW(DL,UL,PL,DR,UR,PR,
     &                      PM,UM,RIL,RIR,
     &                      ALPHA,BETA,PREF,gam,
     &                      TOLPRE,NRITER)

*     
*-----------------------------------------------------------------*
*                                                                 *
C     Exact Riemann Solver for the Time-Dependent                 *
C     One Dimensional Euler Equations                             *
C     For more info, look at JCP 204(2005) pp 193-221		  *
C     by Liu, Khoo and Wang 		                          *
C     as well as Int. Journal For Num. Meth. In Fluids 28 (1998)  *
C     pp 395-418 by Ivings, Causon and Toro			  *
*                                                                 *
*-----------------------------------------------------------------*
*     
      IMPLICIT NONE
     
C     Declaration of variables:
     
      INTEGER i, NRITER, NATURE

      REAL*8 DL,UL,PL,DR,UR,PR
      REAL*8 PM,UM,RIL,RIR
      REAL*8 ALPHA,BETA,PREF,gam

      REAL*8 G1, G2, G3, G4, G5, G6, G7
      REAL*8 PREF1,PREF3
      REAL*8 KCR,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5
      REAL*8 CL, CR, PLL, PRR
      REAL*8 FL, FR, FLD, FRD, UDIFF
      REAL*8 POLD, RHOOLD, PRATIO, CHANGE
      REAL*8 AK, BK, QRT, PR1, PR2, DRP
      REAL*8 TOLPRE, VTOL
      REAL*8 DEBUG

C     TOLPRE is precision tolerance for convergence of newton-raphson
C     iteration and NRITER is max # of newton-raphson iterations
c     TOLPRE = 1.0E-03
c     NRITER = 100

C      WRITE(*,*) DL,UL,PL,DR,UR,PR
C      WRITE(*,*) ALPHA,BETA,PREF,gam

C     INPUT: LEFT = GAS, RIGHT = TAIT
C
C     NATURE = 1 --> RAREFACTION/SHOCK
C              2 --> SHOCK/RAREFACTION
C              3 --> SHOCK/SHOCK
C              4 --> RAREFACTION/RAREFACTION
C  

*     
C     Constants for PG EOS
*     
      G1 = (gam - 1.0)/(2.0*gam)
      G2 = (gam + 1.0)/(2.0*gam)
      G3 = 2.0/(gam - 1.0)
      G4 = 2.0/(gam + 1.0)
      G5 = (gam - 1.0)/(gam + 1.0)
      G6 = gam - 1.0
      G7 = gam/(gam  -1.0)
*     
C     Constants for Tait EOS
*     
      PREF   = -PREF
      PREF1  = (BETA-1.0)/(2.0*BETA)
      PREF3  = 2.0/(BETA-1.0)
      PLL    = PL + PREF
      PRR    = PR + PREF
*     
C     Compute sound speeds
*     
      CL = DSQRT(gam*PL/DL)
      CR = DSQRT(ALPHA*BETA*DR**(BETA-1.0))
      IF(CL.LE.0.0 .OR. CR.LE.0.0)THEN
         WRITE(*,*) 'CL = ', CL, ' CR = ',CR
      ENDIF
      
C     To determine nature of problem
      NATURE = 0
      KCR    = 1.0 - (PREF/PRR)**PREF1
      PARAM1 = CL*G3*((PR/PL)**G1 - 1.0)
      PARAM2 = DSQRT((1.0-(PLL/PRR)**(-1.0/BETA))*(PL-PR)/DR)
      PARAM3 = CR*PREF3*( (PLL/PRR)**PREF1 - 1.0 )
      PARAM4 = DSQRT(G3*PL/DL)*(PR/PL-1.0)/DSQRT(1+PR/(G5*PL))
      PARAM5 = -(G3*CL+PREF3*CR*KCR)
C      WRITE(*,*) PARAM1, PARAM2, PARAM3, PARAM4, PARAM5 

      UDIFF = UL - UR

      IF(PL.GT.PR .AND. PARAM1.LE.UDIFF .AND. UDIFF.LE.PARAM2)THEN
         NATURE = 1
      ENDIF

      IF(PL.LT.PR .AND. PARAM3.LE.UDIFF .AND. UDIFF.LE.PARAM4)THEN
         NATURE = 2
      ENDIF

      IF(PL.GT.PR .AND. UDIFF.GT.PARAM2)THEN
         NATURE = 3
      ENDIF

      IF(PL.LT.PR .AND. UDIFF.GT.PARAM4)THEN
         NATURE = 3
      ENDIF

      IF(PL.GT.PR .AND. PARAM5.LT.UDIFF .AND. UDIFF.LT.PARAM1)THEN
         NATURE = 4
      ENDIF

      IF(PL.LT.PR .AND. PARAM5.LT.UDIFF .AND. UDIFF.LT.PARAM3)THEN
         NATURE = 4
      ENDIF

      IF(UDIFF.LT.PARAM5)THEN
         NATURE = 4
      ENDIF


      IF(PL.EQ.PR)THEN
         IF(UL.GT.UR)THEN
            NATURE = 3
         ELSE
            NATURE = 4
         ENDIF
      ENDIF

      IF(UDIFF.EQ.0.0)THEN
         IF(PL.GT.PR)THEN
            NATURE = 1
         ELSE
            NATURE = 2
         ENDIF
      ENDIF

      IF(NATURE.EQ.0)THEN
        WRITE(*,*) 'NATURE OF RIEMANN PROBLEM HAS NOT BEEN DETERMINED'
      ENDIF
C      WRITE(*,*) 'NATURE = ', NATURE

C     NATURE OF THE RIEMANN PROBLEM HAS BEEN DETERMINED


C     SOLVER FOR EACH CASE

*     

      UDIFF = -UDIFF


C     INITIALIZE NEWTON ITERATION
      
      IF(NATURE.EQ.1 .OR. NATURE.EQ.2)THEN
         PM = 0.5*(PL+PR)
      ELSEIF(NATURE.EQ.3)THEN
         PM = MAX(PL,PR)
      ELSEIF(NATURE.EQ.4)THEN
         PM = MIN(PL,PR)
      ENDIF
      POLD = PM
         
*     
*     
      DO 10 I = 1, NRITER
*     
C     GAS MEDIUM RELATIONS
*     
         IF(NATURE.EQ.1 .OR. NATURE.EQ.4)THEN
*     
C     Rarefaction wave
*     
            PRATIO = POLD/PL
            FL     = G3*CL*(PRATIO**G1 - 1.0)
            FLD    = PRATIO**(-G2)/(DL*CL)
         ELSE
*     
C     Shock wave
*     
            AK   = G4/DL
            BK   = G5*PL
            QRT  = SQRT(AK/(BK + POLD))
            FL   = (POLD - PL)*QRT
            FLD  = (1.0 - 0.5*(POLD - PL)/(BK + POLD))*QRT
         ENDIF
*     
C     WATER MEDIUM RELATIONS
*     
         IF(NATURE.EQ.2 .OR. NATURE.EQ.4)THEN
*     
C     Rarefaction wave
*     
            PRATIO = (POLD  +PREF)/(PR  +PREF)
            FR   = PREF3*CR*(PRATIO**PREF1 -1.0)
            PR1  = (POLD  +PREF)**(PREF1-1.0)
            PR2  = (PR    +PREF)**PREF1
            FRD  = CR*PR1/(PR2 *BETA)
         ELSE
*     
C     Shock wave
*     
            RHOOLD = ((POLD  +PREF)/ALPHA)**(1/BETA)
            FR  = DSQRT((POLD  -PR)*(RHOOLD-DR)/(DR*RHOOLD))
            DRP = RHOOLD**(1.0-BETA)/(ALPHA*BETA)
            IF(FR.NE.0.0)THEN
              FRD = (RHOOLD*(RHOOLD  -DR) +
     &              (POLD  -PR)*DR*DRP)/(DR*RHOOLD**2)
              FRD = FRD*0.5/FR
            ELSE
              FRD = 1.0/(DR*CR)
            ENDIF

         ENDIF

         PM     = POLD - (FL + FR + UDIFF)/(FLD + FRD)
         CHANGE = 2.0*ABS((PM- POLD)/(PM+ POLD))
         IF(CHANGE.LE.TOLPRE)GOTO 20
         IF(PM.LT.0.0)PM = TOLPRE
         POLD  = PM
*     
 10   CONTINUE
*     
      WRITE(*,*) ' *** Warning: '
      WRITE(*,*) 'Newton-Raphson reached max num. iterations ', NRITER
      WRITE(*,*) 'without converging to the desired tolerance', TOLPRE
      WRITE(6,*) PM, PL, PR
      WRITE(*,*) 'INPUT RIEMANN (DL,UL,PL,DR,UR,PR)'
      WRITE(*,*) DL,UL,PL,DR,UR,PR
      WRITE(*,*) 'OUTPUT (PM,UM,RIL,RIR)'
      WRITE(*,*) PM,UM,RIL,RIR
      WRITE(*,*) 'ALPHA,PREF = ', ALPHA,PREF
      WRITE(*,*) ' *** '
*     
 20   CONTINUE

      UM = 0.5*(UL + UR + FR - FL)

      IF (NATURE.EQ.1 .OR. NATURE.EQ.4) THEN
         RIL  = DL*(PM/PL)**(1.0/gam)
      ELSE
         RIL  = DL*(G7*PM  - 0.5*(PM  -PL))/(G7*PL  +0.5*(PM -PL))
      ENDIF

      IF (NATURE.EQ.2 .OR. NATURE.EQ.4) THEN
         RIR  = DR*((PM  +PREF)/(PR  +PREF))**(1/BETA)
      ELSE
         RIR  = DR*((PM  +PREF)/(PR  +PREF))**(1/BETA)
      ENDIF




C     CHECK SOLUTIONS

      IF(NATURE.EQ.1)THEN
c         IF(UM.GE.MAX(UL,UR) .AND. PR.LE.PM .AND. PM.LE.PL)THEN
         IF(PR.LE.PM .AND. PM.LE.PL)THEN
         ELSE
            WRITE(6,*) 'ERROR TYPE 1'
         ENDIF
      ENDIF

      IF(NATURE.EQ.2)THEN
c         IF(UM.LE.MIN(UL,UR) .AND. PL.LE.PM .AND. PM.LE.PR)THEN
         IF(PL.LE.PM .AND. PM.LE.PR)THEN
         ELSE
            WRITE(6,*) 'ERROR TYPE 2'
         ENDIF
      ENDIF

      IF(NATURE.EQ.3)THEN
c         IF(PM.GE.MAX(PL,PR) .AND. UR.LE.UM .AND. UM.LE.UL)THEN
         IF(PM.GE.MAX(PL,PR))THEN
         ELSE
            WRITE(6,*) 'ERROR TYPE 3'
         ENDIF
      ENDIF

      IF(NATURE.EQ.4)THEN
c         IF(PM.LE.MIN(PL,PR) .AND. UL.LE.UM .AND. UM.LE.UR)THEN
         IF(PM.LE.MIN(PL,PR))THEN
         ELSE
            WRITE(6,*) 'ERROR TYPE 4'
         ENDIF
      ENDIF

      END
