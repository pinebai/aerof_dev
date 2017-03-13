      SUBROUTINE ERIEMANNWW(DL,UL,PL,DR,UR,PR,
     &                      PM,UM,RIL,RIR,
     &                      ALPHAL,BETAL,PREFL,
     &                      ALPHAR,BETAR,PREFR,errcod,
     &                      pcutl,pcutr,rcutl,rcutr,
     &                      TOLPRE,NRITER)
*
*----------------------------------------------------------------------*
*                                                                      *
C     Exact Riemann Solver for the Time-Dependent                      *
C     One Dimensional Euler Equations                                  *
*                                                                      *
*----------------------------------------------------------------------*
*
      IMPLICIT NONE
*
C     Declaration of variables:
*
      INTEGER I, NRITER,errcod
      REAL*8  TOLPRE
*
      REAL*8  DL, UL, PL, DR, UR, PR
      REAL*8  ALPHAL,BETAL,PREFL
      REAL*8  ALPHAR,BETAR,PREFR
      REAL*8  RIL,RIR,UM,PM

      REAL*8  CL,CR
      REAL*8  GL1,GL2,GL3,GR1,GR2,GR3
      REAL*8  PLL,PRL,PLR,PRR
      REAL*8  PRATIO,POLD,UDIFF,POLDL,POLDR
      REAL*8  FL, FLD, DLP, FR, FRD, DRP
      REAL*8  AAPR1, AAPR2, RHOOLD, CHANGE
      REAL*8 pcutl,rcutl,pcutr,rcutr
      REAL*8 pcut

      errcod = 0
*********************************************************************

      PL = PREFL + ALPHAL*DL**BETAL
      PR = PREFR + ALPHAR*DR**BETAR

      PCUT = MAX(pcutl,pcutr)
*********************************************************************
      
c     TOLPRE  = 1.E-3
c     NRITER  = 100

      IF (PL.LT.0.0 .OR. PR.LT.0.0) WRITE(*,*) 'NEGATIVE PRESSURES',
     &                              DL,UL,PL,DR,UR,PR

C     Compute constants of the fluids
      PREFL = -PREFL
      PREFR = -PREFR

      GL1  = (BETAL-1.0)/(2.0*BETAL)
      GL2  = GL1 - 1.0
      GL3  = 2.0/(BETAL-1.0)
      GR1  = (BETAR-1.0)/(2.0*BETAR)
      GR2  = GR1 - 1.0
      GR3  = 2.0/(BETAR-1.0)

      PLL    = PL + PREFL
      PRL    = PR + PREFL
      PRR    = PR + PREFR
      PLR    = PL + PREFR

C     Compute sound speeds
      CL = DSQRT(ALPHAL*BETAL*DL**(BETAL-1.0))
      CR = DSQRT(ALPHAR*BETAR*DR**(BETAR-1.0))

C     Initialization of Newton-Raphson iterations
      PM    = 0.5*(PR  +PL)
      POLD  = PM
      UDIFF = UR - UL

      DO 10 I = 1, NRITER

        POLDL = PREFL + POLD
        POLDR = PREFR + POLD

C  LEFT FLUID MEDIUM RELATIONS
        IF(POLD.LE.PL)THEN
C      Rarefaction wave
          PRATIO = POLDL/PLL
          FL     = GL3*CL*(PRATIO**GL1 -1.0)
          AAPR1  = POLDL**GL2
          AAPR2  = PLL**GL1
          FLD    = CL*AAPR1/(AAPR2 *BETAL)
           
        ELSE
C      Shock wave
          RHOOLD = (POLDL/ALPHAL)**(1/BETAL)
          FL  = DSQRT((POLD  -PL)*(RHOOLD-DL)/(DL*RHOOLD))
          DLP = RHOOLD**(1.0-BETAL)/(ALPHAL*BETAL)
          IF(FL.NE.0.0)THEN
            FLD = (RHOOLD*(RHOOLD  -DL) +
     &            (POLD-PL)*DL*DLP)/(DL*RHOOLD**2)
            FLD = FLD*0.5/FL
          ELSE
            FLD = 1.0/(DL*CL)
          ENDIF

        ENDIF

C  RIGHT FLUID MEDIUM RELATIONS

        IF(POLD.LE.PR)THEN
C      Rarefaction wave
          PRATIO = POLDR/PRR
          FR     = GR3*CR*(PRATIO**GR1 -1.0)
          AAPR1  = POLDR**GR2
          AAPR2  = PRR**GR1
          FRD    = CR*AAPR1/(AAPR2 *BETAR)
           
        ELSE
C      Shock wave
          RHOOLD = (POLDR/ALPHAR)**(1/BETAR)
          FR  = DSQRT((POLD  -PR)*(RHOOLD-DR)/(DR*RHOOLD))
          DRP = RHOOLD**(1.0-BETAR)/(ALPHAR*BETAR)
          IF(FR.NE.0.0)THEN
            FRD = (RHOOLD*(RHOOLD  -DR) +
     &            (POLD-PR)*DR*DRP)/(DR*RHOOLD**2)
            FRD = FRD*0.5/FR
          ELSE
            FRD = 1.0/(DR*CR)
          ENDIF

        ENDIF


        PM     = POLD - (FL + FR + UDIFF)/(FLD + FRD)
        IF (PM.LT.PCUT) PM = PCUT
        CHANGE = 2.0*ABS((PM- POLD)/(PM+ POLD))
        IF(CHANGE.LE.TOLPRE)GOTO 20
        IF(PM.LT.0.0)PM = TOLPRE
        POLD  = PM
*
 10   CONTINUE
*
      WRITE(*,*) 'Newton-Raphson reached max num. iterations ', NRITER
      WRITE(*,*) 'without converging to the desired tolerance', TOLPRE
      WRITE(*,*) PM, PL, PR
      WRITE(*,*) 'INPUT RIEMANN (DL,UL,PL,DR,UR,PR)'
      WRITE(*,*) DL,UL,PL,DR,UR,PR
      WRITE(*,*) 'OUTPUT (PM,UM,RIL,RIR)'
      WRITE(*,*) PM,UM,RIL,RIR
*
        errcod = 1
 20   CONTINUE

      UM = 0.5*(UL + UR + FR - FL)

      IF (PM.LE.PL) THEN
        RIL  = DL*((PM  +PREFL)/(PL  +PREFL))**(1/BETAL)
      ELSE
        RIL  = DL*((PM  +PREFL)/(PL  +PREFL))**(1/BETAL)
      ENDIF

      IF (PM.LE.PR) THEN
        RIR  = DR*((PM  +PREFR)/(PR  +PREFR))**(1/BETAR)
      ELSE
        RIR  = DR*((PM  +PREFR)/(PR  +PREFR))**(1/BETAR)
      ENDIF

      IF (RIL.LT.rcutl) RIL = rcutl
      IF (RIR.LT.rcutr) RIR = rcutr
      

      END
