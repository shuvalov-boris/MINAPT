C   IMSL ROUTINE NAME   - ICSSG
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - DGC/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL SUBROUTINE ICSSCV
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - ICSSF,ICSSH
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ICSSG  (Z,Y,L,V,N,H,RO,AWK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
      REAL               Z(1),Y(1),L(1),V(1),H,RO,AWK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            ITER
      REAL               DELTA,ERR,GF1,GF2,GF3,GF4,R1,R2,R3,R4,TAU,TOL
C                                  TOL = MACHINE PRECISION PARAMETER
      DATA               TOL/0.9537E-06/
      DATA               TAU/1.618034/
C                                  CALCULATE THE MINIMUM OF THE
C                                    CROSS VALIDATION FUNCTION.
C                                  FIRST EXECUTABLE STATEMENT
      R1 = 1.0
      R2 = 10.0
      CALL ICSSH (Z,Y,L,V,N,H,R2,GF2,AWK)
    5 CALL ICSSH (Z,Y,L,V,N,H,R1,GF1,AWK)
      IF (GF1.GT.GF2) GO TO 10
      R2 = R1
      GF2 = GF1
      R1 = R1*.1
      GO TO 5
   10 R3 = 10.*R2
   15 CALL ICSSH (Z,Y,L,V,N,H,R3,GF3,AWK)
      IF (GF3.GT.GF2) GO TO 20
      R2 = R3
      GF2 = GF3
      R3 = 10.*R3
      GO TO 15
   20 R2 = R3
      GF2 = GF3
      ITER = 0
      DELTA = (R2-R1)/TAU
      R4 = R1+DELTA
      R3 = R2-DELTA
      CALL ICSSH (Z,Y,L,V,N,H,R3,GF3,AWK)
      CALL ICSSH (Z,Y,L,V,N,H,R4,GF4,AWK)
   25 IF (GF3-GF4) 30,30,35
   30 R2 = R4
      GF2 = GF4
      R4 = R3
      GF4 = GF3
      DELTA = DELTA/TAU
      R3 = R2-DELTA
      CALL ICSSH (Z,Y,L,V,N,H,R3,GF3,AWK)
      ITER = ITER+1
      GO TO 40
   35 R1 = R3
      GF1 = GF3
      R3 = R4
      GF3 = GF4
      DELTA = DELTA/TAU
      R4 = R1+DELTA
      CALL ICSSH (Z,Y,L,V,N,H,R4,GF4,AWK)
      ITER = ITER+1
   40 ERR = (R2-R1)/(R1+R2)
      IF (ERR.GT.100.0*TOL) GO TO 25
      RO = (R1+R2)*.5
      RETURN
      END
