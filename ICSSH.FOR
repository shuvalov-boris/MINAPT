C   IMSL ROUTINE NAME   - ICSSH
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
C   REQD. IMSL ROUTINES - ICSSF
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
      SUBROUTINE ICSSH  (Z,Y,L,V,N,H,RO,CV,AWK)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N
      REAL               Z(1),Y(1),L(1),V(1),H,RO,CV,AWK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I,N2
      REAL               CUO,SN,SINFP,SQINF,TAR
C!!!      DATA               SINFP/0.7237E+76/
C                                  SQINF = SQRT(LARGEST POSITIVE REAL)
C                                  FIRST EXECUTABLE STATEMENT
c*       заменил DATA на вычисление. слишком большой поядок числа !!!!!
      a=1.
      sinfp = huge(a)
      SQINF = SQRT(SINFP)
C                                  CALCULATE THE VALUE OF THE CROSS
C                                    VALIDATION FUNCTION
      DO 5 I=1,N
    5 Y(I) = Z(I)
      CALL ICSSF (N,RO,H,Y,L,AWK)
      SN = 0.0
      DO 10 I=1,N
         TAR = Y(I)-Z(I)
         SN = SN+TAR*TAR
   10 CONTINUE
      SN = SN/N
      N2 = N-2
      TAR = 2.0
      CUO = 1.0/RO
      DO 15 I=1,N2
         TAR = TAR+(1.0/(1.0+V(I)*CUO))
   15 CONTINUE
      CV = 1.0-TAR/N
      CV = CV*CV
      IF (CV.EQ.0.0) GO TO 20
      CV = SN/CV
      RETURN
   20 CV = SQINF
      RETURN
      END
