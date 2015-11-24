
C**********************************************************************C
C*                                                                    *C
C*                                                                    *C
C*        R-AËÃOPÈTM MÈHÈMÈÇAÖÈÈ ÖEËEBOÉ ÔÓHKÖÈÈ C PACTßÆEHÈEM        *C
C*    ÏPOCTPAHCTBA B HAÏPABËEHÈÈ PAÇHOCTÈ ÄBÓX ÏOCËEÄOBATEËÜHÛX       *C
C*                    ÃPAÄÈEHTOB                                      *C
C*                                                                    *C
C* ÏAPAMETPÛ: KOD(6) - MACCÈB ÖEËO×ÈCËEHHÛX ÏAPAMETPOB,ÈCÏOËÜÇÓÞTCß:  *C
C*            KOD(1) - MAKCÈMAËÜHOE ×ÈCËO XPAHÈMÛX BEKTOPOB B MACCÈBE X
C*            KOD(2) - TEKÓÙEE ×ÈCËO BEKTOPOB B MACCÈBE X             *C
C*            KOD(5)=MITER=KR - MAK ÈMAËÜHOE ×ÈCËO ÈTEPAÖÈÉ           *C
C*            KOD(4)=MCALC - MAKCÈMAËÜHOE ×ÈCËO BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ    *C
C*            KOD(6)=BEGMET - ÏPÈÇHAK BOÇOÁHOBËEHÈß METOÄA            *C
C*            X - MACCÈB OÏTÈMÈÇÈPÓEMÛX ÏEPEMEHHÛX PAÇMEPHOCTÈ N      *C
C*            FX - ÇHA×EHÈE ÖEËEBOÉ ÔÓHKÖÈÈ                           *C
C*            G - ÃPAÄÈEHT ÔÓHKÖÈÈ                                    *C
C*            GM=Z1 - ÄOÏÓCTÈMÛE MÈHÈMAËÜHÛE BEËÈ×ÈHÛ KOMÏOHEHT       *C
C*                               ÃPAÄÈEHTA                            *C
C*            DS - HE ÈCÏOËÜÇÓETCß                                    *C
C*            PARM(8) - MACCÈB ÄEÓCTBÈTEËÜHÛX ÏAPAMETPOB,ÈCÏOËÜÇÓÞTCß:*C
C*            PARM(2)=STPDIR=H - HA×AËÜHÛÉ ØAÃ ÏO HAÏPABËEHÈÞ         *C
C*            PARM(3)=MAXSTP=Z3    - MAKCÈMAËÜHÛÉ ØAÃ                 *C
C*            PARM(4)=ACCDX=MINH - MÈHÈMAËÜHÛÉ ØAÃ ÏO HAÏPABËEHÈÞ     *C
C*            PARM(5)=CONSAD=ALPHA-KOÝÔÔÈÖÈEHT PACTßÆEHÈß ÏPOCTPAHCTBA*C
C*            PARM(7)=WAITER=Z2 - TO×HOCTÜ ÏPEÄCTABËEHÈß ×ÈCËA        *C
C*            ECËÈ TEKÓÙEE ÇHA×EHÈE PAÇHOCTÈ ÃPAÄÈEHTOB MEHÜØE PARM(7)*C
C*            ÈËÈ H>PARM(3), TO ÏPOÈCXOÄÈT PECTAPT AËÃOPÈTMA          *C
C*            F - ÖEËEBAß ÔÓHKÖÈß                                     *C
C*            MIN - HE ÈCÏOËÜÇÓETCß                                   *C
C*            OUT - ÏPOÃPAMMA PACC×ETA KPÈTEPÈEB OKOH×AHÈß            *C
C*            A - PAÁO×ÈÉ MACCÈB PAÇMEPHOCTÜÞ 5N                      *C
C*            B1 - MACCÈB PAÇMEPHOCTÜÞ N*N                            *C
C* ÇAME×AHÈE: ÏPOÃPAMMA OÄHOMEPHOÉ MÈHÈMÈÇAÖÈÈ BCTPOEHA B AËÃOPÈTM    *C
C*            IRC - KOÄ ÇABEPØEHÈß AËÃOPÈTMA:                         *C
C*            IRC=0 - ÄOCTÈÃHÓT MÈHÈMAËÜHÛÉ ØAÃ ÏO HAÏPABËEHÈÞ        *C
C*            IRC=2 - ÈC×EPÏAH PECÓPC KOËËÈ×ECTBA ÈTEPAÖÈÉ            *C
C*            IRC=3 - ÈC×EPÏAH PECÓPC KOËËÈ×ECTBA BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ  *C
C*            IRC=4 - KOMÏOHEHTA ÃPAÄÈEHTA ÄOCTÈÃËA CBOEÃO            *C
C*                        MÈHÈMAËÜHOÃO ÇHA×EHÈß                       *C
C*            IRC=7 - AËÃOPÈTM HE MOÆET HAÉTÈ MÈHÈMÓM C ÇAÄAHHOÉ      *C
C*                    TO×HOCTÜÞ ÈÇ-ÇA BËÈßHÈß OØÈÁOK OKPÓÃËEHÈß.      *C
C*                    PEKOMEHÄÓETCß BOCÏOËÜÇOBATÜCß AËÃOÏÈTMOM DRALG  *C
C*            B XOÄE PAÁOTÛ AËÃOPÈTMA MEHßÞTCß ÏAPAMETPÛ:             *C
C*            KOD(4),KOD(5),PARM(2)                                   *C
C*                                                                    *C
C*            BÛÇÛBAÞTCß ÏPOÃPAMMÛ: F,OUT                             *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      SUBROUTINE RALG(KOD,N,X,FX,G,GM,DS,PARM,F,MIN,OUT,A,B1,IRC)
      DIMENSION KOD(6),X(N,1),FX(1),G(N),GM(N),DS(N),PARM(8),A(N,5),
     *B1(N,N),FXX(2)
      REAL NH/3./
      K=1
      ALPHA=PARM(5)
      K5=0
      KR=KOD(5)
      H=PARM(2)
      Z2=PARM(7)
      Z3=PARM(3)
      Z1=GM(1)
      B2=1./ALPHA-1
      B2O=B2
C*       KOD(6)=1 - ÏPOÄOËÆEHÈE ÏPEBAHOÉ PAÁOTÛ
      IF(KOD(6).EQ.1) GO TO 1111
  100 DO 60 I=1,N
      DO 61 J=1,N
   61 B1(I,J)=0.0
   60 B1(I,I)=1.0
C
      CALL  F(X,FX,G,3)
      F1=FX(1)
      KOD(4)=KOD(4)-1
C
      DD=0
      DO 62 I=1,N
   62 DD=DD+G(I)**2
      DD=SQRT(DD)
      IF(DD.LT.Z1) IRC=4
      IF (DD.LT.Z1) GO TO 200
      DO 23 I=1,N
      A(I,3)=G(I)
   23 A(I,5)=X(I,1)
      K2=0
      K1=1
  104 D=PARM(2)/DD
      DO 10 I=1,N
   10 X(I,1)=X(I,1)-G(I)*D
C
      CALL  F(X,FX,G,1)
C*       ÏOÄC×ET ×ÈCËA BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ
      KOD(4)=KOD(4)-1
      F2=FX(1)
C
      IF (F2.GE.F1) GO TO 11
      K1=0
      F1=F2
      K2=K2+1
      IF (K2.LT.5) GO TO 104
      PARM(2)=PARM(2)*10.
      K2=0
      GO TO 104
   11 IF (K1.NE.1) GO TO 12
      DO 13 I=1,N
   13 X(I,1)=A(I,5)
      PARM(2)=PARM(2)/2.
      GO TO 104
   12 D=0
      DO 25 I=1,N
   25 D=D+(X(I,1)-A(I,5))**2
      D=SQRT(D)
      PARM(2)=D/NH
C*        ÏOÄC×ET ×ÈCËA ÈTEPAÖÈÉ
      FXX(1)=F2
      FXX(2)=F1
      DO 1000 I=1,N
      A(I,2)=A(I,5)
1000  A(I,1)=X(I,1)
      CALL  OUT(KOD,N,A,FXX,DS,PARM,IRC)
      IF(IRC.EQ.2)    GO TO 200
1111  IRC=0
      IF(PARM(2).LT.PARM(4)) GO TO 200
C
      CALL  F(X,FX,G,3)
      F1=FX(1)
      KOD(4)=KOD(4)-1
      IF(KOD(4).LE.0) IRC=3
      IF(KOD(4).LE.0) GO TO 200
      DD=0
      DO 26 I=1,N
   26 DD=DD+G(I)**2
      DD=SQRT(DD)
      IF (DD.LT.Z1) IRC=4
      IF (DD.LT.Z1) GO TO 200
  101 K=K+1
      IFLG=0
      IF(A(1,3).EQ.375.) IFLG=1
C     IF(IFLG.EQ.1) GO TO 300
      DO 31 I=1,N
      A(I,4)=0.
      DO 31 J=1,N
   31 A(I,4)=A(I,4)+B1(J,I)*G(J)
  300 CONTINUE
      IF(IFLG.EQ.1) B2=A(2,3)
      DO 32 I=1,N
      A(I,3)=A(I,3)-A(I,4)
      IF(IFLG.EQ.1) A(I,3)=G(N+I)
   32 CONTINUE
      D=0
      DO 33 I=1,N
   33 D=D+A(I,3)**2
      D=SQRT(D)
      IF (D.GE.Z2) GO TO 34
      DO 35 I=1,N
   35 A(I,3)=A(I,4)
      GO TO 102
34    III=1
      DO 36 I=1,N
   36 A(I,3)=A(I,3)/D
      DO 37 I=1,N
      D=0
      DO 38 J=1,N
   38 D=D+B1(I,J)*A(J ,3)
      D=D*B2
      DO 39 J=1,N
   39 B1(I,J)=B1(I,J)+A(J,3)*D
   37 CONTINUE
      D=0
      DO 40 I=1,N
   40 D=D+A(I,3)*A(I,4)
      D=B2*D
      DO 41 I=1,N
   41 A(I,3)=A(I,3)*D+A(I,4)
  102 D=0
      DO 42 I=1,N
   42 D=D+A(I,3)**2
      D=SQRT(D)
      DO 43 I=1,N
      G(I)=0
      DO 44 J=1,N
   44 G(I)=G(I)+B1(I,J)*A(J,3)
   43 CONTINUE
      IF (D.GE.Z2) GO TO 45
      D=0
      DO 46 I=1,N
   46 D=D+(X(I,1)-A(I,5))**2
      D=SQRT(D)
      PARM(2)=D/NH
      GO TO 100
45    III=2
      DO 47 I=1,N
   47 G(I)=G(I)/D
      DO 48 I=1,N
   48 A(I,5)=X(I,1)
      DK1=0.
      TIMEN=1.
      PH=PARM(2)
      K1=1
103   IF(KOD(4).LE.0) IRC=3
      IF(KOD(4).LE.0) GO TO 200
      DO 49 I=1,N
   49 X(I,1)=X(I,1)-PH*G(I)
      K1=K1+1
      CALL  F(X,FX,G,1)
      F2=FX(1)
C
C*        ÏOÄC×ET KOËËÈ×ECTBA BÛ×ÈCËEHÈÉ ÖEËEBOÉ ÔÓHKÖÈÈ
      KOD(4)=KOD(4)-1
C
      DK1=DK1+1.
      IF (F2.LT.F1) GO TO 491
      IF((KOD(1).LE.1).OR.(K1.EQ.2)) GO TO 50
C***    ÇAHECEHÈE ËÓ×ØÈX OÏTÈMÈÇÈPÓEMÛX ÏEPEMEHHÛX B MACCÈB X
      CALL ROUTE(KOD,N,X,FX)
      FX(2)=F1
      DO 490 I=1,N
490   X(I,2)=X(I,1)+PH*G(I)
      GO TO 50
491   F1=F2
      IF((K1-K1/3*3).NE.0) GO TO 103
C     PH=PH*2.
C     DK1=DK1*2.
C     TIMEN=TIMEN+TIMEN
      GO TO 103
   50 PARM(2)=DK1*PARM(2)/NH
C
      K1=DK1
C     CALL    F(X,FX,G,2)
C     B2=B2O
C     KOD(4)=KOD(4)-N
      F1=F2
C
C*          ÏOÄC×ET ×ÈCËA BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ
      IF(KOD(4).LE.0) IRC=3
      IF(KOD(4).LE.0) GO TO 200
C
C     DD=0
C     DO 51 I=1,N
C  51 DD=DD+G(I)**2
C     DD=SQRT(DD)
C     IF(DD.LT.Z1) IRC=4
C     IF (DD.LT.Z1) GO TO 200
      D=0
      DO 52 I=1,N
   52 D=D+(X(I,1)-A(I,5))**2
      D=SQRT(D)
      IF (PARM(2).LE.Z3)GOTO53
      PARM(2)=D/NH
      B2=B2O
      GO TO 100
C*          ÏOÄC×ET ×ÈCËA ÈTEPAÖÈÉ
53    DO 1001 I=1,N
      A(I,2)=A(I,5)
1001  A(I,1)=X(I,1)
      FXX(1)=F2
      FXX(2)=F1
      CALL  OUT(KOD,N,A,FXX,DS,PARM,IRC)
      IF(IRC.EQ.2) GO TO 200
      IRC=0
      CALL  F(X,FX,G,2)
      B2=B2O
      DD=0.
      DO 51 I=1,N
51    DD=DD+G(I)*G(I)
      DD=SQRT(DD)
      IF(DD.LT.Z1) IRC=4
      IF(DD.LT.Z1) GO TO 200
      IF(PARM(2).LT.PARM(4)) GO TO 200
      IF(PARM(2).LE.D*1000.) GO TO 101
C*       ÏPOBEPKA HA OØÈÁKÈ OKPÓÃËEHÈß
      IRC=7
      IF(D.LT.1.E-30) GO TO 200
C     HMIN=AMAX1(1.E-30,PARM(4))
C     IF(D.LT.HMIN) D=HMIN
      D=PARM(2)/D
      DO 54 I=1,N
      DO 55 J=1,N
   55 B1(I,J)=D*B1(I,J)
   54 CONTINUE
C     III=4
      PARM(2)=PARM(2)/D
      GO TO 101
  200 RETURN
C     DEBUG INIT(B1,B2)
      END
