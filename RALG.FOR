
C**********************************************************************C
C*                                                                    *C
C*                                                                    *C
C*        R-A��OP�TM M�H�M��A��� �E�EBO� ��HK��� C PACT��EH�EM        *C
C*    �POCTPAHCTBA B HA�PAB�EH�� PA�HOCT� �B�X �OC�E�OBATE��H�X       *C
C*                    �PA��EHTOB                                      *C
C*                                                                    *C
C* �APAMETP�: KOD(6) - MACC�B �E�O��C�EHH�X �APAMETPOB,�C�O�����TC�:  *C
C*            KOD(1) - MAKC�MA��HOE ��C�O XPAH�M�X BEKTOPOB B MACC�BE X
C*            KOD(2) - TEK��EE ��C�O BEKTOPOB B MACC�BE X             *C
C*            KOD(5)=MITER=KR - MAK �MA��HOE ��C�O �TEPA���           *C
C*            KOD(4)=MCALC - MAKC�MA��HOE ��C�O B���C�EH�� ��HK���    *C
C*            KOD(6)=BEGMET - �P��HAK BO�O�HOB�EH�� METO�A            *C
C*            X - MACC�B O�T�M���P�EM�X �EPEMEHH�X PA�MEPHOCT� N      *C
C*            FX - �HA�EH�E �E�EBO� ��HK���                           *C
C*            G - �PA��EHT ��HK���                                    *C
C*            GM=Z1 - �O��CT�M�E M�H�MA��H�E BE����H� KOM�OHEHT       *C
C*                               �PA��EHTA                            *C
C*            DS - HE �C�O����ETC�                                    *C
C*            PARM(8) - MACC�B �E�CTB�TE��H�X �APAMETPOB,�C�O�����TC�:*C
C*            PARM(2)=STPDIR=H - HA�A��H�� �A� �O HA�PAB�EH��         *C
C*            PARM(3)=MAXSTP=Z3    - MAKC�MA��H�� �A�                 *C
C*            PARM(4)=ACCDX=MINH - M�H�MA��H�� �A� �O HA�PAB�EH��     *C
C*            PARM(5)=CONSAD=ALPHA-KO������EHT PACT��EH�� �POCTPAHCTBA*C
C*            PARM(7)=WAITER=Z2 - TO�HOCT� �PE�CTAB�EH�� ��C�A        *C
C*            EC�� TEK��EE �HA�EH�E PA�HOCT� �PA��EHTOB MEH��E PARM(7)*C
C*            ��� H>PARM(3), TO �PO�CXO��T PECTAPT A��OP�TMA          *C
C*            F - �E�EBA� ��HK���                                     *C
C*            MIN - HE �C�O����ETC�                                   *C
C*            OUT - �PO�PAMMA PACC�ETA KP�TEP�EB OKOH�AH��            *C
C*            A - PA�O��� MACC�B PA�MEPHOCT�� 5N                      *C
C*            B1 - MACC�B PA�MEPHOCT�� N*N                            *C
C* �AME�AH�E: �PO�PAMMA O�HOMEPHO� M�H�M��A��� BCTPOEHA B A��OP�TM    *C
C*            IRC - KO� �ABEP�EH�� A��OP�TMA:                         *C
C*            IRC=0 - �OCT��H�T M�H�MA��H�� �A� �O HA�PAB�EH��        *C
C*            IRC=2 - �C�EP�AH PEC�PC KO����ECTBA �TEPA���            *C
C*            IRC=3 - �C�EP�AH PEC�PC KO����ECTBA B���C�EH�� ��HK���  *C
C*            IRC=4 - KOM�OHEHTA �PA��EHTA �OCT���A CBOE�O            *C
C*                        M�H�MA��HO�O �HA�EH��                       *C
C*            IRC=7 - A��OP�TM HE MO�ET HA�T� M�H�M�M C �A�AHHO�      *C
C*                    TO�HOCT�� ��-�A B���H�� O���OK OKP���EH��.      *C
C*                    PEKOMEH��ETC� BOC�O���OBAT�C� A��O��TMOM DRALG  *C
C*            B XO�E PA�OT� A��OP�TMA MEH��TC� �APAMETP�:             *C
C*            KOD(4),KOD(5),PARM(2)                                   *C
C*                                                                    *C
C*            B���BA�TC� �PO�PAMM�: F,OUT                             *C
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
C*       KOD(6)=1 - �PO�O��EH�E �PEBAHO� PA�OT�
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
C*       �O�C�ET ��C�A B���C�EH�� ��HK���
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
C*        �O�C�ET ��C�A �TEPA���
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
C*        �O�C�ET KO����ECTBA B���C�EH�� �E�EBO� ��HK���
      KOD(4)=KOD(4)-1
C
      DK1=DK1+1.
      IF (F2.LT.F1) GO TO 491
      IF((KOD(1).LE.1).OR.(K1.EQ.2)) GO TO 50
C***    �AHECEH�E �����X O�T�M���P�EM�X �EPEMEHH�X B MACC�B X
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
C*          �O�C�ET ��C�A B���C�EH�� ��HK���
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
C*          �O�C�ET ��C�A �TEPA���
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
C*       �POBEPKA HA O���K� OKP���EH��
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
