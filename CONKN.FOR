C**********************************************************************C
C*                                                                    *C
C*        ��OP��O��BAH�E ���OB � CO��TCTB����X �M XAPAKTEP�CT�K       *C
C*  B �OP��KE PAC�O�O�EH�� ���OB B�O�� �PA�MEHTA KOHT�PA,             *C
C*                 �PEO�PA�OBAH�E KOOP��HAT                           *C
C*  ���OB B �APAMETP��ECKOE �PE�CTAB�EH�E (PAC�ET �X �O�O�EH�� �O     *C
C*                KOOP��HATE T - MACC�B TK)                           *C
C*                                                                    *C
C*                                                                    *C
C* �APAMETP�: X(N),Y(N),T(N) - X,Y � T - KOOP��HAT� TO�EK KOHT�PA     *C
C*            XK(NKNOT),YK(NKNOT),TK(NKNOT) - TO�E ��� ���OB          *C
C*            NBEG,NEND - HA�A�O � KOHE� �PA�MEHTA KOHT�PA, KOTOP��   *C
C*  �O T-KOOP��HATE �E��T B ��A�O�OHE, �OKP�BA��EMC� ���AM�           *C
C*            DIST(I)-M�H�MA��HOE PACCTO�H�E ���A IABS(NPOS(I)) �O    *C
C*  KOHT�PA (B�XO�)                                                   *C
C*            NPOS(I) - HOMEP ���A, �O�AB�E�O HA I-OE MECTO �OC�E     *C
C*  ��OP��O��BAH�� �O T-KOOP��HATE(MACC�B TK) (B�XO�). EC�� NPOS(I)<0,*C
C*  TO �AHH�� ��E� �E��T BHE �PE�E�OB KOHT�PA � HE �B��ETC� ����A���M *C
C*  K O�HOM� �� KOH�OB KOHT�PA.                                       *C
C*            NUSL,NUSR - ��C�O ���OB C OTP��ATE��H�M� BE����HAM�     *C
C*  NPOS C�EBA � C�PABA OT KPAEB �PA�MEHTA COOTBETCTBEHHO             *C
C*            NDOUBL - ��C�O �AP ���OB, �E�A��X B �PE�E�AX �PA�MEHTA �*C
C*  �ME���X O��HAKOB�� T-KOOP��HAT� �POEK���. EC�� KK ���OB �ME�T     *C
C*  O��HAKOB�� �POEK��� TK, TO K NKNOT �O�AB��ETC� HE 1, A KK-1       *C
C*            IFLG: 0 - ��OP��O��BAT� MACC�B� NPOS,DIST,A KOOP��HAT�  *C
C*  ���OB XK,YK HE ��OP��O��BAT�. 1 - ��OP��O��BAT� BCE �T� MACC�B�   *C
C*                                                                    *C
C* A��OP�TM:  PACC��T�BA�TC� TK(I) - PACCTO�H�� B�O�� �OMAHO�,        *C
C*  �POXO���E� �EPE� TO�K� (X(J),Y(J)),OT TO�K� (X(1),Y(1)) �O TO�K�, *C
C*  HAXO���E�C� HA M�H�MA��HOM PACCTO�H�� OT ���A (XK(I),YK(I))       *C
C*            �� ���OB B�XO����X �A KAKO�-���O KPA� �PA�MEHTA KOHT�PA *C
C*  B���PAETC� ��E� C M�H�MA��H�M PACCTO�H�EM �O KPA�HE� TO�K� (�PAH��*C
C*  H�� ��E�). BE����HA TK ��� HE�O ���ET T(1)-(PACCTO�H�E OT ���A �O  *
C*  TO�K� (X(1),Y(1)), EC�� �TOT KPA� - HA�A�O � T(N)+ (PACCTO�H�E �O *C
C*  TO�K� (X(N),Y(N)), EC�� KOHE�. �P���E ����, �E�A��E BHE KOHT�PA,  *C
C*  B�HOC�T�C� BHA�A�O MACC�BA C COOTBETCTB����M� XAPAKTEP�CT�KAM�,   *C
C*  ���O B KOHE� B TO� �OC�E�OBATE��HOCT�, B KOTOPO� OH� BCTPE�A�TC�  *C
C*  ��HA�A��HO. EC�� TAK�X ���OB HE OKA�A�OC� ��� KAKO�O-���O KPA�,TO *C
C*  PACC��T�BA�TC� BE����H� NBEG �/��� NEND.                          *C
C*                                                                    *C
C* �P�ME�AH��: �TO�� ���E�AT� �O���EH�� �O�H�X �POEK��� ��� �AMKH�TO�O*C
C*  KOHT�PA, HO HE�AMKH�TO�O �PA�MEHTA, C�E��ET BBO��T� TO��KO ����,  *C
C*  COOTBETCTB����E �TOM� �PA�MEHT�.                                  *C
C*             �P� PAC�ETE TK(I) � ��OP��O��BAH�� ���OB �PO�CXO��T    *C
C*  �O�H�� �EPE�OP �E� BC�K�X �X��PEH��                               *C
C*             C�E��ET �OMH�T�, �TO B O��EM C���AE �PAB���HOCT�       *C
C*  �OCTPOEH�� �POEK��� ���OB HE �APAHT�P�ETC�. TPE��ETC� ��A�O�OBA�  *C
C*  �PO�E��PA KOHTPO��                                                *C
C*                                                                    *C
C*     �PO�PAMMA HE �AB�C�T OT C�E����K� TPAHC��TOPA � �BM            *C
C*                                                                    *C
C* B���BAEM�E �O��PO�PAMM�: INCLT - �P�C�TCTB�ET B �TOM �E �A��E      *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      SUBROUTINE CONKN(X,Y,T,N,XK,YK,TK,NKNOT,NBEG,NEND,DIST,NPOS,NUSL,
     *NUSR,NDOUBL,IFLG)
      REAL X(N),Y(N),T(N),XK(1),YK(1),TK(1),XY1(2),XY2(2),DIST(NKNOT)
      INTEGER*2   NPOS(NKNOT)
C     TKMIN=T(N)
C     TKMAX=T(1)
      ICFC=0
      IF(X(1).EQ.X(N).AND.Y(1).EQ.Y(N)) ICFC=1
      DO 100 I=1,NKNOT
      SS=(XK(I)-X(1))**2+(YK(I)-Y(1))**2
      TS=0.
      DO 80 J=2,N
C*        O�PE�E�EH�E M�H�MA��HO�O PACCTO�H�� �O OTPE�KA
C*             (X(I-1),Y(I-1)) - (X(I),Y(I))
       XY1(1)=XK(I)-X(J-1)
       XY1(2)=YK(I)-Y(J-1)
       XY2(1)=X(J)-X(J-1)
       XY2(2)=Y(J)-Y(J-1)
       AL=XY1(1)*XY2(1)+XY1(2)*XY2(2)
       IF(J.EQ.2) AL1=AL
       IF(AL.LE.0.) GO TO 80
C*       TO�KA (XK(I),YK(I)) �E��T "�PABEE" �EBO�O KOH�A OTPE�KA
       INDX=2
       IF(AL.GT.XY2(1)**2+XY2(2)**2) INDX=1
       GO TO (52,54), INDX
C*        M�H�MA��HOE PACCTO�H�E ��� OTPE�KA - �O TO�K� (X(J),Y(J))
52     S2=(XK(I)-X(J))**2+(YK(I)-Y(J))**2
       IF(S2.GT.SS) GO TO 80
       IF(S2.EQ.SS.AND.I.EQ.1) GO TO 80
C*       �PO��O��O �MEH��EH�E TEK��E�O M�H�MA��HO�O PACCTO�H��
       SS=S2
       TS=T(J)
       AL1=1.
       GO TO 80
C*        HEO�XO��MO PACC��T�BAT� TO�K� - �POEK��� ���A I HA
C*          TEK���� OTPE�OK
54    DX=XY2(1)
      DY=XY2(2)
      IND=3
      IF(DY.EQ.0.) IND=2
      IF(DX.EQ.0.) IND=1
      GO TO (55,60,65), IND
C*          DX=0.
55    XC=X(J-1)
      YC=-YK(I)
      GO TO 70
C*          DY=0.
60    XC=XK(I)
      YC=Y(J-1)
      GO TO 70
C*        C��TAEM �O O��E� �OPM��E
65    DIV=1./(DX*DX)+1./(DY*DY)
      XC=((XK(I)/DY+YK(I)/DX)/DY-(Y(J-1)/DY-X(J-1)/DX)/DX)/DIV
      YC=((Y(J-1)/DY-X(J-1)/DX)/DY+(XK(I)/DY+YK(I)/DX)/DX)/DIV
C*       KOHTPO�� HA MEH��EE PACCTO�H�E
70    S2=(XK(I)-XC)**2+(YK(I)-YC)**2
      IF(S2.GT.SS) GO TO 80
      IF(S2.EQ.SS.AND.I.EQ.1) GO TO 80
      AL1=1.
      SS=S2
      TS=T(J-1)+SQRT((XC-X(J-1))**2+(YC-Y(J-1))**2)
80    CONTINUE
      TK(I)=TS
      NPOS(I)=I
      DIST(I)=SS
C     IF(TKMIN.GT.TS) TKMIN=TS
C     IF(TKMAX.LT.TS) TKMAX=TS
C*       �OME�AEM ����, KOTOP�E �E�AT BHE KOHT�PA
      IF((TS.EQ.T(N).AND.INDX.EQ.1).OR.(TS.EQ.0..AND.AL1.LT.0.))
     *NPOS(I)=-NPOS(I)
100   CONTINUE
C
C*     ��OP��O��BAH�E ���OB
C
      NK=NKNOT
      NKNOT=0
      NDOUBL=0
      DO 110 I=1,NK
C*       ��OP��O��BAEM  ����
      NKNOT=NKNOT+1
      TIN=TK(I)
      DIN=DIST(I)
      NIN=NPOS(I)
      XIN=XK(I)
      YIN=YK(I)
      CALL  INCLT(TK,XK,YK,TIN,XIN,YIN,NKNOT,DIST,NPOS,DIN,NIN,NDOUBL,
     *IFLG,NK)
110   CONTINUE
C
C*        PAC�ET NBEG � NEND
C
      NBEG=1
      NEND=N
      NUSL=0
      NUSR=0
      IF(ICFC.EQ.1) GO TO 200
C*         PAC�ET NBEG
      DO 130 I=1,N
      IF(T(I).GE.TK(1)) GO TO 135
130   CONTINUE
135   NBEG=I
140   IF(TK(NKNOT).EQ.T(N)) GO TO 160
C*         PAC�ET NEND
      DO 150 I=1,N
      IF(T(N-I+1).LE.TK(NKNOT)) GO TO 155
150   CONTINUE
155   NEND=N-I+1
C
C*      PAC�ET NUSL,NUSR � TK ��� ���OB BHE KOHT�PA, HO ����A���X
C*          K O�HOM� �� E�O KPAEB
C
160   IF(NPOS(1).GT.0) GO TO 180
C*      �EBEE "�EBO�O" KOH�A KOHT�PA ECT� ����
C*        ��EM ����A����
      NKMIN=1
      DO 170 I=2,NKNOT
      IF(TK(I).GE.0..AND.NPOS(I).GT.0) GO TO 175
      IF(DIST(I).GT.DIST(NKMIN)) GO TO 170
      NKMIN=I
170   CONTINUE
175   NUS=I-1
      NUSL=NUS-1
      IF(TK(I).EQ.0) NUSL=NUSL+1
      IF(TK(I).EQ.0) GO TO 180
C*       �EPEC�ET TK ��� ����A��E�O K HA�A�� ���A � �EPECTAHOBKA
      NPS=IABS(NPOS(NKMIN))
      TKS=-DIST(NKMIN)
      DISTS=DIST(NKMIN)
      NPOS(NKMIN)=NPOS(NUS)
      NPOS(NUS)=NPS
      TK(NKMIN)=TK(NUS)
      TK(NUS)=TKS
      DIST(NKMIN)=DIST(NUS)
      DIST(NUS)=DISTS
      IF(IFLG.EQ.0) GO TO 180
C*      ��OP��O��BAH�E � KOOP��HAT ���OB
      XKS=XK(NKMIN)
      YKS=YK(NKMIN)
      XK(NKMIN)=XK(NUS)
      XK(NUS)=XKS
      YK(NKMIN)=YK(NUS)
      YK(NUS)=YKS
C
180   IF(NPOS(NKNOT).GT.0) GO TO 200
C*      �PABEE "�PABO�O" KOH�A KOHT�PA ECT� ����
C*        ��EM ����A����
      NKMIN=NKNOT
      NK1=NKNOT-1
      DO 190 I=1,NK1
      IF(TK(NKNOT-I).LE.T(N).AND.NPOS(NKNOT-I).GT.0) GO TO 195
      IF(DIST(NKNOT-I).GT.DIST(NKMIN)) GO TO 190
      NKMIN=NKNOT-I
190   CONTINUE
195   NUS=NKNOT-I+1
      NUSR=I-1
      IF(TK(NKNOT-I).EQ.T(N)) NUSR=NUSR+1
      IF(TK(NKNOT-I).EQ.T(N)) GO TO 200
C*       �EPEC�ET TK ��� ����A��E�O K KOH�� ���A � �EPECTAHOBKA
      NPS=IABS(NPOS(NKMIN))
      TKS=T(N)+DIST(NKMIN)
      DISTS=DIST(NKMIN)
      NPOS(NKMIN)=NPOS(NUS)
      NPOS(NUS)=NPS
      TK(NKMIN)=TK(NUS)
      TK(NUS)=TKS
      DIST(NKMIN)=DIST(NUS)
      DIST(NUS)=DISTS
      IF(IFLG.EQ.0) GO TO 200
C*      ��OP��O��BAH�E � KOOP��HAT ���OB
      XKS=XK(NKMIN)
      YKS=YK(NKMIN)
      XK(NKMIN)=XK(NUS)
      XK(NUS)=XKS
      YK(NKMIN)=YK(NUS)
      YK(NUS)=YKS
200   RETURN
      END
      SUBROUTINE INCLT(TK,XK,YK,TIN,XIN,YIN,NKNOT,DIST,NPOS,DIN,NIN,
     *NDOUBL,IFLG,NKLAST)
C*********************************************************************
C*
C*    B MACC�B� TK,XK,YK,DIST,NPOS PA�MEPHOCT�� NKNOT-1 �O�AB��ETC� �O
C*  O�HOM� ��EMEHT�. CO����AETC� MOHOTOHHOCT� �BE���EH�� TK C POCTOM I
C*
C*
C*********************************************************************
      REAL XK(1),YK(1),TK(1),DIST(1)
      INTEGER*2 NPOS(1)
      I1=NKNOT-1
      NK=NKNOT
      NDS=NDOUBL
      DO  4 J=1,NK
      IF(NKNOT.EQ.1.OR.J.EQ.NKNOT) GO TO 5
C*      �O�C�ET ��C�A O��HAKOB� �POEK��� ���OB
      IF(NIN.GT.0.AND.NPOS(J).GT.0.AND.TK(J).EQ.TIN.AND.NDS.EQ.NDOUBL)
     *NDOUBL=NDOUBL+1
      IF(TIN.GT.TK(J)) GO  TO 4
      IF(TIN.EQ.TK(J).AND.(IABS(NPOS(J)).NE.1.OR.IABS(NIN).NE.NKLAST))
     *GO TO 4
C*       BBO��M ��EMEHT� B MACC�B�
      DO 1 K=J,I1
      TK(I1-K+J+1)=TK(I1-K+J)
      DIST(I1-K+J+1)=DIST(I1-K+J)
      NPOS(I1-K+J+1)=NPOS(I1-K+J)
      IF(IFLG.EQ.1) XK(I1-K+J+1)=XK(I1-K+J)
      IF(IFLG.EQ.1) YK(I1-K+J+1)=YK(I1-K+J)
1     CONTINUE
      GO TO 5
4     CONTINUE
5     TK(J)=TIN
      DIST(J)=DIN
      NPOS(J)=NIN
      IF(IFLG.EQ.1) XK(J)=XIN
      IF(IFLG.EQ.1) YK(J)=YIN
      RETURN
      END
