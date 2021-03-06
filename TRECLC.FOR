
C**********************************************************************C
C*                                                                    *C
C*                                                                    *C
C*        �EPEC�ET T-KOOP��HAT �PA�MEHTA KOHT�PA B COOTBETCTB�� CO    *C
C*  �HA�EH��M� T-KOOP��HAT �POEK��� ���OB HA KOHT�P � �X �CT�HH�X     *C
C*  BE����H T                                                         *C
C*                                                                    *C
C* �APAMETP�: X(N),Y(N)-KOOP��HAT� TO�EK KOHT�PA (BXO�). HA B�XO�E -  *C
C*       �EPECOPT�POBAHH�E B COOTBETCTB�� C HOB�M� T-�HA�EH��M�       *C
C*            T(I) - PACCTO�H�E OT TO�K� I �O �EPBO� TO�K� KOHT�PA    *C
C*       �O ���E KOHT�PA, A��POKC�M�P�EMO�O �OMAHO�, �POXO���E�       *C
C*       �EPE� TO�K� (XY(I,1),XY(I,2)) (BXO�).B�XO� - �EPEC��TAHH�E   *C
C*       BE����H�                                                     *C
C*            N - MO�ET �MEH���T�C� HA 1, EC�� KOHT�P �AMKH�T � T-    *C
C*       �POEK��� �EPBO�O � �OC�E�HE�O ���A (���� COB�A�A�T) HE       *C
C*       �P�XO��TC� HA TO�K� HA�A�A KOHT�PA (COB�A�A���� C KOHE�HO�   *C
C*       TO�KO� �AMKH�TO�O KOHT�PA
C*            TK(NKNOT) - ��OP��O�EH�E �O BO�PACTAH�� T KOOP��HAT�*C
C*       ���OB C��A�HA - �POEK��� ���OB HA �PA�MEHT KOHT�PA (BXO�)    *C
C*            NKNOT - ��C�O ���OB.                                    *C
C*            TKOLD(I) - �ECTB�TE��HA� T-KOOP��HATA ���A I (��OP��O�.)*C
C*            KFIRST - HOMEP ���A, �ME��E�O �POEK��� HA KOHT�P TK(1)  *C
C*            WK(N) - PA�O��� MACC�B                                  *C
C*                                                                    *C
C* A��OP�TM: BE����HE T(I), �O�AB�E� B �HTEPBA� (TK(J+1)-TK(J)),      *C
C*     �P�CBA�BAETC� �HA�EH�E TKOLD(KFIRST+J-1)+(T(I)-TK(I))/         *C
C*     /(TK(J+1)-TK(J))*(TKOLD(KFIRST+J)-TKOLD(KFIST+J-1)).           *C
C*     BE����H� X,Y COPT�P��TC� B �OP��KE �BE���EH�� �O����B�E�C�     *C
C*     T-KOOP��HAT�. EC�� KOHT�P �AMKH�T, TO N MO�ET ��T� �P�         *C
C*     O�PE�E�EHH�X �C�OB��X �MEH��EHA HA 1. �PE��O�A�AETC�, �TO:     *C
C*   1.�OP��K� C�E�OBAH�� ���OB ��� TK � TKOLD OT���A�TC� ����        *C
C*     ��K���ECK�M C�B��OM.                                           *C
C*   2.B MACC�BAX T,TK � TKOLD HET �B�X COB�A�A���X �HA�EH��,EC��     *C
C*     KOHT�P HE �AMKH�T. EC�� �AMKH�T, TO MO��T COB�A�AT� TO��KO     *C
C*     �EPBA� � �OC�E�H�� KOOP��HAT� � MACC�BOB T � TK.               *C
C*                                                                    *C
C* �P�ME�AH��:  T-KOOP��HAT� ���A I, �O���EH�E  �PO�PAMMO� TPREP �    *C
C*     �ME��E�O  �HA�EH�E IWK(I)<0 HE �O��H� �P�C�TCTBOBAT� B TK �    *C
C*     TKOLD. TO�K� KOHT�PA - TO�K� B ��A�O�OHE �H�EKCOB NBEG � NEND  *C
C*     H�KAK�X �POBEPOK HA KOPPEKTHOCT� �AHH�X HE �E�A�. TO�EK        *C
C*     �PA�MEHTA ME��� ���AM� NKNOT � 1 HE �O��HO ��T�              *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      SUBROUTINE TRECLC(X,Y,T,N,TK,TKOLD,NKNOT,KFIRST,WK)
      REAL T(N),X(N),Y(N),TK(1),TKOLD(1),WK(N)
C*      ��EM ��E� KB, ����A���� K T(1): TK(KB)<=T(1)
      NKNOT1=NKNOT-1
      DO 10 KB=1,NKNOT1
      IF(TK(KB+1).GT.T(1).AND.TK(KB).LE.T(1)) GO TO 15
10    CONTINUE
      KB=NKNOT
15    CONTINUE
      KB1=KB+1
      IF(KB1.GT.NKNOT) KB1=1
C*        B��PAC�BAH�E �OC�E�HE� TO�K� KOHT�PA, EC�� �TO HEO�XO��MO
      N1=N
      IF(X(1).NE.X(N).OR.Y(1).NE.Y(N)) GO TO 20
      IF(TK(KB).EQ.T(1)) GO TO 20
      N1=N-1
20    CONTINUE
C*       �EPEC�ET T-KOOP��HAT� �PA�MEHTA KOHT�PA
      NBT=0
      T1=T(1)
      DTK=TK(KB1)-TK(KB)
      IF(TK(KB1).LT.TK(KB)) DTK=T(N)-TK(KB)+TK(KB1)-T1
      TKS=TK(KB)
      IF(TK(KB1).LT.TK(KB)) TKS=TK(KB)-T(N)+T1
      KBO=KB+KFIRST-1
      IF(KBO.GT.NKNOT) KBO=KBO-NKNOT
      KBO1=KBO+1
      DTKO=TKOLD(KBO1)-TKOLD(KBO)
      DO 50 I=1,N1
      IF(KBO.EQ.1.AND.NBT.EQ.0) NBT=I
      T(I)=TKOLD(KBO)+(T(I)-TKS)/DTK*DTKO
      IF(I.EQ.N1) GO TO 50
      IF(T(I+1).LE.TKS+DTK) GO TO 50
C*       �EPECKO���� HA C�E������ OTPE�OK KOHT�PA
      KB=KB1
      KBO=KBO1
      KB1=KB+1
      IF(KB1.GT.NKNOT) KB1=1
      DTK=TK(KB1)-TK(KB)
      IF(DTK.NE.0.) GO TO 55
C*      �EPECKO���� C O�HO�O KPA� �AMKH�TO�O KOHT�PA HA �P��O�
      KB=KB1
      KB1=KB+1
      IF(KB1.GT.NKNOT) KB1=1
      DTK=TK(KB1)-TK(KB)
55    IF(TK(KB1).LT.TK(KB)) DTK=T(N)-TK(KB)+TK(KB1)-T1
      TKS=TK(KB)
C     IF(TK(KB1).LT.TK(KB)) TKS=TK(KB)-T(N)+T1
      IF(KBO.EQ.NKNOT) KBO=1
      KBO1=KBO+1
      DTKO=TKOLD(KBO1)-TKOLD(KBO)
      IF(KBO.EQ.1) NBT=0
50    CONTINUE
C
C*      ��OP��O��BAH�E TO�EK KOHT�PA CO��ACHO HOB�M �HA�EH��M T
C
C*        ��OP��O��BAEM X-KOOP��HAT�
      N=N1
      IF(NBT.EQ.1) GO TO 100
      NP=NBT
      DO 60 I=1,N
      WK(I)=X(NP)
      NP=NP+1
      IF(NP.GT.N) NP=1
60    CONTINUE
      DO 65 I=1,N
65    X(I)=WK(I)
C*      ��OP��O��BAEM Y-KOOP��HAT�
      NP=NBT
      DO 70 I=1,N
      WK(I)=Y(NP)
      NP=NP+1
      IF(NP.GT.N) NP=1
70    CONTINUE
      DO 75 I=1,N
75    Y(I)=WK(I)
C*      ��OP��O��BAEM T-KOOP��HAT�
      NP=NBT
      DO 80 I=1,N
      WK(I)=T(NP)
      NP=NP+1
      IF(NP.GT.N) NP=1
80    CONTINUE
      DO 85 I=1,N
85    T(I)=WK(I)
C
100   RETURN
      END