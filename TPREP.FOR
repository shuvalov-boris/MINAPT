C**********************************************************************C
C*                                                                    *C
C*                                                                    *C
C*        �CTAHOBKA COOTBETCTB�� ���OB � �PA�MEHTA KOHT�PA �O         *C
C*     T-KOOP��HATE ��� �APAMETP��ECKO�O C��A�H-�PE�CTAB�EH�� KOHT�PA *C
C*                                                                    *C
C* �APAMETP�: XY(N,2) - KOOP��HAT� TO�EK KOHT�PA (BXO�). HA B�XO�E -  *C
C*       CO C��A�EHH�M� "B��POCAM�", EC�� C��A��BAH�E �AKA��BA�OC�.   *C
C*            T(I) - PACCTO�H�E OT TO�K� I �O �EPBO� TO�K� KOHT�PA    *C
C*       �O ���E KOHT�PA, A��POKC�M�P�EMO�O �OMAHO�, �POXO���E�       *C
C*       �EPE� TO�K� (XY(I,1),XY(I,2)) (B�XO�)                        *C
C*            XK(NKNOT),YK(NKNOT)           - X,Y     KOOP��HAT�      *C
C*       ���OB C��A�HA (BXO�).                                        *C
C*            TK(I) - T-KOOP��HATA ���A I,      PACC��T�BAETC�        *C
C*       �OCPE�CTBOM �POEKT�POBAH�� ���OB HA �PA�MEHT KOHT�PA.        *C
C*            NBEG,NEND - HOMEPA TO�EK, �KA��BA���E HA HA�A�O � KOHE� *C
C*       MAKC�MA��HO�O OTPE�KA �PA�MEHTA,KOTOP�� �E��T BH�TP� OTPE�KA,*C
C*       BK���A��E�O �POEK��� ���OB (�O T-KOOP��HATE) (B�XO�)         *C
C*            NUSL,NUSR - ��C�O ���OB, KOTOP�E �OC�E �X ��OP��O��BAH��*C
C*       �O T-KOOP��HATE, �O����B�E�C� �P� �POEKT�POBAH�� ���OB HA    *C
C*       KOHT�P, HE �ME�T TO�EK �PA�MEHTA, �E�A��X ME��� TAK�M ���OM  *C
C*       � E�O COCE�OM, � HAXO��TC� "�EBEE" ��� "�PABEE" �PA�MEHTA    *C
C*            EPS - KP�TEP�� ��� PGM C��A��BAH�� B��POCOB.            *C
C*            NDOUBL - ��C�O ���OB �ME���X ��EHT��H�� T-KOOP��HAT�    *C
C*       �POEK��� �X HA �PA�MEHT KOHT�PA. B C���AE, KO��A �EPB�� �    *C
C*       �OC�E�H�� ���� COB�A�A�T C HA�A��HO� TO�KO� �AMKH�TO�O       *C
C*       KOHT�PA (PABHO� KOHE�HO�), TO T-KOOP��HAT� �T�X ���OB PA����-*C
C*       H� � NDOUBL HE ��KC�P�ET �TO �����POBAH�E ���OB              *C
C*            IFLG=AB  -  ��PAB������ �APAMETP :                      *C
C*              B=1 - C��A��BAT� KOOP��HAT� TO�EK �PA�MEHTA HEO�XO��MO*C
C*              B=0 - TAKO� HEO�XO��MOCT� HET                         *C
C*              A=1 - XK � YK HEO�XO��MO ��OP��O��T� B COOTBETCTB�� C *C
C*       �O����B�E�C� T-KOOP��HATO�    �POEK��� �TO�O ���A            *C
C*              A=0 - ��OP��O��BAT� HE HA�O                           *C
C*              WK(MAX(N*3,2*NKNOT))-PA�O��� MACC�B.HA B�XO�E -       *C
C*       KBA�PAT� M�H�MA��H�X PACCTO�H�� OT                           *C
C*               ���OB �O �PA�MEHTA KOHT�PA                           *C
C*            IWK(I)     - �E�O��C�EHH�� MACC�B. HA B�XO�E - HOMEP�   *C
C*       �����, ������������� �� ����������� TK. PA�MEPHOCT� - NKNOT  *C
C*       EC�� IWK(I)<0, TO �TO ��E� (IABS(IWK(I)))                    *C
C*        BHE �PA�MEHTA KOHT�PA (CM O��CAH�E                          *C
C*       �APAMETPOB NUSL,NUSR)                                        *C
C*            IRC - KO� BO�BPATA :                                    *C
C*                0 - HOPMA��HOE �ABEP�EH�E                           *C
C*                1 - �AHH�E HE AHA����POBA��C� HA HA����E B��POCOB   *C
C*      �P� �A�AH�� TAKO�O AHA���A, T. K. ��C�O B��PACOB > 25%        *C
C*                2 - �OP��OK C�E�OBAH�� �POEK��� ���OB(TK) ��K���ECK�*C
C*      HE COB�A�AET C �OP��KOM C�E�OBAH�� ���OB                      *C
C*                3 - �KB�BA�EHTHO IRC=1&IRC=2                        *C
C*                4 - MA�O N (<7 �P� B=1 ��� <2 �P� B=0)              *C
C*                5 - NKNOT < 2                                       *C
C*                6 - HEBEPHO �O��OTOB�EH� BXO�H�E MACC�B� ���        *C
C*      �PO�PAMM� C��A��BAH�� B��POCOB ICSMOU                         *C
C*                7 - ���� ���������� �����������, �� � TK ���        *C
C*      ����������� �����                                             *C
C*                                                                    *C
C* A��OP�TM: 1.C��A��BA�TC�      B��POC� (EC�� �TO �AKA�AHO)          *C
C*           2.EC�� ��C�O B��POCOB > 25% - C��A��BAH�E HE �PO��BO��TC�*C
C*           3.KOOP��HATA TO�K� I C��A��BAETC�, EC�� OHA OTK�OH�ETC�  *C
C*     OT C��A�H-KP�BO� >EPS*(T(I+3)-T(I-3))/6.                       *C
C*           4.�POEKT�P��M ���� HA �PA�MEHT KOHT�PA, PACC��T�BA� �    *C
C*     ��OP��O��BA� TK,IWK,WK � XK,YK (EC�� TPE��ETC�). C��TAEM NBEG, *C
C*     NEND,NUSL,NUSR � NDOUBL.                                       *C
C*                                                                    *C
C* �P�ME�AH�E: EC�� �P� �POEKT�POBAH�� ��OP��O�EHH�X PAHEE ���OB HET  *C
C*             �BEPEHHOCT� B ��K���ECKOM COB�A�EH�� �OP��KA C�E�OBAH��*C
C*             �X �POEK���, TO PEKOMEH��ETC� BBO��T� B �PA�MEHT       *C
C*             KOHT�PA KOOP��HAT� ���OB.                              *C
C*             TO�K� �PA�MEHTA HE �O��H� COB�A�AT�, �A �CK���EH�EM    *C
C*             KPA�H�X TO�EK �P� �AMKH�TOM KOHT�PE                    *C
C*             �PE��O�A�AETC� �P� C��A��BAH��, �TO �EPB�E � �OC�E�H�E *C
C*             3 TO�K� �PA�MEHTA HAXO��TC� HA CBOEM "�CT�HHOM" MECTE  *C
C*             (O�PAH�E�EH�E �PO�PAMM� ICSMOU)                        *C
C*             CMOTP�TE TAK�E KOMMEHTAP�� �PO�AMM� CONKN              *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      SUBROUTINE TPREP(XY,T,N,XK,YK,TK,NKNOT,NBEG,NEND,NUSL,NUSR,NDOUBL,
     *EPS,WK,IWK,IFLG,IRC)
      REAL XY(N,2),T(N),XK(1),YK(1),TK(1),WK(N,3)
      INTEGER*2 IWK(NKNOT)
C		����� ���������� ��� ��������
      write(*,*) 'TPREP'
      open(18, FILE='tprep.log', action='write')
	write(18, *) '-------TPREP.LOG---------- '
	write(18, *) 'N =', N, ', NKNOT =', NKNOT
	write(18, *) 'XY T:'
	write(18, 1) (XY(I, 1), XY(I, 2), T(I), I = 1, N)
    1 format(3F12.4)
      close(18)
C
C*          �POBEPKA KOPPEKTHOCT� BXO�H�X �APAMETPOB
      IRC=4
      NMAX=2
      IF(IFLG-IFLG/2*2.EQ.1) NMAX=7
      IF(N.LT.NMAX) GO TO 300
      IRC=5
      IF(NKNOT.LT.2) GO TO 300
C*          ������� �� ������ (� ������ �� �����) ?
      ICN=0
      IF(XK(1).EQ.XK(NKNOT).AND.YK(1).EQ.YK(NKNOT)) ICN=1
C*          PAC�ET PACCTO�H�� �O TO�EK �OMAHO� B�O�� KOHTYPA
      T(1)=0.
      DO 10 I=2,N
      T(I)=T(I-1)+SQRT((XY(I,1)-XY(I-1,1))**2+(XY(I,2)-XY(I-1,2))**2)
10    CONTINUE
      IRC=0
      IF(IFLG-IFLG/2*2.EQ.0) GO TO 100
C
C*        AHA��� �AHH�X HA B��POC�
C
      DO 20 I=1,N
      WK(I,1)=XY(I,1)
      WK(I,2)=XY(I,2)
20    WK(I,3)=T(I)
C*          C��A��BAEM  TO�K� �O KOOP��HATE X
      CALL  ICSMOU(WK(1,3),WK,N,1.,EPS,N/4+1,WK(1,4),IER)
      IF(IER.NE.0) IRC=6
      IF(IRC.EQ.6) GO TO 300
C*          C��A��BAEM  TO�K� �O KOOP��HATE Y
      CALL  ICSMOU(WK(1,3),WK(1,2),N,1.,EPS,N/4+1,WK(1,4),IER)
      IF(IER.NE.0) IRC=6
      IF(IRC.EQ.6) GO TO 300
C*            NVILD - ��C�O B��POCOB
      NVILD=0
      DO 30 I=1,N
      IF(WK(I,1).EQ.XY(I,1).AND.WK(I,2).EQ.XY(I,2)) GO TO 30
      NVILD=NVILD+1
30    CONTINUE
C
      IF(NVILD*4.GT.N) IRC=1
      IF(NVILD*4.GT.N.OR.NVILD.EQ.0) GO TO 100
C*          �EPEC�ET T
      DO 60 I=1,N
      DO 60 J=1,2
60    XY(I,J)=WK(I,J)
      T(1)=0.
      DO 70 I=2,N
70    T(I)=T(I-1)+SQRT((XY(I,1)-XY(I-1,1))**2+(XY(I,2)-XY(I-1,2))**2)
C
C*       PAC�ET TK,NBEG,NEND, � T�.,  Y�OP��O��BAH�E Y��OB B �OP��KE
C*   BO�PACTAH�� TK
C
100   CALL  CONKN(XY,XY(1,2),T,N,XK,YK,TK,NKNOT,NBEG,NEND,WK,IWK,NUSL,
     *NUSR,NDOUBL,IFLG/2)
      IF(ICN.EQ.1) NDOUBL=NDOUBL-1
      IF(NDOUBL.LT.0) NDOUBL=0
C*          �POBEPKA HA ��K���ECK�� ��OP��O�EHHOCT� �POEK��� ���OB
      ICU=-1
      I1=NKNOT
      NKBEG=0
      DO 110 I=1,NKNOT
      IF(IABS(IWK(I)).EQ.1) NKBEG=I-1
      IF(IABS(IWK(I)).NE.IABS(IWK(I1)+1).AND.IABS(IWK(I)).NE.1)
     * GO TO 120
      IF(IABS(IWK(I)).EQ.1.AND.IABS(IWK(I1)).NE.NKNOT) GO TO 120
      I1=I
110   CONTINUE
      ICU=0
      GO TO 130
C*          �������� �� ���������� �������� ������� ���������� �����
120   I1=NKNOT
      ICU=-1
      NKBEG=0
      DO 125 I=1,NKNOT
      IF(IABS(IWK(I)).EQ.1) NKBEG=I-1
      IF(IABS(IWK(I)).NE.IABS(IWK(I1)-1).AND.IABS(IWK(I)).NE.NKNOT)
     * GO TO 127
      IF(IABS(IWK(I)).EQ.NKNOT.AND.IABS(IWK(I1)).NE.1) GO TO 127
      I1=I
125   CONTINUE
      ICU=1
      GO TO 130
127   IF(ICN.EQ.0) GO TO 130
C*           ������ ������������� �������. ������ ������� ����������
C*    ������� � ���������� ����� � ����� ��������� �� ����������� 
C*                     ���������������
      DO I=1,NKNOT
        IF(IABS(IWK(I)).EQ.1) NK1=I
        IF(IABS(IWK(I)).EQ.NKNOT) NKLAST=I
      END DO
      IWKS=IWK(NK1)
      IWK(NK1)=IWK(NKLAST)
      IWK(NKLAST)=IWKS
C*       ���������: ���� ���������� �����������?
      ICU=-1
      I1=NKNOT
      NKBEG=0
      DO I=1,NKNOT
        IF(IABS(IWK(I)).EQ.1) NKBEG=I-1
        IF(IABS(IWK(I)).NE.IABS(IWK(I1)+1).AND.IABS(IWK(I)).NE.1)
     *   GO TO 128
        IF(IABS(IWK(I)).EQ.1.AND.IABS(IWK(I1)).NE.NKNOT) GO TO 128
      I1=I
      END DO
      ICU=0
      GO TO 130
C*          �������� �� ���������� �������� ������� ���������� �����
128   I1=NKNOT
      ICU=-1
      NKBEG=0
      DO I=1,NKNOT
         IF(IABS(IWK(I)).EQ.1) NKBEG=I-1
         IF(IABS(IWK(I)).NE.IABS(IWK(I1)-1).AND.IABS(IWK(I)).NE.NKNOT)
     *   GO TO 130
         IF(IABS(IWK(I)).EQ.NKNOT.AND.IABS(IWK(I1)).NE.1) GO TO 130
         I1=I
      END DO
      ICU=1
C
130   IF(IFLG/2.EQ.1.AND.ICU.GE.0) GO TO 300
      IF(IFLG/2.EQ.1.AND.ICU.LT.0) GO TO 290
      IF(IABS(IWK(1)).EQ.1.AND.ICU.EQ.0) GO TO 300
C*          ��� ����� - ���� �� �������������. �������� � ������������
C*    � �������� ���������� ����� ������� TK,IWK,WK
C      
C*          ���������� ����������
      DO I=1,NKNOT
        NU=IABS(IWK(I))
        WK(NU+NKNOT,1)=WK(I,1)
      END DO
C
      DO I=1,NKNOT
        WK(I,1)=WK(I+NKNOT,1)
      END DO
C*          ������������ ������ TK
      TKMAX=0.
      DO I=1,NKNOT
        NU=IABS(IWK(I))
          IF(IWK(I).GT.0.AND.TKMAX.LT.TK(I)) TKMAX=TK(I)
        WK(NU+NKNOT,1)=TK(I)
      END DO
C
      DO I=1,NKNOT
        TK(I)=WK(I+NKNOT,1)
      END DO
C*           
      IF(ICU.LT.0) GO TO 290
      IF(ICU.EQ.0) GO TO 140
C*          ���� "���������������". ����������� ������ ��
      DO I=1,NKNOT
        IF(IWK(I).GT.0) TK(I)=TKMAX-TK(I)
      END DO
      IRC=0
      IF(IABS(IWK(1)).EQ.1.AND.ICU.EQ.0) GO TO 300
C
140      IRC=7
      IF(XY(1,1).NE.XY(N,1).OR.XY(1,2).NE.XY(N,2)) GO TO 300
      IRC=0
C*          �������� ������� �������. ������������� TK
      IF(ICU.EQ.1) GO TO  150
      DO I=1,NKBEG
        TK(NKNOT-I+1)=TK(NKNOT-I+1)+T(N)
      END DO
      GO TO 160
150   NKA=NKNOT-NKBEG-1
      IF(NKA.EQ.0) GO TO 160
C*          �������� �� ��� ����������������� �����
      DO I=NKBEG+2,NKNOT
        TK(I)=T(N)+TK(I)
      END DO
C*          TK(1) -> 0
160   DO I=2,NKNOT
        TK(I)=TK(I)-TK(1)
      END  DO
      TK(1)=0.
      GO TO 300
C*      HET ��K���ECKO� ��OP��O�EHHOCT�
290   IRC=IRC+2
300   RETURN
      END
