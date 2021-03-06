C**********************************************************************C
C*                                                                    *C
C*                                                                    *C
C*        C��A�H O��CAH�E      �PA�MEHTA KOHT�PA (�APAMETP��ECKOE     *C
C*         �PE�CTAB�EH�E) C PAC�ETOM TO�HOCT� A��POKC�MA���           *C
C*                                                                    *C
C* �APAMETP�: XY(N,2) - KOOP��HAT� TO�EK KOHT�PA (BXO�), �O��OTOB�EH- *C
C*       H�E �PO�PAMMAM� TPREP � TRECLC.                              *C
C*            T(I) - �APAMETP��ECKA� KOOP��HATA TO�EK KOHT�PA         *C
C*            C(NKNOT,3,2) - �APAMETP� C��A�HA B  �OPMATE             *C
C*       ��� IMSL. �OC�E�H�� �H�EKC - HOMEP KOOP��HAT�:1-X,2-Y(B�XO�) *C
C*            IC=NKNOT-1                                              *C
C*            XK(NKNOT),YK(NKNOT),TK(NKNOT) - X,Y � T KOOP��HAT�      *C
C*       ���OB C��A�HA (BXO�)                                         *C
C*            EPS -  �O����B�A�C� TO�HOCT� A��POKC�MA��� (��������.   *�
C*		����������)					      *C
C*            IFLG=AB  -  ��PAB������ �APAMETP :                      *C
C*              B=0 - �POBO��T� A��POKC�MA��� KOHT�PA                 *C
C*              B=1 - PACC��TAT� TO��KO C��A�H KO������EHT� � EPS     *C
C*              A=0 - A��POKC�MA��� C O�T�MA��HO� PACCTAHOBKO� ���OB  *C
C*              A=1 - �E�                                             *C
C*            WK((NKNOT+6)*N*2+3*NKNOT) - PA�O��� MACC�B              *C
C*            IRC - KO� BO�BPATA :                                    *C
C*                0 - HOPMA��HOE �ABEP�EH�E                           *C
C*                1 - O���KA B �HA�EH��X N ��� NKNOT                  *C
C*                2 - BXO�H�E MACC�B� HEBEPHO �O��OTOB�EH�            *C
C*                3 - HE �OCT��H�TA H��HA� TO�HOCT� O�T�M��A��� �P�   *C
C*      O�T�MA��HO� PACCTAHOBKE ���OB. HA B�XO�E - HA�����A� KOH��-   *C
C*      ��PA��� ���OB                                                 *C
C*                                                                    *C
C* A��OP�TM: B �AB�C�MOCT� OT �AKA�A �POBO��TC� C��A�H A��POKC�MA���  *C
C*     ��� �POCTO PAC�ET C��A�H KO������EHTOB �O MACC�BAM XK,YK � TK  *C
C*     C PAC�ETOM KP�TEP�� A��POKC�MA��� �PA�MEHTA KOHT�PA.           *C
C*     EC�� KOHT�P �AMKH�T, TO �OC�E �POBE�EH�� A��POKC�MA��� OH      *C
C*     "�AM�KAETC�" �PO�PAMMO� ICSPLN, �TO �P�BO��T K HEKOTOPOM�      *C
C*     �X���EH�� A��POKC�MA��� (EC��, KOHE�HO, OHA �AKA�AHA). EPS     *C
C*     �P� �TOM �EPEC��T�BAETC� �OBTOPHO.                             *C
C*           KOHT�P C��TAETC� �AMKH�T�M, EC�� XK(1)=XK(NKNOT) �       *C
C*     YK(1)=YK(NKNOT)                                                *C
C*                                                                    *C
C* �P�ME�AH�E: �P� �POBE�EH�� A��POKC�MA��� �O��HO B��O�H�T�C�        *C
C*             O�PAH��EH�E: NKNOT HE �O�EE 28 - O�PAH��EH�E �PO�PAMM  *C
C*             ICSFKU � ICSFKV � TO��KO.                              *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      SUBROUTINE SPAPP(XY,T,N,XK,YK,TK,NKNOT,C,IC,EPS,WK,IFLG,IRC)
      common/nksave/ nks
      REAL XY(N,2),T(N),XK(1),YK(1),TK(1),WK(N,1),C(IC,3,2)
	DATA NKMAX/28/
	write(*, *) 'SPAPP: N =', N, ', NKNOT =', NKNOT
C*          �CTAHOBKA ��A�KA: �AMKH�T KOHT�P ��� HET?
      ICZ=0
      IF(XK(1).EQ.XK(NKNOT).AND.YK(1).EQ.YK(NKNOT)) ICZ=1
C*          �POBEPKA KOPPEKTHOCT� BXO�H�X �APAMETPOB
      IRC=0
      IF(IFLG-IFLG/2*2.EQ.1.AND.N.GE.4.AND.NKNOT.GT.1) GO TO 185
      IRC=1
      IF(N.LT.4.OR.NKNOT.LT.2.OR.NKNOT.GT.NKMAX) GO TO 300
      IRC=0
      nks=nknot
      IF(IFLG/2-IFLG/4*2.EQ.1) GO TO 150
C
C*       O�T�MA��H�� �O��OP �O�O�EH�� ���OB
C
	CALL  ICSVKU(T,XY,N,TK,NKNOT,WK,C,ic,EPS,WK(NKNOT+NKNOT+1,1),
     *IER)
      IF(IER.EQ.33) IRC=3
      IF(IER.GT.33) IRC=2
      IF(IRC.EQ.2) GO TO 300
      GO TO 180
C
C*       A��POKC�MA��� �P� ��KC�POBAHH�X �O�O�EH��X Y��OB
C
150   CALL  ICSFKU(T,XY,N,0,TK,NKNOT,WK,C,ic,EPS,WK(NKNOT+NKNOT+1
     *,1),IER)
      IF(IER.NE.0) IRC=2
      IF(IRC.EQ.2) GO TO 300
C*       �EPEKA��BAEM KOOP��HAT� ���OB
180   DO 183 I=1,ic
      XK(I)=WK(I,1)
183   YK(I)=WK(NKNOT+I,1)

C*       PACC��T�BAEM KOOP��HAT� �OC�E�HE�O ���A
      CALL  ICSEVU(TK,XK,NKNOT,C,ic,TK(NKNOT),XK(NKNOT),1,IER)
      CALL  ICSEVU(TK,YK,NKNOT,C(1,1,2),ic,TK(NKNOT),YK(NKNOT),1,IER)
      IF(ICZ.EQ.0) GO TO 187
C*       KOHT�P �AMKH�T. �E�AEM �AM�KAH�E E�O A��POKC�MA���
      XK(NKNOT)=(XK(1)+XK(NKNOT))/2.
      YK(NKNOT)=(YK(1)+YK(NKNOT))/2.
      XK(1)=XK(NKNOT)
      YK(1)=YK(NKNOT)
C*        PAC�ET C��A�H KO������EHTOB �O MACC�BAM XK,YK � TK
185   IF(ICZ.EQ.1) CALL  ICSPLN(TK,XK,NKNOT,C,ic,WK,IER)
      IF(ICZ.EQ.0) CALL  ICSCCU(TK,XK,NKNOT,C,ic,IER)
	IF(IER.NE.0) IRC=2
      IF(IRC.EQ.2) GO TO 300
      IF(ICZ.EQ.1) CALL  ICSPLN(TK,YK,NKNOT,C(1,1,2),ic,WK,IER)
	IF(ICZ.EQ.0) CALL  ICSCCU(TK,YK,NKNOT,C(1,1,2),ic,IER)
	IF(IER.NE.0) IRC=2
      IF(IRC.EQ.2) GO TO 300
C*         PAC�ET TO�HOCT� A��POKC�MA���
187   EPS=0.
      CALL  ICSEVU(TK,XK,NKNOT,C,ic,T,WK,N,IER)
      CALL  ICSEVU(TK,YK,NKNOT,C(1,1,2),ic,T,WK(1,2),N,IER)
c      dt=t(2)-t(1)
      DO 190 I=1,N
c      if(i.ne.n.and.i.ne.1) dt=t(i+1)-t(i-1)
c      if(i.eq.n) dt=t(n)-t(n-1)
190    EPS=EPS+((XY(I,1)-WK(I,1))**2+(XY(I,2)-WK(I,2))**2)
      eps=sqrt(eps/(n-1)) 
c190   eps=sqrt(eps/(t(n)-t(1))/2.)
300   RETURN
      END