C**********************************************************************
C*
C*
C*     C��T�BAH�E ��������������� �H�OPMA��� B �OPMATE C�CTEM�
C*  SEAFLOW (�EO�PA���ECK�E KOOP��HAT�) � �PEO�PA�OBAH�E EE B �EKAPTOB�
C*  KOOP��HAT�. ��������
C*   �������������� ���������� ��� ������ � ������� �� ������� ��������
C*   ���������� (AVHRR-���������� 1 ��)
C*
C* �APAMETP�: NAMFIL - �M� �A��A C �AHH�M�
C*            FI,LAMBDA - TO�KA HA�A�A KOOP��HAT (�EO�PA���ECK�X)
C*	      idate0(7) - ���� ������ ������� �������: ���, �����, 
C*		�����, ����, ������.
C*            XY(N,2) - �������������� KOOP��HAT� ��������� �������
C*            V(N,2) - ����� ���������� ������� (� �����) � ������ ����-
C*             �������� ���������� (� ������)
C*            METPOS(n+1) - ����� ������� � ������ �������
C*            WK - PA�O��� MACC�B PA�MEPHOCT�� N*2
C*            IFLG:0 - �A�� NAMFIL HE OTKP�T
C*                 1 - OTKP�T
C*                 2,3 - �PO��CT�T� O��H ��OK �AHH�X (2-�A�� HE OTKP�T,
C*                       3-OTKP�T)
C*	      ITYPCN - ��� �������������� ������:
C*		0 - xy � 10�� �� ������, 1 - xy � �������� � 
C*		�������������� ������� ���������
C*	      MAXREC - ������������ ����� ����������� ����� � �����
C*            IRC - KO� BO�BPATA
C*                 0 - C��TAH ��OK �H�OPMA���
C*                 2 - TO�E, �TO � 0 + �OCT��H�T KOHE� �A��A �
C*                       �A�� �AKP�T
C*                 4 - �A�� HE HA��EH �P� OTKP�T��
C*                 5 - O���KA BBO�A-B�BO�A,��� ������ �� ���������������
C*                     ��� �A�� ��CT
C*		   6 - ����� ����� ����� ������ ����� ������ maxrec
C*
C* A��OP�TM :  �PO��CKAETC� ��� C��T�BAETC� ��OK �AHH�X. C��TAHH�E
C*             �AHH�E �EPEBO��T�C� B �EKAPTOB�� C�CTEM� KOOP��HAT,
C*	       ��� �������� � �������������� � ������� � 
C*	       �����(fi,lambda).
C*             TO�HOCT� �PEO�PA�OBAH�� O��CAHA B �PO�PAMME DGKM
C*
C* �P�ME�AH�E.B��PAHHA� �EKAPTOBA� C�CTEMA KOOP��HAT B OT����� OT
C*   MEPKATOPOBCKO� �POEK��� �MEET PABHOMEPH�� �KA�� KOOP��HAT. �O�TOM�
C*   �P� O�PE�E�EH�� CME�EH�� ME��� �B�M� ��KCE�AM�, OTCTO���M� �O
C*   ��POTE �P�� OT �P��A HA L KM � HA PSI �PA��COB �O �O��OTE �O�B�T�C�
C*   O���KA B X-KOOP��HATE ERR=PI*L*PSI*SIN(LAMBDA)/180. ��� ERR=1KM
C*   �P� L=57KM PSI=1/SIN(LAMBDA).
C* PEKOMEH�A���.�EHTPOM KOOP��HAT ����E BCE�O B��PAT� �EHTP �CC�E��EMO�O
C*   PA�OHA OKEAHA.C�E��ET ���E�AT� PACC�ETOB �O����X PACCTO�H�� B�O��
C*   MEP���AHA B�A�� OT �EHTPA KOOP��HAT.
C* ����������. ������� ��� (idate0 � ���� ������� idate) �� ������ ����
C*   ������ ����.
C*
C*
C***********************************************************************
      SUBROUTINE RDGINF(NAMFIL,FI,LAMBDA,idate0,XY,V,N,METPOS,
     *WK,IFLG,itypcn,maxrec,IRC)
      REAL V(1),XY(2,1),WL1(2),WL2(2),RL1(2),RL2(2),WK(1),
     *lambda,tm(2)
      character*(*) NAMFIL
      INTEGER*2 INTG1(2),INTG2(2),METPOS(1),idate0(5),idate(5,2),mn(12)
      data mn/31,28,31,30,31,30,31,31,30,31,30,31/
      IRC=0
      IF(IFLG-IFLG/2*2.EQ.1) GO TO 10
C*       OTKP�BAEM �A��
      OPEN (8,ERR=300,FILE=NAMFIL,STATUS='OLD')
10    READ (8,20,END=290,ERR=290) I1,INBLOK,INTG1,RL1,INTG2,RL2
20    FORMAT(I1,I2/,2I5,2F9.3,/,2i5,2f9.3)
      IF(INBLOK.NE.2.or.i1.ne.1) GO TO 290
      read(8,21,err=290,end=290) nrec,metpos(1)
21    format(2i3)
      if(nrec.gt.maxrec) irc=6
      if(nrec.gt.maxrec) go to 310
      IF(IFLG.LE.1) GO TO 35
C*      �PO��CK ��OKA �AHH�X
      DO 30 I=1,NREC
      READ (8,40,ERR=290,END=290)
30    CONTINUE
      GO TO 254
C
C*             BBO� ������
C
35    n=nrec
c*		������ �������������� ������� time 
      je=1
      id0=0
c*		������ ������ ��� � ����. ���� idate0
      me=idate0(2)-1
      if(me.eq.0) go to 23
        do 22 i=1,me
22      id0=id0+mn(i)
23    id0=id0+idate0(3)
      idy0=0
      if(idate0(1)-idate0(1)/4*4.eq.0) idy0=1
      if(me.ge.2) id0=id0+idy0
c*		�������� ������� ��� � ����. ���� idate0
      h0=idate0(4)+idate0(5)/60.
c
      DO 60 I=1,NREC
      READ(8,40,ERR=290,END=290) INTG1,WL1,INTG2,WL2,ns1,(idate(j,1),
     *j=1,5),v(i+n),ns2,(idate(j,2),j=1,5),v(i+1+n)
40    FORMAT(2I5,2F9.3,2I5,2F9.3,2(i5,5i3,f9.3))
      XY(1,I)=WL1(2)
      XY(2,I)=WL1(1)
      metpos(i+1)=ns1
c*		�������� ������� ��� ��� idate
      do 26 j=1,je
        id=0
c*		������ ������ ��� � ����. ���� idate
        me=idate(2,j)-1
        if(me.eq.0) go to 25
          do 24 k=1,me
24        id=id+mn(k)
25      id=id+idate(3,j)
        idy=0
        if(idate(1,j)-idate(1,j)/4*4.eq.0) idy=1
        if(me.ge.2) id=id+idy
c*		�������� ������� ��� � ����. ���� idate
        h=idate(4,j)+idate(5,j)/60.
        idd=id-id0
        ieq=0
	if(idate0(1).eq.idate(1,j)-1.or.idate(1,j)-idate0(1).lt.-1) 
     *  ieq=1
        if(ieq.eq.1) idd=idd+365+idy0
        ieq=0
	if(idate0(1).eq.idate(1,j)+1.or.idate0(1)-idate(1,j).lt.-1) 
     *  ieq=1
        if(ieq.eq.1) idd=idd-365-idy
	tm(j)=idd*24.+h-h0
26	continue
	v(i)=tm(1)
60    CONTINUE
C*      �EPEC�ET �AHH�X 
      CALL  CONVGD(V,XY,WK,N,FI,LAMBDA,dt,2,itypcn)
C*      �PEO�PA�OBAH�E MACC�BA XY B �OPMAT XY(N,2)
      DO 251 I=1,N
      WK(I)=XY(1,I)
251   WK(I+N)=XY(2,I)
      DO 252 I=1,N
      XY(I,1)=WK(I)
252   XY(N+I,1)=WK(N+I)
C
C*      �POBEPKA: ECT� �� E�E ��OK� �AHH�X B �A��E ?
254   READ (8,255,ERR=270,END=270) I
255   FORMAT(I1)
      BACK SPACE 8
      if(i.ne.1) go to 270
      GO TO 320
C*      KOHE� �A��A
270   IRC=IRC+2
      GO TO 310
C*      O���KA BBO�A-B�BO�A
290   IRC=5
      GO TO 310
C*      �A�� HE HA��EH
300   IRC=4
      go to 320
310   close(8,err=320)
320   RETURN
      END