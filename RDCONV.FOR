C**********************************************************************
C*
C*
C*     C��T�BAH�E KOHT�PHO-BEKTOPHO� �H�OPMA��� B �OPMATE C�CTEM�
C*  SEAFLOW (�EO�PA���ECK�E KOOP��HAT�) � �PEO�PA�OBAH�E EE B �EKAPTOB�
C*  KOOP��HAT�: � �������� BEKTOPOB CKOPOCT� - B CM/CEK. ��������
C*   �������������� ���������� ��� ������ �� �������� �������� �������-
C*   ��� (AVHRR-���������� 1 ��)
C*
C* �APAMETP�: NAMFIL - �M� �A��A C �AHH�M�
C*          paspim(4,2) - ������� ����������� - ���������� ��������
C*            ������ ���� ����������� (paspim(1-4,1)) � ���������
C*          �������: x,y-� ��������, s,d-������,������� � ��������.
C*          paspim(1-4,2) - ��������� ������� ������� ����.(�����) 
C*            FI,LAMBDA - TO�KA HA�A�A KOOP��HAT (�EO�PA���ECK�X)
C*          idate0(7) - ���� ������ ������� �������: ���, �����, 
C*          �����, ����, ������, �������, ����.
C*            XY(N,2) - KOOP��HAT� TO�EK KOHT�PA ��� KOOP��HAT� TO�EK,
C*             ��E ��MEP��AC� CKOPOCT�. PA�MEPHOCT� N*2
C*            V(N,2) - CKOPOCT� TE�EH�� B TO�KAX XY, EC�� C��T�BA�AC�
C*                     BEKTOPHA� �H�OPMA���
C*            METPOS(METPOS(1)+1) - �O�O�EH�� B MACC�BE XY TO�EK -
C*                      HA�A�       �PA�MEHTOB KOHT�POB. METPOS(1) -
C*                 ��C�O �PA�MEHTOB .
C*            TIME(2) - 1:BPEM� ������ ������, �� �������� ��������� ��������
C*                 ���������� (������������ idate0) � ����� (���� �
C*             � ����� ���� ����, ����� - ��� ����� .123)
C*                      2:�������� ����� �������� ��� ���������� 
C*                 ����� ��������� (� �������)
C*            WK - PA�O��� MACC�B PA�MEPHOCT�� N*2
C*            IFLG:0 - �A�� NAMFIL HE OTKP�T
C*                 1 - OTKP�T
C*                 2,3 - �PO��CT�T� O��H ��OK �AHH�X (2-�A�� HE OTKP�T,
C*                       3-OTKP�T)
C*          ITYPCN - ��� �������������� ������:
C*          0 - xy � 10�� �� ������, 1 - xy � �������� � 
C*          �������������� ������� ���������
C*          MAXREC - ������������ ����� ����������� ����� � �����
C*            IRC - KO� BO�BPATA
C*                 0 - C��TAH ��OK KOHT�PHO� �H�OPMA���
C*                 1 - ��OK BEKTOPOB
C*                 2,3 - TO�E, �TO � 0,1 + �OCT��H�T KOHE� �A��A �
C*                       �A�� �AKP�T
C*                 4 - �A�� HE HA��EH �P� OTKP�T��
C*                 5 - O���KA BBO�A-B�BO�A ,��� ��������� ��� ������,
C*                     ��� �A�� ��CT
C*             6 - ����� ����� ����� ������ ����� ������ maxrec
C*             7 - � ����� ����������� ����
C*
C* A��OP�TM :  �PO��CKAETC� ��� C��T�BAETC� ��OK �AHH�X. C��TAHH�E
C*             �AHH�E �EPEBO��T�C� B �EKAPTOB�� C�CTEM� KOOP��HAT,
C*           ��� �������� � �������������� � ������� � 
C*           �����(fi,lambda).
C*             TO�HOCT� �PEO�PA�OBAH�� O��CAHA B �PO�PAMME DGKM
C*             �P� C��T�BAH�� KOHT�PHO� �H�OPMA��� TO�K�-COCE��,
C*             �����P����E �P�� �P��A, B��PAC�BA�TC�.
C*
C* �P�ME�AH�E.B��PAHHA� �EKAPTOBA� C�CTEMA KOOP��HAT B OT����� OT
C*   MEPKATOPOBCKO� �POEK��� �MEET PABHOMEPH�� �KA�� KOOP��HAT. �O�TOM�
C*   �P� O�PE�E�EH�� CME�EH�� ME��� �B�M� ��KCE�AM�, OTCTO���M� �O
C*   ��POTE �P�� OT �P��A HA L KM � HA PSI �PA��COB �O �O��OTE �O�B�T�C�
C*   O���KA B X-KOOP��HATE ERR=PI*L*PSI*SIN(LAMBDA)/180. ��� ERR=1KM
C*   �P� L=57KM PSI=1/SIN(LAMBDA). ��� PAC�ETA CKOPOCTE� TO�HOCT� B�O�HE
C*   �OCTATO�HA. �P� �POBE�EH�� PA�OT C KOHT�PAM� HEO�XO��MO ��T�
C*   BH�MATE��H�M. (�������������� ��� p������ L �p���������:
C*   L=2*r*sin((lam1+lam2)/2), ��� r-p����� �����; lam1,lam2 - ��p���
C*   �����-������).
C*            �������� ����� �������� ������� ����� �������� ��������
C*   ������ ��� dt=12 ����� � �������� <=3m/���. ERR ~ V*dt
C* PEKOMEH�A���.�EHTPOM KOOP��HAT ����E BCE�O B��PAT� �EHTP �CC�E��EMO�O
C*   PA�OHA OKEAHA.C�E��ET ���E�AT� PACC�ETOB �O����X PACCTO�H�� B�O��
C*   MEP���AHA B�A�� OT �EHTPA KOOP��HAT.
C* ����������. ������� ��� (idate0 � ���� ������ idate) �� ������ ����
C*   ������ ����.
C*
C*
C***********************************************************************
	SUBROUTINE RDCONV(NAMFIL,paspim,FI,LAMBDA,idate0,XY,V,N,METPOS,
     *TIME,WK,IFLG,itypcn,maxrec,IRC)
	REAL V(1),XY(2,1),WL1(2),WL2(2),RL1(2),RL2(2),WK(1),paspim(4,2),
     *lambda,tm(2),time(2)
	character*(*) NAMFIL,line80*80
	INTEGER*2 INTG1(2),INTG2(2),METPOS(1),idate0(7),idate(7,2),mn(12)
	data mn/31,28,31,30,31,30,31,31,30,31,30,31/
	IRC=0
	IF(IFLG-IFLG/2*2.EQ.1) GO TO 10
C*       OTKP�BAEM �A��
	OPEN (8,ERR=300,FILE=NAMFIL,STATUS='OLD')
10    READ (8,20,END=290,ERR=290) I1,INBLOK,INTG1,RL1,INTG2,RL2
20    FORMAT(I1,I2/,2I5,2F9.3/,2I5,2F9.3)
	IF(INBLOK.LT.0.OR.INBLOK.GT.1.or.i1.ne.1) GO TO 290
	irdt=0
	time(1)=.123
c*          ��p�������, ���� �� � ����� ����
	read (8,201,err=290,end=290) line80
201   format(a80)
	back space 8
	lenlin=LEN_TRIM(line80)
	if((lenlin.eq.24.and.inblok.eq.1).or.(lenlin.eq.45.and.inblok.eq.
     *0)) go to 2000
c*          ���� ���
	if(inblok.eq.0) irc=7
	if(inblok.eq.0) go to 320
	irdt=1
	read(8,21,err=290,end=290) nrec
	go to 211
c
2000  if(inblok.eq.0) read(8,21,err=290,end=290) nrec,idate
	if(inblok.eq.1) read(8,21,err=290,end=290) nrec,(idate(i,1),i=1,7)
21    format(15i3)
211   if(nrec.gt.maxrec) irc=6
	if(nrec.gt.maxrec) go to 310
	if(irdt.eq.1) go to 28
c*          ������ �������������� ������� time(1) � ��������� ���� 
c*          ��� ������� ��������� dt
	dt=0.
	je=1
	if(inblok.eq.0) je=2
	id0=0
c*          ������ ������ ��� � ����. ���� idate0
	me=idate0(2)-1
	if(me.eq.0) go to 23
	  do 22 i=1,me
22      id0=id0+mn(i)
23    id0=id0+idate0(3)
	idy0=0
	if(idate0(1)-idate0(1)/4*4.eq.0) idy0=1
	if(me.ge.2) id0=id0+idy0
c*          �������� ������� ��� � ����. ���� idate0
	h0=idate0(4)+idate0(5)/60.+idate0(6)/3600.+idate0(7)/360000.
c*          �������� ������� ��� ��� idate
	do 26 j=1,je
	  id=0
c*          ������ ������ ��� � ����. ���� idate
	  me=idate(2,j)-1
	  if(me.eq.0) go to 25
	    do 24 i=1,me
24        id=id+mn(i)
25      id=id+idate(3,j)
	  idy=0
	  if(idate(1,j)-idate(1,j)/4*4.eq.0) idy=1
	  if(me.ge.2) id=id+idy
c*          �������� ������� ��� � ����. ���� idate
	  h=idate(4,j)+idate(5,j)/60.+idate(6,j)/3600.+idate(7,j)/360000.
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
26    continue
	time(1)=tm(1)
	if(inblok.eq.0) dt=(tm(2)-tm(1))*3600.
	if(inblok.eq.0) time(2)=dt/60.
c*          ���������� �������� �����������
28    do 29 i=1,2
	paspim(i,1)=intg1(i)      
	paspim(i,2)=intg2(i)      
	paspim(2+i,1)=rl1(i)      
29    paspim(2+i,2)=rl2(i)      
	IF(INBLOK.EQ.0) IRC=1
	IF(IFLG.LE.1) GO TO 35
C*      �PO��CK ��OKA �AHH�X
	DO 30 I=1,NREC
	READ (8,40,ERR=290,END=290)
30    CONTINUE
	GO TO 254
C
35    IF(INBLOK.EQ.0) GO TO 200
C
C*             BBO� KOHT�PA
C
	IND=1
	METPOS(1)=1
	METPOS(2)=1
	DO 60 I=1,NREC
	READ(8,40,ERR=290,END=290) INTG1,WL1,INTG2,WL2
40    FORMAT(2I5,2F9.3,2I5,2F9.3) 
	XY(1,IND)=WL1(2)
	XY(2,IND)=WL1(1)
	IF(IND.EQ.1) GO TO 50
	IF(WL1(2).EQ.XS.AND.WL1(1).EQ.YS) GO TO 50
C*       KOHTPO�� OKOH�AH�� �PA�MEHTA KOHT�PA.
	XY(1,IND)=XS
	XY(2,IND)=YS
	IND=IND+1
	METPOS(1)=METPOS(1)+1
	METPOS(METPOS(1)+1)=IND
	XY(1,IND)=WL1(2)
	XY(2,IND)=WL1(1)
50    IF(IND.EQ.1) GO TO 55
C*        KOHTPO�� �����POBAH�� TO�EK
	IF(XY(1,IND).EQ.XY(1,IND-1).AND.XY(2,IND).EQ.XY(2,IND-1))
     * IND=IND-1
55    IND=IND+1
	XS=WL2(2)
	YS=WL2(1)
60    CONTINUE
C*      BK���EH�E B MACC�B �OC�E�HE�O ��EMEHTA KOHT�PA
	XY(1,IND)=XS
	XY(2,IND)=YS
	IF(IND.GT.1.AND.XY(1,IND-1).EQ.XS.AND.XY(2,IND-1).EQ.YS)
     * IND=IND-1
	N=IND
	GO TO 250
C
C*        BBO� BEKTOPOB
C
200   IND=1
	DO 240 I=1,NREC
	READ(8,40,ERR=290,END=290) INTG1,WL1,INTG2,WL2
	XY(1,IND)=WL1(2)
	XY(2,IND)=WL1(1)
C	write(*,41)	XY(1, ind), xy(2, ind)
	WK((IND-1)*2+1)=WL2(2)
	WK(IND*2)=WL2(1)
	IND=IND+1
240   CONTINUE
C
	N=IND-1
C	write(*,41)	(XY(1, i), xy(2, i), i=1,n)
C*      �EPEC�ET �AHH�X 
250   CALL  CONVGD(V,XY,WK,N,FI,LAMBDA,dt,INBLOK+1,itypcn)
C*      �PEO�PA�OBAH�E MACC�BA XY B �OPMAT XY(N,2)
	DO 251 I=1,N
	WK(I)=XY(1,I)
251   WK(I+N)=XY(2,I)
	DO 252 I=1,N
	XY(I,1)=WK(I)
252   XY(N+I,1)=WK(N+I)
C
C*      �POBEPKA: ECT� �� E�E ��OK� �AHH�X B �A��E ???
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