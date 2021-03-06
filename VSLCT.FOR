C**********************************************************************C
C*                                                                    *C
C*    ����� ��������� �� nbo-������ �����, ��� �� ��������            *C
C*   �� ������� �������� ������ ����� nbc �������� �����������        *C
C*   ������� ������� � ������� �������� ������� �� ������� ��         *C
C*   ���������� �� ����� dstmax                                       *C
C*                                                                    *C
C* �APAMETP�: dat,idat,lentyp - ������� � ��������� inpsfd            *C
C*            nbg - ������ ����� � ���� ������ ���������              *C
C*            x,y.vx.vy - ���������� � ���������� ���������� �������� *C
C*            nv - ����� ���������� �������� ��������                 *C
C*            dstmax - ������������ ���������� �� ������� �� �������  *C
C*                                                                    *C
C* ����������. �������� �������� ����� ���� ���������� ����� �������� *C
C*          ��� ���� �������.                                         *C
C*                                                                    *C
C**********************************************************************C
	SUBROUTINE vslct(dat,idat,lentyp,nbc,nbo,nbg,x,y,vx,vy,nv,dstmax)
	REAL X(1),Y(1),vx(1),vy(1),dat(1),xy1(2),xy2(2)
	INTEGER   idat(1),ibgi(2)
	integer*4 lentyp(2,1),ibgc,ibgc2,nbg(2,nbo),ibg
C*          ��������� ����� � ��������
	ibgc=lentyp(1,nbc)
	npnt=idat(ibgc+3)
	nf=0
	if(lentyp(2,nbc).eq.1) nf=idat(ibgc+4)
	ibgc2=(ibgc+nf)/2+4
	nv=0
	ic=1
	if(dat(ibgc2).ne.dat(ibgc2+npnt-1).or.dat(ibgc2+npnt).ne.
     *dat(ibgc2+npnt*2-1)) ic=0
C
	do 100 i=1,nbo
c*          ��������� ����� ���������
	ibg=nbg(1,i)
	nvi=nbg(2,i)
	  do 10 j=2,3
10      ibgi(j-1)=nvi*j
c
		  DO 100 j=1,nvi
	    SS=(dat(ibg+j-1)-dat(ibgc2))**2+(dat(ibg+nvi+j-1)-
     *    dat(ibgc2+npnt))**2
	  inout=0
		DO 80 k=1,npnt-1
C*        O�PE�E�EH�E M�H�MA��HO�O PACCTO�H�� �O OTPE�KA
		XY1(1)=dat(ibg+j-1)-dat(ibgc2+k-1)
		XY1(2)=dat(ibg+nvi+j-1)-dat(ibgc2+npnt+k-1)
		XY2(1)=dat(ibgc2+k)-dat(ibgc2+k-1)
		XY2(2)=dat(ibgc2+npnt+k)-dat(ibgc2+npnt+k-1)
		AL=XY1(1)*XY2(1)+XY1(2)*XY2(2)
		IF(AL.LE.0.) GO TO 80
C*       TO�KA �E��T "�PABEE" �EBO�O KOH�A OTPE�KA
		INDX=2
		IF(AL.GT.XY2(1)**2+XY2(2)**2) INDX=1
		GO TO (52,54), INDX
C*        M�H�MA��HOE PACCTO�H�E - �O TO�K� ������� ����� �������
52          S2=(dat(ibg+j-1)-dat(ibgc2+k))**2+(dat(ibg+nvi+j-1)-
     *      dat(ibgc2+npnt+k))**2
		IF(S2.GT.SS) GO TO 80
C*       �PO��O��O �MEH��EH�E TEK��E�O M�H�MA��HO�O PACCTO�H��
		SS=S2
		inout=1
		if(k.eq.npnt-1.and.ic.eq.0) inout=0
		GO TO 80
C*        HEO�XO��MO PACC��T�BAT� TO�K� - �POEK��� ������� HA
C*          TEK���� OTPE�OK
54          DX=XY2(1)
		DY=XY2(2)
		INDX=3
		IF(DY.EQ.0.) INDX=2
		IF(DX.EQ.0.) INDX=1
		GO TO (55,60,65), INDX
C*          DX=0.
55          XC=dat(ibgc2+k-1)
		YC=-dat(ibg+nvi+j-1)
		GO TO 70
C*          DY=0.
60          XC=dat(ibg+j-1)
		YC=dat(ibgc2+npnt+k-1)
		GO TO 70
C*        C��TAEM �O O��E� �OPM��E
65          DIV=1./(DX*DX)+1./(DY*DY)
		XC=((dat(ibg+j-1)/DY+dat(ibg+nvi+j-1)/DX)/DY-
     *      (dat(ibgc2+npnt+k-1)/DY-dat(ibgc2+k-1)/DX)/DX)/DIV
		YC=((dat(ibgc2+npnt+k-1)/DY-dat(ibgc2+k-1)/DX)/DY+
     *      (dat(ibg+j-1)/DY+dat(ibg+nvi+j-1)/DX)/DX)/DIV
C*       KOHTPO�� HA MEH��EE PACCTO�H�E
70          S2=(dat(ibg+j-1)-XC)**2+(dat(ibg+nvi+j-1)-YC)**2
		IF(S2.GT.SS) GO TO 80
		IF(S2.EQ.SS.AND.k.EQ.1) GO TO 80
		SS=S2
		inout=1
80          CONTINUE
	if(inout.eq.0.or.ss.gt.dstmax**2) go to 100
c*          ��������� �������� ������� �������� ���������� ��
c*           ���������� ����� ������� � ��������� �������� �� �������
	nv=nv+1
	x(nv)=dat(ibg+j-1)
	y(nv)=dat(ibg+nvi+j-1)
	vx(nv)=dat(ibg+ibgi(1)+j-1)
	vy(nv)=dat(ibg+ibgi(2)+j-1)
100         CONTINUE
	RETURN
	END