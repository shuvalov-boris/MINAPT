c*********************************************************************
c*
c*  ����� �������� ������� fun �������� ����������� �����������
c*
c* ���������: method - ����� ������ �����������
c*          kod,n,x,fx,d,dm,ds,parm,minl,out,a,b - ��������� ������ -
c*                ������� � ������� 17 - ���������� ssplib
c*          numopx - ������ � �������� ������������� ����������.
c*                ���� numopx(i)<0, �� ��� ��������, ��� ���������� �
c*              �������� � ��������� numopx(i-1):iabs(numopx(i)) -
c*              �� ��������������. numopx(i)=0 - ����� ������.
c*          itp - ������ ���������� �������:
c*              itp<=0 - ������ �� ����������
c*              itp=1(11) - ���������� ��������� � �������� ��������
c*                      ����������, ������� (x,fx),��� �������� � 
c*                    ����� ������.
c*              itp=2(12) - �������� ���, ����� �������������
c*                    �������� ����������
c*              itp=3(13) - �������� ���
c*             itp>=10 - ���������� �������� ��������� savmui �����
c*                itp/10 ��������
c*          istp - ����� ��������, ����� ������� ����������
c*               ������������� ����������
c*          irc - ��� ��������:
c*              irc=0 - ���������� ����������
c*              irc>0 - ��� �������� ������ �����������
c*              irc=-1,(-2) - ����������� ������� �������� �� 
c*                ��������� �������� (����������)
c*              irc=-10 - ������� ��������� ������ numopx, ��� 
c*                      �������������� ���������� <1
c*
c*
c* ��������: �������������� ������� ������������ fun(x,fx). ���� �����
c*           ����������� �� ��������� � ������������ �� �����������
c*         ��������� ������� ���������(grad) � �������� ������������
c*         ������� �����(gesse), �� ��� �������������� ��������.
c* �����������. ����� ���������� X �� ������ ���� >100. ��� ������
c*         ����� ����������� ���������� ����� ������� minfun �������
c*         common-������� /mfcom/ �������� 5*n ��������������� ����.
c* ���������. 1.��� ���������� ������ ������ � ����������� �������� ���
c*         ���������� ������� ��������� �/��� ������� ����� 
c*           ���������� ������ �p�������� ngrad=-1 � ��������� minfun
c*          2.��� itp>10 - ���������� ���������� ������������:
c*          savmui(xbest,fxbest)
c*          3.��� ��������������� �������� �������� ����������
c*          ������������� ��� ������� ��������� ������:
c*          f1 - �������� ����������
c*          f2 - �������� �� ��������� ��������
c*          f3 - ������� ���������� ��������� savmui
c*          ��������� savmui ����� ���� ��� ����������������, ��� 
c*          � ���������� ������. ��� ���������� ��������� ����-
c*          �������� ���������� ���������� ��������� ������� �����-
c*          ��� ������� � ��������� � ����� optxfx.sav . ��� 
c*          ���������� ����������� �������� ������������� ����������-
c*          �� ���������� resxfx(xbest,fxbest,irc) � ����� ��������
c*          irc=0 ��� ���������� ���������� � irc=1 ��� ������� �����
c*
c*
c*********************************************************************
	SUBROUTINE minfun(METHOD,KOD,N,X,FX,D,DM,DS,PARM,MINL,OUT,A,B,IRC,
     *numopx,itp,ISTP)
	external subf
	REAL*4 X(1),FX(1),D(1),DM(1),DS(1),PARM(8),A(1),B(1),FXOPT,XOPT
	REAL*4 NGRD,mlt
	integer*4 ibinar(100),two
	INTEGER KOD(8),NUMopx(1)
	COMMON/mfcmin/KOD1(8),IW,nnow,NN,M,IS,MAXIT,MAXFX,NGRAD,FXOPT,
     *NGRD,ityp1,istps,kd4,kd5
     */mfcom/XOPT(500)
     */grcom/ig,ing,ifst,nggr,mlt,ibinar
	data two/2/
C		����� ������� ����������
      open(unit=15, file='minfun.log', action='write')
	write(15, *) 'N = ', N
	write(15, *) 'FI: ', X(1)
	write(15, *) 'LAMBDA: ', X(2)
	write(15, *) 'NODES:'
	nknot = (n - 5) / 2
	do i = 3, 3 + nknot - 1
		write(15, *) X(I), X(I+nknot)
	enddo
	write(15, *) 'A0 = ', X(2 + 2*nknot + 1)
	write(15, *) 'Rm = ', X(2 + 2*nknot + 2)
	write(15, *) 'R0 = ', X(2 + 2*nknot + 3)
c*          ��������� �p��p���� grad
	ifst=1
	nggr=15
	mlt=1.0
	ig=0
	ing=0
	do 166 i=1,100
166   ibinar(i)=0
	itype=itp-itp/10*10
	istps=itp/10
	if(istps.lt.0) istps=0
c
c*         ������������ ������� �������������� ����������
c
	if(n.gt.100.and.itype.gt.0) print 1
1     format(1x,'��������������: ����� ���������� ������ �������������')
	do 2 i=1, 100
2       ibinar(i)=0
	irc=-10
	if(n.lt.1) return
	do 4 i=1,n
4       xopt(n+i)=x(i)
	nnow=n
	if(numopx(1).eq.0) go to 122
c*         �� ��� ���������� ������� fun ��������������. ���������
c*           ������ �������-�������� (ibinar) � ������ ��������������
c*         ���������� X
	nxo=0
c*           �������� ������������ ������ numopx
	n1=n+1
	if(numopx(1).lt.0) go to 300
	do 5 i=2,n1
	  if(numopx(i).eq.0) go to 6
	  if(iabs(numopx(i)).gt.n) go to 300
	  if(iabs(numopx(i)).le.iabs(numopx(i-1))) go to 300
5     continue
6     nnow=0
	ie=i-1
	do 10 i=1,ie
	  nx=numopx(i)
	  ibg=nx
	  ien=nx
	  if(nx.gt.0) go to 7
c*           ���� ���� ����������� ����������
	  ibg=nxo+1
	  ien=iabs(nx)
7       do 8 j=ibg,ien
	    nw=j/31
		if(nw*31.ne.j) nw=nw+1
		idn=(nw)*31-j
8         ibinar(nw)=ibinar(nw)+two**idn
	  if(nxo+1.eq.ibg) go to 10
C*      ����p ��������p����� ��p�������
	  iend=ibg-1
	  ibg=nxo+1
	  do 9 j=ibg,iend
		nnow=nnow+1
9         x(nnow)=xopt(n+j)
10      nxo=ien
	if(ien.eq.n) go to 12
c*      ����p ��������p����� ��p�������
	ibg=ien+1
	do 11 i=ibg,n
	  nnow=nnow+1
11      x(nnow)=xopt(n+i)
12    if(nnow.lt.1) go to 300
c---------------------------------------------------------------------
122   MAXFX=KOD(4)
	MAXIT=KOD(5)
	IS=ISTP
	ityp1=itype
	nn=n
	DO 13 I=1,8
   13 KOD1(I)=KOD(I)
	NGRD=1.E30
	IW=-1
	write (15, *) 'MINFUN: call FUN'
	call fun(xopt(n+1),fx)
	write (15, *) 'MINFUN: FUN is completed'
	if (fx(1).eq.0.) then 
	  write(33, *) 'FX is NULL!'
	  write(15, *) 'FX is NULL! TERMINATED!'
	  RETURN
	endif
c*          ��������� ���������� ����������� � ������ xopt
	do 30 i=1,n
30    xopt(i)=xopt(n+i)
	fxopt=fx(1)
	kd4=0
	kd5=0
c
	NGRAD=-1
	if(itype.ge.1) PRINT 32,N,FX(1),(XOPT(N+I),I=1,N)
   32 FORMAT (3X,'N=',I3,3X,999(2X,'FX&X=',2X,6D10.3,/))
	if(itype.ge.2) PRINT 34,KOD,D(1),DM(1),DS(1)
   34 FORMAT (3X,'KOD=',8I8,/,3X,' D=',D10.3,' DM=',D10.3,' DS=',D10.3)
	if(itype.ge.2) PRINT 35,PARM
   35 FORMAT (3X,'PARM=',4d10.3)
	CALL gettim(ihr,imin,isec,i100th)
	timeb=3600.*ihr+60.*imin+isec+0.01*i100th
	CALL METHOD(KOD1,N,X,FX,D,DM,DS,PARM,SUBF,MINL,OUT,A,B,IRC)
	CALL gettim(ihr,imin,isec,i100th)
	timend=3600.*ihr+60.*imin+isec+0.01*i100th
	time=timend-timeb
	if(kd5.ne.0) irc=-1
	if(kd4.ne.0) irc=-2
	MX=MAXFX-KOD1(4)-kd4
	MT=MAXIT-KOD1(5)-kd5
	if(itype.ge.1) WRITE(6,40) MT,MX,NGRAD,FXOPT,NGRD
   40 FORMAT(5X,I5,5X,I5,7X,I5,4X,D18.10,2X,E18.10)
	if(itype.ge.1) PRINT 43,(XOPT(I),I=1,N)
   43 FORMAT(1X,'X=',2X,5D14.7)
	if(itype.ge.1) WRITE(6,37) IRC,time
37    FORMAT(3X,'KO� BO�BPATA',2X,I4,' BPEM� PA�OT� ',F10.2,' CEK.')
	if(irc.lt.0) print *,'  ������� ���������� ����������'
	do 50 i=1,n
50    x(i)=xopt(i)
	fx(1)=fxopt
	if(itp.ge.10) then
		call  savmui(x,fx)
		write(15,*) 'MINFUN: Line185: call SAVMUI'
	endif
300   close(15)
      return
	END

c*************************************************************
c*
c*          SUBF
c*
c*
c*************************************************************
	SUBROUTINE subf(X,FX,GR,M)
	REAL*4 X(1),FX(1),GR(1),FXSAVE,ACC,NGRD,mlt
	REAL*8 NGNOW
	COMMON/mfcmin/KOD1(8),IW,NNOW,N,L,ISTP,MAXIT,MAXFX,NGRAD,fxopt,
     *NGRD,itype,istps,kd4,kd5
     */mfcom/ XOPT(500)
     */grcom/ ig,ing,ifst,nggr,mlt,ibinar(100)
	INTEGER NM(2),idold/0/,isav/0/,f1/15104/,f2/15360/,f3/15616/
	integer*4 ibinar,maska,s2(31)/16#40000000,16#20000000,16#10000000,
     *16#08000000,16#04000000,16#02000000,16#01000000,16#00800000,
     *16#00400000, 
     *16#00200000,16#00100000,16#00080000,16#00040000,16#00020000,
     *16#00010000,
     *16#00008000,16#00004000,16#00002000,16#00001000,16#00000800,
     *16#00000400,
     *16#00000200,16#00000100,16#00000080,16#00000040,16#00000020,
     *16#00000010,
     *16#00000008,16#00000004,16#00000002,16#00000001/,maskas
	 write(15, *) 'SUBF'
	write(*, *) 'SUBF'
	IDELT=MAXIT-KOD1(5)
	IF(IW.EQ.-1.and.itype.ge.1) WRITE(6,660)
660   FORMAT(3X,'  ��C�O ',2X,'  ��C�O  ',2X,'   ��C�O  ',7X,'�HA�EH�E',
     *10X,' HOPMA')
	IF(IW.EQ.-1.and.itype.ge.1) WRITE(6,662)
662   FORMAT(3X,'�TEPA���',2X,'B��. �-��',2X,'�PA��EHTOB',7X,'  �-��  ',
     *10X,' �PA�-TA')
	IW=0
	IF(NGRAD.NE.-1) GO TO 100
c*          acc - �������� ��������� �����
	acc=epsilon(acc)
	DO 99 I=1,n
	  xopt(n*2+I)=100.e0*ACC
	  xopt(n*3+i)=0.
99      xopt(n*4+I)=5.
	ngrad=0
100   CONTINUE
c
c*          ����� �� ������� ����������
c
C101   call  key(1,icodsm)
C	if(icodsm.eq.0) go to 106
c*          ��������� ����������. ��������� ��� �������
C	call  key(0,icodsm)
c*          ���������� ��� 
C	if(icodsm.eq.f3) isav=1
C	if(icodsm.eq.f1) kd4=kod1(4)
C	if(icodsm.eq.f1) kod1(4)=0
C	if(icodsm.eq.f2) kd5=kod1(5)
C	if(icodsm.eq.f2) kod1(5)=0
C	go to 101
c
106   IF(IDOLD.EQ.IDELT) GO TO 1
	ID=IDELT-IDELT/ISTP*ISTP
	IF(ID.NE.0)GO TO 1
	IF(IDELT.EQ.0)GOTO 1
	MFX=MAXFX-KOD1(4)
	if(itype.ge.2) PRINT 90 ,IDELT,MFX,NGRAD,FXSAVE,NGRD
  90  FORMAT(5X,I5,5X,I5,7X,I5,4X,D18.10,2X,E18.10)
	if(itype.ge.3) PRINT 92,(XOPT(N+I),I=1,N)
  92  FORMAT(10X,'�HA�EH�� �EPEMEHH�X',1X,6E14.7)
1     continue
c*          ����� ��������� savmui, ���� ���������
	IF(IDOLD.EQ.IDELT.and.isav.eq.0) GO TO 922
	if(istps.eq.0.and.isav.eq.0) go to 922
	if(istps.ne.0) id=idelt-idelt/istps*istps
	if(id.ne.0.and.isav.eq.0) go to 922
	if(idelt.eq.0.and.isav.eq.0) go to 922
	call  savmui(xopt,fx)
      write(15,*) 'MINFUN: Line256: call SAVMUI'
	isav=0
c
922   IDOLD=IDELT
	if(n.eq.nnow) go to 97
c*            ������������ ������� ������� ����������
	nwd=1
	nxw=1
	nxopt=0
	maskas=0
	do 95 i=1,n
	maska=maskas+s2(nxw)
	if(ibinar(nwd).ge.maska) go to 93
c*            ���������� � ������� i - ��������������
	nxopt=nxopt+1
	xopt(n+i)=x(nxopt)
	go to 94
c*            ���������� � ������� i �������������
93    maskas=maska
c*            �������� ������ �������� ������� ibinar
94    nxw=nxw+1
	if(nxw.le.31) go to 95
c*            ������ ����� ��������
	nxw=1
	nwd=nwd+1
	maskas=0
95    continue
	go to 108
c
c*        ���������� ������� ��p�������
97    do 102 i=1,n
102   xopt(n+i)=x(i)
c
108   NM(2)=M
	NM(1)=N
	IF(M-M/2*2.EQ.0) GO TO 2
	write(15, *) 'SUBF: call FUN'
	CALL FUN(XOPT(N+1),FX)
	write(15, *) 'SUBF: FUN returned FX =', FX(1)
	fxsave=fx(1)
	IF(FXOPT.LE.FX(1)) GO TO 2
	FXOPT=FX(1)
	DO 111 I=1,N
111     XOPT(I)=XOPT(N+I)
    2 IF(M/2-M/4*2.EQ.0) GO TO 3
      write (15,*) 'MINFUN: Line295: call GRAD'
	CALL GRAD(XOPT(N+1),FX,GR,NM,xopt(n*2+1),xopt(n*3+1),xopt(n*4+1))
C*     PAC�ET HOPM� �PA��EHTA
	NGNOW=0.D0
	DO 20 I=1,N
20    NGNOW=NGNOW+GR(I)*GR(I)
	NGNOW=DSQRT(NGNOW)
	GN=NGNOW
	IF(GN.LT.NGRD) NGRD=GN
C
	NGRAD=NGRAD+1
3     IF(M.GT.3) then
		CALL GESSE(XOPT(N+1),FX,GR,NM,xopt(n*2+1),xopt(n*3+1),
	*xopt(n*4+1))	
		write(*,*) 'MINFUN: line312: call GESSE'
	ENDIF
	ig=0
	RETURN
	END
