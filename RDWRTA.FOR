$NOTSTRICT
c*********************************************************************
c*
c*
c*          ������ � ������ ������ � �����. ����� ����������� ��� 
c*    ������������ ���������.
c*
c* ���������:   name - ��� ������ (������ � �����������)
c*          icninf - ������ ���������� �� ���������� ������
c*          istinf - ������ ������������ ����������
c*          lentyp - ���� � ������ ������ � ������� idat
c*          idat - ����� ������
c*          iflg - ������:
c*                0 - ������
c*                1 - ������
c*                2 - ������ � ���������� ������
c*          irc - ��� ��������:
c*                0 - ���������� ����������
c*                1 - ��� iflg=2 ����� � ������ name ���
c*                 ����������. ������ �� ������������
c*                2 - ��� iflg<2 ������ ��� ��� ������ ��� ��� 
c*                ��������
c*                3 - ������ ������ ��� ������. ������ 
c*                � ������� ������������
c*                4 - ��� ����� ��� ������ �� ������� ��������
c*
c*
c*********************************************************************
	subroutine rdwrta(name,icnin,istin,lentyp,idat,iflg,irc)

	character where*38
	real fi, lambda
      integer*2 idat(100000),idate0(7),icninf(130), istinf(40), nknot
     *  ,icnin(130), istin(40), nb, itypcn,iflgsp
	integer*4 lentyp(2,101)

	integer*2 bcio(2)
	integer*4 len,nbrw
	character*(*) name
	data maxlen/16384/
	equivalence 
     *  (icninf(5),where),(istinf(3),fi),(istinf(5),lambda),
     *  (istinf(8),idate0),(istinf(31),itypcn),
     *  (istinf(36),nknot),(istinf(39),iflgsp)
c
c*          �������� ������
c
C	call  bopen(bcio,name,0,ir)
    8 format(A, A)
      open(unit=12, file=name, err=200, action='read', status='old')
C	if(iflg.ne.2) go to 200
	if(iflg.eq.0) go to 10
	goto 190
C	if(iflg.eq.2) go to 210
C	open(unit=12, file=name, err=200, action='write', status='old')
C	go to 100
c
c*          ������ ������
10    read(12, *, ERR=91) where
2     format(A38)
      read(12, *, ERR=91) fi, lambda
      read(12, *, ERR=91) idate0
      read(12, *, ERR=91) itypcn, iflgsp
      read(12, *, ERR=91) nknot

	nb=istinf(7)
	len=(lentyp(1,nb+1)-1)*2
	irc=0
	if(len.eq.0) go to 190
c*        �������� ����� ������ ��� ������ � ������ ����� maxlen ����
	nbrw=len/maxlen
	lenl=len-nbrw*maxlen
	if(nbrw.le.0) go to 190
C	do 20 i=1,nbrw
C	call  bread(bcio,idat,maxlen,irc)
C	if(irc.ne.0) go to 220
C20    continue        
C50    if(lenl.le.0) go to 190
C	call  bread(bcio,idat,lenl,irc)
C	if(irc.ne.0) go to 220
C	go to 190
c
c*          ������ � �����
c
100   call  bwrite(bcio,icninf,icninf(1)*2,ir)
	if(ir.eq.-1) go to 230
	if(ir.ne.0) go to 220
	call  bwrite(bcio,istinf,istinf(1)*2,ir)
	if(ir.eq.-1) go to 230
	if(ir.ne.0) go to 220
C	call  bwrite(bcio,lentyp,(istinf(2)+1)*8,ir)
	if(ir.eq.-1) go to 230
	if(ir.ne.0) go to 220
	nb=istinf(7)
	len=(lentyp(1,nb+1)-1)*2
	irc=0
	if(len.eq.0) go to 190
c*        �������� ����� ������ ��� ������ � ������ ����� maxlen ����
	nbrw=len/maxlen
	lenl=len-nbrw*maxlen
	if(nbrw.le.0) go to 150
	do 120 i=1,nbrw
	call  bwrite(bcio,idat,maxlen,irc)
	if(irc.eq.-1) go to 230
	if(irc.ne.0) go to 220
120   continue        
150   if(lenl.le.0) go to 190
	call  bwrite(bcio,idat,lenl,irc)
	if(irc.eq.-1) go to 230
	if(irc.ne.0) go to 220
C   call  bclose(bcio,ir)
190     close(12, err=92)
      istin=istinf
	icnin=icninf
	return
200   irc=2
	return
210   irc=1
	return
220   irc=3
	return
230   irc=4
	return
   91 write(*, *) 'RDWRTA: ������ ������ ������'
      return
   92 write(*, *)	'RDWRTA: ������ �������� �����'
	end