c*********************************************************************
c*
c*
c*          ������ � ������� ����� dat(idat)
c*
c*
c* ���������:   dat,idat,lentyp,lenmax,nb,nbmax, - ������� � inpsfd
c*          nbin - ����� �������� �����
c*          x,y,t - ������� ������ (���� ��� ������� ���� �����)
c*          ibgwd - ���������� ����� ������ �������� ������ � 
c*             ������� ������ dat ����� nbin
c*          n - ����� �������� ��������� ������ (��� iwk=4 - �����
c*             �������������� ����������) 
c*          iwk - ��� ������:
c*             1 - ������� ��� ����� � 1 �� 2
c*             2 - ������� ����
c*             3 - �������� (n>0) ��� ������� (n<0) ���� ��� ��������
c*             4 - �������� ��� ������������ ���������������� �����
c*             5 - �������� ����� �����
c*             6 - �������� ����� '���������������' ��������� ���
c*                ���������� � �����
c*             7 - �������� ����� ���������� �������
c*			 8 - �������� ����� ����� �� ����� �������
c*             9 - ��������� ����� ����� �������
c*          iflg - ������ ����������
c*          irc - ��� ��������:
c*             0 - ���������� ����������
c*             1 - ��� ����� � ������� dat(idat)
c*             2 - ��� ��������� ������
c*
c* ����������. ���� ������� ���������� ��� ������ ����� ����� ������,
c*    ������� ��� ������ ����������� ����� ������, ����������� 
c*    ��������� ������.
c*
c*
c*********************************************************************
	subroutine datwrk(dat,idat,lentyp,lenmax,nb,nbmax,nbin,x,y,t,
     *  ibgwd,n,iwk,iflg,irc)
	real*4 dat(1),x(1),y(1),t(1),rl
	integer*2 idat(1),intg(2), nb, nbmax, nbin, nknot, nc
	integer*4 lentyp(2,1),lenmax,ibgwk4,ibeg,ibg,ibg2,len,
     *  ibgfr,ib
	equivalence (rl,intg)
	data lencdt/53/
	go to (10,20,30,60,90,100,140, 200, 210),iwk
c
c*          ����� ���� ����� (����� ������� -> ���� �������)
c
10    irc=1
	ibeg=lentyp(1,nbin)+6
	nknot=idat(ibeg-3)
	nwd=idat(ibeg-2)
	if(nwd-nwd/2*2.eq.1) nwd=nwd+1
	dl=(nknot-1)*14+2-nwd
C	write(6,1) lenmax, nb, dl
C    1 format(' lenmax = ', I9, ', nb = ', I4, ', dl = ', I7)
	if(lenmax.lt.lentyp(1,nb+1)+dl-1) return
c*          ����� � ������� idat(dat) ����
	nwd=-nwd
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibeg,nwd,irc)
	ibeg=ibeg+nknot*4
	nwd=dl-nwd
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibeg,nwd,irc)
c*          �������� ���������� � ������� �����
	ibeg=lentyp(1,nbin)
	idat(ibeg)=0
	idat(ibeg+1)=0
	idat(ibeg+4)=0
	idat(ibeg+5)=0
	lentyp(2,nbin)=2
	return
c
c*          �������� �����
c
20    ibeg=lentyp(1,nbin)
	nwd=ibeg-lentyp(1,nbin+1)
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibeg,nwd,irc)
	return
c
c*          ����������/�������� ����� ��� ��������� � �����
c*     � ������ ���������� ����� ���������� x,y,t - ���������� �������
c*     � ������ �� ���������� - x(n,2) - ����������, y(n,2) - 
c*     ��������, t(n,2) - ������ ������ � �������� (������ ��� 
c*     "���������������" ���������)
c
30    nknot=idat(lentyp(1,nbin)+3)
	itpblk=lentyp(2,nbin)
	if(n.le.0) go to 35
c*          ���� ���������
	irc=1
	dl=n*8
	if(idat(lentyp(1,nbin)+5).eq.2) dl=dl+n*4
	if(itpblk.eq.2) dl=n*18
	if(lentyp(1,nb+1)+dl-1.gt.lenmax) return
c*          ����� � ������� dat(idat) ��� ����
35    ibeg=lentyp(1,nbin)+6+(ibgwd-1)*2
	n2=n*2
	ntimes=3
	if(itpblk.eq.0) ntimes=4
	if(itpblk.eq.0.and.dl.eq.n*12) ntimes=6
	do 50 i=1,ntimes
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibeg,n2,irc)
	if(n.lt.0) go to 50
	if(itpblk.eq.0) go to 45
c*          ��������� ����
	  do 40 j=1,n
	  if(i.eq.1) dat(ibeg/2+j)=x(j)
	  if(i.eq.2) dat(ibeg/2+j)=y(j)
	  if(i.eq.3) dat(ibeg/2+j)=t(j)
40      continue
	  go to 50
c*          ��������� ��������
45      continue
	  do 47 j=1,n
	  if(i.eq.1) dat(ibeg/2+j)=x(j)
	  if(i.eq.2) dat(ibeg/2+j)=x(j+n)
	  if(i.eq.3) dat(ibeg/2+j)=y(j)
	  if(i.eq.4) dat(ibeg/2+j)=y(j+n)
	  if(i.eq.5) dat(ibeg/2+j)=t(j)
	  if(i.eq.6) dat(ibeg/2+j)=t(j+n)
47      continue
50    ibeg=ibeg+n2+nknot*2
	if(itpblk.eq.0) go to 55
c
c*          ������ ������ ������� ��� ������������ �������
	ibeg=ibeg-(ibgwd-1)*2
	nwd=n*12
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibeg,nwd,irc)
c*          ������������ ���� �������
55    idat(lentyp(1,nbin)+3)=nknot+n
	if(itpblk.eq.2) idat(lentyp(1,nbin)+5)=0
	return
c
c*          ������� ��� ������������� ��������������� ����
c
c*  �� �����: x(n) - �������������� ����������,
c*          y(n) - �� ������� (�����, integer*2)
c*          t(52) - ��������� �����������
c*          iflg - ���� ����� ������� �������
c*          ibgwd - ����� �����, ��� ��������� ��������������
c
c*          �������� ����������� ������ � �����
60    irc=2
	if(nbin.eq.0.and.nb.eq.nbmax) return
	irc=1
	leninf=lencdt*2+3*n+6
	lenold=0
	if(nbin.ne.0) lenold=lencdt*2+idat(lentyp(1,nbin)+3)*3+6
	idlen=leninf-lenold
	if(lenmax.lt.lentyp(1,nb+1)+idlen-1) return
	irc=0
c
	if(lenold.eq.0) nb=nb+1
	if(lenold.eq.0) nbin=nb
	if(lenold.eq.0) lentyp(1,nb+1)=lentyp(1,nb)+leninf
	if(lenold.eq.0) go to 65
c*          ���������� �������� ������ � ������
	ibgwk4=112+lentyp(1,nbin)
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibgwk4,idlen
     *  ,irc)
c
c*          ������� ������ � ����
c
65    lentyp(2,nbin)=4
c*          ������������ ���������� �����
	ibg=lentyp(1,nbin)+6
	ibg2=ibg/2
	dat(ibg2-2)=0.123
	idat(ibg-3)=n
	idat(ibg-4)=iflg
	idat(ibg-2)=ibgwd
	idat(ibg-1)=0
c*          ���� ���������� �����������
	do 70 i=1,lencdt
70    dat(ibg2+i)=t(i)
c*          ���� ���������� � ������� fxbest-����������
	ibg2=ibg2+lencdt
	dat(ibg2)=1.e30
	do 75 i=1,n
75    dat(ibg2+i)=x(i)
c*          ���� �������� ����������
	n2=n/2
	ibg2=ibg2+n
	if(n2.eq.0) go to 85
	do 80 i=1,n2
80    dat(ibg2+i)=y(i)
85    if(n-n2*2.eq.0) return
c*          n - �������� �����
	rl=y(n2+1)
	idat(ibg+lencdt*2-1+n*3)=intg(2)
	return
c
c*          �������� ����� �����
c
90    irc=2
	if(nb.ge.nbmax) return
	ibg=lentyp(1,nbin)
	len=lentyp(1,nbin+1)-lentyp(1,nbin)
	irc=1
	ibgfr=lentyp(1,nb+1)
	if(ibgfr+len-1.gt.lenmax) return
c*          � ������ ���� ����� ��� �����������
	irc=0
	do 95 i=0,len-1
95    idat(ibgfr+i)=idat(ibg+i)
	nb=nb+1
	lentyp(1,nb+1)=lentyp(1,nb)+len
	lentyp(2,nb)=lentyp(2,nbin)
	return
c
c*          ��������� ����� ��������� ��� ���������� � �����
c
c*
c* ���������: ��������� �������� x(ibgwd,4),y(ibgwd,2):
c*    � x - x � y ���������� �������� � ��� ���������� ��������
c*    � y - ������ � �������� ������� � ������, ��������� � ������
c*    ��������� ��������
c*    n - ����� ��������, n<=ibgwd
c*    t - �����, �������� � ������� �����
c*    iflg - �������� 6-�� �������� �������
c*    nbin - <=0:
c*          ��� ��������� (itypcn=iabs(nbin))
c*          >0:
c*          ����� ����� � �������� ����������� ������
c*
c*
100   continue
c*          �������� ������� ��������� ������ � �����
	irc=2
	if(nb.ge.nbmax.and.nbin.le.0) return
	ibg=lentyp(1,nb+1)
	if(nbin.gt.0) ibg=lentyp(1,nbin)
	len=n*8
	if(nbin.le.0) len=len+6
	if(nbin.gt.0) iflg=idat(ibg+5)
	if(iflg.eq.2) len=len+n*4
	irc=1
	if(lentyp(1,nb+1)+len-1.gt.lenmax) return
	irc=0
	no=0
	if(nbin.gt.0) no=idat(ibg+3)
	nn=no+n
	if(nbin.gt.0) go to 105
c*          ��������� ������� �����
	dat(ibg/2+1)=t(1)
	idat(ibg+2)=iabs(nbin)
	idat(ibg+3)=n
	idat(ibg+4)=0
	idat(ibg+5)=iflg
	go to 109
c*          ���������� ���� nbin
105   ntime=4
	if(iflg.eq.2) ntime=6
	ib=no+no+ibg+6
	do 107 i=1,ntime
	call  movdat(idat,lentyp,lenmax,nb,nbin,ib,n+n,irc)
107   ib=ib+nn+nn
c*          ��������� �����
109   iy=ibgwd
	iv1=iy+ibgwd
	iv2=iv1+ibgwd
	ibg2=ibg/2+3+no
	n2=nn+nn
	n3=n2+nn
	n4=n3+nn
	n5=n4+nn
	do 110 i=1,n
	dat(ibg2+i)=x(i)
	dat(ibg2+nn+i)=x(iy+i)
	dat(ibg2+n2+i)=x(iv1+i)
110   dat(ibg2+n3+i)=x(iv2+i)
c
	if(iflg.ne.2) go to 130
c*          ������� ������ � ��������
	do 120 i=1,n
	dat(ibg2+n4+i)=y(i)
120   dat(ibg2+n5+i)=y(ibgwd+i)
c
130   if(nbin.gt.0) idat(ibg+3)=nn
	if(nbin.gt.0) return
	nb=nb+1
	lentyp(2,nb)=0
	lentyp(1,nb+1)=lentyp(1,nb)+len
	return
c*
c*          ��������� ����� ���������� �������
c*
c* ���������: x(n),y(n) - ���������� ����� �������
c*          t(1) - �����, �������� � ������� �����
c*          iflg - �������� 6-�� �������� �������
c*          nbin - ��� ��������� (itypcn)
c*          
c* ����������. ����������� ��������� ���� ��������������� ������
c*    � ������ ����� �������� ��������� ���������. ���� �� ����� 
c*    ����� �������������
c*
140   continue
c*          �������� ������� ��������� ������ � ����� � ������
	irc=2
	if(nb.ge.nbmax) return
	irc=1
	ibg=lentyp(1,nb+1)
	len=n*4+6+2
	if(ibg+len-1.gt.lenmax) return
c*          � ����� ������ �����
	irc=0
	np=1
	nf=1
	idat(ibg+6)=1
	if(n.le.1) go to 170
c
c*          ��������� ����� ������� �� ���������
c
	nplook=0
150   nplook=nplook+1
	if(nplook.eq.n) go to 170
	np=np+1
	if(x(np-1).ne.x(np).or.y(np-1).ne.y(np)) go to 150
c*          ����� ���������
	if(nplook.eq.n-1) np=np-1
	if(nplook.eq.n-1) go to 170
c*          ���������� ��������� ������ ������ ���������
	nf=nf+1
	idat(ibg+5+nf)=np
c*          ����� �������� ���������
	ie=n-nplook-nf+1
	do 160 i=1,ie
	x(np+i-1)=x(np+i)
160   y(np+i-1)=y(np+i)
c
c*          ��������� ��������� ����� ����� � ����� ������������� 
c*                      �����
	nplook=nplook+1
	go to 150
c
170   continue
c
c*          ��������� ������� �����
c
	dat(ibg/2+1)=t(1)
	idat(ibg+2)=nbin
	idat(ibg+3)=np
	idat(ibg+4)=nf
	idat(ibg+5)=iflg
c
c*          �������� ��������� � ����
	ibg=(ibg+nf+6)/2
	do 180 i=1,np
	dat(ibg+i)=x(i)
180   dat(ibg+i+np)=y(i)
c*          ��������� ������� ������
	nb=nb+1
	lentyp(1,nb+1)=lentyp(1,nb)+np*4+6+nf+nf-nf/2*2
	lentyp(2,nb)=1
	return
C
C       �������� ����� ����� �������� n, ��������������
C	  ���������� (�� ������) �� ����� ������� (���� nbin)
C
200   nknot = n + 1
      irc=3
	if (n.le.1) goto 222
      irc=2
	if(nb+1.ge.nbmax) goto 222
	ibg=lentyp(1,nbin)
	irc=1
	ibgfr=lentyp(1,nb+1)
	len=6+nknot*6+(nknot-1)*12
	if(ibgfr+len-1.gt.lenmax) goto 222
c*          � ������ ���� ����� ��� �����������
	irc=0
	NC=idat(ibg+3)
	idat(ibgfr)  =0
	idat(ibgfr+1)=0
	idat(ibgfr+2)=idat(ibg+2)
	idat(ibgfr+3)=nknot
	idat(ibgfr+4)=0
	idat(ibgfr+5)=0

       nf=idat(ibg+4)
	 nf1=nf
	 if(nf-nf/2*2.eq.1) nf1=nf1+1
	ibg=ibg+6+nf1
	ibgfr=ibgfr+6

	do i = 1, nknot - 1
	  idx = ((NC - 1) * (i - 1) / (nknot-1)) + 1
	  dat(ibgfr/2+i) = dat(ibg/2+idx)
	  dat(ibgfr/2+i+nknot) = dat(ibg/2+idx+NC)
	enddo
	dat(ibgfr/2+nknot) = dat(ibgfr/2 + 1)
	dat(ibgfr/2+nknot*2) = dat(ibgfr/2 + nknot + 1)	
	call Parametrize(nknot, dat(ibgfr/2+1), dat(ibgfr/2+nknot+1),
     *  dat(ibgfr/2+nknot*2+1))
	write(*, *) 'parametrize', dat(ibgfr/2+1), dat(ibgfr/2+nknot+1),
     *  dat(ibgfr/2+nknot*2+1) 
      nb=nb+1
      lentyp(1,nb+1)=lentyp(1,nb)+6+nknot*6+(nknot-1)*12
      lentyp(2,nb)=2
	nbin=nb
      return
C
C	��������� ����� ����� � ��� ��������������
C*		x(n),y(n) - ���������� �����
c*          nbin - ��� ��������� (itypcn)
C
  210	  nknot = n
  	  write(*, *) 'DATWRK: nb =', nb, ', n =', n, ', nbin =', nbin,
     *    ', len =', len, ', LT(1) =', lentyp(1, nb),
     *    ', LT(2) =', lentyp(2, nb)
        irc = 3
	  if (n.le.1.or.nbin.le.0) goto 222
c*          � ������ ���� ����� ��� �����������
        irc=0
	  ibg=lentyp(1,nbin)/2+3
        write(*, *) 'datwrk: ��������� ������'
	  do i = 1, nknot
	    dat(ibg+i) = x(i)
	    dat(ibg+i+nknot+1) = y(i)
	  enddo
	  dat(ibg+nknot+1) = x(1)
	  dat(ibg+(nknot+1)*2) = y(1)
	  nknot=nknot+1
	  call Parametrize(nknot, dat(ibg+1), 
     *    dat(ibg+nknot+1), dat(ibg+nknot*2+1))
        write(*, *) 'DW DAT: ' 
	  write(*, '(3F12.6)') (dat(ibg+i), dat(ibg+nknot+i)
     *    , dat(ibg+nknot*2+i), i = 1, nknot)
        return

  222 end


C**********************************************************************
C*
C*						Parametrize
C*
C*		�������������� ������� ����� � ������� ������ � �������
C*
C*	���������:
C*		N - ���������� �����
C*		X(N) - x-���������� �����
C*		Y(N) - y-���������� �����
C*		T(N) - t-���������� ����� (�����)
C*
C**********************************************************************
      SUBROUTINE Parametrize(N, X, Y, T)
	  IMPLICIT NONE
	  INTEGER*2, INTENT(IN) :: N
	  REAL, INTENT(IN)    :: X(N), Y(N)
	  REAL, INTENT(OUT)   :: T(N)

	  INTEGER i

	  T(1) = 0.
	  DO i = 2, N
	    T(i) = T(i-1) + sqrt(
     *	  (x(i) - x(i - 1))**2 + (y(i) - y(i - 1))**2)
	  END DO

	END