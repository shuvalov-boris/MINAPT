c*********************************************************************
c*
c*
c*		����� ����� ������� ������ idat
c*
c* ���������:   idat,lentyp,lenmax,nb - ������� � ��������� inpsfd
c*		nbin - ����, � ������� ������������ ���������
c*		ibeg - ������ ����������� ���������
c*		nwd - ����� ���� (<0 - ���������,>0 - �����������)
c*		irc : 0 - ���������� ����������
c*		      1 - ��� ����� ��� ����������� ������
c*
c* ����������. ���� ����� ����, �� ������� ���������� ����������
c*	�������, ����� ��� �����, �� ���� ��������� ��������� �
c*	����� ������ ���������� �� 1 ������
c*
c*
c*********************************************************************
	subroutine movdat(idat,lentyp,lenmax,nb,nbin,ibeg,nwd,irc)
	integer*2 idat(1), nb, nbin, irc
	integer*4 lentyp(2,1),ibeg,lastwd,ibg,lenmax, nwd
	irc=0
	lastwd=lentyp(1,nb+1)
	if(nwd.gt.0) go to 50
c
c*		������� ������ �� �������
c
	nbo=nb
	if(lentyp(1,nbin)-lentyp(1,nbin+1).ne.nwd) go to 20
c*		������� ���� ����
	nb=nb-1
	if(nb.lt.nbin) return
	do 10 i=nbin,nb+1
	lentyp(2,i)=lentyp(2,i+1)
10	lentyp(1,i)=lentyp(1,i+1)+nwd
c*		�������� ������ ������
20	ibg=ibeg-nwd
	do 30 i=ibg,lastwd-1
30	idat(i+nwd)=idat(i)
	if(nbo.ne.nb) return
	go to 70
c
c*		���������� ������
c
50	if(lastwd+nwd-1.gt.lenmax) irc=1
	if(lastwd+nwd-1.gt.lenmax) return
	if(ibeg.eq.lastwd) go to 70 
	do 60 i=1,lastwd-ibeg
60	idat(lastwd+nwd-i)=idat(lastwd-i)
c*		��������� ���� ������
70	do 80 i=nbin+1,nb+1
80	lentyp(1,i)=lentyp(1,i)+nwd
	return
	end