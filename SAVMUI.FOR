c*********************************************************************
c*
c*
c*	���������� ����������� �������� ������� ������� �� 
c*		������� �������� � ���������� �����
c*
c* ��p����p�:   irc-��� ����p���: 0-��p������� ����p�����,
c*			1-������ �p� ���p����, ��� �p� ������
c*
c* ����������. ����������� �������� ���������� �� common-�������� 
c*		���������������� ���������
c*
c*
c*********************************************************************
	subroutine savmui(x, fx)
	common/mfcmin/ kod1(8),iw,nnow,n,m,is,maxit,maxfx,ngrad,fxopt,
     *  ngrd,ityp1,istps,kd4,kd5
	common /mfcom/ xopt(500)
	irc=0
	open (19,FILE='optxfx.sav',ERR=30)
	write(19,10,ERR=30) n,fxopt
10	format(i6,e13.6)
	write(19,20,ERR=30) (xopt(i),i=1,n)
20	format(6e13.6)
	close(19)
	return
30	irc=1
	return
	end
