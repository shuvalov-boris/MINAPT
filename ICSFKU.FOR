C************************************************************************
C*
C*
C*	���������-���������� ��� ��������������� ������-�������������
C*	����� ���������� ������-�������������
C*
C* ���������: �� ��, ��� � � ��������� ICSFKU ��� IMSL, �����:
C*	      XY(NX,2) - X � Y ���������� ����������������� �������
C*	      Y(NKNOT,2) -  X � Y ���������� �����
C*
C* ��������: ���������� ��������� ���������� ������������� ������
C*   ���������� �������. ��������� ICSFK1 � ICSFK2 ���������� �� 
C*   ��������� ICSFKU ��� IMSL ������ ���������� COMMON-��������.
C*   ������ ������������� �������� �����.
C*
C*
C************************************************************************
	SUBROUTINE ICSFKU(T,XY,NX,MODE,TK,NKNOT,Y,C,IC,EPS,WK,IER)
	common /nksave/ nks
	REAL T(1),XY(NX,2),TK(1),Y(1),C(IC,3,2),WK(1)

	CALL ICSFK1(T,XY,NX,MODE,TK,NKNOT,Y,C,IC,EPS1,WK,IER1)
	CALL ICSFK2(T,XY(1,2),NX,MODE,TK,NKNOT,Y(nks+1),C(1,1,2),IC,EPS2,
     *  WK(NX*(NKS+6)+1),IER2)

	IER=IER1
	IF(IER2.NE.0) IER=IER2
c	print 4,ier1,ier2,ier
c4	format(1x,'ier1,ier2,ier=',3i5)
	IF(IER.NE.0) GO TO 10
C*          ���������� ���������� ��������
	EPS=SQRT(EPS1*EPS1+EPS2*EPS2)
c	print 6,EPS1,EPS2,EPS
c6	FORMAT(1X,'EPS1,EPS2,EPS=',3E15.5)
10	RETURN
	END