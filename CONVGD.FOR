      SUBROUTINE CONVGD(V,XY,WK,N,FI,LAMBDA,TIME,KINF,ITYPCN)
C***********
C*
C*     �EPEBO� �AHH�X �� �EO�PA���ECKO� C�CTEM� KOOP��HAT B �EKAPTOB��
C*      KOM�OHEHT� CKOPOCT� V B CM/CEK, KOOP��HAT� TO�EK XY B 10 KM
C*
C*
C* �APAMETP�: V(N,2) - BEKTOPA CKOPOCT� (KINF=2) - B�XO�
C*            XY(2,N) - MACC�B KOOP��HAT B �PA��CAX(BXO�).HA B�XO�E
C*                       B 10KM.
C*            WK(2,N) - TO�K� B KOTOP�E �EPEMECT���C� MAPKEP� CKOPOCT�
C*            FI,LAMBDA - HA�A�O KOOP��HAT (�EO�PA���ECK�X)
C*            TIME - BPEM� �EPEME�EH�� MAPKEPOB CKOPOCT� �� TO�EK
C*                 XY(1,I),XY(2,I) B TO�K� WK(1,I),WK(2,I)
C*            KINF - T�� TPE��EMO�O �PEO�PA�OBAH��:
C*              2 - PACC��TAT� KOOP��HAT� TO�EK KOHT�PA
C*              1 - PACC��TAT� KOOP��HAT� �O�O�EH�� MAPKEPOB CKOPOCT�
C*                  � CAM� CKOPOCT�
C*            itypcn - ��� �������������� ( ������ � pgm rdconv)
C*
C***********************************************************************
      REAL V(N,1),XY(2,1),WK(2,1),LAMS,LAMBDA
      GO TO (140,10),KINF
C
C*       KINF=2
C
C
C*            �PEO�PA�OBAH�E KOOP��HAT KOHT�PA B KM
C
10    if(itypcn.eq.0) go to 15
C*  	�������������� � ������������� �������������� ����������
      do 5 i=1,n
      xy(1,n)=xy(1,n)-fi			 ! index: n or i ???
5     xy(2,n)=xy(2,n)-lambda
      return
C*		�������������� � ��������� ����������       
15    CALL  DGKM(X0,Y0,FI,LAMBDA)
      DO 130 I=1,N
C       WRITE(6,126) XY(1,I),XY(2,I)
C26     FORMAT(3X,'XY=',2E10.3)
        LAMS=XY(2,I)
        FII=XY(1,I)
        CALL  DGKM(XY(1,I),XY(2,I),FII,LAMS)
        CALL  DGKM(X0,Z,FI,LAMS)
        XY(1,I)=(XY(1,I)-X0)/10.
        XY(2,I)=(XY(2,I)-Y0)/10.
130   CONTINUE
C     WRITE(6,26) ((XY(J,I),J=1,2),I=1,N)
C6    FORMAT(3X,10F7.2)
      close(25)
      RETURN
C             KINF=1
C
C*               PAC�ET CKOPOCTE�
C
140   DO 150 I=1,N
      CALL  DGKM(x1,y1,xy(1,i),xy(2,i))
      call  dgkm(x1,z,xy(1,i),wk(2,i))
      CALL  DGKM(x2,y2,wk(1,i),wk(2,i))
      V(I,1)=(x2-x1)/TIME*100000.
      V(I,2)=(y2-y1)/TIME*100000.
150   CONTINUE
C     WRITE(6,26) ((V(J,I),I=1,2),J=1,N)
C     WRITE(6,26) ((XY(J,I),J=1,2),I=1,N)
      go to 10
      end