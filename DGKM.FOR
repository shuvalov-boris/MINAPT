      SUBROUTINE DGKM(X,Y,FI,LAMBDA)
C***************************************************************
C*
C*
C* �PO�PAMMA �EPEBO�A KOOP��HAT �� �PA��COB B K��OMETP�
C*       HA�A�O OTC�ETA - �P�HB�� HA �KBATOPE
C*
C* �APAMTEP�: X,Y - �O��OTA � ��POTA TO�K� B KM
C*            FI,LAMBDA - TO �E B �PA��CAX
C*
C* �P�ME�AH�E. TO�HOCT� PAC�ETA �OCTATO�HA ��� �AHH�X, CH�T� CO
C*             C��TH�KOB�X CH�MKOB B�COKO�O PA�PE�EH�� (1 KM):
C*                                        �O MEP��. �O ��POTE
C* TO�HOCT� MO�E�� �EO��A: �O��C -        -290M       0M
C*      HA 1 �PA��C        45 �PA� ��POT� +80M       23*PI*FI/360M
C* (�EO�� - �����CO��)     �KBATOP -      +110M       0M
C*                   L+ERR=L(�CT�HHA�)
C*
C***************************************************************
      REAL LAMBDA
      DATA R0/6378.245/,PI/3.1415926/,ONEDG/111.24/
      REARTH=R0*(1.-0.00335*SIN(LAMBDA*PI/180.)**2)
      RLITL=REARTH*COS(LAMBDA*PI/180.)
      X=PI*RLITL*FI/180.
      IF(X.LT.0.) X=2.*PI*RLITL+X
      Y=ONEDG*LAMBDA
      RETURN
      END