C**********************************************************************C
C*                                                                    *C
C*                                                                    *C
C*        ÏPOÃPAMMA TECTÈPOBAHÈß ÏOÄÏPOÃPAMM CÏËAÉ-OÏÈCAHÈß KOHTÓPOB  *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      REAL XY(10,2),T(20),XK(5),YK(5),TK(5),WK(100,3),TKOLD(5),V(100,2)
     *,C(4,3,2)
      INTEGER IWK(5),METPOS(20)
      DATA XY/0.,2.,3.,4.,5.,6.,6.,3.,-2.,0.,0.,2.,2.5,3.,2.,0.,-2.,-4.,
     *-5.,0. /,XK/0.,5.02,5.99,3.11,3.28/,
     *YK/0.,2.18,-0.05,-3.98,-6.74/,IFLG/0/,
     *N/9/,NKNOT/5/,TKOLD/0.,1.,2.,3.,4./,KFIRST/1/
      DO 10 I=1,N
10    XY(N+I,1)=XY(I,2)
	T(1)=0.
	DO 5 I=2,N
5	T(I)=T(I-1)+SQRT((XY(I,1)-XY(I-1,1))**2+(XY(N+I,1)-XY(N+I-1,1))
     *  **2)
	TK(1)=0.
	TK(2)=6.479
	TK(3)=8.71
	TK(4)=14.25
	TK(5)=21.375
c	tk(6)=t(n)
	CALL  SPAPP(XY,T,N,XK,YK,TK,NKNOT,C,NKNOT-1,EPS,WK,IFLG,IRC)
	print 66, IRC,XK,YK,C,eps
66	FORMAT(2X,'RC,XK,YK,C,e=',I3,10F6.2,5(/,3X,10F6.2))
c	DO 20 I=1,20
c20	T(I)=26.935/30.*I
	CALL  ICSEVU(TK,XK,NKNOT,C,NKNOT-1,T,WK,n,IRC)
	CALL  ICSEVU(TK,YK,NKNOT,C(1,1,2),NKNOT-1,T,WK(1,2),n,IRC)
	print 73,tk
73	format('tk=',10f6.2)
	print 76, (WK(I,1),I=1,20),T
76	FORMAT(3X,10F6.2)
	print 76, (WK(I,2),I=1,20)
C      DO 15 I=1,5
C      IFLG=0
C      CALL  RDCONV(NAMFIL,141.,34.,XY,V,N,METPOS,42000.,WK,IFLG,IRC)
C      WRITE(6,16) N,IRC,(METPOS(J),J=1,5)
C16    FORMAT(3X,'N,IRC,METPOS(5)=',7I5)
C      IF(IRC.GT.3) RETURN
C15    CONTINUE
C     CALL TPREP(XY,T,N,XK(1),YK(1),TK,NKNOT,NBEG,NEND,NUSL,NUSR,NDOUBL,
C    *EPS,WK,IWK,IFLG,IRC)
C     WRITE(6,56) NKNOT
C6    FORMAT(3X,'NKNOT=',I5)
C     WRITE(6,66) T,(TK(I),I=1,5    ),(WK(I,1),I=1,5    ),(IWK(I),I=1,
C    *5    ),XK,YK,NBEG,NEND,NUSL,NUSR,NDOUBL,IRC
C6    FORMAT(3X,'T,TK',10F6.2,5F6.2,/,3X,'DIST,NPOS',5F6.2,5I4,/,
C    *3X,'XK,YK',10F6.2,/,3X,'NBEG,NEND,NUSL,NUSR,NDOUBL,IRC',6I5)
C     CALL  TRECLC(XY,XY(1,2),T,N,TK,TKOLD,NKNOT,IWK(1),WK)
C     WRITE(6,76) N
C6    FORMAT(3X,'N    =',I5)
C     WRITE(6,66) T,(TK(I),I=1,5    ),(WK(I,1),I=1,5),(IWK(I),I=1,
C    *5    ),XK,YK,NBEG,NEND,NUSL,NUSR,NDOUBL,IRC
C     WRITE(6,86) T,XY
C6    FORMAT(3X,'T,X,Y',15F6.2)
C     RETURN
C     END
C     SUBROUTINE ICSMOU
C     PRINT 1
C     FORMAT(1X,'ICSMOU')
      RETURN
      END