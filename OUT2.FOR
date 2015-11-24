      SUBROUTINE OUT2(KOD,N,X,FX,DS,EST,IER)
        DIMENSION KOD(6),EST(6)
        DIMENSION X(N,2),FX(2),DS(N)
      
	  DF=EST(6)
        KOD(5)=KOD(5)-1
        IER=-1

        IF(KOD(5).LE.0) THEN
	    IER = 2
          GOTO 2
	  END IF
        
	  R=ABS(FX(1)-FX(2))
        IF(R-DF)1,1,2

1       DO 5 I=1,N
          R=ABS(X(I,1)-X(I,2))
          IF(R-DS(I))5,5,2
5       CONTINUE
        IER=0

2       RETURN
      END