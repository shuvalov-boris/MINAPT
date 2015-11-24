
      SUBROUTINE ROUTE (KOD,N,X,FX)
      DIMENSION X(N,1),FX(1),KOD(1)
      IF(KOD(2).EQ.KOD(1))KOD(2)=KOD(2)-1
      J=KOD(2)
    1 IF(J) 4,4,2
    2 DO 3 I=1,N
    3 X(I,J+1)=X(I,J)
      FX(J+1)=FX(J)
      J=J-1
      GO TO 1
    4 KOD(2)=KOD(2)+1
      RETURN
      END
