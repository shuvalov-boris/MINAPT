      SUBROUTINE MINL5(KOD,N,X,FX,S,F,EST,Y,IER)
      DIMENSION X(N),Y(N),S(N)
      DIMENSION EST(1),KOD(1)
      DATA C/1.618034/
      DATA C1/0.381966/
      DATA C2/0.618034/
C      BÛÁOP ÄËÈHÛ OTPEÇKA
      LIM=MIN0(KOD(3),KOD(4))
      L1=LIM
      H=0.
      DO 66 I=1,N
   66 H=H+S(I)*S(I)
      GN=SQRT(H)
      H=EST(2)*EST(5)/GN
      DL=EST(4)/GN
      A=0.
      B=H
      H1=0.
      K=LIM-1
    5 DO 2 I=1,N
    2 Y(I)=X(I)+B*S(I)
      CALL F(Y,FY,1,1)
      LIM=LIM-1
      IF(LIM)50,50,55
   55 IF(FY-FX)3,4,4
    3 A=A+H1
      ALFA1=A+H
      H1=H
      H=C*H
      FX=FY
      B=ALFA1+H
      GO TO 5
    4 RAZN=B-A
      IF(K-LIM)101,6,21
    6 FZ=FX
   19 ALFA1=A+C1*RAZN
      DO 8 I=1,N
    8 Y(I)=X(I)+ALFA1*S(I)
      CALL F(Y,FX,1,1)
      LIM=LIM-1
      IF(LIM.LE.0.AND.FX.LT.FZ)GO TO 53
      IF(LIM)17,17,18
   18 IF(FX-FZ)21,16,16
   16 B=ALFA1
      RAZN=B-A
      IF(RAZN-DL)17,17,19
   17 FX=FZ
      GO TO 44
   21 ALFA2=A+C2*RAZN
      DO 10I=1,N
   10 Y(I)=X(I)+ALFA2*S(I)
      CALL F(Y,FY,1,1)
      LIM=LIM-1
      IF(LIM)40,40,7
C      ÏOÈCK MÈHÈMÓMA ÔÓHKÖÈÈ METOÄOM ÇOËOTOÃO CE×EHÈß
    7 IF(FX-FY)12,11,11
   11 IF(RAZN-DL)40,40,20
   20 RAZN=C2*RAZN
      A=ALFA1
      ALFA1=ALFA2
      FX=FY
      GO TO 21
   12 IF(RAZN-DL)53,53,15
   15 RAZN=C2*RAZN
      B=ALFA2
      ALFA2=ALFA1
      FY=FX
      ALFA1=A+C1*RAZN
      DO 13 I=1,N
   13 Y(I)=X(I)+ALFA1*S(I)
      CALL F(Y,FX,1,1)
      LIM=LIM-1
      IF(LIM)53,53,7
   40 D0=ALFA2
  100 FX=FY
   58 DO 41 I=1,N
   41 X(I)=X(I)+D0*S(I)
  103 EST(2)=D0*GN
  101 KOD(4)=KOD(4)-L1+LIM
      IER=0
      RETURN
C      BÛXOÄ ÏO ËÈMÈTÈPOBAHHOMÛ KOËÈ×ECTBÓ BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ
   50 IF(K-LIM)101,51,45
   45 IF(FX-FY)53,57,57
   57 D0=B
      GO TO 100
   53 D0=ALFA1
      GO TO 58
   51 IF(FX-FY)44,42,42
   42 D0=H
      GO TO 100
   44 EST(2)=0.
      GO TO 101
      END
