      SUBROUTINE  MIN03(KOD,N,X,FX,D,DM,DS,EST,F,MINLN,OUTN,A,B,IER)
      REAL*4    D(N),DS(N),DM(N),X(N,1),FX(1),A(N,1), RK, EST(1),B(1),
     *B3,R,R1,D2,D1,S,FX2,FX3,RB,S1,R2,TWO/2.0/
      INTEGER KOD(1)
      DATA R2/1.4142132/
C     R2=DSQRT(0.2D1)
      N1=N+1
      M1=0
      L=KOD(1)
      RK=0.0
      DO 989 J=1,N
      B3=DM(J)
  989 RK=RK+B3*B3
      IF(KOD(6).NE.1)GO TO 665
  777 IF(KOD(2)-N1)665,137,137
C     COPTÈPOBKA ÏO BOÇPACTAHÈÞ
  390 IER=3
  650 RETURN
  137 DO 230 J=1,N1
      DO 231 I=1,N
      A(I,J)=X(I,J)
  231 CONTINUE
  230 B(J)=FX(J)
      GO TO 135
  665 R=N
      DO 743 I=1,N
  743 A(I,1)=X(I,1)
      B(1)=FX(1)
C     ÏOCTPOEHÈE HA×AËÜHOÃO CÈMÏËEKCA
      R1= SQRT(R+1.0)-1.0
      D2=EST(2)/(R*R2)
      D1=D2*(R1+R)
      D2=D2*R1
      DO 303 J=2,N1
      DO 302 I=1,N
      A(I,J)=D2+A(I,1)
  302 CONTINUE
  303 A(J-1,J)=D1+A(J-1,1)
      KOD(2)=N1
      KM=KOD(2)
  346 DO 382 I=2,N1
      CALL F(A(1,I),B(I),1,1)
      KOD(4)=KOD(4)-1
      IF(KOD(4))390,390,382
  382 CONTINUE
 2087 IF(M1-1)135,2088,2088
 2088 CALL OUTN(KOD,N,X,FX,DS,EST,IER)
      IF(IER)135,135,650
C HAXOÆÄEHÈE MÈHÈMAËÜHOÃO ×ÈCËA
  135 J=1
      DO 601 I=2,N1
      IF(B(J)-B(I))601,601,602
  602 J=I
  601 CONTINUE
      IL=J
      IF(IL-1)234,100,234
  234 DO 441 I=1,N
      B3=A(I,1)
      A(I,1)=A(I,IL)
  441 A(I,IL)=B3
      B3=B(1)
      B(1)=B(IL)
      B(IL)=B3
C  HAXOÆÄEHÈE MAKCÈMAËÜHOÃO ×ÈCËA
  100 J=1
      DO 10 I=2,N1
      IF(B(J)-B(I))12,10,10
   12 J=I
   10 CONTINUE
      M=J
      IF(M-N1)138,39,138
  138 DO 563 I=1,N
      B3=A(I,N1)
      A(I,N1)=A(I,M)
  563 A(I,M)=B3
      B3=B(N1)
      B(N1)=B(M)
      B(M)=B3
C HAXOÆÄEHÈE ÖEHTPA TßÆECTÈ
   39 DO 401 I=1,N
      S=0.0
      DO 19 J=1,N
   19 S=S+A(I,J)
  401 A(I,L+1)=S/N
C OTOÁPAÆEHÈE
      IL=1
   20 DO 21 I=1,N
   21 A(I,L+2)=A(I,L+1)+1.0 *(A(I,L+1)-A(I,N1))
      CALL F(A(1,L+2),FX2,1,1)
      KOD(4)=KOD(4)-1
      IF(KOD(4))390,390,22
   22 IF(FX2-B(1))25,31,31
C PACTßÆEHÈE
   25 DO 26 I=1,N
      A(I,L+3)=A(I,L+1)+TWO *(A(I,L+2)-A(I,L+1))
26    CONTINUE
      CALL F(A(1,L+3),FX3,1,1)
      KOD(4)=KOD(4)-1
      IF(KOD(4))390,390,27
   27 IF(FX3-B(1))28,30,30
   28 CALL  ROUTE(KOD,N,X,FX)
  596 IF(KM.EQ.KOD(1)) KM=KM-1
      J=KM
    1 IF(J)4,4,2
    2 DO 3 I=1,N
    3 A(I,J+1)=A(I,J)
      B(J+1)=B(J)
      J=J-1
      GO TO 1
    4 KM=KM+1
      GO TO (597,598,599,600),IL
  597 DO 29 I=1,N
      A(I,1)=A(I,L+3)
   29 X(I,1)=A(I,L+3)
      FX(1)=FX3
      B(1)=FX3
  101 CALL OUTN(KOD,N,X,FX,DS,EST,IER)
      IF(IER)100,100,650
   30 CALL  ROUTE(KOD,N,X,FX)
      IL=2
      GO TO 596
  598 DO 301 I=1,N
      A(I,1)=A(I,L+2)
  301 X(I,1)=A(I,L+2)
      FX(1)=FX2
      B(1)=FX2
      GO TO 101
   31 DO 34 I=2,N
      IF(FX2-B(I))901,34,34
   34 CONTINUE
   33 IF(FX2-B(N1))36,1000,1000
  901 IL=3
      GO TO 596
  599 DO 45  I=1,N
   45 A(I,2)=A(I,L+2)
      B(2)=FX2
      GO TO 100
C CÆATÈE
   36 DO 111 I=1,N
  111 A(I,N1)=A(I,L+2)
      B(N1)=FX2
 1000 DO 801 I=1,N
  801 A(I,L+4)=A(I,L+1)+0.5*(A(I,N1)-A(I,L+1))
      CALL F(A(1,L+4),FX4,1,1)
      KOD(4)=KOD(4)-1
      IF(KOD(4))390,390,51
   51 IF(FX4-B(N1))38,139,139
   38 IL=4
      GO TO 596
  600 DO 41 I=1,N
   41 A(I,2)=A(I,L+4)
      B(2)=FX4
      GO TO 100
  139 RB=0.0
      DO 833 J=1,N1
      S=0.0
      DO 899 I=1,N1
      S1=A(I,J)-A(I,L+1)
  899 S=S+S1*S1
      IF(RB-S)855,833,833
  855 RB=S
  833 CONTINUE
      IF(RB.GT.RK)GO TO 211
      IER=4
      GO TO 650
  211 M1=M1+1
      DO 200 J=2,N1
      DO 200 I=1,N
C PEÄÓKÖÈß
  200 A(I,J)=A(I,1)+0.5  *(A(I,J)-A(I,1))
      DO 55 I=2,N1
      CALL F(A(1,I),B(I),1,1)
C     IF(B(I)-B(1))52,53,53
C  53 CONTINUE
C     GO TO 135
      CALL  ROUTE(KOD,N,X,FX)
      DO 54 J=1,N
   54 X(J,1)=A(J,I)
   55 FX(1)=B(I)
C     B(1)=B(I)
      GO TO 2088
C     DEBUG INIT
      END