      SUBROUTINE MIN05(KOD,N,X,FXM,D,DM,DS,EST,F,MINLN,OUTN,A,B,IER)
      DIMENSION X(N,1),A(N,N),D(N),DM(N),DS(N),B(N,4),FXM(1)
      DIMENSION KOD(6),EST(6)
      NM1=N-1
      IF(KOD(6).NE.0) GO TO 3
C     ÇAÄAHÈE HA×AËÜHOÉ ÏPßMOÉ MATPÈÖÛ
      DO 2 J=1,N
      DO 1 I=1,N
    1 A(I,J)=0.
    2 A(J,J)=1.
      DO 33 I=1,N
   33 B(I,1)=D(I)
C     HA×AËO ÈTEPAÖÈOHHOÃO ÏPOÖECCA
    3 F1=FXM(1)
      DO 4 J=1,N
    4 B(J,2)=X(J,1)
      SUM=0.
      F2=FXM(1)
      IPR=0
C     HA×AËO ÏOÈCKA HAÏPABËEHÈß ÓÁÛBAHÈß ÔÓHKÖÈÈ
      DO 7 J=1,N
   41 DO 5 I=1,N
    5 B(I,4)=B(I,2)+A(I,J)*B(J,1)
      CALL F(B(1,4),FX,1,1)
      KOD(4)=KOD(4)-1
      IF(KOD(4))321,321,51
C     BÛXOÄ ÏO ËÈMÈTÈPOBAHHOMÓ KOËÈ×ECTBÓ BÛ×ÈCËEHÈÉ ÔÓHKÖÈÈ
  321 F2=FX
   32 IER=3
   22 IF(F2-F1)222,220,220
  222 CALL ROUTE(KOD,N,X,FXM)
      FXM(1)=F2
      DO 223 I=1,N
  223 X(I,1)=B(I,2)
      GO TO 220
   51 IF(FX-F2)15,52,52
   52 DO 50 I=1,N
   50 B(I,4)=B(I,2)-A(I,J)*B(J,1)
      CALL F(B(1,4),FX,1,1)
      KOD(4)=KOD(4)-1
      IF(KOD(4))321,321,59
   59 IF(FX-F2)53,54,54
   53 DO 55 I=1,N
   55 A(I,J)=-A(I,J)
      GO TO 15
C     KOHEÖ ÏOÈCKA HAÏPABËEHÈß ÓÁÛBAHÈß ÔÓHKÖÈÈ
C     KPÈTEPÈÉ ÏPEÄEËÜHOÃO ÄPOÁËEHÈß ØAÃA
   54 IF(B(J,1).GT.DM(J)) GO TO 57
      IPR=IPR+1
      GO TO 7
   57 B(J,1)=B(J,1)/10.
      GO TO 41
C     ÏOÈCK MÈHÈMÓMA HA J-OM HAÏPABËEHÈÈ
   15 FX=F2
      EST(2)=B(J,1)
      CALL MINLN(KOD,N,B(1,2),F2,A(1,J),F,EST,B(1,4),IER)
      IF(IER.EQ.1)GO TO 22
      IF (KOD(4).LE.0) GO TO 32
      IF(EST(2)) 48,47,48
   47 B(J,1)=B(J,1)/2
      GO TO 58
   48 IF(B(J,1).GT.EST(2)) B(J,1)=EST(2)
C     OÏPEÄEËEHÈE HAÏPABËEHÈß MCOXP C MAKCÈMAËÜHÛM ÈÇMEHEHÈEM ÔÓHKÖÈÈ
   58 R=FX-F2
      IF(R-SUM)7,7,6
    6 MCOXP=J
      SUM=R
    7 CONTINUE
      IF (IPR.NE.N) GO TO 8
C     BÛXOÄ ÏO ÏPEÄEËÜHOMÓ ÄPOÁËEHÈÞ ØAÃA
      IER=4
      GO TO 22
    8 DO 10 J=1,N
   10 B(J,3)=2.0*B(J,2)-X(J,1)
      CALL F(B(1,3),F3,1,1)
      KOD(4)=KOD(4)-1
      IF(KOD(4))32,32,9
    9 R=F3-F1
      CALL ROUTE(KOD,N,X,FXM)
C     ÈCCËEÄOBAHÈE HAÏPABËEHÈß MCOXP
      IF (R)11,60,60
   11 R=2.0*(F1-2.0*F2+F3)*((F1-F2-SUM)/R)**2
      IF(R-SUM)12,60,60
   12 IF(MCOXP-N)13,115,115
   13 DO 14 J=MCOXP,NM1
      J1=J+1
      DO 14 I=1,N
   14 A(I,J)=A(I,J1)
  115 R=0.
      DO 16 I=1,N
      A(I,N)=B(I,2)-X(I,1)
   16 R=A(I,N)**2+R
      R1=SQRT(R)
      EST(2)=R1
C     ÏOÈCK MÈHÈMÓMA B HAÉÄEHHOM HAÏPABËEHÈÈ
      CALL MINLN(KOD,N,B(1,2),F2,A(1,N),F,EST,B(1,4),IER)
      IF(EST(2))71,72,71
   71 B(N,1)=EST(2)
      GO TO 61
   72 B(N,1)=R1/2
   61 FXM(1)=F2
      DO 62 I=1,N
   62 X(I,1)=B(I,2)
   73 IF(IER.EQ.1)GO TO 221
      IF(KOD(4).GT.0)GO TO 19
      IER=3
      GO TO 221
   60 R=F2-F3
      IF(R) 61,61,63
   63 FXM(1)=F3
      DO 64 I=1,N
   64 X(I,1)=B(I,3)
C     ÏPOBEPKA KPÈTEPÈß OKOH×AHÈß C×ETA ÏO TO×HOCTÈ
   19 CALL OUTN(KOD,N,X,FXM,DS,EST,IER)
  220 IF(KOD(2).EQ.1)GO TO 225
C     HAXOÆÄEHÈE ÏPÈPAÙEHÈß ÏO APÃÓMEHTAM È BEËÈ×ÈHÛ CÄBÈÃA
  221 D0=0.
      DO 25 I=1,N
      R=ABS(X(I,1)-X(I,2))
      D(I)=R
      R=R*R
   25 D0=D0+R
      EST(2)=SQRT(D0)
  225 IF(IER.LT.0) GO TO 3
      RETURN
      END
