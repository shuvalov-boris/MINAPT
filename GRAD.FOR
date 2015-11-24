C**********************************************************************C
C*                                                                    *C
C*                                                                    *C
C*        ÏPOÃPAMMA ×ÈCËEHHOÃO PAC×ETA ÃPAÄÈEHTA È                    *C
C*      BEPXHEÃO TPEÓÃOËÜHÈKA MATPÈÖÛ ÃECCE (ÏO CTOËÁÖAM)             *C
C*                                                                    *C
C* ÏAPAMETPÛ:  X(N) - MACCÈB OÏTÈMÈÇÈPÓEMÛX ÏEPEMEHHÛX PAÇMEPHÃCTÈ N  *C
C*        FX(L) - MACCÈB ÇHA×EHÈÉ ÖEËEBOÉ ÔÓHKÖÈÈ                     *C
C*             GR - MACCÈB CO CËEÄÓÞÙEÉ CTPÓKTÓPOÉ:                   *C
C*                N ÝËEMEHTOB - ÄËß ÃPAÄÈEHTA                         *C
C*                N(N+1)/2 ÝËEMEHTOB ÄËß MATPÈÖÛ ÃECCE                *C
C*           N(N+1)/2 ÝËEMEHTOB - PAÁO×ÈÉ MACCÈB                      *C
C*           NM(2): NM(1)=N NM(2)=M (CM. ÏPOÃPAMMÓ SUBF)              *C
C*           STP - MACCÈB ,ÈMEÞÙÈÉ CËEÄÓÞÙÓÞ CTPÓKTÓPÓ:               *C
C*           N ÝËEMEHTOB - ØAÃ PAC×ETA ÃPAÄÈEHTA                      *C
C*        DGR(N) - OÖEHKÈ ÄÈAÃOHAËÜHÛX ÝËEMEHTOB MATPÈÖÛ ÃECCE        *C
C*        MULT(N) - BECOBÛE KOÝÔÔÈÖÈEHTÛ ÄËß ØAÃA: ÏPÈ PAC×ETE        *C
C*        MATPÈÖÛ ÃECCE ØAÃ ÓMHOÆAETCß HA COOTBETCTBÓÞÙÈÉ KOÝÔÔÈÖÈEHT *C
C*                                                                    *C
C* AËÃOPÈTM PAC×ETA ÃPAÄÈEHTA :    ÄËß KAÆÄOÉ ÏEPEMEHHOÉ              *C
C*        OCÓÙECTBËßETCß PAC×ET TO×HOCTÈ ÃPAÄÈEHTA, OÁÓCËOBËEHHÛÉ     *C
C*        TO×HOCTÜÞ ÏPEÄCTABËEHÈß ×ÈCËA (OHA PABEHA 100%/BET1) È      *C
C*        TO×HOCTÈ ËÈHEÉHOÉ AÏÏPOKCÈMAÖÈÈ ÖEËEBOÉ ÔÓHKÖÈÈ (OHA PABHA  *C
C*        100%/BET2). OÖEHKA BTOPOÉ ÏPOÈÇBOÄHOÉ ÏPOÈCXOÄÈT ×EPEÇ      *C
C*        MAXGR BÛ×ÈCËEHÈÉ ÃPAÄÈEHTA. ØAÃ PAÇHOCTHOÉ CXEMÛ ÏOÄÁÈPAETCß*C
C*        TAK, ×TOÁÛ BÛÏOËHßËOCÜ ÓCËOBÈE: BET2>BET1 È BET2<BET1*3.    *C
C*        HAÈËÓ×ØAß TO×HOCTÜ PAC×ETA ÁÓÄET ÏPÈ TAKOM ØAÃE, KOÃÄA      *C
C*        BET1=BET2(T.E. ÏOÃPEØHOCTÜ,OÁÓCËOBËEHHAß TO×HOCTÜÞ ÏPEÄ-    *C
C*        CTABËEHÈß ×ÈCËA È ÏOÃPEØHOCTÜ ËÈHEÉHOÉ CXEMÛ PAC×ETA        *C
C*        KOMÏOHEHTÛ ÃAÄÈEHTA COBÏAÄAÞT).TO×HOCTÜ PAC×ETA ÃPAÄÈEHTA   *C
C*        PABHA MAX(100%/BET1 È 100%/BET2). ECËÈ OHA PABHA 100%/BET2  *C
C*        È > 10% ,TO ÃPAÄÈEHT C×ÈTAETCß ÏO CXEME BTOPOÃO ÏOPßÄKA     *C
C*        TO×HOCTÈ                                                    *C
C* AËÃOPÈTM PAC×ETA BEPXEHÃO TPEÓÃOËÜHÈKA MATPÈÖÛ ÃECCE:              *C
C*        ÄËß KAÆÄOÃO ÝËEMEHTA MATPÈÖÛ C×ÈTAETCß TO×HOCTÜ,OÁÓCËOBËE-  *C
C*        HHAß TO×HOCTÜÞ ÏPEÄCTABËEHÈß ×ÈCËA. ECËÈ ÝTA TO×HOCTÜ       *C
C*        BÛØE 0.1%,TO ØAÃ ÄEËÈTCß HA ×ÈCËO MLT, ECËÈ HÈÆE 1%,TO      *C
C*        ÓMHOÆAETCß HA ÝTO ×ÈCËO. ECËÈ ÏOËÓ×EHHAß MATPÈÖA HE ßBËß-   *C
C*        ETCß ÏOËOÆÈTEËÜHO OÏPEÄEËEHHOÉ, TO K HEÉ ÄOÁABËßETCß MAT-   *C
C*        PÈÖA EI, ÃÄE I - EÄÈHÈ×HAß, A E - MOÄÓËÜ MÈHÈMAËÜHOÃO COÁ-  *C
C*        CTBEHHOÃO ×ÈCËA ÈCXOÄHOÉ MATPÈÖÛ.                           *C
C*                                                                    *C
C* OÃPAHÈ×EHÈE: HE ÁOËEE 3100 ÏEPEMEHHÛX                              *C
C*                                                                    *C
C* ÏPÈME×AHÈE : ECËÈ INDXD=1, TO ÏOËÜÇOBATEËÜ CAM C×ÈTAET             *C
C*        MATPÈÖÓ ÃECCE                                               *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      SUBROUTINE GRAD(X,FX,GR,NM,STP,DGR,MULT)
      COMMON /INDXDR/ INDXD
     */grcom/ig,ing,fst,ngrad,mlt,ibinar(100) 
      DIMENSION BINAR(100)
      REAL*4 X(1),FX,DFX,STP2,STEP,FXR,FXL,FXMAX,BET1,BET2,STPMIN,GRI,
     *GR(1),STP(1),ACC,MULT(1),FXLR,MLT,DVD,GRJ,STPI,
     *GRMAX,GESII,
     *GESJJ,EIGN,DGR(1),FXMAX1
      INTEGER NM(2),NGRAD,MAXGR/15/,FST
      integer*4 ibinar,intg
      EQUIVALENCE (INTG,RL),(BINAR,IBINAR)
      JK(J,K)=J*(J-1)/2+K
      M=NM(2)
      M2=M/2
      FST=FST+1
      IF((INDXD.EQ.1).AND.(M2-M2/2*2.EQ.0)) RETURN
      IF(INDXD.EQ.1) M=M-M/4*4
      N=NM(1)
      NDER=0
      IF(MLT.GT.1.E0) GO TO 1
      ig=1
      AN=1.E0/N
      acc=epsilon(acc)
      MLT=3.E0**AN
      DVD=1.E0/MLT
      NN12=N*(N+1)/2
      NNS=NN12+N+1
1     NGRAD=NGRAD+1
c
      NWN=0
      NNOW=0
      NWD=1
C***           PAC×ET ÃPAÄÈEHTA
C
      DO 101 I=1,N
C*             HAÄO ËÈ C×ÈTATÜ ÄAHHÓÞ KOMÏOHEHTÓ?
      NUMBWD=I/31
      IDELT=I-NUMBWD*31
      IF(IDELT.NE.0) NUMBWD=NUMBWD+1
      IF(IBINAR(NUMBWD).EQ.0) GO TO 11
      IF(IDELT.EQ.0) IDELT=31
      INDX=IBINAR(NUMBWD)/2**(IDELT-1)
      INDX=INDX-INDX/2*2
      IF(INDX.EQ.1) GO TO 101
C
11    STPSAV=STP(I)
      NNOW=NNOW+1
12    STEP=STP(I)
      STPI=STP(I)
      IF(M.GT.3.AND.INDXD.NE.1) STEP=STEP*MULT(I)
      DGRI=DGR(I)
      STPMIN=AMAX1(ABS(X(I))*ACC*10.,1.E-30)
      IF(STEP.LT.STPMIN) STEP=STPMIN
      XSV=X(I)
      X(I)=X(I)+STEP
      CALL FUN(X,FXR)
      IF(FX.EQ.FXR.AND.FST.EQ.1.AND.STPSAV*2.E3.GT.STP(I)) STP(I)=
     *STP(I)*10.
      IF(FX.EQ.FXR.AND.FST.EQ.1.AND.STPSAV*2.E3.GT.STP(I)) X(I)=XSV
      IF(FX.EQ.FXR.AND.FST.EQ.1.AND.STPSAV*200.0.GT.STP(I)) GO TO 12
      IF((M.LT.4).AND.((NGRAD.NE.0).AND.(NGRAD.LT.MAXGR)))GO TO 20
C
C*             PAC×ET  GESSE(I,I)
C
      STP2=STEP+STEP
      X(I)=X(I)-STP2
      NGRAD=0
      CALL FUN(X,FXL)
      DGRI=((FXR-FX)+(FXL-FX))/STP2
      X(I)=X(I)+STP2
      IF(ABS(DGRI).GT.1.E-30) GO TO 19
      IF(DGRI.LT.0) DGRI=-1.E-30
      IF(DGRI.GE.0.) DGRI=1.E-30
19    IF(M.EQ.4.OR.M.EQ.5) GO TO 55
      GO TO 53
20    CONTINUE
C
      DFX=ABS(FXR-FX)
      FXMAX=AMAX1(ABS(FX),ABS(FXR))
      IF(FXMAX.LT.1.E-30) FXMAX=1.E-30
      BET1=AMIN1(DFX/(FXMAX*ACC),STEP/(ACC*AMAX1(ABS(X(I)),ABS(X(I)-
     *STEP))))
      IF(BET1.GT.1.E+7) BET1=1.E+7
      GRI=(FXR-FX)/STEP
C
      BET2=ABS(GRI/DGRI)
      IF(BET1-BET2) 30,30,40
C              BET2>=BET1  . ËÈHEÉHAß AÏÏPOKCÈMAÖÈß ÓÄOBËETBOPÈTEËÜHA
30    IF((BET2.LT.BET1*3.).AND.(BET1.LT.1.E+5)) GO TO 50
      STP(I)=STEP*2.
      IF(BET1.GT.1.E+5) STP(I)=STEP/3.
      GO TO 50
C
C         BET1>BET2  ËÈHEÉHAß AÏÏPOKCÈMAÖÈß HE ÓÄOBËETBOPÈTEËÜHA
40    IF((BET1.LT.BET2*3.).AND.(BET2.LT.1.E+5)) GO TO 50
      STP(I)=STP(I)/3.
      IF(BET2.GT.1.E+5) STP(I)=STEP*4.
50    CONTINUE
C
      IF(STP(I).LT.STPMIN) STP(I)=STPMIN
      IF((STP(I).EQ.STPMIN).OR.(AMAX1(BET1,BET2).LT.20.))GO TO 52
      IF(((BET2.GT.30.E0).OR.(BET2.GT.BET1/3.)).AND.(NGRAD.NE.0))GOTO55
C
C*     HAXOÄÈMCß BÁËÈÇÈ TO×KÈ MÈHÈMÓMA ÈËÈ TO×KÈ ÏEPEÃÈÁA.
C*    KOMÏOHEHTÛ ÃPAÄÈEHTA C×ÈTAEM ÏO CXEME BTOPOÃO ÏOPßÄKA TO×HOCTÈ
C
52    IF(NGRAD.EQ.0) GO TO 53
      STP2=STEP+STEP
      X(I)=X(I)-STP2
      NDER=NDER+1
      AII=BET1
      AII=BET2
      CALL FUN(X,FXL)
      DGRI=((FXR-FX)+(FXL-FX))/STP2
      X(I)=X(I)+STP2
      IF(ABS(DGRI).GT.1.E-30) GO TO 53
      DGRI=1.E-30
      IF(DGRI.LT.0.E0) DGRI=-1.E-30
53    GRI=(FXR-FXL)/STP2
      FXMAX1=AMAX1(ABS(FXL),ABS(FXR))
      BET1=AMIN1(ABS(FXR-FXL)/(FXMAX1*ACC),STP2/(ACC*AMAX1(ABS(X(I))
     *,ABS(X(I)-STP2))))
      IF(BET1.LT.10.E0) STP(I)=STEP*1.5
      IF(BET1.GT.200.E0) STP(I)=STEP/2.
C
55    IF(GRI.EQ.0..AND.(M.EQ.2.OR.M.EQ.3.OR.M.GT.5)) STP(I)=STP(I)*3.
      MULT(I)=MULT(I)/(STP(I)/STPI)
      IF(M.EQ.2.OR.M.EQ.3.OR.M.GT.5) GR(I)=GRI
      IF(M.GT.3) GR(NN12+I)=FXR
      IF(M.LE.3) DGR(I)=DGRI
      IF(M.GT.3) DGR(I)=STEP
      IF((M.LE.3).OR.(INDXD.EQ.1))GO TO 100
C
C*        TPEÁÓETCß PAC×ÈTÛBATÜ MATPÈÖÓ GESSE
C
      GESII=((FXR-FX)+(FXL-FX))/(STEP*STEP)
      JKI=JK(I,I)
      GR(N+JKI)=GESII
      IF(I.EQ.1) GO TO 95
      I1=I-1
      DO 90 J=1,I1
      NUMBWD=J/31
      NWN=NWN+1
      IDELT=J-NUMBWD*31
      IF(IDELT.NE.0) NUMBWD=NUMBWD+1
      IF(IBINAR(NUMBWD).EQ.0) GO TO  57
      IF(IDELT.EQ.0) IDELT=31
      INDX=IBINAR(NUMBWD)/2**(IDELT-1)
      INDX=INDX-INDX/2*2
      IF(INDX.EQ.1) GO TO  90
57    X(J)=X(J)+DGR(J)
      CALL FUN(X,FXLR)
      X(J)=X(J)-DGR(J)
      K=NN12+J
      FXRR=GR(K)
      DFX=(FXLR-FXR)+(FX-FXRR)
      GR(N+NWN)=DFX/(DGR(J)*STEP)
      GR(NNS+NWD)=GR(N+NWN)
      NWD=NWD+1
      FXMAX=AMAX1(ABS(FX),ABS(FXRR),ABS(FXL),ABS(FXLR),1.E-30)
      BET1=ABS(DFX)/(FXMAX*ACC)
      GRJ=GR(J)
      IF(M.NE.2.AND.M.NE.3.AND.M.LE.5) GRJ=(FXRR -FX)/DGR(J)
      IF(BET1.GT.20.E0) GO TO 60
C      TO×HOCTÜ PAC×ETA GESSE ÏO F>1%
      IF(((DFX.LE.0.).AND.(GRI.GT.0.)).OR.((DFX.GT.0.).AND.(GRI.LE.
     *0.)))MULT(I)=MULT(I)*MLT
      IF(((DFX.LE.0.).AND.(GRJ.GT.0.)).OR.((DFX.GT.0.).AND.(GRJ.LE.0)))
     *MULT(J)=MULT(J)*MLT
      GO TO 70
C     TO×HOCTÜ PAC×ETA GESSE(I,J)<1%
60    IF(BET1.LT.200.E0) GO TO 70
      MULT(I)=MULT(I)*DVD
      MULT(J)=MULT(J)*DVD
70    CONTINUE
C         ÏPOBEPKA TO×HOCTÈ GESSE(I,J) ÏO ÃPAÄÈEHTÓ
      GRMAX=AMAX1(ABS(GRI),ABS(GRJ),1.E-30)
      BET2=ABS(DFX)/(GRMAX*ACC)
      GESJJ=GR(N+JK(J,J))
      IF(BET2.GT.20.E0) GO TO 80
C         TO×HOCTÜ ÏO ÃPAÄÈEHTÓ >1%
      IF(((DFX.LE.0.).AND.(GESII.LE.0.)).OR.((DFX.GT.0.).AND.(GESII.GT.
     *0.))) MULT(I)=MULT(I)*MLT
      IF(((DFX.LE.0.).AND.(GESJJ.LE.0.)).OR.((DFX.GT.0.).AND.(GESJJ.GT.
     *0.))) MULT(J)=MULT(J)*MLT
      GO TO 85
80    IF(BET2.LE.200.E0) GO TO 85
C         TO×HOCTÜ PAC×ETA GESSE(I,J) < 0.1%
      MULT(I)=MULT(I)*DVD
      MULT(J)=MULT(J)*DVD
85    CONTINUE
      IF(MULT(J).GT.1000.) MULT(J)=1000.
      IF(MULT(J).LT.0.001) MULT(J)=0.001
90    CONTINUE
C         ÏOÈCK KOÝÔÔÈÖÈEHTA  MULT ÄËß GESSE(I,I)
95    FXMAX1=AMAX1(ABS(FXL),FXMAX)
      BET1=ABS(GESII)*STEP*STEP/(FXMAX1*ACC)
      IF(BET1.GT.20.E0) GO TO 96
      MULT(I)=MULT(I)*MLT
      GO TO 97
96    IF(BET1.LT.200.E0) GO TO 97
      MULT(I)=MULT(I)*DVD
97    BET2=ABS(GESII*STEP*STEP/(GRI*ACC))
      IF(BET2.GT.20.E0) GO TO 98
      MULT(I)=MULT(I)*MLT
      GO TO 99
98    IF(BET2.LT.200.E0) GO TO 99
      MULT(I)=MULT(I)*DVD
99    CONTINUE
      GR(NNS+NWD)=GESII
      NWD=NWD+1
      IF(MULT(I).GT.1000.E0) MULT(I)=1000.E0
      IF(MULT(I).LT.0.001) MULT(I)=0.001
100   X(I)=XSV
      NWN=NWN+1
101   CONTINUE
      IF(NDER.EQ.NNOW) NGRAD=0
C     CALL  GRD(X,FX,G1,NGR)
C     DO 222 I=1,N
C22   G1(I)=ABS((G1(I)-GR(I))/G1(I))*100.
C     WRITE(6,666) (G1(I),I=1,N)
C66   FORMAT(3X,'DGR=',9(E11.3,'%'),/)
      IF(M.LT.4) RETURN
C
C               BOCCTAHOBËEHÈE MACCÈBA DGR
C
      DO 110 I=1,N
110   DGR(I)=ABS(GR(N+JK(I,I))*DGR(I))
      IF(INDXD.EQ.1) RETURN
C               ÏOÄÃOTOBKA PAÁO×EÃO MACCÈBA
C
C*       BOCCTAHOBËEHÈE ÏOËOÆÈTEËÜHOÉ OÏPEÄEËEHHOCTÈ : GESSE(I,I)=
C*       GESSE(I,I)-EIGN ,ÃÄE EIGN - MÈHÈMAËÜHOE OTPÈÖATEËÜHOE
C*       COÁCTBEHHOE ×ÈCËO GESSE
C
C     CALL  GESS(X,FX,G1,NM,STP,DGR,MLT)
      CALL EIGEN(GR(NNS+1),R,NNOW,1)
C     DO 1000 I=1,NN12
C000  G2(I)=ABS((GR(N+I)-G1(I))/G1(I))*100.
C     WRITE(6,603) (G2(I),I=1,NN12),(G1(I),I=1,3)
C03   FORMAT(3X,'GESSE=',9(E11.3,'%'),/)
C     CALL  EIGEN(G1,R,N,1)
C     DO 1100 I=1,N
C100  G2(I)=ABS((GR(NNS+JK(I,I))-G1(JK(I,I)))/G1(JK(I,I)))*100.
C     WRITE(6,604) (G2(I),I=1,N) ,G1(1),G1(3)
C04   FORMAT(3X,'EIGEN=',9(E11.3,'%'),/)
      EIGN=GR(NNOW*(NNOW+1)/2+NNS)
      IF(EIGN.GE.0.E0) RETURN
      DO 130 I=1,N
130   GR(N+JK(I,I))=GR(N+JK(I,I))-EIGN
C
      RETURN
      END
      BLOCK DATA
      COMMON /INDXDR/ IN
     */grcom/ ig,ing,fst,nggr,mlt,ibinar(100)
      real*4 mlt
      integer fst
      integer*4 ibinar
      DATA IN/1/
      data ig/0/,ing/0/,fst/1/,nggr/15/,mlt/1.0/,ibinar/100*0/
      END
