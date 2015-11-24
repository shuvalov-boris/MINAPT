C**********************************************************************C
C*                                                                    *C
C*        ÓÏOPßÄO×ÈBAHÈE ÓÇËOB È COÏÓTCTBÓÞÙÈX ÈM XAPAKTEPÈCTÈK       *C
C*  B ÏOPßÄKE PACÏOËOÆEHÈß ÓÇËOB BÄOËÜ ÔPAÃMEHTA KOHTÓPA,             *C
C*                 ÏPEOÁPAÇOBAHÈE KOOPÄÈHAT                           *C
C*  ÓÇËOB B ÏAPAMETPÈ×ECKOE ÏPEÄCTABËEHÈE (PAC×ET ÈX ÏOËOÆEHÈÉ ÏO     *C
C*                KOOPÄÈHATE T - MACCÈB TK)                           *C
C*                                                                    *C
C*                                                                    *C
C* ÏAPAMETPÛ: X(N),Y(N),T(N) - X,Y È T - KOOPÄÈHATÛ TO×EK KOHTÓPA     *C
C*            XK(NKNOT),YK(NKNOT),TK(NKNOT) - TOÆE ÄËß ÓÇËOB          *C
C*            NBEG,NEND - HA×AËO È KOHEÖ ÔPAÃMEHTA KOHTÓPA, KOTOPÛÉ   *C
C*  ÏO T-KOOPÄÈHATE ËEÆÈT B ÄÈAÏOÇOHE, ÏOKPÛBAÞÙEMCß ÓÇËAMÈ           *C
C*            DIST(I)-MÈHÈMAËÜHOE PACCTOßHÈE ÓÇËA IABS(NPOS(I)) ÄO    *C
C*  KOHTÓPA (BÛXOÄ)                                                   *C
C*            NPOS(I) - HOMEP ÓÇËA, ÏOÏABØEÃO HA I-OE MECTO ÏOCËE     *C
C*  ÓÏOPßÄO×ÈBAHÈß ÏO T-KOOPÄÈHATE(MACCÈB TK) (BÛXOÄ). ECËÈ NPOS(I)<0,*C
C*  TO ÄAHHÛÉ ÓÇEË ËEÆÈT BHE ÏPEÄEËOB KOHTÓPA È HE ßBËßETCß ÁËÈÆAÉØÈM *C
C*  K OÄHOMÓ ÈÇ KOHÖOB KOHTÓPA.                                       *C
C*            NUSL,NUSR - ×ÈCËO ÓÇËOB C OTPÈÖATEËÜHÛMÈ BEËÈ×ÈHAMÈ     *C
C*  NPOS CËEBA È CÏPABA OT KPAEB ÔPAÃMEHTA COOTBETCTBEHHO             *C
C*            NDOUBL - ×ÈCËO ÏAP ÓÇËOB, ËEÆAÙÈX B ÏPEÄEËAX ÔPAÃMEHTA È*C
C*  ÈMEÞÙÈX OÄÈHAKOBÓÞ T-KOOPÄÈHATÓ ÏPOEKÖÈÈ. ECËÈ KK ÓÇËOB ÈMEÞT     *C
C*  OÄÈHAKOBÓÞ ÏPOEKÖÈÞ TK, TO K NKNOT ÄOÁABËßETCß HE 1, A KK-1       *C
C*            IFLG: 0 - ÓÏOPßÄO×ÈBATÜ MACCÈBÛ NPOS,DIST,A KOOPÄÈHATÛ  *C
C*  ÓÇËOB XK,YK HE ÓÏOPßÄO×ÈBATÜ. 1 - ÓÏOPßÄO×ÈBATÜ BCE ÝTÈ MACCÈBÛ   *C
C*                                                                    *C
C* AËÃOPÈTM:  PACC×ÈTÛBAÞTCß TK(I) - PACCTOßHÈß BÄOËÜ ËOMAHOÉ,        *C
C*  ÏPOXOÄßÙEÉ ×EPEÇ TO×KÈ (X(J),Y(J)),OT TO×KÈ (X(1),Y(1)) ÄO TO×KÈ, *C
C*  HAXOÄßÙEÉCß HA MÈHÈMAËÜHOM PACCTOßHÈÈ OT ÓÇËA (XK(I),YK(I))       *C
C*            ÈÇ ÓÇËOB BÛXOÄßÙÈX ÇA KAKOÉ-ËÈÁO KPAÉ ÔPAÃMEHTA KOHTÓPA *C
C*  BÛÁÈPAETCß ÓÇEË C MÈHÈMAËÜHÛM PACCTOßHÈEM ÄO KPAÉHEÉ TO×KÈ (ÃPAHÈ×*C
C*  HÛÉ ÓÇEË). BEËÈ×ÈHA TK ÄËß HEÃO ÁÓÄET T(1)-(PACCTOßHÈE OT ÓÇËA ÄO  *
C*  TO×KÈ (X(1),Y(1)), ECËÈ ÝTOT KPAÉ - HA×AËO È T(N)+ (PACCTOßHÈE ÄO *C
C*  TO×KÈ (X(N),Y(N)), ECËÈ KOHEÖ. ÄPÓÃÈE ÓÇËÛ, ËEÆAÙÈE BHE KOHTÓPA,  *C
C*  BÛHOCßTÜCß BHA×AËO MACCÈBA C COOTBETCTBÓÞÙÈMÈ XAPAKTEPÈCTÈKAMÈ,   *C
C*  ËÈÁO B KOHEÖ B TOÉ ÏOCËEÄOBATEËÜHOCTÈ, B KOTOPOÉ OHÈ BCTPE×AÞTCß  *C
C*  ÈÇHA×AËÜHO. ECËÈ TAKÈX ÓÇËOB HE OKAÇAËOCÜ ÄËß KAKOÃO-ËÈÁO KPAß,TO *C
C*  PACC×ÈTÛBAÞTCß BEËÈ×ÈHÛ NBEG È/ÈËÈ NEND.                          *C
C*                                                                    *C
C* ÏPÈME×AHÈß: ×TOÁÛ ÈÇÁEÆATÜ ÏOËÓ×EHÈß ËOÆHÛX ÏPOEKÖÈÉ ÄËß ÇAMKHÓTOÃO*C
C*  KOHTÓPA, HO HEÇAMKHÓTOÃO ÔPAÃMEHTA, CËEÄÓET BBOÄÈTÜ TOËÜKO ÓÇËÛ,  *C
C*  COOTBETCTBÓÞÙÈE ÝTOMÓ ÔPAÃMEHTÓ.                                  *C
C*             ÏPÈ PAC×ETE TK(I) È ÓÏOPßÄO×ÈBAHÈÈ ÓÇËOB ÏPOÈCXOÄÈT    *C
C*  ÏOËHÛÉ ÏEPEÁOP ÁEÇ BCßKÈX ÓXÈÙPEHÈÉ                               *C
C*             CËEÄÓET ÏOMHÈTÜ, ×TO B OÁÙEM CËÓ×AE ÏPABÈËÜHOCTÜ       *C
C*  ÏOCTPOEHÈß ÏPOEKÖÈÉ ÓÇËOB HE ÃAPAHTÈPÓETCß. TPEÁÓETCß ÄÈAËOÃOBAß  *C
C*  ÏPOÖEÄÓPA KOHTPOËß                                                *C
C*                                                                    *C
C*     ÏPOÃPAMMA HE ÇABÈCÈT OT CÏEÖÈÔÈKÈ TPAHCËßTOPA È ÝBM            *C
C*                                                                    *C
C* BÛÇÛBAEMÛE ÏOÄÏPOÃPAMMÛ: INCLT - ÏPÈCÓTCTBÓET B ÝTOM ÆE ÔAÉËE      *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      SUBROUTINE CONKN(X,Y,T,N,XK,YK,TK,NKNOT,NBEG,NEND,DIST,NPOS,NUSL,
     *NUSR,NDOUBL,IFLG)
      REAL X(N),Y(N),T(N),XK(1),YK(1),TK(1),XY1(2),XY2(2),DIST(NKNOT)
      INTEGER*2   NPOS(NKNOT)
C     TKMIN=T(N)
C     TKMAX=T(1)
      ICFC=0
      IF(X(1).EQ.X(N).AND.Y(1).EQ.Y(N)) ICFC=1
      DO 100 I=1,NKNOT
      SS=(XK(I)-X(1))**2+(YK(I)-Y(1))**2
      TS=0.
      DO 80 J=2,N
C*        OÏPEÄEËEHÈE MÈHÈMAËÜHOÃO PACCTOßHÈß ÄO OTPEÇKA
C*             (X(I-1),Y(I-1)) - (X(I),Y(I))
       XY1(1)=XK(I)-X(J-1)
       XY1(2)=YK(I)-Y(J-1)
       XY2(1)=X(J)-X(J-1)
       XY2(2)=Y(J)-Y(J-1)
       AL=XY1(1)*XY2(1)+XY1(2)*XY2(2)
       IF(J.EQ.2) AL1=AL
       IF(AL.LE.0.) GO TO 80
C*       TO×KA (XK(I),YK(I)) ËEÆÈT "ÏPABEE" ËEBOÃO KOHÖA OTPEÇKA
       INDX=2
       IF(AL.GT.XY2(1)**2+XY2(2)**2) INDX=1
       GO TO (52,54), INDX
C*        MÈHÈMAËÜHOE PACCTOßHÈE ÄËß OTPEÇKA - ÄO TO×KÈ (X(J),Y(J))
52     S2=(XK(I)-X(J))**2+(YK(I)-Y(J))**2
       IF(S2.GT.SS) GO TO 80
       IF(S2.EQ.SS.AND.I.EQ.1) GO TO 80
C*       ÏPOÈÇOØËO ÓMEHÜØEHÈE TEKÓÙEÃO MÈHÈMAËÜHOÃO PACCTOßHÈß
       SS=S2
       TS=T(J)
       AL1=1.
       GO TO 80
C*        HEOÁXOÄÈMO PACC×ÈTÛBATÜ TO×KÓ - ÏPOEKÖÈÞ ÓÇËA I HA
C*          TEKÓÙÈÉ OTPEÇOK
54    DX=XY2(1)
      DY=XY2(2)
      IND=3
      IF(DY.EQ.0.) IND=2
      IF(DX.EQ.0.) IND=1
      GO TO (55,60,65), IND
C*          DX=0.
55    XC=X(J-1)
      YC=-YK(I)
      GO TO 70
C*          DY=0.
60    XC=XK(I)
      YC=Y(J-1)
      GO TO 70
C*        C×ÈTAEM ÏO OÁÙEÉ ÔOPMÓËE
65    DIV=1./(DX*DX)+1./(DY*DY)
      XC=((XK(I)/DY+YK(I)/DX)/DY-(Y(J-1)/DY-X(J-1)/DX)/DX)/DIV
      YC=((Y(J-1)/DY-X(J-1)/DX)/DY+(XK(I)/DY+YK(I)/DX)/DX)/DIV
C*       KOHTPOËÜ HA MEHÜØEE PACCTOßHÈE
70    S2=(XK(I)-XC)**2+(YK(I)-YC)**2
      IF(S2.GT.SS) GO TO 80
      IF(S2.EQ.SS.AND.I.EQ.1) GO TO 80
      AL1=1.
      SS=S2
      TS=T(J-1)+SQRT((XC-X(J-1))**2+(YC-Y(J-1))**2)
80    CONTINUE
      TK(I)=TS
      NPOS(I)=I
      DIST(I)=SS
C     IF(TKMIN.GT.TS) TKMIN=TS
C     IF(TKMAX.LT.TS) TKMAX=TS
C*       ÏOME×AEM ÓÇËÛ, KOTOPÛE ËEÆAT BHE KOHTÓPA
      IF((TS.EQ.T(N).AND.INDX.EQ.1).OR.(TS.EQ.0..AND.AL1.LT.0.))
     *NPOS(I)=-NPOS(I)
100   CONTINUE
C
C*     ÓÏOPßÄO×ÈBAHÈE ÓÇËOB
C
      NK=NKNOT
      NKNOT=0
      NDOUBL=0
      DO 110 I=1,NK
C*       ÓÏOPßÄO×ÈBAEM  ÓÇËÛ
      NKNOT=NKNOT+1
      TIN=TK(I)
      DIN=DIST(I)
      NIN=NPOS(I)
      XIN=XK(I)
      YIN=YK(I)
      CALL  INCLT(TK,XK,YK,TIN,XIN,YIN,NKNOT,DIST,NPOS,DIN,NIN,NDOUBL,
     *IFLG,NK)
110   CONTINUE
C
C*        PAC×ET NBEG È NEND
C
      NBEG=1
      NEND=N
      NUSL=0
      NUSR=0
      IF(ICFC.EQ.1) GO TO 200
C*         PAC×ET NBEG
      DO 130 I=1,N
      IF(T(I).GE.TK(1)) GO TO 135
130   CONTINUE
135   NBEG=I
140   IF(TK(NKNOT).EQ.T(N)) GO TO 160
C*         PAC×ET NEND
      DO 150 I=1,N
      IF(T(N-I+1).LE.TK(NKNOT)) GO TO 155
150   CONTINUE
155   NEND=N-I+1
C
C*      PAC×ET NUSL,NUSR È TK ÄËß ÓÇËOB BHE KOHTÓPA, HO ÁËÈÆAÉØÈX
C*          K OÄHOMÓ ÈÇ EÃO KPAEB
C
160   IF(NPOS(1).GT.0) GO TO 180
C*      ËEBEE "ËEBOÃO" KOHÖA KOHTÓPA ECTÜ ÓÇËÛ
C*        ÈÙEM ÁËÈÆAÉØÈÉ
      NKMIN=1
      DO 170 I=2,NKNOT
      IF(TK(I).GE.0..AND.NPOS(I).GT.0) GO TO 175
      IF(DIST(I).GT.DIST(NKMIN)) GO TO 170
      NKMIN=I
170   CONTINUE
175   NUS=I-1
      NUSL=NUS-1
      IF(TK(I).EQ.0) NUSL=NUSL+1
      IF(TK(I).EQ.0) GO TO 180
C*       ÏEPEC×ET TK ÄËß ÁËÈÆAÉØEÃO K HA×AËÓ ÓÇËA È ÏEPECTAHOBKA
      NPS=IABS(NPOS(NKMIN))
      TKS=-DIST(NKMIN)
      DISTS=DIST(NKMIN)
      NPOS(NKMIN)=NPOS(NUS)
      NPOS(NUS)=NPS
      TK(NKMIN)=TK(NUS)
      TK(NUS)=TKS
      DIST(NKMIN)=DIST(NUS)
      DIST(NUS)=DISTS
      IF(IFLG.EQ.0) GO TO 180
C*      ÓÏOPßÄO×ÈBAHÈE È KOOPÄÈHAT ÓÇËOB
      XKS=XK(NKMIN)
      YKS=YK(NKMIN)
      XK(NKMIN)=XK(NUS)
      XK(NUS)=XKS
      YK(NKMIN)=YK(NUS)
      YK(NUS)=YKS
C
180   IF(NPOS(NKNOT).GT.0) GO TO 200
C*      ÏPABEE "ÏPABOÃO" KOHÖA KOHTÓPA ECTÜ ÓÇËÛ
C*        ÈÙEM ÁËÈÆAÉØÈÉ
      NKMIN=NKNOT
      NK1=NKNOT-1
      DO 190 I=1,NK1
      IF(TK(NKNOT-I).LE.T(N).AND.NPOS(NKNOT-I).GT.0) GO TO 195
      IF(DIST(NKNOT-I).GT.DIST(NKMIN)) GO TO 190
      NKMIN=NKNOT-I
190   CONTINUE
195   NUS=NKNOT-I+1
      NUSR=I-1
      IF(TK(NKNOT-I).EQ.T(N)) NUSR=NUSR+1
      IF(TK(NKNOT-I).EQ.T(N)) GO TO 200
C*       ÏEPEC×ET TK ÄËß ÁËÈÆAÉØEÃO K KOHÖÓ ÓÇËA È ÏEPECTAHOBKA
      NPS=IABS(NPOS(NKMIN))
      TKS=T(N)+DIST(NKMIN)
      DISTS=DIST(NKMIN)
      NPOS(NKMIN)=NPOS(NUS)
      NPOS(NUS)=NPS
      TK(NKMIN)=TK(NUS)
      TK(NUS)=TKS
      DIST(NKMIN)=DIST(NUS)
      DIST(NUS)=DISTS
      IF(IFLG.EQ.0) GO TO 200
C*      ÓÏOPßÄO×ÈBAHÈE È KOOPÄÈHAT ÓÇËOB
      XKS=XK(NKMIN)
      YKS=YK(NKMIN)
      XK(NKMIN)=XK(NUS)
      XK(NUS)=XKS
      YK(NKMIN)=YK(NUS)
      YK(NUS)=YKS
200   RETURN
      END
      SUBROUTINE INCLT(TK,XK,YK,TIN,XIN,YIN,NKNOT,DIST,NPOS,DIN,NIN,
     *NDOUBL,IFLG,NKLAST)
C*********************************************************************
C*
C*    B MACCÈBÛ TK,XK,YK,DIST,NPOS PAÇMEPHOCTÜÞ NKNOT-1 ÄOÁABËßETCß ÏO
C*  OÄHOMÓ ÝËEMEHTÓ. COÁËÞÄAETCß MOHOTOHHOCTÜ ÓBEËÈ×EHÈß TK C POCTOM I
C*
C*
C*********************************************************************
      REAL XK(1),YK(1),TK(1),DIST(1)
      INTEGER*2 NPOS(1)
      I1=NKNOT-1
      NK=NKNOT
      NDS=NDOUBL
      DO  4 J=1,NK
      IF(NKNOT.EQ.1.OR.J.EQ.NKNOT) GO TO 5
C*      ÏOÄC×ET ×ÈCËA OÄÈHAKOBÛ ÏPOEKÖÈÉ ÓÇËOB
      IF(NIN.GT.0.AND.NPOS(J).GT.0.AND.TK(J).EQ.TIN.AND.NDS.EQ.NDOUBL)
     *NDOUBL=NDOUBL+1
      IF(TIN.GT.TK(J)) GO  TO 4
      IF(TIN.EQ.TK(J).AND.(IABS(NPOS(J)).NE.1.OR.IABS(NIN).NE.NKLAST))
     *GO TO 4
C*       BBOÄÈM ÝËEMEHTÛ B MACCÈBÛ
      DO 1 K=J,I1
      TK(I1-K+J+1)=TK(I1-K+J)
      DIST(I1-K+J+1)=DIST(I1-K+J)
      NPOS(I1-K+J+1)=NPOS(I1-K+J)
      IF(IFLG.EQ.1) XK(I1-K+J+1)=XK(I1-K+J)
      IF(IFLG.EQ.1) YK(I1-K+J+1)=YK(I1-K+J)
1     CONTINUE
      GO TO 5
4     CONTINUE
5     TK(J)=TIN
      DIST(J)=DIN
      NPOS(J)=NIN
      IF(IFLG.EQ.1) XK(J)=XIN
      IF(IFLG.EQ.1) YK(J)=YIN
      RETURN
      END
