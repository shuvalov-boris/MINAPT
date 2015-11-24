
C**********************************************************************C
C*                                                                    *C
C*                                                                    *C
C*        ÏEPEC×ET T-KOOPÄÈHAT ÔPAÃMEHTA KOHTÓPA B COOTBETCTBÈÈ CO    *C
C*  ÇHA×EHÈßMÈ T-KOOPÄÈHAT ÏPOEKÖÈÈ ÓÇËOB HA KOHTÓP È ÈX ÈCTÈHHÛX     *C
C*  BEËÈ×ÈH T                                                         *C
C*                                                                    *C
C* ÏAPAMETPÛ: X(N),Y(N)-KOOPÄÈHATÛ TO×EK KOHTÓPA (BXOÄ). HA BÛXOÄE -  *C
C*       ÏEPECOPTÈPOBAHHÛE B COOTBETCTBÈÈ C HOBÛMÈ T-ÇHA×EHÈßMÈ       *C
C*            T(I) - PACCTOßHÈE OT TO×KÈ I ÄO ÏEPBOÉ TO×KÈ KOHTÓPA    *C
C*       ÏO ÄÓÃE KOHTÓPA, AÏÏPOKCÈMÈPÓEMOÃO ËOMAHOÉ, ÏPOXOÄßÙEÉ       *C
C*       ×EPEÇ TO×KÈ (XY(I,1),XY(I,2)) (BXOÄ).BÛXOÄ - ÏEPEC×ÈTAHHÛE   *C
C*       BEËÈ×ÈHÛ                                                     *C
C*            N - MOÆET ÓMEHÜØÈTÜCß HA 1, ECËÈ KOHTÓP ÇAMKHÓT È T-    *C
C*       ÏPOEKÖÈß ÏEPBOÃO È ÏOCËEÄHEÃO ÓÇËA (ÓÇËÛ COBÏAÄAÞT) HE       *C
C*       ÏPÈXOÄÈTCß HA TO×KÓ HA×AËA KOHTÓPA (COBÏAÄAÞÙÓÞ C KOHE×HOÉ   *C
C*       TO×KOÉ ÇAMKHÓTOÃO KOHTÓPA
C*            TK(NKNOT) - ÓÏOPßÄO×EHÛE ÏO BOÇPACTAHÈÞ T KOOPÄÈHATÛ*C
C*       ÓÇËOB CÏËAÉHA - ÏPOEKÖÈÈ ÓÇËOB HA ÔPAÃMEHT KOHTÓPA (BXOÄ)    *C
C*            NKNOT - ×ÈCËO ÓÇËOB.                                    *C
C*            TKOLD(I) - ÄECTBÈTEËÜHAß T-KOOPÄÈHATA ÓÇËA I (ÓÏOPßÄO×.)*C
C*            KFIRST - HOMEP ÓÇËA, ÈMEÞÙEÃO ÏPOEKÖÈÞ HA KOHTÓP TK(1)  *C
C*            WK(N) - PAÁO×ÈÉ MACCÈB                                  *C
C*                                                                    *C
C* AËÃOPÈTM: BEËÈ×ÈHE T(I), ÏOÏABØEÉ B ÈHTEPBAË (TK(J+1)-TK(J)),      *C
C*     ÏPÈCBAÈBAETCß ÇHA×EHÈE TKOLD(KFIRST+J-1)+(T(I)-TK(I))/         *C
C*     /(TK(J+1)-TK(J))*(TKOLD(KFIRST+J)-TKOLD(KFIST+J-1)).           *C
C*     BEËÈ×ÈHÛ X,Y COPTÈPÓÞTCß B ÏOPßÄKE ÓBEËÈ×EHÈß ÏOËÓ×ÈBØEÉCß     *C
C*     T-KOOPÄÈHATÛ. ECËÈ KOHTÓP ÇAMKHÓT, TO N MOÆET ÁÛTÜ ÏPÈ         *C
C*     OÏPEÄEËEHHÛX ÓCËOBÈßX ÓMEHÜØEHA HA 1. ÏPEÄÏOËAÃAETCß, ×TO:     *C
C*   1.ÏOPßÄKÈ CËEÄOBAHÈß ÓÇËOB ÄËß TK È TKOLD OTËÈ×AÞTCß ËÈØÜ        *C
C*     ÖÈKËÈ×ECKÈM CÄBÈÃOM.                                           *C
C*   2.B MACCÈBAX T,TK È TKOLD HET ÄBÓX COBÏAÄAÞÙÈX ÇHA×EHÈÉ,ECËÈ     *C
C*     KOHTÓP HE ÇAMKHÓT. ECËÈ ÇAMKHÓT, TO MOÃÓT COBÏAÄATÜ TOËÜKO     *C
C*     ÏEPBAß È ÏOCËEÄHßß KOOPÄÈHATÛ Ó MACCÈBOB T È TK.               *C
C*                                                                    *C
C* ÏPÈME×AHÈß:  T-KOOPÄÈHATÛ ÓÇËA I, ÏOËÓ×EHÛE  ÏPOÃPAMMOÉ TPREP È    *C
C*     ÈMEÞÙEÃO  ÇHA×EHÈE IWK(I)<0 HE ÄOËÆHÛ ÏPÈCÓTCTBOBATÜ B TK È    *C
C*     TKOLD. TO×KÈ KOHTÓPA - TO×KÈ B ÄÈAÏOÇOHE ÈHÄEKCOB NBEG È NEND  *C
C*     HÈKAKÈX ÏPOBEPOK HA KOPPEKTHOCTÜ ÄAHHÛX HE ÄEËAÞ. TO×EK        *C
C*     ÔPAÃMEHTA MEÆÄÓ ÓÇËAMÈ NKNOT È 1 HE ÄOËÆHO ÁÛTÜ              *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      SUBROUTINE TRECLC(X,Y,T,N,TK,TKOLD,NKNOT,KFIRST,WK)
      REAL T(N),X(N),Y(N),TK(1),TKOLD(1),WK(N)
C*      ÈÙEM ÓÇEË KB, ÁËÈÆAÉØÈÉ K T(1): TK(KB)<=T(1)
      NKNOT1=NKNOT-1
      DO 10 KB=1,NKNOT1
      IF(TK(KB+1).GT.T(1).AND.TK(KB).LE.T(1)) GO TO 15
10    CONTINUE
      KB=NKNOT
15    CONTINUE
      KB1=KB+1
      IF(KB1.GT.NKNOT) KB1=1
C*        BÛÁPACÛBAHÈE ÏOCËEÄHEÉ TO×KÈ KOHTÓPA, ECËÈ ÝTO HEOÁXOÄÈMO
      N1=N
      IF(X(1).NE.X(N).OR.Y(1).NE.Y(N)) GO TO 20
      IF(TK(KB).EQ.T(1)) GO TO 20
      N1=N-1
20    CONTINUE
C*       ÏEPEC×ET T-KOOPÄÈHATÛ ÔPAÃMEHTA KOHTÓPA
      NBT=0
      T1=T(1)
      DTK=TK(KB1)-TK(KB)
      IF(TK(KB1).LT.TK(KB)) DTK=T(N)-TK(KB)+TK(KB1)-T1
      TKS=TK(KB)
      IF(TK(KB1).LT.TK(KB)) TKS=TK(KB)-T(N)+T1
      KBO=KB+KFIRST-1
      IF(KBO.GT.NKNOT) KBO=KBO-NKNOT
      KBO1=KBO+1
      DTKO=TKOLD(KBO1)-TKOLD(KBO)
      DO 50 I=1,N1
      IF(KBO.EQ.1.AND.NBT.EQ.0) NBT=I
      T(I)=TKOLD(KBO)+(T(I)-TKS)/DTK*DTKO
      IF(I.EQ.N1) GO TO 50
      IF(T(I+1).LE.TKS+DTK) GO TO 50
C*       ÏEPECKO×ÈËÈ HA CËEÄÓÞÙÈÉ OTPEÇOK KOHTÓPA
      KB=KB1
      KBO=KBO1
      KB1=KB+1
      IF(KB1.GT.NKNOT) KB1=1
      DTK=TK(KB1)-TK(KB)
      IF(DTK.NE.0.) GO TO 55
C*      ÏEPECKO×ÈËÈ C OÄHOÃO KPAß ÇAMKHÓTOÃO KOHTÓPA HA ÄPÓÃOÉ
      KB=KB1
      KB1=KB+1
      IF(KB1.GT.NKNOT) KB1=1
      DTK=TK(KB1)-TK(KB)
55    IF(TK(KB1).LT.TK(KB)) DTK=T(N)-TK(KB)+TK(KB1)-T1
      TKS=TK(KB)
C     IF(TK(KB1).LT.TK(KB)) TKS=TK(KB)-T(N)+T1
      IF(KBO.EQ.NKNOT) KBO=1
      KBO1=KBO+1
      DTKO=TKOLD(KBO1)-TKOLD(KBO)
      IF(KBO.EQ.1) NBT=0
50    CONTINUE
C
C*      ÓÏOPßÄO×ÈBAHÈE TO×EK KOHTÓPA COÃËACHO HOBÛM ÇHA×EHÈßM T
C
C*        ÓÏOPßÄO×ÈBAEM X-KOOPÄÈHATÓ
      N=N1
      IF(NBT.EQ.1) GO TO 100
      NP=NBT
      DO 60 I=1,N
      WK(I)=X(NP)
      NP=NP+1
      IF(NP.GT.N) NP=1
60    CONTINUE
      DO 65 I=1,N
65    X(I)=WK(I)
C*      ÓÏOPßÄO×ÈBAEM Y-KOOPÄÈHATÓ
      NP=NBT
      DO 70 I=1,N
      WK(I)=Y(NP)
      NP=NP+1
      IF(NP.GT.N) NP=1
70    CONTINUE
      DO 75 I=1,N
75    Y(I)=WK(I)
C*      ÓÏOPßÄO×ÈBAEM T-KOOPÄÈHATÓ
      NP=NBT
      DO 80 I=1,N
      WK(I)=T(NP)
      NP=NP+1
      IF(NP.GT.N) NP=1
80    CONTINUE
      DO 85 I=1,N
85    T(I)=WK(I)
C
100   RETURN
      END