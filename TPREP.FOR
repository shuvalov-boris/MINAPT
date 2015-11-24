C**********************************************************************C
C*                                                                    *C
C*                                                                    *C
C*        ÓCTAHOBKA COOTBETCTBÈß ÓÇËOB È ÔPAÃMEHTA KOHTÓPA ÏO         *C
C*     T-KOOPÄÈHATE ÄËß ÏAPAMETPÈ×ECKOÃO CÏËAÉH-ÏPEÄCTABËEHÈß KOHTÓPA *C
C*                                                                    *C
C* ÏAPAMETPÛ: XY(N,2) - KOOPÄÈHATÛ TO×EK KOHTÓPA (BXOÄ). HA BÛXOÄE -  *C
C*       CO CÃËAÆEHHÛMÈ "BÛÁPOCAMÈ", ECËÈ CÃËAÆÈBAHÈE ÇAKAÇÛBAËOCÜ.   *C
C*            T(I) - PACCTOßHÈE OT TO×KÈ I ÄO ÏEPBOÉ TO×KÈ KOHTÓPA    *C
C*       ÏO ÄÓÃE KOHTÓPA, AÏÏPOKCÈMÈPÓEMOÃO ËOMAHOÉ, ÏPOXOÄßÙEÉ       *C
C*       ×EPEÇ TO×KÈ (XY(I,1),XY(I,2)) (BÛXOÄ)                        *C
C*            XK(NKNOT),YK(NKNOT)           - X,Y     KOOPÄÈHATÛ      *C
C*       ÓÇËOB CÏËAÉHA (BXOÄ).                                        *C
C*            TK(I) - T-KOOPÄÈHATA ÓÇËA I,      PACC×ÈTÛBAETCß        *C
C*       ÏOCPEÄCTBOM ÏPOEKTÈPOBAHÈß ÓÇËOB HA ÔPAÃMEHT KOHTÓPA.        *C
C*            NBEG,NEND - HOMEPA TO×EK, ÓKAÇÛBAÞÙÈE HA HA×AËO È KOHEÖ *C
C*       MAKCÈMAËÜHOÃO OTPEÇKA ÔPAÃMEHTA,KOTOPÛÉ ËEÆÈT BHÓTPÈ OTPEÇKA,*C
C*       BKËÞ×AÞÙEÃO ÏPOEKÖÈÈ ÓÇËOB (ÏO T-KOOPÄÈHATE) (BÛXOÄ)         *C
C*            NUSL,NUSR - ×ÈCËO ÓÇËOB, KOTOPÛE ÏOCËE ÈX ÓÏOPßÄO×ÈBAHÈß*C
C*       ÏO T-KOOPÄÈHATE, ÏOËÓ×ÈBØEÉCß ÏPÈ ÏPOEKTÈPOBAHÈÈ ÓÇËOB HA    *C
C*       KOHTÓP, HE ÈMEÞT TO×EK ÔPAÃMEHTA, ËEÆAÙÈX MEÆÄÓ TAKÈM ÓÇËOM  *C
C*       È EÃO COCEÄOM, È HAXOÄßTCß "ËEBEE" ÈËÈ "ÏPABEE" ÔPAÃMEHTA    *C
C*            EPS - KPÈTEPÈÉ ÄËß PGM CÃËAÆÈBAHÈß BÛÁPOCOB.            *C
C*            NDOUBL - ×ÈCËO ÓÇËOB ÈMEÞÙÈX ÈÄEHTÈ×HÓÞ T-KOOPÄÈHATÓ    *C
C*       ÏPOEKÖÈÈ ÈX HA ÔPAÃMEHT KOHTÓPA. B CËÓ×AE, KOÃÄA ÏEPBÛÉ È    *C
C*       ÏOCËEÄHÈÉ ÓÇËÛ COBÏAÄAÞT C HA×AËÜHOÉ TO×KOÉ ÇAMKHÓTOÃO       *C
C*       KOHTÓPA (PABHOÉ KOHE×HOÉ), TO T-KOOPÄÈHATÛ ÝTÈX ÓÇËOB PAÇËÈ×-*C
C*       HÛ È NDOUBL HE ÔÈKCÈPÓET ÝTO ÄÓÁËÈPOBAHÈE ÓÇËOB              *C
C*            IFLG=AB  -  ÓÏPABËßÞÙÈÉ ÏAPAMETP :                      *C
C*              B=1 - CÃËAÆÈBATÜ KOOPÄÈHATÛ TO×EK ÔPAÃMEHTA HEOÁXOÄÈMO*C
C*              B=0 - TAKOÉ HEOÁXOÄÈMOCTÈ HET                         *C
C*              A=1 - XK È YK HEOÁXOÄÈMO ÓÏOPßÄO×ÈTÜ B COOTBETCTBÈÈ C *C
C*       ÏOËÓ×ÈBØEÉCß T-KOOPÄÈHATOÉ    ÏPOEKÖÈÉ ÝTOÃO ÓÇËA            *C
C*              A=0 - ÓÏOPßÄO×ÈBATÜ HE HAÄO                           *C
C*              WK(MAX(N*3,2*NKNOT))-PAÁO×ÈÉ MACCÈB.HA BÛXOÄE -       *C
C*       KBAÄPATÛ MÈHÈMAËÜHÛX PACCTOßHÈÉ OT                           *C
C*               ÓÇËOB ÄO ÔPAÃMEHTA KOHTÓPA                           *C
C*            IWK(I)     - ÖEËO×ÈCËEHHÛÉ MACCÈB. HA BÛXOÄE - HOMEPÀ   *C
C*       ÓÇËÎÂ, ÓÏÎÐßÄÎ×ÅÍÍÛÅ ÏÎ ÂÎÇÐÀÑÒÀÍÈÞ TK. PAÇMEPHOCTÜ - NKNOT  *C
C*       ECËÈ IWK(I)<0, TO ÝTO ÓÇEË (IABS(IWK(I)))                    *C
C*        BHE ÔPAÃMEHTA KOHTÓPA (CM OÏÈCAHÈE                          *C
C*       ÏAPAMETPOB NUSL,NUSR)                                        *C
C*            IRC - KOÄ BOÇBPATA :                                    *C
C*                0 - HOPMAËÜHOE ÇABEPØEHÈE                           *C
C*                1 - ÄAHHÛE HE AHAËÈÇÈPOBAËÈCÜ HA HAËÈ×ÈE BÛÁPOCOB   *C
C*      ÏPÈ ÇAÄAHÈÈ TAKOÃO AHAËÈÇA, T. K. ×ÈCËO BÛÁPACOB > 25%        *C
C*                2 - ÏOPßÄOK CËEÄOBAHÈß ÏPOEKÖÈÉ ÓÇËOB(TK) ÖÈKËÈ×ECKÈ*C
C*      HE COBÏAÄAET C ÏOPßÄKOM CËEÄOBAHÈß ÓÇËOB                      *C
C*                3 - ÝKBÈBAËEHTHO IRC=1&IRC=2                        *C
C*                4 - MAËO N (<7 ÏPÈ B=1 ÈËÈ <2 ÏPÈ B=0)              *C
C*                5 - NKNOT < 2                                       *C
C*                6 - HEBEPHO ÏOÄÃOTOBËEHÛ BXOÄHÛE MACCÈBÛ ÄËß        *C
C*      ÏPOÃPAMMÛ CÃËAÆÈBAHÈß BÛÁPOCOB ICSMOU                         *C
C*                7 - ÓÇËÛ ÖÈÊËÈ×ÅÑÊÈ ÓÏÎÐßÄÎ×ÅÍÛ, ÍÎ Ó TK ÍÅÒ        *C
C*      ÌÎÍÎÒÎÍÍÎÃÎ ÐÎÑÒÀ                                             *C
C*                                                                    *C
C* AËÃOPÈTM: 1.CÃËAÆÈBAÞTCß      BÛÁPOCÛ (ECËÈ ÝTO ÇAKAÇAHO)          *C
C*           2.ECËÈ ×ÈCËO BÛÁPOCOB > 25% - CÃËAÆÈBAHÈE HE ÏPOÈÇBOÄÈTCß*C
C*           3.KOOPÄÈHATA TO×KÈ I CÃËAÆÈBAETCß, ECËÈ OHA OTKËOHßETCß  *C
C*     OT CÏËAÉH-KPÈBOÉ >EPS*(T(I+3)-T(I-3))/6.                       *C
C*           4.ÏPOEKTÈPÓÅM ÓÇËÛ HA ÔPAÃMEHT KOHTÓPA, PACC×ÈTÛBAß È    *C
C*     ÓÏOPßÄO×ÈBAß TK,IWK,WK È XK,YK (ECËÈ TPEÁÓETCß). C×ÈTAEM NBEG, *C
C*     NEND,NUSL,NUSR È NDOUBL.                                       *C
C*                                                                    *C
C* ÏPÈME×AHÈE: ECËÈ ÏPÈ ÏPOEKTÈPOBAHÈÈ ÓÏOPßÄO×EHHÛX PAHEE ÓÇËOB HET  *C
C*             ÓBEPEHHOCTÈ B ÖÈKËÈ×ECKOM COBÏAÄEHÈÈ ÏOPßÄKA CËEÄOBAHÈß*C
C*             ÈX ÏPOEKÖÈÉ, TO PEKOMEHÄÓETCß BBOÄÈTÜ B ÔPAÃMEHT       *C
C*             KOHTÓPA KOOPÄÈHATÛ ÓÇËOB.                              *C
C*             TO×KÈ ÔPAÃMEHTA HE ÄOËÆHÛ COBÏAÄATÜ, ÇA ÈCKËÞ×EHÈEM    *C
C*             KPAÉHÈX TO×EK ÏPÈ ÇAMKHÓTOM KOHTÓPE                    *C
C*             ÏPEÄÏOËAÃAETCß ÏPÈ CÃËAÆÈBAHÈÈ, ×TO ÏEPBÛE È ÏOCËEÄHÈE *C
C*             3 TO×KÈ ÔPAÃMEHTA HAXOÄßTCß HA CBOEM "ÈCTÈHHOM" MECTE  *C
C*             (OÃPAHÈE×EHÈE ÏPOÃPAMMÛ ICSMOU)                        *C
C*             CMOTPÈTE TAKÆE KOMMEHTAPÈÈ ÏPOÃAMMÛ CONKN              *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      SUBROUTINE TPREP(XY,T,N,XK,YK,TK,NKNOT,NBEG,NEND,NUSL,NUSR,NDOUBL,
     *EPS,WK,IWK,IFLG,IRC)
      REAL XY(N,2),T(N),XK(1),YK(1),TK(1),WK(N,3)
      INTEGER*2 IWK(NKNOT)
C		Âûâîä ïàðàìåòðîâ äëÿ ïðîàåðêè
      write(*,*) 'TPREP'
      open(18, FILE='tprep.log', action='write')
	write(18, *) '-------TPREP.LOG---------- '
	write(18, *) 'N =', N, ', NKNOT =', NKNOT
	write(18, *) 'XY T:'
	write(18, 1) (XY(I, 1), XY(I, 2), T(I), I = 1, N)
    1 format(3F12.4)
      close(18)
C
C*          ÏPOBEPKA KOPPEKTHOCTÈ BXOÄHÛX ÏAPAMETPOB
      IRC=4
      NMAX=2
      IF(IFLG-IFLG/2*2.EQ.1) NMAX=7
      IF(N.LT.NMAX) GO TO 300
      IRC=5
      IF(NKNOT.LT.2) GO TO 300
C*          ÇÀÌÊÍÓÒ ËÈ ÊÎÍÒÓÐ (Â ÑÌÛÑËÅ ÏÎ ÓÇËÀÌ) ?
      ICN=0
      IF(XK(1).EQ.XK(NKNOT).AND.YK(1).EQ.YK(NKNOT)) ICN=1
C*          PAC×ET PACCTOßHÈÉ ÄO TO×EK ËOMAHOÉ BÄOËÜ KOHTYPA
      T(1)=0.
      DO 10 I=2,N
      T(I)=T(I-1)+SQRT((XY(I,1)-XY(I-1,1))**2+(XY(I,2)-XY(I-1,2))**2)
10    CONTINUE
      IRC=0
      IF(IFLG-IFLG/2*2.EQ.0) GO TO 100
C
C*        AHAËÈÇ ÄAHHÛX HA BÛÁPOCÛ
C
      DO 20 I=1,N
      WK(I,1)=XY(I,1)
      WK(I,2)=XY(I,2)
20    WK(I,3)=T(I)
C*          CÃËAÆÈBAEM  TO×KÈ ÏO KOOPÄÈHATE X
      CALL  ICSMOU(WK(1,3),WK,N,1.,EPS,N/4+1,WK(1,4),IER)
      IF(IER.NE.0) IRC=6
      IF(IRC.EQ.6) GO TO 300
C*          CÃËAÆÈBAEM  TO×KÈ ÏO KOOPÄÈHATE Y
      CALL  ICSMOU(WK(1,3),WK(1,2),N,1.,EPS,N/4+1,WK(1,4),IER)
      IF(IER.NE.0) IRC=6
      IF(IRC.EQ.6) GO TO 300
C*            NVILD - ×ÈCËO BÛÁPOCOB
      NVILD=0
      DO 30 I=1,N
      IF(WK(I,1).EQ.XY(I,1).AND.WK(I,2).EQ.XY(I,2)) GO TO 30
      NVILD=NVILD+1
30    CONTINUE
C
      IF(NVILD*4.GT.N) IRC=1
      IF(NVILD*4.GT.N.OR.NVILD.EQ.0) GO TO 100
C*          ÏEPEC×ET T
      DO 60 I=1,N
      DO 60 J=1,2
60    XY(I,J)=WK(I,J)
      T(1)=0.
      DO 70 I=2,N
70    T(I)=T(I-1)+SQRT((XY(I,1)-XY(I-1,1))**2+(XY(I,2)-XY(I-1,2))**2)
C
C*       PAC×ET TK,NBEG,NEND, È TÄ.,  YÏOPßÄO×ÈBAHÈE YÇËOB B ÏOPßÄKE
C*   BOÇPACTAHÈß TK
C
100   CALL  CONKN(XY,XY(1,2),T,N,XK,YK,TK,NKNOT,NBEG,NEND,WK,IWK,NUSL,
     *NUSR,NDOUBL,IFLG/2)
      IF(ICN.EQ.1) NDOUBL=NDOUBL-1
      IF(NDOUBL.LT.0) NDOUBL=0
C*          ÏPOBEPKA HA ÖÈKËÈ×ECKÓÞ ÓÏOPßÄO×EHHOCTÜ ÏPOEKÖÈÉ ÓÇËOB
      ICU=-1
      I1=NKNOT
      NKBEG=0
      DO 110 I=1,NKNOT
      IF(IABS(IWK(I)).EQ.1) NKBEG=I-1
      IF(IABS(IWK(I)).NE.IABS(IWK(I1)+1).AND.IABS(IWK(I)).NE.1)
     * GO TO 120
      IF(IABS(IWK(I)).EQ.1.AND.IABS(IWK(I1)).NE.NKNOT) GO TO 120
      I1=I
110   CONTINUE
      ICU=0
      GO TO 130
C*          ÏÐÎÂÅÐÊÀ ÍÀ ÖÈÊËÈ×ÅÑÊÈ ÎÁÐÀÒÍÛÉ ÏÎÐßÄÎÊ ÑËÅÄÎÂÀÍÈß ÓÇËÎÂ
120   I1=NKNOT
      ICU=-1
      NKBEG=0
      DO 125 I=1,NKNOT
      IF(IABS(IWK(I)).EQ.1) NKBEG=I-1
      IF(IABS(IWK(I)).NE.IABS(IWK(I1)-1).AND.IABS(IWK(I)).NE.NKNOT)
     * GO TO 127
      IF(IABS(IWK(I)).EQ.NKNOT.AND.IABS(IWK(I1)).NE.1) GO TO 127
      I1=I
125   CONTINUE
      ICU=1
      GO TO 130
127   IF(ICN.EQ.0) GO TO 130
C*           ÑÏËÀÉÍ ÀÏÏÐÎÊÑÈÌÀÖÈÈ ÇÀÌÊÍÓÒ. ÌÅÍßÅÌ ÏÎÐßÄÎÊ ÑËÅÄÎÂÀÍÈß
C*    ÏÅÐÂÎÃÎ È ÏÎÑËÅÄÍÅÃÎ ÓÇËÎÂ È ÂÍÎÂÜ ÏÐÎÂÅÐßÅÌ ÍÀ ÖÈÊËÈ×ÅÑÊÓÞ 
C*                     ÓÏÎÐßÄÎ×ÅÍÍÎÑÒÜ
      DO I=1,NKNOT
        IF(IABS(IWK(I)).EQ.1) NK1=I
        IF(IABS(IWK(I)).EQ.NKNOT) NKLAST=I
      END DO
      IWKS=IWK(NK1)
      IWK(NK1)=IWK(NKLAST)
      IWK(NKLAST)=IWKS
C*       ÏÐÎÂÅÐßÅÌ: ÓÇËÛ ÖÈÊËÈ×ÅÑÊÈ ÓÏÎÐßÄÎ×ÅÍÛ?
      ICU=-1
      I1=NKNOT
      NKBEG=0
      DO I=1,NKNOT
        IF(IABS(IWK(I)).EQ.1) NKBEG=I-1
        IF(IABS(IWK(I)).NE.IABS(IWK(I1)+1).AND.IABS(IWK(I)).NE.1)
     *   GO TO 128
        IF(IABS(IWK(I)).EQ.1.AND.IABS(IWK(I1)).NE.NKNOT) GO TO 128
      I1=I
      END DO
      ICU=0
      GO TO 130
C*          ÏÐÎÂÅÐÊÀ ÍÀ ÖÈÊËÈ×ÅÑÊÈ ÎÁÐÀÒÍÛÉ ÏÎÐßÄÎÊ ÑËÅÄÎÂÀÍÈß ÓÇËÎÂ
128   I1=NKNOT
      ICU=-1
      NKBEG=0
      DO I=1,NKNOT
         IF(IABS(IWK(I)).EQ.1) NKBEG=I-1
         IF(IABS(IWK(I)).NE.IABS(IWK(I1)-1).AND.IABS(IWK(I)).NE.NKNOT)
     *   GO TO 130
         IF(IABS(IWK(I)).EQ.NKNOT.AND.IABS(IWK(I1)).NE.1) GO TO 130
         I1=I
      END DO
      ICU=1
C
130   IF(IFLG/2.EQ.1.AND.ICU.GE.0) GO TO 300
      IF(IFLG/2.EQ.1.AND.ICU.LT.0) GO TO 290
      IF(IABS(IWK(1)).EQ.1.AND.ICU.EQ.0) GO TO 300
C*          ÁÛË ÇÀÊÀÇ - ÓÇËÛ ÍÅ ÓÏÎÐßÄÎ×ÈÂÀÒÜ. ÏÐÈÂÎÄÈÌ Â ÑÎÎÒÂÅÒÑÒÂÈÈ
C*    Ñ ÏÎÐßÄÊÎÌ ÑËÅÄÎÂÀÍÈß ÓÇËÎÂ ÌÀÑÑÈÂÛ TK,IWK,WK
C      
C*          ÑÎÐÒÈÐÎÂÊÀ ÐÀÑÑÒÎßÍÈÉ
      DO I=1,NKNOT
        NU=IABS(IWK(I))
        WK(NU+NKNOT,1)=WK(I,1)
      END DO
C
      DO I=1,NKNOT
        WK(I,1)=WK(I+NKNOT,1)
      END DO
C*          ÓÏÎÐßÄÎ×ÈÂÅÌ ÌÀÑÑÈÂ TK
      TKMAX=0.
      DO I=1,NKNOT
        NU=IABS(IWK(I))
          IF(IWK(I).GT.0.AND.TKMAX.LT.TK(I)) TKMAX=TK(I)
        WK(NU+NKNOT,1)=TK(I)
      END DO
C
      DO I=1,NKNOT
        TK(I)=WK(I+NKNOT,1)
      END DO
C*           
      IF(ICU.LT.0) GO TO 290
      IF(ICU.EQ.0) GO TO 140
C*          ÓÇËÛ "ÀÍÒÈÓÏÎÐßÄÎ×ÅÍÛ". ÏÐÅÎÁÐÀÇÓÅÌ ÌÀÑÑÈÂ ÒÊ
      DO I=1,NKNOT
        IF(IWK(I).GT.0) TK(I)=TKMAX-TK(I)
      END DO
      IRC=0
      IF(IABS(IWK(1)).EQ.1.AND.ICU.EQ.0) GO TO 300
C
140      IRC=7
      IF(XY(1,1).NE.XY(N,1).OR.XY(1,2).NE.XY(N,2)) GO TO 300
      IRC=0
C*          ÔÐÀÃÌÅÍÒ ÊÎÍÒÓÐÀ ÇÀÌÊÍÓÒ. ÏÅÐÅÑ×ÈÒÛÂÀÅÌ TK
      IF(ICU.EQ.1) GO TO  150
      DO I=1,NKBEG
        TK(NKNOT-I+1)=TK(NKNOT-I+1)+T(N)
      END DO
      GO TO 160
150   NKA=NKNOT-NKBEG-1
      IF(NKA.EQ.0) GO TO 160
C*          ÏÅÐÅÑ×ÅÒ ÒÊ ÄËß ÀÍÒÈÓÏÎÐßÄÎ×ÅÍÍÛÕ ÓÇËÎÂ
      DO I=NKBEG+2,NKNOT
        TK(I)=T(N)+TK(I)
      END DO
C*          TK(1) -> 0
160   DO I=2,NKNOT
        TK(I)=TK(I)-TK(1)
      END  DO
      TK(1)=0.
      GO TO 300
C*      HET ÖÈKËÈ×ECKOÉ ÓÏOPßÄO×EHHOCTÈ
290   IRC=IRC+2
300   RETURN
      END
