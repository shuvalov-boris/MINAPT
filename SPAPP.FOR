C**********************************************************************C
C*                                                                    *C
C*                                                                    *C
C*        CÏËAÉH OÏÈCAHÈE      ÔPAÃMEHTA KOHTÓPA (ÏAPAMETPÈ×ECKOE     *C
C*         ÏPEÄCTABËEHÈE) C PAC×ETOM TO×HOCTÈ AÏÏPOKCÈMAÖÈÈ           *C
C*                                                                    *C
C* ÏAPAMETPÛ: XY(N,2) - KOOPÄÈHATÛ TO×EK KOHTÓPA (BXOÄ), ÏOÄÃOTOBËEH- *C
C*       HÛE ÏPOÃPAMMAMÈ TPREP È TRECLC.                              *C
C*            T(I) - ÏAPAMETPÈ×ECKAß KOOPÄÈHATA TO×EK KOHTÓPA         *C
C*            C(NKNOT,3,2) - ÏAPAMETPÛ CÏËAÉHA B  ÔOPMATE             *C
C*       ÏÏÏ IMSL. ÏOCËEÄHÈÉ ÈHÄEKC - HOMEP KOOPÄÈHATÛ:1-X,2-Y(BÛXOÄ) *C
C*            IC=NKNOT-1                                              *C
C*            XK(NKNOT),YK(NKNOT),TK(NKNOT) - X,Y È T KOOPÄÈHATÛ      *C
C*       ÓÇËOB CÏËAÉHA (BXOÄ)                                         *C
C*            EPS -  ÏOËÓ×ÈBØAßCß TO×HOCTÜ AÏÏPOKCÈMAÖÈÈ (ÑÐÅÄÍÅÊÂ.   *Ñ
C*		ÎÒÊËÎÍÅÍÈÅ)					      *C
C*            IFLG=AB  -  ÓÏPABËßÞÙÈÉ ÏAPAMETP :                      *C
C*              B=0 - ÏPOBOÄÈTÜ AÏÏPOKCÈMAÖÈÞ KOHTÓPA                 *C
C*              B=1 - PACC×ÈTATÜ TOËÜKO CÏËAÉH KOÝÔÔÈÖÈEHTÛ È EPS     *C
C*              A=0 - AÏÏPOKCÈMAÖÈß C OÏTÈMAËÜHOÉ PACCTAHOBKOÉ ÓÇËOB  *C
C*              A=1 - ÁEÇ                                             *C
C*            WK((NKNOT+6)*N*2+3*NKNOT) - PAÁO×ÈÉ MACCÈB              *C
C*            IRC - KOÄ BOÇBPATA :                                    *C
C*                0 - HOPMAËÜHOE ÇABEPØEHÈE                           *C
C*                1 - OØÈÁKA B ÇHA×EHÈßX N ÈËÈ NKNOT                  *C
C*                2 - BXOÄHÛE MACCÈBÛ HEBEPHO ÏOÄÃOTOBËEHÛ            *C
C*                3 - HE ÄOCTÈÃHÓTA HÓÆHAß TO×HOCTÜ OÏTÈMÈÇAÖÈÈ ÏPÈ   *C
C*      OÏTÈMAËÜHOÉ PACCTAHOBKE ÓÇËOB. HA BÛXOÄE - HAÈËÓ×ØAß KOHÔÈ-   *C
C*      ÃÓPAÖÈß ÓÇËOB                                                 *C
C*                                                                    *C
C* AËÃOPÈTM: B ÇABÈCÈMOCTÈ OT ÇAKAÇA ÏPOBOÄÈTCß CÏËAÉH AÏÏPOKCÈMAÖÈß  *C
C*     ÈËÈ ÏPOCTO PAC×ET CÏËAÉH KOÝÔÔÈÖÈEHTOB ÏO MACCÈBAM XK,YK È TK  *C
C*     C PAC×ETOM KPÈTEPÈß AÏÏPOKCÈMAÖÈÈ ÔPAÃMEHTA KOHTÓPA.           *C
C*     ECËÈ KOHTÓP ÇAMKHÓT, TO ÏOCËE ÏPOBEÄEHÈß AÏÏPOKCÈMAÖÈÈ OH      *C
C*     "ÇAMÛKAETCß" ÏPOÃPAMMOÉ ICSPLN, ×TO ÏPÈBOÄÈT K HEKOTOPOMÓ      *C
C*     ÓXÓÄØEHÈÞ AÏÏPOKCÈMAÖÈÈ (ECËÈ, KOHE×HO, OHA ÇAKAÇAHA). EPS     *C
C*     ÏPÈ ÝTOM ÏEPEC×ÈTÛBAETCß ÏOBTOPHO.                             *C
C*           KOHTÓP C×ÈTAETCß ÇAMKHÓTÛM, ECËÈ XK(1)=XK(NKNOT) È       *C
C*     YK(1)=YK(NKNOT)                                                *C
C*                                                                    *C
C* ÏPÈME×AHÈE: ÏPÈ ÏPOBEÄEHÈÈ AÏÏPOKCÈMAÖÈÈ ÄOËÆHO BÛÏOËHßTÜCß        *C
C*             OÃPAHÈ×EHÈE: NKNOT HE ÁOËEE 28 - OÃPAHÈ×EHÈE ÏPOÃPAMM  *C
C*             ICSFKU È ICSFKV È TOËÜKO.                              *C
C*                                                                    *C
C*                                                                    *C
C**********************************************************************C
      SUBROUTINE SPAPP(XY,T,N,XK,YK,TK,NKNOT,C,IC,EPS,WK,IFLG,IRC)
      common/nksave/ nks
      REAL XY(N,2),T(N),XK(1),YK(1),TK(1),WK(N,1),C(IC,3,2)
	DATA NKMAX/28/
	write(*, *) 'SPAPP: N =', N, ', NKNOT =', NKNOT
C*          ÓCTAHOBKA ÔËAÆKA: ÇAMKHÓT KOHTÓP ÈËÈ HET?
      ICZ=0
      IF(XK(1).EQ.XK(NKNOT).AND.YK(1).EQ.YK(NKNOT)) ICZ=1
C*          ÏPOBEPKA KOPPEKTHOCTÈ BXOÄHÛX ÏAPAMETPOB
      IRC=0
      IF(IFLG-IFLG/2*2.EQ.1.AND.N.GE.4.AND.NKNOT.GT.1) GO TO 185
      IRC=1
      IF(N.LT.4.OR.NKNOT.LT.2.OR.NKNOT.GT.NKMAX) GO TO 300
      IRC=0
      nks=nknot
      IF(IFLG/2-IFLG/4*2.EQ.1) GO TO 150
C
C*       OÏTÈMAËÜHÛÉ ÏOÄÁOP ÏOËOÆEHÈÉ ÓÇËOB
C
	CALL  ICSVKU(T,XY,N,TK,NKNOT,WK,C,ic,EPS,WK(NKNOT+NKNOT+1,1),
     *IER)
      IF(IER.EQ.33) IRC=3
      IF(IER.GT.33) IRC=2
      IF(IRC.EQ.2) GO TO 300
      GO TO 180
C
C*       AÏÏPOKCÈMAÖÈß ÏPÈ ÔÈKCÈPOBAHHÛX ÏOËOÆEHÈßX YÇËOB
C
150   CALL  ICSFKU(T,XY,N,0,TK,NKNOT,WK,C,ic,EPS,WK(NKNOT+NKNOT+1
     *,1),IER)
      IF(IER.NE.0) IRC=2
      IF(IRC.EQ.2) GO TO 300
C*       ÏEPEKA×ÈBAEM KOOPÄÈHATÛ ÓÇËOB
180   DO 183 I=1,ic
      XK(I)=WK(I,1)
183   YK(I)=WK(NKNOT+I,1)

C*       PACC×ÈTÛBAEM KOOPÄÈHATÛ ÏOCËEÄHEÃO ÓÇËA
      CALL  ICSEVU(TK,XK,NKNOT,C,ic,TK(NKNOT),XK(NKNOT),1,IER)
      CALL  ICSEVU(TK,YK,NKNOT,C(1,1,2),ic,TK(NKNOT),YK(NKNOT),1,IER)
      IF(ICZ.EQ.0) GO TO 187
C*       KOHTÓP ÇAMKHÓT. ÄEËAEM ÇAMÛKAHÈE EÃO AÏÏPOKCÈMAÖÈÈ
      XK(NKNOT)=(XK(1)+XK(NKNOT))/2.
      YK(NKNOT)=(YK(1)+YK(NKNOT))/2.
      XK(1)=XK(NKNOT)
      YK(1)=YK(NKNOT)
C*        PAC×ET CÏËAÉH KOÝÔÔÈÖÈEHTOB ÏO MACCÈBAM XK,YK È TK
185   IF(ICZ.EQ.1) CALL  ICSPLN(TK,XK,NKNOT,C,ic,WK,IER)
      IF(ICZ.EQ.0) CALL  ICSCCU(TK,XK,NKNOT,C,ic,IER)
	IF(IER.NE.0) IRC=2
      IF(IRC.EQ.2) GO TO 300
      IF(ICZ.EQ.1) CALL  ICSPLN(TK,YK,NKNOT,C(1,1,2),ic,WK,IER)
	IF(ICZ.EQ.0) CALL  ICSCCU(TK,YK,NKNOT,C(1,1,2),ic,IER)
	IF(IER.NE.0) IRC=2
      IF(IRC.EQ.2) GO TO 300
C*         PAC×ET TO×HOCTÈ AÏÏPOKCÈMAÖÈÈ
187   EPS=0.
      CALL  ICSEVU(TK,XK,NKNOT,C,ic,T,WK,N,IER)
      CALL  ICSEVU(TK,YK,NKNOT,C(1,1,2),ic,T,WK(1,2),N,IER)
c      dt=t(2)-t(1)
      DO 190 I=1,N
c      if(i.ne.n.and.i.ne.1) dt=t(i+1)-t(i-1)
c      if(i.eq.n) dt=t(n)-t(n-1)
190    EPS=EPS+((XY(I,1)-WK(I,1))**2+(XY(I,2)-WK(I,2))**2)
      eps=sqrt(eps/(n-1)) 
c190   eps=sqrt(eps/(t(n)-t(1))/2.)
300   RETURN
      END