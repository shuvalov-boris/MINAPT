      SUBROUTINE DGKM(X,Y,FI,LAMBDA)
C***************************************************************
C*
C*
C* ÏPOÃPAMMA ÏEPEBOÄA KOOPÄÈHAT ÈÇ ÃPAÄÓCOB B KÈËOMETPÛ
C*       HA×AËO OTC×ETA - ÃPÈHBÈ× HA ÝKBATOPE
C*
C* ÏAPAMTEPÛ: X,Y - ÄOËÃOTA È ØÈPOTA TO×KÈ B KM
C*            FI,LAMBDA - TO ÆE B ÃPAÄÓCAX
C*
C* ÏPÈME×AHÈE. TO×HOCTÜ PAC×ETA ÄOCTATO×HA ÄËß ÄAHHÛX, CHßTÛ CO
C*             CÏÓTHÈKOBÛX CHÈMKOB BÛCOKOÃO PAÇPEØEHÈß (1 KM):
C*                                        ÏO MEPÈÄ. ÏO ØÈPOTE
C* TO×HOCTÜ MOÄEËÈ ÃEOÈÄA: ÏOËÞC -        -290M       0M
C*      HA 1 ÃPAÄÓC        45 ÃPAÄ ØÈPOTÛ +80M       23*PI*FI/360M
C* (ÃEOÈÄ - ÝËËÈÏCOÈÄ)     ÝKBATOP -      +110M       0M
C*                   L+ERR=L(ÈCTÈHHAß)
C*
C***************************************************************
      REAL LAMBDA
      DATA R0/6378.245/,PI/3.1415926/,ONEDG/111.24/
      REARTH=R0*(1.-0.00335*SIN(LAMBDA*PI/180.)**2)
      RLITL=REARTH*COS(LAMBDA*PI/180.)
      X=PI*RLITL*FI/180.
      IF(X.LT.0.) X=2.*PI*RLITL+X
      Y=ONEDG*LAMBDA
      RETURN
      END