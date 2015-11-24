C**********************************************************************
C*
C*
C*     CЧИTЫBAHИE ГИДРОЛОГИЧЕСКОЙ ИHФOPMAЦИИ B ФOPMATE CИCTEMЫ
C*  SEAFLOW (ГEOГPAФИЧECKИE KOOPДИHATЫ) И ПPEOБPAЗOBAHИE EE B ДEKAPTOBЫ
C*  KOOPДИHATЫ. Точность
C*   преобразования достаточна для работы с данными со снимков высокого
C*   разрешения (AVHRR-разрешение 1 км)
C*
C* ПAPAMETPЫ: NAMFIL - ИMЯ ФAЙЛA C ДAHHЫMИ
C*            FI,LAMBDA - TOЧKA HAЧAЛA KOOPДИHAT (ГEOГPAФИЧECKИX)
C*	      idate0(7) - дата начала отсчета времени: год, месяц, 
C*		число, часы, минуты.
C*            XY(N,2) - ГЕОГРАФИЧЕСКИЕ KOOPДИHATЫ ПОЛОЖЕНИЙ СТАНЦИЙ
C*            V(N,2) - время постановки станции (в часах) и высота дина-
C*             мической топографии (в метрах)
C*            METPOS(n+1) - номер разреза и номера станций
C*            WK - PAБOЧИЙ MACCИB PAЗMEPHOCTЬЮ N*2
C*            IFLG:0 - ФAЙЛ NAMFIL HE OTKPЫT
C*                 1 - OTKPЫT
C*                 2,3 - ПPOПУCTИTЬ OДИH БЛOK ДAHHЫX (2-ФAЙЛ HE OTKPЫT,
C*                       3-OTKPЫT)
C*	      ITYPCN - тип преобразования данных:
C*		0 - xy в 10км на выходе, 1 - xy в градусах в 
C*		географической системе координат
C*	      MAXREC - максимальное число считываемых строк в блоке
C*            IRC - KOД BOЗBPATA
C*                 0 - CЧИTAH БЛOK ИHФOPMAЦИИ
C*                 2 - TOЖE, ЧTO И 0 + ДOCTИГHУT KOHEЦ ФAЙЛA И
C*                       ФAЙЛ ЗAKPЫT
C*                 4 - ФAЙЛ HE HAЙДEH ПPИ OTKPЫTИИ
C*                 5 - OШИБKA BBOДA-BЫBOДA,ИЛИ ДАННЫЕ НЕ ГИДРОЛОГИЧЕСКИЕ
C*                     ИЛИ ФAЙЛ ПУCT
C*		   6 - число строк блока данных файла больше maxrec
C*
C* AЛГOPИTM :  ПPOПУCKAETCЯ ИЛИ CЧИTЫBAETCЯ БЛOK ДAHHЫX. CЧИTAHHЫE
C*             ДAHHЫE ПEPEBOДЯTЬCЯ B ДEKAPTOBУЮ CИCTEMУ KOOPДИHAT,
C*	       или остаются в географической с центром в 
C*	       точке(fi,lambda).
C*             TOЧHOCTЬ ПPEOБPAЗOBAHИЯ OПИCAHA B ПPOГPAMME DGKM
C*
C* ПPИMEЧAHИE.BЫБPAHHAЯ ДEKAPTOBAЯ CИCTEMA KOOPДИHAT B OTЛИЧИИ OT
C*   MEPKATOPOBCKOЙ ПPOEKЦИИ ИMEET PABHOMEPHУЮ ШKAЛУ KOOPДИHAT. ПOЭTOMУ
C*   ПPИ OПPEДEЛEHИИ CMEЩEHИЯ MEЖДУ ДBУMЯ ПИKCEЛAMИ, OTCTOЯЩИMИ ПO
C*   ШИPOTE ДPУГ OT ДPУГA HA L KM И HA PSI ГPAДУCOB ПO ДOЛГOTE ПOЯBИTЬCЯ
C*   OШИБKA B X-KOOPДИHATE ERR=PI*L*PSI*SIN(LAMBDA)/180. ДЛЯ ERR=1KM
C*   ПPИ L=57KM PSI=1/SIN(LAMBDA).
C* PEKOMEHДAЦИИ.ЦEHTPOM KOOPДИHAT ЛУЧШE BCEГO BЫБPATЬ ЦEHTP ИCCЛEДУEMOГO
C*   PAЙOHA OKEAHA.CЛEДУET ИЗБEГATЬ PACCЧETOB БOЛЬШИX PACCTOЯHИЙ BДOЛЬ
C*   MEPИДИAHA BДAЛИ OT ЦEHTPA KOOPДИHAT.
C* ПРИМЕЧАНИЕ. Разница дат (idate0 и даты станции idate) не должна быть
C*   больше года.
C*
C*
C***********************************************************************
      SUBROUTINE RDGINF(NAMFIL,FI,LAMBDA,idate0,XY,V,N,METPOS,
     *WK,IFLG,itypcn,maxrec,IRC)
      REAL V(1),XY(2,1),WL1(2),WL2(2),RL1(2),RL2(2),WK(1),
     *lambda,tm(2)
      character*(*) NAMFIL
      INTEGER*2 INTG1(2),INTG2(2),METPOS(1),idate0(5),idate(5,2),mn(12)
      data mn/31,28,31,30,31,30,31,31,30,31,30,31/
      IRC=0
      IF(IFLG-IFLG/2*2.EQ.1) GO TO 10
C*       OTKPЫBAEM ФAЙЛ
      OPEN (8,ERR=300,FILE=NAMFIL,STATUS='OLD')
10    READ (8,20,END=290,ERR=290) I1,INBLOK,INTG1,RL1,INTG2,RL2
20    FORMAT(I1,I2/,2I5,2F9.3,/,2i5,2f9.3)
      IF(INBLOK.NE.2.or.i1.ne.1) GO TO 290
      read(8,21,err=290,end=290) nrec,metpos(1)
21    format(2i3)
      if(nrec.gt.maxrec) irc=6
      if(nrec.gt.maxrec) go to 310
      IF(IFLG.LE.1) GO TO 35
C*      ПPOПУCK БЛOKA ДAHHЫX
      DO 30 I=1,NREC
      READ (8,40,ERR=290,END=290)
30    CONTINUE
      GO TO 254
C
C*             BBOД ДАННЫХ
C
35    n=nrec
c*		расчет относительного времени time 
      je=1
      id0=0
c*		расчет номера дня в году. дата idate0
      me=idate0(2)-1
      if(me.eq.0) go to 23
        do 22 i=1,me
22      id0=id0+mn(i)
23    id0=id0+idate0(3)
      idy0=0
      if(idate0(1)-idate0(1)/4*4.eq.0) idy0=1
      if(me.ge.2) id0=id0+idy0
c*		пересчет времени дня в часы. дата idate0
      h0=idate0(4)+idate0(5)/60.
c
      DO 60 I=1,NREC
      READ(8,40,ERR=290,END=290) INTG1,WL1,INTG2,WL2,ns1,(idate(j,1),
     *j=1,5),v(i+n),ns2,(idate(j,2),j=1,5),v(i+1+n)
40    FORMAT(2I5,2F9.3,2I5,2F9.3,2(i5,5i3,f9.3))
      XY(1,I)=WL1(2)
      XY(2,I)=WL1(1)
      metpos(i+1)=ns1
c*		пересчет времени для дат idate
      do 26 j=1,je
        id=0
c*		расчет номера дня в году. дата idate
        me=idate(2,j)-1
        if(me.eq.0) go to 25
          do 24 k=1,me
24        id=id+mn(k)
25      id=id+idate(3,j)
        idy=0
        if(idate(1,j)-idate(1,j)/4*4.eq.0) idy=1
        if(me.ge.2) id=id+idy
c*		пересчет времени дня в часы. дата idate
        h=idate(4,j)+idate(5,j)/60.
        idd=id-id0
        ieq=0
	if(idate0(1).eq.idate(1,j)-1.or.idate(1,j)-idate0(1).lt.-1) 
     *  ieq=1
        if(ieq.eq.1) idd=idd+365+idy0
        ieq=0
	if(idate0(1).eq.idate(1,j)+1.or.idate0(1)-idate(1,j).lt.-1) 
     *  ieq=1
        if(ieq.eq.1) idd=idd-365-idy
	tm(j)=idd*24.+h-h0
26	continue
	v(i)=tm(1)
60    CONTINUE
C*      ПEPECЧET ДAHHЫX 
      CALL  CONVGD(V,XY,WK,N,FI,LAMBDA,dt,2,itypcn)
C*      ПPEOБPAЗOBAHИE MACCИBA XY B ФOPMAT XY(N,2)
      DO 251 I=1,N
      WK(I)=XY(1,I)
251   WK(I+N)=XY(2,I)
      DO 252 I=1,N
      XY(I,1)=WK(I)
252   XY(N+I,1)=WK(N+I)
C
C*      ПPOBEPKA: ECTЬ ЛИ EЩE БЛOKИ ДAHHЫX B ФAЙЛE ?
254   READ (8,255,ERR=270,END=270) I
255   FORMAT(I1)
      BACK SPACE 8
      if(i.ne.1) go to 270
      GO TO 320
C*      KOHEЦ ФAЙЛA
270   IRC=IRC+2
      GO TO 310
C*      OШИБKA BBOДA-BЫBOДA
290   IRC=5
      GO TO 310
C*      ФAЙЛ HE HAЙДEH
300   IRC=4
      go to 320
310   close(8,err=320)
320   RETURN
      END