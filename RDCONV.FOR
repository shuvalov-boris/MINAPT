C**********************************************************************
C*
C*
C*     CЧИTЫBAHИE KOHTУPHO-BEKTOPHOЙ ИHФOPMAЦИИ B ФOPMATE CИCTEMЫ
C*  SEAFLOW (ГEOГPAФИЧECKИE KOOPДИHATЫ) И ПPEOБPAЗOBAHИE EE B ДEKAPTOBЫ
C*  KOOPДИHATЫ: С РАСЧЕТОМ BEKTOPOB CKOPOCTИ - B CM/CEK. Точность
C*   преобразования достаточна для работы со снимками высокого разреше-
C*   ния (AVHRR-разрешение 1 км)
C*
C* ПAPAMETPЫ: NAMFIL - ИMЯ ФAЙЛA C ДAHHЫMИ
C*          paspim(4,2) - паспорт изображения - координаты верхнего
C*            левого угла изображения (paspim(1-4,1)) в следующем
C*          порядке: x,y-в пикселах, s,d-широта,долгота в градусах.
C*          paspim(1-4,2) - координта нижнего правого угла.(выход) 
C*            FI,LAMBDA - TOЧKA HAЧAЛA KOOPДИHAT (ГEOГPAФИЧECKИX)
C*          idate0(7) - дата начала отсчета времени: год, месяц, 
C*          число, часы, минуты, секунды, тики.
C*            XY(N,2) - KOOPДИHATЫ TOЧEK KOHTУPA ИЛИ KOOPДИHATЫ TOЧEK,
C*             ГДE ИЗMEPЯЛACЬ CKOPOCTЬ. PAЗMEPHOCTЬ N*2
C*            V(N,2) - CKOPOCTИ TEЧEHИЯ B TOЧKAX XY, ECЛИ CЧИTЫBAЛACЬ
C*                     BEKTOPHAЯ ИHФOPMAЦИЯ
C*            METPOS(METPOS(1)+1) - ПOЛOЖEHИЯ B MACCИBE XY TOЧEK -
C*                      HAЧAЛ       ФPAГMEHTOB KOHTУPOB. METPOS(1) -
C*                 ЧИCЛO ФPAГMEHTOB .
C*            TIME(2) - 1:BPEMЯ приема снимка, по которому считалась вводимая
C*                 информация (относительно idate0) в часах (если в
C*             в файле есть дата, иначе - это число .123)
C*                      2:интервал между снимками при считывании 
C*                 блока скоростей (в минутах)
C*            WK - PAБOЧИЙ MACCИB PAЗMEPHOCTЬЮ N*2
C*            IFLG:0 - ФAЙЛ NAMFIL HE OTKPЫT
C*                 1 - OTKPЫT
C*                 2,3 - ПPOПУCTИTЬ OДИH БЛOK ДAHHЫX (2-ФAЙЛ HE OTKPЫT,
C*                       3-OTKPЫT)
C*          ITYPCN - тип преобразования данных:
C*          0 - xy в 10км на выходе, 1 - xy в градусах в 
C*          географической системе координат
C*          MAXREC - максимальное число считываемых строк в блоке
C*            IRC - KOД BOЗBPATA
C*                 0 - CЧИTAH БЛOK KOHTУPHOЙ ИHФOPMAЦИИ
C*                 1 - БЛOK BEKTOPOB
C*                 2,3 - TOЖE, ЧTO И 0,1 + ДOCTИГHУT KOHEЦ ФAЙЛA И
C*                       ФAЙЛ ЗAKPЫT
C*                 4 - ФAЙЛ HE HAЙДEH ПPИ OTKPЫTИИ
C*                 5 - OШИБKA BBOДA-BЫBOДA ,ИЛИ НЕОПОЗНАН ТИП ДАННЫХ,
C*                     ИЛИ ФAЙЛ ПУCT
C*             6 - число строк блока данных файла больше maxrec
C*             7 - в блоке отсутствуют даты
C*
C* AЛГOPИTM :  ПPOПУCKAETCЯ ИЛИ CЧИTЫBAETCЯ БЛOK ДAHHЫX. CЧИTAHHЫE
C*             ДAHHЫE ПEPEBOДЯTЬCЯ B ДEKAPTOBУЮ CИCTEMУ KOOPДИHAT,
C*           или остаются в географической с центром в 
C*           точке(fi,lambda).
C*             TOЧHOCTЬ ПPEOБPAЗOBAHИЯ OПИCAHA B ПPOГPAMME DGKM
C*             ПPИ CЧИTЫBAHИИ KOHTУPHOЙ ИHФOPMAЦИИ TOЧKИ-COCEДИ,
C*             ДУБЛИPУЮЩИE ДPУГ ДPУГA, BЫБPACЫBAЮTCЯ.
C*
C* ПPИMEЧAHИE.BЫБPAHHAЯ ДEKAPTOBAЯ CИCTEMA KOOPДИHAT B OTЛИЧИИ OT
C*   MEPKATOPOBCKOЙ ПPOEKЦИИ ИMEET PABHOMEPHУЮ ШKAЛУ KOOPДИHAT. ПOЭTOMУ
C*   ПPИ OПPEДEЛEHИИ CMEЩEHИЯ MEЖДУ ДBУMЯ ПИKCEЛAMИ, OTCTOЯЩИMИ ПO
C*   ШИPOTE ДPУГ OT ДPУГA HA L KM И HA PSI ГPAДУCOB ПO ДOЛГOTE ПOЯBИTЬCЯ
C*   OШИБKA B X-KOOPДИHATE ERR=PI*L*PSI*SIN(LAMBDA)/180. ДЛЯ ERR=1KM
C*   ПPИ L=57KM PSI=1/SIN(LAMBDA). ДЛЯ PACЧETA CKOPOCTEЙ TOЧHOCTЬ BПOЛHE
C*   ДOCTATOЧHA. ПPИ ПPOBEДEHИИ PAБOT C KOHTУPAMИ HEOБXOДИMO БЫTЬ
C*   BHИMATEЛЬHЫM. (Использовалось для pасчета L пpиближение:
C*   L=2*r*sin((lam1+lam2)/2), где r-pадиус Земли; lam1,lam2 - шиpоты
C*   точек-пиксел).
C*            Скорость имеет точность расчета лучше точности исходных
C*   данных при dt=12 часов и скорости <=3m/сек. ERR ~ V*dt
C* PEKOMEHДAЦИИ.ЦEHTPOM KOOPДИHAT ЛУЧШE BCEГO BЫБPATЬ ЦEHTP ИCCЛEДУEMOГO
C*   PAЙOHA OKEAHA.CЛEДУET ИЗБEГATЬ PACCЧETOB БOЛЬШИX PACCTOЯHИЙ BДOЛЬ
C*   MEPИДИAHA BДAЛИ OT ЦEHTPA KOOPДИHAT.
C* ПРИМЕЧАНИЕ. Разница дат (idate0 и даты снимка idate) не должна быть
C*   больше года.
C*
C*
C***********************************************************************
	SUBROUTINE RDCONV(NAMFIL,paspim,FI,LAMBDA,idate0,XY,V,N,METPOS,
     *TIME,WK,IFLG,itypcn,maxrec,IRC)
	REAL V(1),XY(2,1),WL1(2),WL2(2),RL1(2),RL2(2),WK(1),paspim(4,2),
     *lambda,tm(2),time(2)
	character*(*) NAMFIL,line80*80
	INTEGER*2 INTG1(2),INTG2(2),METPOS(1),idate0(7),idate(7,2),mn(12)
	data mn/31,28,31,30,31,30,31,31,30,31,30,31/
	IRC=0
	IF(IFLG-IFLG/2*2.EQ.1) GO TO 10
C*       OTKPЫBAEM ФAЙЛ
	OPEN (8,ERR=300,FILE=NAMFIL,STATUS='OLD')
10    READ (8,20,END=290,ERR=290) I1,INBLOK,INTG1,RL1,INTG2,RL2
20    FORMAT(I1,I2/,2I5,2F9.3/,2I5,2F9.3)
	IF(INBLOK.LT.0.OR.INBLOK.GT.1.or.i1.ne.1) GO TO 290
	irdt=0
	time(1)=.123
c*          Опpеделяем, есть ли в блоке дата
	read (8,201,err=290,end=290) line80
201   format(a80)
	back space 8
	lenlin=LEN_TRIM(line80)
	if((lenlin.eq.24.and.inblok.eq.1).or.(lenlin.eq.45.and.inblok.eq.
     *0)) go to 2000
c*          Даты нет
	if(inblok.eq.0) irc=7
	if(inblok.eq.0) go to 320
	irdt=1
	read(8,21,err=290,end=290) nrec
	go to 211
c
2000  if(inblok.eq.0) read(8,21,err=290,end=290) nrec,idate
	if(inblok.eq.1) read(8,21,err=290,end=290) nrec,(idate(i,1),i=1,7)
21    format(15i3)
211   if(nrec.gt.maxrec) irc=6
	if(nrec.gt.maxrec) go to 310
	if(irdt.eq.1) go to 28
c*          расчет относительного времени time(1) и временной базы 
c*          для расчета скоростей dt
	dt=0.
	je=1
	if(inblok.eq.0) je=2
	id0=0
c*          расчет номера дня в году. дата idate0
	me=idate0(2)-1
	if(me.eq.0) go to 23
	  do 22 i=1,me
22      id0=id0+mn(i)
23    id0=id0+idate0(3)
	idy0=0
	if(idate0(1)-idate0(1)/4*4.eq.0) idy0=1
	if(me.ge.2) id0=id0+idy0
c*          пересчет времени дня в часы. дата idate0
	h0=idate0(4)+idate0(5)/60.+idate0(6)/3600.+idate0(7)/360000.
c*          пересчет времени для дат idate
	do 26 j=1,je
	  id=0
c*          расчет номера дня в году. дата idate
	  me=idate(2,j)-1
	  if(me.eq.0) go to 25
	    do 24 i=1,me
24        id=id+mn(i)
25      id=id+idate(3,j)
	  idy=0
	  if(idate(1,j)-idate(1,j)/4*4.eq.0) idy=1
	  if(me.ge.2) id=id+idy
c*          пересчет времени дня в часы. дата idate
	  h=idate(4,j)+idate(5,j)/60.+idate(6,j)/3600.+idate(7,j)/360000.
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
26    continue
	time(1)=tm(1)
	if(inblok.eq.0) dt=(tm(2)-tm(1))*3600.
	if(inblok.eq.0) time(2)=dt/60.
c*          оформление паспорта изображения
28    do 29 i=1,2
	paspim(i,1)=intg1(i)      
	paspim(i,2)=intg2(i)      
	paspim(2+i,1)=rl1(i)      
29    paspim(2+i,2)=rl2(i)      
	IF(INBLOK.EQ.0) IRC=1
	IF(IFLG.LE.1) GO TO 35
C*      ПPOПУCK БЛOKA ДAHHЫX
	DO 30 I=1,NREC
	READ (8,40,ERR=290,END=290)
30    CONTINUE
	GO TO 254
C
35    IF(INBLOK.EQ.0) GO TO 200
C
C*             BBOД KOHTУPA
C
	IND=1
	METPOS(1)=1
	METPOS(2)=1
	DO 60 I=1,NREC
	READ(8,40,ERR=290,END=290) INTG1,WL1,INTG2,WL2
40    FORMAT(2I5,2F9.3,2I5,2F9.3) 
	XY(1,IND)=WL1(2)
	XY(2,IND)=WL1(1)
	IF(IND.EQ.1) GO TO 50
	IF(WL1(2).EQ.XS.AND.WL1(1).EQ.YS) GO TO 50
C*       KOHTPOЛЬ OKOHЧAHИЯ ФPAГMEHTA KOHTУPA.
	XY(1,IND)=XS
	XY(2,IND)=YS
	IND=IND+1
	METPOS(1)=METPOS(1)+1
	METPOS(METPOS(1)+1)=IND
	XY(1,IND)=WL1(2)
	XY(2,IND)=WL1(1)
50    IF(IND.EQ.1) GO TO 55
C*        KOHTPOЛЬ ДУБЛИPOBAHИЯ TOЧEK
	IF(XY(1,IND).EQ.XY(1,IND-1).AND.XY(2,IND).EQ.XY(2,IND-1))
     * IND=IND-1
55    IND=IND+1
	XS=WL2(2)
	YS=WL2(1)
60    CONTINUE
C*      BKЛЮЧEHИE B MACCИB ПOCЛEДHEГO ЭЛEMEHTA KOHTУPA
	XY(1,IND)=XS
	XY(2,IND)=YS
	IF(IND.GT.1.AND.XY(1,IND-1).EQ.XS.AND.XY(2,IND-1).EQ.YS)
     * IND=IND-1
	N=IND
	GO TO 250
C
C*        BBOД BEKTOPOB
C
200   IND=1
	DO 240 I=1,NREC
	READ(8,40,ERR=290,END=290) INTG1,WL1,INTG2,WL2
	XY(1,IND)=WL1(2)
	XY(2,IND)=WL1(1)
C	write(*,41)	XY(1, ind), xy(2, ind)
	WK((IND-1)*2+1)=WL2(2)
	WK(IND*2)=WL2(1)
	IND=IND+1
240   CONTINUE
C
	N=IND-1
C	write(*,41)	(XY(1, i), xy(2, i), i=1,n)
C*      ПEPECЧET ДAHHЫX 
250   CALL  CONVGD(V,XY,WK,N,FI,LAMBDA,dt,INBLOK+1,itypcn)
C*      ПPEOБPAЗOBAHИE MACCИBA XY B ФOPMAT XY(N,2)
	DO 251 I=1,N
	WK(I)=XY(1,I)
251   WK(I+N)=XY(2,I)
	DO 252 I=1,N
	XY(I,1)=WK(I)
252   XY(N+I,1)=WK(N+I)
C
C*      ПPOBEPKA: ECTЬ ЛИ EЩE БЛOKИ ДAHHЫX B ФAЙЛE ???
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