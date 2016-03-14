c*********************************************************************
c*
c*  поиск минимума функции fun методами безусловной минимизации
c*
c* параметры: method - адрес метода минимизации
c*          kod,n,x,fx,d,dm,ds,parm,minl,out,a,b - параметры метода -
c*                описаны в выпуске 17 - библиотека ssplib
c*          numopx - массив с номерами фиксированных переменных.
c*                если numopx(i)<0, то это означает, что переменные с
c*              номерами в диапазоне numopx(i-1):iabs(numopx(i)) -
c*              не оптимизируемые. numopx(i)=0 - конец списка.
c*          itp - индекс управления печатью:
c*              itp<=0 - ничего не печатается
c*              itp=1(11) - печатаются начальные и конечные значения
c*                      переменных, функции (x,fx),код возврата и 
c*                    время работы.
c*              itp=2(12) - печатать все, кроме промежуточных
c*                    значений переменных
c*              itp=3(13) - печатать все
c*             itp>=10 - необходимо вызывать программу savmui через
c*                itp/10 итераций
c*          istp - число итераций, после которых печатаются
c*               промежуточные результаты
c*          irc - код возврата:
c*              irc=0 - нормальное завершение
c*              irc>0 - код возврата метода минимизации
c*              irc=-1,(-2) - минимизацию прервал оператор по 
c*                окончании итерации (немедленно)
c*              irc=-10 - неверно составлен массив numopx, или 
c*                      оптимизируемых переменных <1
c*
c*
c* алгоритм: минимизируется функция пользователя fun(x,fx). если метод
c*           минимизации не поисковый и пользователь не предоставил
c*         программы расчета градиента(grad) и верхнего треугольника
c*         матрицы гессе(gesse), то они рассчитываются численно.
c* ограничение. число переменных X не должно быть >100. для снятия
c*         этого ограничения необходимо перед вызовом minfun завести
c*         common-область /mfcom/ размером 5*n четырехбайтовых слов.
c* замечание. 1.Для повторного вызова метода с сохранением констант для
c*         численного расчета градиента и/или матрицы Гессе 
c*           необходимо убрать пpисвоение ngrad=-1 в программе minfun
c*          2.При itp>10 - вызывается программма пользователя:
c*          savmui(xbest,fxbest)
c*          3.При оптимизационных расчетах возможно диалоговое
c*          вмешательство при нажатии следующих клавиш:
c*          f1 - прервать немедленно
c*          f2 - прервать по окончании итерации
c*          f3 - вызвать немедленно программу savmui
c*          Программа savmui может быть как пользовательской, так 
c*          и программой пакета. При отсутствии программы поль-
c*          зователя происходит сохранение наилучших текущих значе-
c*          ний функции и аргумента в файле optxfx.sav . Для 
c*          считывания оптимальных значений рекомендуется пользовать-
c*          ся программой resxfx(xbest,fxbest,irc) с кодом возврата
c*          irc=0 при нормальном завершении и irc=1 при ошибках ввода
c*
c*
c*********************************************************************
	SUBROUTINE minfun(METHOD,KOD,N,X,FX,D,DM,DS,PARM,MINL,OUT,A,B,IRC,
     *numopx,itp,ISTP)
	external subf
	REAL*4 X(1),FX(1),D(1),DM(1),DS(1),PARM(8),A(1),B(1),FXOPT,XOPT
	REAL*4 NGRD,mlt
	integer*4 ibinar(100),two
	INTEGER KOD(8),NUMopx(1)
	COMMON/mfcmin/KOD1(8),IW,nnow,NN,M,IS,MAXIT,MAXFX,NGRAD,FXOPT,
     *NGRD,ityp1,istps,kd4,kd5
     */mfcom/XOPT(500)
     */grcom/ig,ing,ifst,nggr,mlt,ibinar
	data two/2/
C		Вывод входных параметров
      open(unit=15, file='minfun.log', action='write')
	write(15, *) 'N = ', N
	write(15, *) 'FI: ', X(1)
	write(15, *) 'LAMBDA: ', X(2)
	write(15, *) 'NODES:'
	nknot = (n - 5) / 2
	do i = 3, 3 + nknot - 1
		write(15, *) X(I), X(I+nknot)
	enddo
	write(15, *) 'A0 = ', X(2 + 2*nknot + 1)
	write(15, *) 'Rm = ', X(2 + 2*nknot + 2)
	write(15, *) 'R0 = ', X(2 + 2*nknot + 3)
c*          константы пpогpаммы grad
	ifst=1
	nggr=15
	mlt=1.0
	ig=0
	ing=0
	do 166 i=1,100
166   ibinar(i)=0
	itype=itp-itp/10*10
	istps=itp/10
	if(istps.lt.0) istps=0
c
c*         формирование массива оптимизируемых переменных
c
	if(n.gt.100.and.itype.gt.0) print 1
1     format(1x,'предупреждение: число переменных больше максимального')
	do 2 i=1, 100
2       ibinar(i)=0
	irc=-10
	if(n.lt.1) return
	do 4 i=1,n
4       xopt(n+i)=x(i)
	nnow=n
	if(numopx(1).eq.0) go to 122
c*         не все переменные функции fun оптимизируемые. формируем
c*           массив флажков-статусов (ibinar) и вектор оптимизируемых
c*         переменных X
	nxo=0
c*           проверка корректности списка numopx
	n1=n+1
	if(numopx(1).lt.0) go to 300
	do 5 i=2,n1
	  if(numopx(i).eq.0) go to 6
	  if(iabs(numopx(i)).gt.n) go to 300
	  if(iabs(numopx(i)).le.iabs(numopx(i-1))) go to 300
5     continue
6     nnow=0
	ie=i-1
	do 10 i=1,ie
	  nx=numopx(i)
	  ibg=nx
	  ien=nx
	  if(nx.gt.0) go to 7
c*           идет блок фиксируемых переменных
	  ibg=nxo+1
	  ien=iabs(nx)
7       do 8 j=ibg,ien
	    nw=j/31
		if(nw*31.ne.j) nw=nw+1
		idn=(nw)*31-j
8         ibinar(nw)=ibinar(nw)+two**idn
	  if(nxo+1.eq.ibg) go to 10
C*      отбоp оптимизиpуемых пеpеменных
	  iend=ibg-1
	  ibg=nxo+1
	  do 9 j=ibg,iend
		nnow=nnow+1
9         x(nnow)=xopt(n+j)
10      nxo=ien
	if(ien.eq.n) go to 12
c*      отбоp оптимизиpуемых пеpеменных
	ibg=ien+1
	do 11 i=ibg,n
	  nnow=nnow+1
11      x(nnow)=xopt(n+i)
12    if(nnow.lt.1) go to 300
c---------------------------------------------------------------------
122   MAXFX=KOD(4)
	MAXIT=KOD(5)
	IS=ISTP
	ityp1=itype
	nn=n
	DO 13 I=1,8
   13 KOD1(I)=KOD(I)
	NGRD=1.E30
	IW=-1
	write (15, *) 'MINFUN: call FUN'
	call fun(xopt(n+1),fx)
	write (15, *) 'MINFUN: FUN is completed'
	if (fx(1).eq.0.) then 
	  write(33, *) 'FX is NULL!'
	  write(15, *) 'FX is NULL! TERMINATED!'
	  RETURN
	endif
c*          Занесение начального приближения в массив xopt
	do 30 i=1,n
30    xopt(i)=xopt(n+i)
	fxopt=fx(1)
	kd4=0
	kd5=0
c
	NGRAD=-1
	if(itype.ge.1) PRINT 32,N,FX(1),(XOPT(N+I),I=1,N)
   32 FORMAT (3X,'N=',I3,3X,999(2X,'FX&X=',2X,6D10.3,/))
	if(itype.ge.2) PRINT 34,KOD,D(1),DM(1),DS(1)
   34 FORMAT (3X,'KOD=',8I8,/,3X,' D=',D10.3,' DM=',D10.3,' DS=',D10.3)
	if(itype.ge.2) PRINT 35,PARM
   35 FORMAT (3X,'PARM=',4d10.3)
	CALL gettim(ihr,imin,isec,i100th)
	timeb=3600.*ihr+60.*imin+isec+0.01*i100th
	CALL METHOD(KOD1,N,X,FX,D,DM,DS,PARM,SUBF,MINL,OUT,A,B,IRC)
	CALL gettim(ihr,imin,isec,i100th)
	timend=3600.*ihr+60.*imin+isec+0.01*i100th
	time=timend-timeb
	if(kd5.ne.0) irc=-1
	if(kd4.ne.0) irc=-2
	MX=MAXFX-KOD1(4)-kd4
	MT=MAXIT-KOD1(5)-kd5
	if(itype.ge.1) WRITE(6,40) MT,MX,NGRAD,FXOPT,NGRD
   40 FORMAT(5X,I5,5X,I5,7X,I5,4X,D18.10,2X,E18.10)
	if(itype.ge.1) PRINT 43,(XOPT(I),I=1,N)
   43 FORMAT(1X,'X=',2X,5D14.7)
	if(itype.ge.1) WRITE(6,37) IRC,time
37    FORMAT(3X,'KOД BOЗBPATA',2X,I4,' BPEMЯ PAБOTЫ ',F10.2,' CEK.')
	if(irc.lt.0) print *,'  Задание прекращено оператором'
	do 50 i=1,n
50    x(i)=xopt(i)
	fx(1)=fxopt
	if(itp.ge.10) then
		call  savmui(x,fx)
		write(15,*) 'MINFUN: Line185: call SAVMUI'
	endif
300   close(15)
      return
	END

c*************************************************************
c*
c*          SUBF
c*
c*
c*************************************************************
	SUBROUTINE subf(X,FX,GR,M)
	REAL*4 X(1),FX(1),GR(1),FXSAVE,ACC,NGRD,mlt
	REAL*8 NGNOW
	COMMON/mfcmin/KOD1(8),IW,NNOW,N,L,ISTP,MAXIT,MAXFX,NGRAD,fxopt,
     *NGRD,itype,istps,kd4,kd5
     */mfcom/ XOPT(500)
     */grcom/ ig,ing,ifst,nggr,mlt,ibinar(100)
	INTEGER NM(2),idold/0/,isav/0/,f1/15104/,f2/15360/,f3/15616/
	integer*4 ibinar,maska,s2(31)/16#40000000,16#20000000,16#10000000,
     *16#08000000,16#04000000,16#02000000,16#01000000,16#00800000,
     *16#00400000, 
     *16#00200000,16#00100000,16#00080000,16#00040000,16#00020000,
     *16#00010000,
     *16#00008000,16#00004000,16#00002000,16#00001000,16#00000800,
     *16#00000400,
     *16#00000200,16#00000100,16#00000080,16#00000040,16#00000020,
     *16#00000010,
     *16#00000008,16#00000004,16#00000002,16#00000001/,maskas
	 write(15, *) 'SUBF'
	write(*, *) 'SUBF'
	IDELT=MAXIT-KOD1(5)
	IF(IW.EQ.-1.and.itype.ge.1) WRITE(6,660)
660   FORMAT(3X,'  ЧИCЛO ',2X,'  ЧИCЛO  ',2X,'   ЧИCЛO  ',7X,'ЗHAЧEHИE',
     *10X,' HOPMA')
	IF(IW.EQ.-1.and.itype.ge.1) WRITE(6,662)
662   FORMAT(3X,'ИTEPAЦИЙ',2X,'BЫЧ. Ф-ИЙ',2X,'ГPAДИEHTOB',7X,'  Ф-ИИ  ',
     *10X,' ГPAД-TA')
	IW=0
	IF(NGRAD.NE.-1) GO TO 100
c*          acc - точность реального числа
	acc=epsilon(acc)
	DO 99 I=1,n
	  xopt(n*2+I)=100.e0*ACC
	  xopt(n*3+i)=0.
99      xopt(n*4+I)=5.
	ngrad=0
100   CONTINUE
c
c*          Опрос на наличие прерываний
c
C101   call  key(1,icodsm)
C	if(icodsm.eq.0) go to 106
c*          Произошло прерывание. Считываем код клавиши
C	call  key(0,icodsm)
c*          Распознаем код 
C	if(icodsm.eq.f3) isav=1
C	if(icodsm.eq.f1) kd4=kod1(4)
C	if(icodsm.eq.f1) kod1(4)=0
C	if(icodsm.eq.f2) kd5=kod1(5)
C	if(icodsm.eq.f2) kod1(5)=0
C	go to 101
c
106   IF(IDOLD.EQ.IDELT) GO TO 1
	ID=IDELT-IDELT/ISTP*ISTP
	IF(ID.NE.0)GO TO 1
	IF(IDELT.EQ.0)GOTO 1
	MFX=MAXFX-KOD1(4)
	if(itype.ge.2) PRINT 90 ,IDELT,MFX,NGRAD,FXSAVE,NGRD
  90  FORMAT(5X,I5,5X,I5,7X,I5,4X,D18.10,2X,E18.10)
	if(itype.ge.3) PRINT 92,(XOPT(N+I),I=1,N)
  92  FORMAT(10X,'ЗHAЧEHИЯ ПEPEMEHHЫX',1X,6E14.7)
1     continue
c*          Вызов программы savmui, если требуется
	IF(IDOLD.EQ.IDELT.and.isav.eq.0) GO TO 922
	if(istps.eq.0.and.isav.eq.0) go to 922
	if(istps.ne.0) id=idelt-idelt/istps*istps
	if(id.ne.0.and.isav.eq.0) go to 922
	if(idelt.eq.0.and.isav.eq.0) go to 922
	call  savmui(xopt,fx)
      write(15,*) 'MINFUN: Line256: call SAVMUI'
	isav=0
c
922   IDOLD=IDELT
	if(n.eq.nnow) go to 97
c*            формирование полного вектора переменных
	nwd=1
	nxw=1
	nxopt=0
	maskas=0
	do 95 i=1,n
	maska=maskas+s2(nxw)
	if(ibinar(nwd).ge.maska) go to 93
c*            переменная с номером i - оптимизируемая
	nxopt=nxopt+1
	xopt(n+i)=x(nxopt)
	go to 94
c*            переменная с номером i фиксированная
93    maskas=maska
c*            контроль номера элемента массива ibinar
94    nxw=nxw+1
	if(nxw.le.31) go to 95
c*            меняем номер элемента
	nxw=1
	nwd=nwd+1
	maskas=0
95    continue
	go to 108
c
c*        запоминаем текущие пеpеменные
97    do 102 i=1,n
102   xopt(n+i)=x(i)
c
108   NM(2)=M
	NM(1)=N
	IF(M-M/2*2.EQ.0) GO TO 2
	write(15, *) 'SUBF: call FUN'
	CALL FUN(XOPT(N+1),FX)
	write(15, *) 'SUBF: FUN returned FX =', FX(1)
	fxsave=fx(1)
	IF(FXOPT.LE.FX(1)) GO TO 2
	FXOPT=FX(1)
	DO 111 I=1,N
111     XOPT(I)=XOPT(N+I)
    2 IF(M/2-M/4*2.EQ.0) GO TO 3
      write (15,*) 'MINFUN: Line295: call GRAD'
	CALL GRAD(XOPT(N+1),FX,GR,NM,xopt(n*2+1),xopt(n*3+1),xopt(n*4+1))
C*     PACЧET HOPMЫ ГPAДИEHTA
	NGNOW=0.D0
	DO 20 I=1,N
20    NGNOW=NGNOW+GR(I)*GR(I)
	NGNOW=DSQRT(NGNOW)
	GN=NGNOW
	IF(GN.LT.NGRD) NGRD=GN
C
	NGRAD=NGRAD+1
3     IF(M.GT.3) then
		CALL GESSE(XOPT(N+1),FX,GR,NM,xopt(n*2+1),xopt(n*3+1),
	*xopt(n*4+1))	
		write(*,*) 'MINFUN: line312: call GESSE'
	ENDIF
	ig=0
	RETURN
	END
