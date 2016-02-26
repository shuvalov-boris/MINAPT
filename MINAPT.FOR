!*********************************************************************
!*
!*
!*              M I N A P T
!*
!*       Расчет динамической топографии к востоку от Японии
!*                (синоптический подход)
!*
!*    nbobj(nobj) - список блоков-обьектов, участвующих в оптимизации
!*          если nbobj(i)<0, то параметры блока не оптимизируются и
!*      в список переменных не включаются
!*    nblopt(nopt) - список блоков с данными для восстановления через
!*          минимизацию параметров обьектов
!*    ifl(1)=1 - вычесть из полей скорости и/или динамической топографии
!*        вклады в них блоков не участвующих в оптимизации (требуется
!*      для ускорения счета)
!*    ifl(2)=1 - в конце оптимизации печатаются рассогласования полей
!*
!*
!*********************************************************************
	common/gddat/ icninf,istinf, nknot, iflgtp, iflgvr
     1	, iflgtd, x0, xr, epsil, lentyp, nb, idat, nbs, nbvel
	common /mindat/ifl(10),nobj,nbobj(20),nopt,nblopt(20)
	common /optres/ fxopt, iFunTrialNumber
$include:'geodyn.for'
      EXTERNAL ragl, min03, min05, min14, minl5, OUT2
	INTEGER MINFUN_TRIAL_NUMBER, OPT_VAR_NUMBER
	PARAMETER (
     *  MINFUN_TRIAL_NUMBER   = 33*33*10,	!Число обращений к минимизационной функции
     *  OPT_VAR_NUMBER	= 33		!Число возможных оптимизируемых переменных
     *)
      REAL*4 A(OPT_VAR_NUMBER*(MINFUN_TRIAL_NUMBER+4)), 
     *    B(OPT_VAR_NUMBER*MINFUN_TRIAL_NUMBER),
!     *	PARM(8)/0.0, 1.e-02, 300., 1.e-5, 3., 1.e-4, 1.e-5, 1.e-1/,
     *	PARM(8)/0.0, 1., 300., 1., 3., 1., 1., 1./,
     *	D(MINFUN_TRIAL_NUMBER)  /MINFUN_TRIAL_NUMBER * 0.01/,
     *	DM(MINFUN_TRIAL_NUMBER) /MINFUN_TRIAL_NUMBER * 0.e0/,
     *	DS(MINFUN_TRIAL_NUMBER) /MINFUN_TRIAL_NUMBER * 1.e-5/
	INTEGER !KOD(8)/11,0,9,1800,200,3*0/,
     *    KOD(8)/						! min03
     *		OPT_VAR_NUMBER,	!
     *		0,						!
     *		0,
     *		MINFUN_TRIAL_NUMBER,			! кол-во вычислений мин. функции???
     *		200,					! кол-во итераций ?
     *		0,						! точность предсталения числа (0-простая, 1-двойная) ???
     *		2*0/,
     *	NUMF(OPT_VAR_NUMBER+1)/23,-24,0,4,5,6,7,8,9,10,11,12,13,14,15	 !#todo
     *		,16,17,18,19,20,21,22,23,24,25,26,27,28,29,-30,0,0,0,33/
	* !,NUMF(33)/-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
	* !-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0
     * !, 0, 0/
	REAL X(OPT_VAR_NUMBER*MINFUN_TRIAL_NUMBER),
     *     FX(MINFUN_TRIAL_NUMBER*10+1),datblk(100),prm(97)
	integer nprm(20),nspblk(20),ibgobj(20),ibprm(20), np
	integer idatbl(200),ipasp(6),iprm(194), optimizationIsEnabled
	equivalence (datblk,idatbl,ipasp),(ipasp(3),nspbl),
     *(ipasp(4),npr),(ipasp(5),nf),(datblk(4),prm),(idatbl(7),iprm)
	real t(NMAX),wk(NMAX*(NKMAX+6)*8), piw(4,2)
	integer*2 iwk(202), nbcont, nbs
	data MAXREC/1000/, fxopt/999999999.0/
	
	write (*, '(30I5)')  NUMF
    1 format(30I5)
	open(23, FILE='minapt.log', ACTION='write')
	count = nargs()
	write(*, *) 'Количество параметров коммандной строки', count

C			ЧТЕНИЕ ВХОДНЫХ ДАННЫХ
      open(unit=11, file=NAMARC, action='read', status='old')
      read(11, *, ERR=999) where, fi, lambda, idate0							 
	*	, optimizationIsEnabled
	
	open(unit=33, file=RESFILE, action='write')	   
	
	nobj = 1
		
   10 read(11, *, END=999) nknot, itypcn, iflgtp, iflgsp, iflgvr, iflgtd
     1	, x0, xr, eps, epsil
      nb = 0
c*		расчет узлов
c	
      mr1=maxrec+1
	lnam=LEN_TRIM(where)
	call  inpsfd(fi,lambda,idate0,itypcn,where(:lnam),
     *  nb,lentyp,dat,idat,
     *  NBMAX,LENMAX,piw,MAXREC,iwk,wk,wk(mr1*2+1),wk(mr1*4+1),irc)
	write(6,17) irc
   17 format(3x,'программа inpsfd завершилась с irc=',i4)
	ircis=irc	
C			Поиск блока с контуром
	do i = 1, nb
	  if (lentyp(2, i).eq.1) nbcont = i
	  if (lentyp(2, i).eq.0) nbvel = i
      enddo	
C			Вывод координат точек контура в файл для визуадизации
      cibg = lentyp(1, nbcont) + 4
      incnt = idat(lentyp(1, nbcont) + 3)
      call WriteXY(incnt, dat(cibg), dat(cibg+incnt), 0., 0., 
     *  fi, lambda, 'contour.txt', 0, dat(cibg))

   15 format(2F11.3)
C			Расчет радиуса вихря (x0)
      write(*, *) ' Расчет радиуса вихря'
	x0sum = 0.0
	x0min = 100000.0 !
	x0max = -x0min
      do i = 0, incnt - 1
	  x0cur = sqrt((dat(cibg + i))**2 + (dat(cibg + incnt + i))**2)
	  x0sum = x0sum + x0cur
	  x0max = max(x0max, x0cur) 
	  x0min = min(x0min, x0cur)
	enddo
	x0 = x0sum / incnt
	xr = x0 * 2.5
	
      nbs=nbcont					   
      call  datwrk(dat,idat,lentyp,lenmax,nb,nbmax,nbs,x,y,t,ibgwk4,
     *  nknot,8,iflg,irc)	
      write(6, *) ' Подпрограмма DATWRK(8) завершилась с irc =', irc
      if(irc.ne.0) write (*, *) 'DATWRK:  Нет места в архиве'
      indcnt = idat(lentyp(1, nbs)+3)
	indbg = lentyp(1, nbs) / 2 + 4
	write(*, '(3F11.5)') (dat(indbg+i), 
     *  dat(indbg+i+indcnt), dat(indbg + i + indcnt * 2), i=0,indcnt-1)
      lenc=(indcnt-1)*3
	write(*, *) 'MA: lenc =', lenc
C	nknot = indcnt
	call CalculateSplineCoefficients(NKNOT+1, dat(indbg),  
     *  dat(indbg+indcnt), dat(indbg+indcnt*2), dat(indbg+indcnt*3), 
     *  dat(indbg+indcnt*3+lenc), WK,	IRCcsc)
     	idat(lentyp(1, nbs)+5)=1
      call WriteXY(indcnt, dat(indbg), dat(indbg+indcnt), 0., 0.,
     *  fi, lambda, 'nodes0.txt', 2, dat(indbg+indcnt*2))
	call WriteSplineXY(int(indcnt), dat(indbg), dat(indbg+indcnt), 
     *  dat(indbg+indcnt*2), dat(indbg+indcnt*3), fi, lambda,
     *  'spline0.txt')
	
      nbobj(1) = nbs

c		сортировка блока скоростей по возрастанию расстояния до центра вихря
C			Вывод координат начала векторов скоростей в файл для визуадизации
      ivibg = lentyp(1, nbvel) / 2 + 4
	ivcnt = idat(lentyp(1, nbvel) + 3)
	call WriteXY(ivcnt, dat(ivibg), dat(ivibg+ivcnt), 0., 0., 
     *  fi, lambda, 'velbeg.txt', 0, dat(ivibg))
           
	open(31, FILE='distbeforesort.txt', ACTION='write')
      do j = 1, ivcnt
        dr = sqrt(dat(ivibg+j)**2+dat(ivibg+ivcnt+j)**2)
	  write(31, *) dr
	  end do	  
	close(31)

      write(*, *) 'minapt: X(1), Y(1)', dat(ivibg), dat(ivibg+ivcnt)
	call Sort(ivcnt, dat(ivibg), dat(ivibg+ivcnt), dat(ivibg+ivcnt*2), 
     *   dat(ivibg +ivcnt*3), 0., 0., fi, lambda)												  	

	open(32, FILE='distsort.txt', ACTION='write')
      do j = 0, ivcnt-1
        dr = sqrt(dat(ivibg+j)**2+dat(ivibg+ivcnt+j)**2)
	  write(32, *) dr
	end do	  
	close(32)

c
c*               заполнение блока архива (блок не последний)   *******
c
c-----------------------------------------------------------------
c*          константы блоков, описывающих обьекты		????? на что влияет? убрать?
	n=0
	do 8 i=1,nobj
	  ibgobj(i)=lentyp(1,iabs(nbobj(i)))  !начало блока
	  nprm(i)=idat(ibgobj(i)+3)	  !число точек
	  if(lentyp(2,iabs(nbobj(i))).eq.2) then
		nprm(i)=nprm(i)*2
!		write(*, *) 'MINAPT: lentyp(2,iabs(nbobj(i))) = ', 
!     *		lentyp(2,iabs(nbobj(i))), ', nprm = ', nprm(i)
	  endif
	  nspblk(i)=idat(ibgobj(i)+2)
	  nf=idat(ibgobj(i)+4)		!число фрагментов контура?
	  if(nf-nf/2*2.eq.1) nf=nf+1
	  ibgobj(i)=ibgobj(i)/2+nf/2+4
	  ibprm(i)=n+1
	  if(nbobj(i).gt.0) n=n+nprm(i)
!	write(*, *) 'MINAPT: nprm(', i, ') = ', nprm(i)
8     continue
!      write (*, *) 'MINAPT: constant block, n = ', n, ', nb = '
!	*  , nbobj(1), ', nobj = ', nobj
C	, nbobj(1)


	do 110 i=1,nobj
	  if(nbobj(i).le.0) go to 110
	  do 105 j=0,nprm(i)-1
105       x(ibprm(i)+j)=dat(ibgobj(i)+j)
110   continue	   
c
C		Вывод данных блоков для проверки
	open(15, FILE='idat.log', action='write')
	cibg=lentyp(1,nbs)
	 write(15, 111) (idat(i), i=1, 7000)
	close(15)
  111 format(I10)
	open(16, FILE='dat.log', action='write')
	cibg=lentyp(1,nbs)
	 write(16, 112) (dat(i), i=1, 4000)
  112 format(F15.6)
      close(16)

C		Инициализация массива (X) оптимизируемых параметров
      n = 0
      X(1) = 0.0
	X(2) = 0.0
	write(23, *) 'nknot =', nknot
	cibg=lentyp(1,nbs)
      ibg=(cibg+6)/2
      do i = 1, nknot
	  X(i + 2) = dat(ibg + i)
	  X(i + nknot + 2) = dat(ibg + nknot + 1 + i)
	  n = n + 2
	enddo
      X(2+ nknot * 2 + 1) = 13. * 48.5 / 2
	X(2+ nknot * 2 + 2) = 5.
	X(2+ nknot * 2 + 3) = 13

      n = 2 + n + 3

	write(23, *) ' N = ', N

	write(23, *) ' Nodes0:'
	write(23, 16) (dat(ibg + i), dat(ibg + i + nknot), i = 1, nknot+1)
   16 format(2F15.6)
      write(23, *) 'X: ', (X(I), I=1,n)

C      call CheckRm()	

!	open(unit=33, file=RESFILE, action='write')
	
	if (optimizationIsEnabled.eq.1) then 
	  write (*, *) 'MINAPT: call MINFUN'
	  CALL minfun(min03,KOD,n,X,FX,D,DM,DS,PARM,MINL5,OUT2,A,B,IRC,
     *	NUMF,3,1)
	  write (*, *) 'minfun. irc = ', irc
	  print *,'n,x=',n,(x(i),i=1,n)
	endif
C	 call  fun(x,fx)

  999 close(26)
      close(33)
	close(23)	
200	END

c-------------------------------------------------------------
	block data
	common /mindat/ifl(10),nobj,nbobj(20),nopt,nblopt(20)
	data nobj/1/,nbobj/3,6,8,10,16*0/
	data nopt/2/,nblopt/5,13,18*0/,ifl/0,9*0/ 
	end



	SUBROUTINE Fun(X,FX)
c**********************************************************************
c*
c*						FUN
c*
c*          Расчет целевой функции
c*
c*	Параметры:
c*		X - массив с параметрами для расчета целевой функции (вход)
c*		    X(1) = cx				   - x-координата центра вихря
c*		    X(2) = cy				   - y-координата центра вихря
c*		    X(3)..X(3+nknot-1)		   - x-компоненты узлов сплайна
c*		    X(3+nknot)..X(3+nknot*2-1) - y-компоненты узлов сплайна
c*		    X(3+nknot*2)   = A0	- максимальный перепад функции тока 
c*								  (амплитуда)
c*		    X(3+nknot*2+1) = Rm	- точка (радиус) зоны максимальных  
c*								  скоростей до центра вихря
c*		    X(3+nknot*2+2) = R0	- точка (радиус) на переферии вихря,
c*								  где скорости нулевые
c*				nknot - количество узлов сплайна
c*		FX - значение целевой функции (выход)
c*
c*	Используется модель вихря, основанная на радиальной функции тока:
c*                    { A0*(1-r**2/R0*Rm)           0  <  r <= Rm
c*			   F(r)={ A0/(R0(R0-Rm))*(r-R0)**2    Rm <= r <= R0
c*                    { 0							        r >  R0	
c***********************************************************************
	common/gddat/ icninf,istinf, nknot, iflgtp, iflgvr
     1	, iflgtd, x0, xr, epsil, lentyp, nb, idat, nbs, nbvel
C	common /mindat/ifl(10),nobj,nbobj(20),nopt,nblopt(20)
	common /optres/ fxopt, iFunTrialNumber
$include:'geodyn.for'
	real x(1), A0, Rm, R0
	real t(NMAX),wk(NMAX*(NKMAX+6)*8)
	real vg(VMAX), vag(VMAX), dist(VMAX), vr(VMAX),
	*  cx, cy, dRm
	integer*2 nbs

      write(*, *) 'FUN: running, nknot =', nknot
	write(33, *) 'Начало выполнения FUN'

C		Извлечение параметров из массива X
      cx = X(1)
	cy = X(2)  
      A0 = X(2 + nknot * 2 + 1)
      Rm = X(2 + nknot * 2 + 2)
      R0 = X(2 + nknot * 2 + 3)
	if (Rm.eq.R0) then
	  write(*, *) 'Rm = R0 =', R0, Rm
	  print*, 'MINAPT: ERROR: Rm = R0, Terminated'
	  return
	endif
	dRm = R0
c		копирование блока узлов из массива X в массив dat
!      write(*, *) 'nbvel =', nbvel
!      write(*, *) 'lentyp =', (lentyp(i, nbvel), i=1, 2)
!      write(*, *) 'копирование блока узлов в массив dat'      
      call  datwrk(dat,idat,lentyp,lenmax,nb,nbmax,nbs,X(3),X(3+nknot),t
     *  ,ibgwk4,nknot,9,iflg,irc)
	write(*, *) ' Подпрограмма DATWRK(9) завершилась с irc =', irc

C		Вывод данных блоков для проверки
! 	open(25, FILE='idat_.log', action='write')
!	write(25, 111) (idat(i), i=1, 7000)
!	close(25)
!  111 format(I10)
!	open(16, FILE='dat_.log', action='write')
!      write(16, 112) (dat(i), i=1, 4000)
!  112 format(F15.6)
!      close(16) 

C		Вывод значений параметров целевой функции после оптимизации
      write(33, *) 'FUN: X:'
      write(33, 5) X(1), X(2)
    5 format(5x,'X0'10x,'Y0'/2F15.6)
	ibgfr=lentyp(1,nbs)
	write(33, *) ' nodes:    x       y          t'
	write(33, 18) (dat(ibgfr/2+3+i), dat(ibgfr/2+3+nknot+1+i),
     *   dat(ibgfr/2+3+(nknot+1)*2 + i), i = 1, nknot + 1)
   18 format(3F15.6)
      write(33, 7) X(3+nknot*2), X(4+nknot*2), X(5+nknot*2)
    7 format(10x,'A0',13x,'Rm',13x,'R0'/3F15.6)	 

c*          Расчет сплайн-коэффициентов
      cibg=lentyp(1,nbs) / 2 + 4
      lenc=(nknot)*3
	call CalculateSplineCoefficients(NKNOT+1, dat(cibg), 
     *	dat(cibg+nknot+1), dat(cibg + (nknot+1) * 2),  
     *	dat(cibg + (nknot+1) * 3), 
     *    dat(cibg + (nknot+1) * 3 + lenc), WK,	IRCcsc)
      if (IRCcsc.ne.0) then
	  write(33, *) 'Ошибка при рассчете коэффициентов сплайна, irc =',
     *	irc
	  goto 999
	endif		  	  
	idat(lentyp(1, nbs)+5) = 1
!      write(*, *) ' Сплайн-коэффициенты посчитаны'
	
C		Создание файла для визуализации оптимизированных узлов
	call WriteXY(nknot+1, dat(cibg), dat(cibg+nknot+1), 0., 0., fi,  
     *  lambda, 'nodes_opt.txt', 2, dat(cibg+(nknot+1)*2))
C		Создание файла для визуализации соответствующего сплайна
	call WriteSplineXY(int(nknot+1), dat(cibg), dat(cibg+nknot+1), 
     *  dat(cibg+(nknot+1)*2), dat(cibg+(nknot+1)*3), fi, lambda,
     *  'spline_opt.txt')

C*		РАСЧЕТ ТАНГЕНЦИАЛЬНЫХ КОМПОНЕНТ СКОРОСТИ
!      write(*, *) 'Разложение векторов скорости на компоненты'
      call vprep(fi, lambda, dat, idat, lentyp, nbvel, nbs, x0, xr, 
     *  dist,vg,vag,tga,vr,nv,iflgvr,epsil,irc)
	write(*,160) irc
160	format(3x,'программа vprep. irc=',i4)
!!!! ПРОВЕРИТЬ МЕНЯЕТ ЛИ VPREP УЗЛЫ ???									//TODO

C		Создание файла для визуализации соответствующего сплайна
	call WriteSplineXY(int(nknot+1), dat(cibg), dat(cibg+nknot+1), 
     *  dat(cibg+(nknot+1)*2), dat(cibg+(nknot+1)*3), fi, lambda,
     *  'spline_opt1.txt')

c               nbvel (блок скоростей)
      ibg=lentyp(1,nbvel)	 !начало блока скоростей в массиве idat
	nvi=idat(ibg+3)		 !количество маркеров скоростей
	ix=ibg/2+4			 !начало x-координат точек в массиве dat
	iy=ix+nvi			 !начало y-координат точек в массиве dat
	write(*, *) 'nbvel =', nbvel, ', ibg =', ibg, ', nvi =', nvi

c		Вывод в файл расстояний до центра вихря упорядоченных скоростей
!	open(32, FILE='dists.txt', ACTION='write')
!      write(32, *) 'cx&cy =', cx, cy
!      do j = 0, nvi-1
!        dr = sqrt((dat(ix+j)-cx)**2+(dat(iy+j)-cy)**2)
!	  write(32, *) dr
!	end do	  
!	close(32)

C		Расчет A0 и Rm путем перебора интервалов по радиальному профилю...
C		Формулы для A0 и Rm посчитаны путем минимизации целевой функции
C		При первой итерации производится расчет целевой функции с
C			параметрами A0 и Rm по умолчанию
	CALL ApproximateProfile(nvi, dat(ix), dat(iy), int(nknot+1),  
     *	dat(cibg), dat(cibg+nknot+1), dat(cibg + (nknot+1) * 2),  
     *	dat(cibg + (nknot+1) * 3), Vg, Vag, 
     *    Fi, Lambda, CX, CY, Rm, A0, R0)

      write(33, '(3F15.6)') A0, Rm, R0

	fx = ObjectiveFunctionValue(nvi, dat(ix), dat(iy), Vg, Vag, 
     *    Fi, Lambda, CX, CY, A0, Rm, R0, 1)

c*		Вывод pезультов pасчета
C      write(33, 13) ircis, irctp, ircsp, irc, 
C     *  itypcn, iflgtp, iflgsp, iflgvr, iflgtd, nknot, x0min, x0, x0max
C     *  , xr, sperr, nv
C   13 format(4I2, ' ', 5I1, I3, 5F8.3, I5)	 

        write (33, *) 'Thus FX = ', FX
	  write (33, *) 'FXOPT = ', FXopt
	  if (FX < FXopt) then
	    FXopt = FX
	    write(33, *) ' NEW OPTIMUM!$ ', FX
	  end if
	  iFunTrialNumber = iFunTrialNumber + 1
  999 write(*, *) 'RUN: the end iFTN = ', iFunTrialNumber
      
      close(11)
	close(22)

c	nknot = nknot - 1

C	write(33, *) 'FUN: X: after all'
C      write(33, 5) X(1), X(2)
C      write(33, *) 'Nodes:    X             Y'
C	write(33, 6) (X(i), X(nknot+i),  i = 3, 2 + nknot)
C      write(33, 7) X(3+nknot*2), X(4+nknot*2), X(5+nknot*2)

C		проверка на изменение оптимизируемых параметров в функции FUN
	cibg = lentyp(1, nbs) / 2 + 4

	call CheckOnEqualityOptimizationBlock(2+nknot*2+3, int(nknot), X,  
     *    cx, cy, dat(cibg), dat(cibg+nknot+1), R0, irc)
	if (irc > 0) then 
	  write(*, *) 'Ошибка: изменен блок оптимизируемых',
	*    ' переменных. irc =', irc
	  stop 'изменен блок оптимизируемых переменных'
	end if


	write(33, 456)
  456 format(150('$'))
	
	END



c*********************************************************************
c*
c*							Swap
c*
c*		обмен значениями между двумя переменными A и B
c*
c*********************************************************************
	SUBROUTINE Swap(A, B)
	  IMPLICIT NONE
	  REAL, INTENT(INOUT) :: A, B
	
	  REAL tmp

	  tmp = A
	  A = B
	  B = tmp

	END													
	

c*********************************************************************
c*
c*							Swap4
c*
c*		обмен значениями между двумя элементами массивов
c*
c*	Параметры:
c*		X, Y, VX, VY, Z - массивы
c*		i, j - индексы элементов массивов для обмена значениями
c*
c*********************************************************************
	SUBROUTINE Swap4(X, Y, VX, VY, Z, i, j)
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: i, j
	  REAL, INTENT(INOUT) :: X(1), Y(1), VX(1), VY(1), Z(1)

	  CALL Swap( X(i),  X(j))
	  CALL Swap(VX(i), VX(j))
	  CALL Swap( Y(i),  Y(j))
	  CALL Swap(VY(i), VY(j))
	  CALL Swap( Z(i),  Z(j))

	END



c*********************************************************************
c*
c*							Sort
c*
c*		Сортировка скоростей блока по возрастанию расстояния 
c*			   до центра вихря
c*
c*        Используется алгоритм сортировки методом выбора
c*
c*	Параметраы:
c*		N - количество векторов скоростей
c*		X(N)  - x-координаты начальной точки векторов
c*		Y(N)  - y-координаты начальной точки векторов
c*		VX(N) -	x-координаты скоростей
c*		VY(N) - y-координаты скоростей
c*        X0, Y0 - координаты центра вихря
c*		fi, lambda - начало географической системы координат 
c*                     (долгота, шиpота)
c*					   
c*********************************************************************
	SUBROUTINE Sort(N, X, Y, VX, VY, X0, Y0, fi, lambda)
	  IMPLICIT NONE
	  INTEGER, INTENT(IN)    :: N
	  REAL,    INTENT(IN)    :: X0, Y0, fi, lambda
	  REAL,    INTENT(INOUT) :: X(N), Y(N), VX(N), VY(N)

	  REAL dist(N), minimum
	  INTEGER min_idx, i, j

	  DO i = 1, N
	    dist(i) = sqrt((X0 - X(i))**2 + (Y0 - Y(i))**2)
	  END DO

	  DO i = 1, N-1
	    minimum = dist(i)
	    min_idx = i
	    DO j = i+1, N
	      IF (dist(j).lt.minimum) THEN
	        minimum = dist(j)
			min_idx = j
	      ENDIF
	    END DO
		CALL Swap4(X, Y, VX, VY, dist, i, min_idx)
	  END DO
	  CALL WriteXY(N, X, Y, 0., 0., fi, lambda, 'velbegsort.txt', 1, 
     *    dist)
	END



c*********************************************************************
c*
c*						CalculateSplineCoefficients
c*
c*		Расчет коэффициентов сплайна по узлам
c*
c*	Параметры:
c*		NKNOT - количество узлов сплайна (вход)
c*		XK - x-координаты узлов сплайна	 (вход)
c*		YK - y-координаты узлов сплайна	 (вход)
c*		TK - t-координата узлов сплайна	 (вход)
c*		CX - сплайн-коэффициенты для x-координаты (выход)
c*		CY - сплайн-коэффициенты для y-координаты (выход)
c*		IRC - ошибка выполнения подпрограммы
c*********************************************************************
	SUBROUTINE CalculateSplineCoefficients(NKNOT, XK, YK, TK, CX, CY,
     *	WK, IRC)
	  IMPLICIT NONE
	  INTEGER, INTENT(IN)  :: NKNOT
        REAL, INTENT(IN)     :: XK(1), YK(1), TK(1)
	  REAL, INTENT(OUT)    :: CX(1), CY(1)
	  REAL				   :: WK(1)
	  INTEGER, INTENT(OUT) :: IRC

	  INTEGER irc1, irc2, i

	  write(*, *) 'CalculateSplineCoefficients'
	  write(*, *) 'NKNOT =', 	NKNOT
	  write(*, *) 'XK YK TK'
	  write(*, *) (XK(i), YK(i), TK(i), i = 1, NKNOT) 

	  IRC = 0

	  if((XK(1).eq.XK(NKNOT)).and.(YK(1).eq.YK(NKNOT))) goto 5
c*        сплайн не замкнутый
        write(*, *) ' сплайн не замкнутый'
	  call icsccu(TK, XK, NKNOT, CX, NKNOT - 1, irc1)
	  if(irc1.ne.0) print *,'  Ошибка сплайн-аппроксимации irc=',irc1
	  call icsccu(TK, YK, NKNOT, CY, NKNOT - 1, irc2)
	  if(irc2.ne.0) print *,'  Ошибка сплайн-аппроксимации irc=',irc2
	  if (irc1.ne.0.or.irc2.ne.0) IRC = 1
	  RETURN
c*         сплайн замкнутый
5       continue
	  write(*, *) ' сплайн замкнутый'
	  call icspln(TK, XK, NKNOT, CX, NKNOT - 1, WK, irc1)
	  if(irc1.ne.0) print *,'  Ошибка сплайн-аппроксимации irc=',irc1
	  if (irc1.eq.131) write(*, *) 'TK: ', (TK(I), I=1,NKNOT)
	  call icspln(TK, YK, NKNOT, CY, NKNOT - 1, WK, irc2)
	  if(irc2.ne.0) print *,'  Ошибка сплайн-аппроксимации irc=',irc2
	  if (irc2.eq.131) write(*, *) 'TK: ', (TK(I), I=1,NKNOT)
	  if (irc1.ne.0.or.irc2.ne.0) IRC = 1

	END



c*********************************************************************
c*
c*						VelocityTeor
c*
c*  	   Расчет теоретической скорости согласно модели вихря
c*	    
c*	   A0, Rm, R0 описаны в заголовке подрпограммы Fun
c*	   dr - расстояние от точки, в которой необходимо посчитать
c*			теоретическую скорость, до центра вихря
c*
c*		Функция возвращает значение теоретической скорости в точке,
c*			находящейся на расстояии dr от центра вихря	
c*********************************************************************
      REAL FUNCTION VelocityTeor(dr, A0, Rm, R0)
        IMPLICIT NONE
	  REAL, INTENT(IN) :: dr, A0, Rm, R0
	  
	    if (dr.lt.Rm) then 
	      VelocityTeor = - 2 * A0 * dr / (R0 * Rm)
        else 
	      VelocityTeor = 2 * A0 / (R0 * (R0 - Rm)) * (dr - R0) !*(-1) ? см. Ф. тока	  
	  endif
      END



c*********************************************************************
c*
c*					ApproximateProfile
c*			
c*	  	Аппроксимация радиального профиля скорости вихря 
c*		путем расчета амплитуды A0 и радиального значения максимальной скорости
c*
c*	Параметры:
c*		N - количество скоростей
c*		X(N) - x-координаты скоростей
c*		Y(N) - y-координаты скоростей
c*		V(N) - значения скоростей
c*		CX, CY - координаты центра вихря
c*		INDEX - количество скоростей, лежащих между центром вихря и 
c*				предполагаемым Rm
c*		Vl - сумма значений скоростей, лежащих "не правее" от скорости
c*			   с номером INDEX
c*		Ng - количество	скоростей, лежащих "правее" от скорости
c*			   с номером INDEX
c*		...
c*********************************************************************
	SUBROUTINE ApproximateProfile(N, X, Y, NK, XK, YK, TK, SC, Vt, Vn, 
     *  Fi, Lambda, CX, CY, Rm, A0, R0)
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: N, NK
        REAL,    INTENT(IN) :: X(N), Y(N), XK(NK), YK(NK), TK(NK), 
     *                         SC(NK * 3 * 2), Vt(N), Vn(N), 
     *						 Fi, Lambda, CX, CY, R0
	  REAL                :: Rm, A0

        REAL, EXTERNAL :: ObjectiveFunctionValue

	  REAL :: fx, fxopt, A0opt, RMopt

	  REAL    :: Vl, Rl, VRl, RRl, Vg, Rg, VRg, RRg
        INTEGER :: Nl, Ng

	  INTEGER i, irc1, irc2
	  REAL dr1, dr2, dRm, dRm1, dRm2, Rm1, Rm2, A01, A02
	  
	  write(*, *) 'run ApproximateProfile'

	  write(*, *) 'Rm, A0, R0 =', Rm, A0, R0

	  A0opt = A0; RMopt = Rm

	  FXopt = ObjectiveFunctionValue(N, X, Y, Vt, Vn, Fi, Lambda,  
     *                         CX, CY, A0, Rm, R0, 0)
	  write(*, *) 'FXinit = ', FXopt

        DO i = 1, N - 1 !или через каждый километр
	    call CalculateSumParameters(N, X, Y, Vt, CX, CY,   
     *    i, Vl, Rl, VRl, RRl, Nl, Vg, Rg, VRg, RRg, Ng)
C		расчет границ интервала, для которого считаются параметры A0 и Rm
	    dr1 = sqrt((X(i)-cx)**2 + (Y(i)-cy)**2)
	    dr2 = sqrt((X(i+1)-cx)**2 + (Y(i+1)-cy)**2)	
C		расчет Rm и A0 | первая система	  
	    call CalculateRmA0(1, Rm1, A01, R0, dr1, dr2, dRm1, irc1, 
     *      Vl, Rl, VRl, RRl, Nl, Vg, Rg, VRg, RRg, Ng)
	    if (iand(irc1, 2).eq.2.and.dRm1.lt.dRm) then
	      dRm = dRm1
	      Rm = Rm1
	      A0 = -A01
	    endif
C		расчет Rm и A0 | вторая система
	    call CalculateRmA0(2, Rm2, A02, R0, dr1, dr2, dRm2, irc2, 
     *      Vl, Rl, VRl, RRl, Nl, Vg, Rg, VRg, RRg, Ng)  
          if (iand(irc2, 2).eq.2.and.dRm2.lt.dRm) then
	      dRm = dRm2
	      Rm = Rm2
	      A0 = -A02
	    endif
	    if (irc1.gt.0.or.irc2.gt.0) then 
	      !write(*, *) A0, Rm, R0, i	
	      cycle
	    endif

	  FX = ObjectiveFunctionValue(N, X, Y, Vt, Vn, Fi, Lambda,  
     *                         CX, CY, A0, Rm, R0, 0)

	  if (FX < FXopt) then
	    write(*, *) ' Neww MINIMUM                     ', FX
	    write(*, *) ' A0, Rm = ', A0, Rm 
	    !CALL ObjFuncResultsToLog(N, X, Y, NK, XK, YK, TK, 
!     *	SC, Vt, Vn, Fi, Lambda, CX, CY, A0, Rm, R0)
	    FXopt = FX
	    A0opt = A0; RMopt = RMopt
	  end if
 
	  END DO

	  A0 = A0opt; Rm = RMopt

	END


c*********************************************************************
c*
c*					ObjectiveFunctionValue
c*
c*			Расчет целевой функции
c*
c*********************************************************************
	REAL FUNCTION ObjectiveFunctionValue(N, X, Y, Vt, Vn, 
     *    Fi, Lambda, CX, CY, A0, Rm, R0, Mode)
        IMPLICIT NONE
	  INTEGER, INTENT(IN) :: N, Mode
	  REAL,    INTENT(IN) :: X(N), Y(N), Vt(N), Vn(N), Fi, Lambda,
     *                         CX, CY, A0, Rm, R0
	  REAL, EXTERNAL :: VelocityTeor

	  INTEGER :: i
	  REAL    :: fx = 0., dr, xy(2), Vteor, v(1), dx, dy

	  OPEN(19, FILE='func.log', ACTION='write')
        fx = 0.0
        DO i = 1, N
	    dr = sqrt((cx - X(i))**2 + (cy - Y(i))**2)
	    IF (dr > R0) cycle
	    Vteor = VelocityTeor(dr, A0, Rm, R0)
	    FX = fx + (Vt(i) - Vteor)**2
		if (Mode == 1) Fx = fx + Vn(i)**2
	   ! xy(1)=X(i)
         ! xy(2)=Y(i)
         ! CALL cnkmdg(fi,lambda,xy,v,1,dx,dy,0,0)
          !write(19, *) dr, xy(1), xy(2), Vteor, Vt(i)!,dist(i)+x0
        END DO   
	  CLOSE(19)

	  ObjectiveFunctionValue = fx / N

	END
	  
	 
c*********************************************************************
c*
c*				ObjFuncResultsToLog
c*
c*		Вывод оптимального результата
c*
c*********************************************************************	 
      SUBROUTINE ObjFuncResultsToLog(N, X, Y, NK, XK, YK, TK, 
     *	SC, Vt, Vn, Fi, Lambda, CX, CY, A0, Rm, R0)
	  IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, NK
	  REAL,    INTENT(IN) :: X(N), Y(N), XK(NK), YK(NK), TK(NK), 
     *		Vt(N), Vn(N), SC(NK*2*3), Fi, Lambda, CX, CY, A0, Rm, R0

	  REAL, EXTERNAL :: VelocityTeor 

	  REAL :: dr, Vteor
	  INTEGER :: i

		! Вывод радиального профиля тангенциальных скоростей  
	  open(20, FILE='profileR.opt', ACTION='write')
	  open(26, FILE='profileT.opt', ACTION='write')
        do i = 1, N
	    dr = sqrt((CX - X(I))**2 + (CY - Y(i))**2)
	    if (dr.gt.R0) cycle
	    Vteor = VelocityTeor(dr, A0, Rm, R0)
          write(20, 32) dr, Vt(i)
	    write(26, 32) dr, Vteor
        enddo   
        close(20)
	  close(26)
   32   format(3F12.4)
      
C	Создание файла для визуализации оптимизированных узлов
	  CALL WriteXY(NK, XK, YK, CX, CY, Fi, Lambda, 'nodes.opt', 2, TK)
C	Создание файла для визуализации оптимизированного сплайна	  	
        CALL WriteSplineXY(NK, XK, YK, TK, SC, Fi, Lambda, 'spline.opt')
C	Создание файла для визуализации легенды (тангенциальной скорости)
        CALL WriteXY(N, X, Y, CX, CY, Fi, Lambda,'legend_vg.txt', 1, Vt)

	END



c*********************************************************************
c*
c*					ToCircularModel
c*
c*********************************************************************
	SUBROUTINE ToCircularModel(Vt, Vn, R0, NK, XK, YK, TK, SC, 
     *    CX, CY, Fi, Lambda)
!	  IMPLICIT NONE
	  


	END
	  
c*********************************************************************
c*
c*					CalculateSumParameters
c*			
c*	  	Расчет параметров целевой функции, зависящих от выбора Rm
c*			и представляющих собой суммы значений радиусов и скоростей
c*
c*	Параметры:
c*		N - количество скоростей
c*		X(N) - x-координаты скоростей
c*		Y(N) - y-координаты скоростей
c*		V(N) - значения скоростей
c*		CX, CY - координаты центра вихря
c*		INDEX - количество скоростей, лежащих между центром вихря и 
c*				предполагаемым Rm
c*		Vl - сумма значений скоростей, лежащих "не правее" от скорости
c*			   с номером INDEX
c*		Ng - количество	скоростей, лежащих "правее" от скорости
c*			   с номером INDEX
c*		...
c*********************************************************************
	SUBROUTINE CalculateSumParameters(N, X, Y, V, CX, CY, INDEX,
     *    Vl, Rl, VRl, RRl, Nl, Vg, Rg, VRg, RRg, Ng)
	  IMPLICIT NONE
	  INTEGER, INTENT(IN)  :: N, INDEX
        REAL,    INTENT(IN)  :: X(N), Y(N), V(N), CX, CY
	  REAL,    INTENT(OUT) :: Vl, Rl, VRl, RRl, Vg, Rg, VRg, RRg
        INTEGER, INTENT(OUT) :: Nl, Ng

	  INTEGER i
	  REAL dr

	  Vl = 0.; Rl = 0.; VRl = 0.; RRl = 0.
	  Vg = 0.; Rg = 0.; VRg = 0.; RRg = 0.
	  Nl = 0; Ng = 0		 

!	  write(*, *) 'CSP: N =', N, ', index =', INDEX
        do i = 1, N
	    dr = sqrt((CX - X(i))**2 + (CY - Y(i))**2)
	    if (i <= INDEX) then
	      Vl = Vl + V(i)
	      Rl = Rl + dr
		  VRl = VRl + V(i) * dr
	      RRl = RRl + dr * dr
		  Nl = Nl + 1
	    else
	      Vg = Vg + V(i)
	      Rg = Rg + dr
		  VRg = VRg + V(i) * dr
		  RRg = RRg + dr * dr	
		  Ng = Ng + 1
	    endif
	  enddo

	END



c*********************************************************************
c*
c*					CalculateRmA0
c*			
c*	  	Расчет параметров Rm и A0 целевой функции на интервале
c*
c*	Параметры:
c*		INDEX - (=1|2) первую или вторую систему корней
c*		Rm, A0, R0  - см. в заголовке п/п FUN
c*		DR1, DR2 - правая и левая границы интервала соответственно
c*		Vl,...,Ng - см. в заголовке п/п	CalculateSumParameters
c*		IRC = ab(биты) - код возврата
c*		  = -1 - ошибка значения переменной INDEX
c*		  a=1 - 0 <= Rm <= R0
c*		  b=1 - Rm лежит в интервале, для которого и рассчитывался: 
c*				DR1 <= Rm <=DR2
c*********************************************************************
	SUBROUTINE CalculateRmA0(INDEX, Rm, A0, R0, DR1, DR2, DRm, IRC, 
     *    Vl, Rl, VRl, RRl, Nl, Vg, Rg, VRg, RRg, Ng)
	  IMPLICIT NONE
        INTEGER, INTENT( IN) :: INDEX, Nl, Ng
	  REAL,    INTENT(OUT) :: Rm, A0, DRm
	  REAL,    INTENT( IN) :: R0, DR1, DR2,  
     *                          Vl, Rl, VRl, RRl, Vg, Rg, VRg, RRg
	  INTEGER, INTENT(OUT) :: IRC

	  REAL    :: factor
	  INTEGER :: a = 0, b = 0

	  if (INDEX.eq.1) factor = 1.0
	  if (INDEX.eq.2) factor = -1.0
	  if (INDEX.ne.1.and.INDEX.ne.2) then
	    IRC = -1
	    return
	  endif

	  DRm = 100500.

	  Rm = (R0 *(VRl*(RRg-2*R0*Rg+R0*R0*Ng)-2*RRl*(VRg-R0*Vg)) + 
     *        factor *
	*        R0*sqrt((2*RRl*(VRg-R0*Vg)-VRl*(RRg-2*R0*Rg+R0*R0*Ng))**2
     *	    + 4*RRl*(VRg-R0*Vg)*
     *		 (VRl*(RRg-2*R0*Rg+R0*R0*Ng)-RRl*(VRg-R0*Vg)))) / 
     *		 (2*VRl*(RRg-2*R0*Rg+R0*R0*Ng)-RRl*(VRg-R0*Vg))
C		проверка на корректность параметра Rm
	  if (Rm.gt.0.and.Rm.lt.R0) then
C		расчет A0
          a = 1
		A0 = (Rm*(VRg-R0*Vg)-VRl*(R0-Rm)) / 
     *       (Rm*Rm*(RRg-2*R0*Rg+R0*R0*Ng)+RRl*(R0-Rm)**2) *
     *	   (R0*Rm*(R0-Rm)) / 2
	    DRm = min(abs(Rm-DR1), abs(Rm-DR2))
!          write(*, '(A, I1)') ' CORRECT RM ', INDEX
	    if (Rm.ge.DR1.and.Rm.le.DR2) then
		  b = 1
!		  write(*, '(A, I1)') ' TRUE RM', INDEX
	    end if
!	    write(*, *) 'A0 =', A0
	  endif
	  IRC = 2 * a + b
	END



c*********************************************************************
c*
c*						WriteXY
c*
c*		Перевод координат в географическую систему координат и 
c*			вывод их в файл
c*
c*	Параметры:
c*		N - количество точек	
c*		X(N) - x-координаты точек
c*		Y(N) - y-координаты точек
c*		cx, cy - координаты центра вихря
c*		FI, LAMBDA - начало географической системы координат (д. и ш.)
c*		FILENAME - имя файла вывода
c*		IFLG = (ab) - флаги вывода 
c*		   a=1 - выводить координаты cx cy
c*		   b=1 - выводить дополнительный столбец LEGEND(N)
c*
c*********************************************************************
      SUBROUTINE WriteXY(N, X, Y, CX, CY, FI, LAMBDA, FILENAME, IFLG, 
     *    LEGEND)
	  INTEGER,   INTENT(IN) :: N, IFLG
	  REAL,      INTENT(IN) :: X(N), Y(N), FI, LAMBDA, LEGEND(N)
	  CHARACTER, INTENT(IN) :: FILENAME*(*)

	  REAL xxyy(2), v(1)

	  open(14, FILE=FILENAME, ERR=9, ACTION='write')

	  IF (iand(IFLG, 2).eq.2) THEN
	    xxyy(1) = cx
	    xxyy(2) = cy
	    call cnkmdg(FI, LAMBDA, xxyy, v, 1, dx, dy, 0, 0)
	    write(14, 5) xxyy
	  END IF

	  DO i = 1, N
	    xxyy(1) = X(i)
	    xxyy(2) = Y(i)
          call cnkmdg(FI,LAMBDA,xxyy,v,1,dx,dy,0,0)   
	    IF (iand(IFLG, 1).ne.1) THEN
		  write(14, 5) xxyy
		ELSE 
	      write(14, 5) xxyy, LEGEND(i)
	    END IF
	  END DO
  5   format(3F12.6)
	  close(14, ERR=9)
!	  write(*, *) 'WXY ', FILENAME, ' OK'
	  RETURN

    9   write(*, *) 'WriteXY: ERROR file opening!'
 !       write(*, *) 'WriteXY: n =', N
  !      write(*, *) 'WriteXY: file', FILENAME
        RETURN
   10   write(*, *) 'WriteXY: ERROR file closing!', FILENAME
	END	



c*********************************************************************
c*
c*						WriteSplineXY
c*
c*		Расчет промежуточных точек (между узлами) сплайна
c*			в географической системе координат и вывод их в файл
c*
c*	Параметры:
c*		NK - количество узлов сплайна	
c*		XK(NK) - x-координаты узлов сплайна
c*		YK(NK) - y-координаты узлов сплайна
c*		TK(NK) - t-координаты узлов сплайна
c*		C(NK,3,2) - параметры сплйна
c*		FI, LAMBDA - начало географической системы координат (д. и ш.)
c*		XC, YC - координаты центра вихря
c*		FILENAME - имя файла вывода
c*									 
c*********************************************************************
	SUBROUTINE WriteSplineXY(NK, XK, YK, TK, C, FI, LAMBDA,	FILENAME)
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: NK
	  REAL, INTENT(IN) :: XK(NK), YK(NK), TK(NK), c((nk-1)*3*2),
     *    FI, LAMBDA
	  CHARACTER(LEN = *), INTENT(IN) :: FILENAME

	  INTEGER :: N, IFLGW = 0, i, ier, lenc
	  PARAMETER (N = 200)
	  REAL x(N+1), y(N+1), legend(1), ts
	  CHARACTER :: sp_coef_filename*50
	  CHARACTER(LEN = *), PARAMETER :: scoef = "coef_"

	  sp_coef_filename = scoef // FILENAME
	  write(*,*) 'new sp_coef filename is ', sp_coef_filename
	  open(17, FILE=sp_coef_filename, ACTION='write')
	  write(17, *) (C(i), i=1,(nk-1)*6)
	  close(17)	
        lenc=(nk-1)*3
        write(*, *) 'WSXY: lenc =', lenc
        DO i = 0, N 
	    ts = TK(1) + (TK(NK) - TK(1)) * i / N 
		call  icsevu(tk,xk,nk,c,        nk-1,ts,x(i+1),1,ier)
	    call  icsevu(tk,yk,nk,c(lenc+1),nk-1,ts,y(i+1),1,ier)
	  END DO

	  CALL WriteXY(N+1, x, y, 0., 0., FI, LAMBDA, FILENAME, IFLGW, 
     *    legend)
	END



c*********************************************************************
c*
c*						CheckOnEqualityOptimizationBlock
c*
c*		Проверка на изменение переменных оптимизируемого блока
c*
c*	Параметры:
c*		N - размер блока оптимизируемых переменных
c*		NK - количество оптимизируемых узлов (незамкнутых)
c*		X(N) - блок оптимизируемых переменных
c*		X0, Y0 - координаты центра вихря
c*		XK, YK - x,y-координаты узлов сплайна
c*		R0 - см. описание в заголовка п/п FUN
c*		IRC	- код возврата
c*			0 - нормальное завершение
c*			>0 - найдено различие в значении переменной(ых),
c*				IRC указывает на первое расхождение в блоке X
c*
c*********************************************************************
	SUBROUTINE CheckOnEqualityOptimizationBlock(N, NK, X, X0, Y0, 
     *    XK, YK, R0, IRC)
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: N, NK
	  REAL, INTENT(IN) :: X(N), X0, Y0, XK(NK+1), YK(NK+1), R0
	  INTEGER, INTENT(OUT) :: IRC

	  INTEGER i

	  IRC = 0

	  IF (X0.ne.X(1)) IRC = 1
	  IF (Y0.ne.X(2)) IRC = 2
	  DO i = 1, nk
	    IF (X(2+i).ne.XK(i)) IRC = 2+i
		IF (X(2+i+nk).ne.YK(i)) IRC = 2+i+nk
        END DO
	
        IF (X(3).ne.XK(NK+1)) IRC = 3
	  IF (X(3+nk).ne.YK(NK+1)) IRC = 3+nk

	  IF (R0.ne.X(2+nk*2+3)) IRC = 2+nk*2+3

	  IF (IRC > 0) THEN
	    write(*, ('(2F12.6)')) X(1), X0, X(2), Y0
	    DO i = 1, NK
	      write(*, *) X(2+i), XK(i)
	    END DO
	    write(*, *) X(3), XK(NK+1)
	    DO i = 1, NK
	      write(*, *) X(2+NK+i), YK	(i)
	    END DO
	    write(*, *) X(3+NK), YK(NK+1)
	    write(*, *)	X(5+NK*2), R0
!	    ... // TODO
	  END IF
	END



c*********************************************************************
c*
c*						CheckRm
c*
c*		Проверка корректности формул для Rm и A0
c*
c*********************************************************************
      SUBROUTINE CheckRm()
        REAL x(4), y(4), vg(4), cx, cy
	  INTEGER idx, N
	  DATA x /0.25, 2.0, 10, 12.2/, y /4 * 0./, vg /1, 8, 20, 5/, 
     *    R0/12.933333/, cx, cy /2 * 0.0/
	  DATA N /4/, idx/2/

	  write(*, *) 'Проверка формул расчета A0 и Rm'

	  CALL CalculateSumParameters(N, x, y, vg, cx, cy,   
     *    idx, Vl, Rl, VRl, RRl, Nl, sVg, Rg, VRg, RRg, Ng)

        write(*, *) ' SumVl =', Vl
	  write(*, *) ' SumRl =', Rl
	  write(*, *) ' SumVRl =', VRl
	  write(*, *) ' SumRRl =', RRl
        write(*, *) ' SumVg =', sVg
	  write(*, *) ' SumRg =', Rg
	  write(*, *) ' SumVRg =', VRg
	  write(*, *) ' SumRRg =', RRg
	  write(*, *) ' NumLess =', Nl
	  write(*, *) ' NumGreater =', Ng
	  write(*, *) ' R0 =', R0

	  rm1 = (R0 *(VRl*(RRg-2*R0*Rg+R0*R0*Ng)-2*RRl*(VRg-R0*sVg)) + 
	*        R0*sqrt((2*RRl*(VRg-R0*sVg)-VRl*(RRg-2*R0*Rg+R0*R0*Ng))**2
     *	    + 4*RRl*(VRg-R0*sVg)*
     *		 (VRl*(RRg-2*R0*Rg+R0*R0*Ng)-RRl*(VRg-R0*sVg)))) / 
     *		 (2*(VRl*(RRg-2*R0*Rg+R0*R0*Ng)-RRl*(VRg-R0*sVg)))

	  rm2 = (R0 *(VRl*(RRg-2*R0*Rg+R0*R0*Ng)-2*RRl*(VRg-R0*sVg)) - 
     *       R0*sqrt((2*RRl*(VRg-R0*sVg)-VRl*(RRg-2*R0*Rg+R0*R0*Ng))**2
     *	    + 4*RRl*(VRg-R0*sVg)*
     *		 (VRl*(RRg-2*R0*Rg+R0*R0*Ng)-RRl*(VRg-R0*sVg)))) / 
     *		 (2*(VRl*(RRg-2*R0*Rg+R0*R0*Ng)-RRl*(VRg-R0*sVg)))
	  write(*, *) 'Rm =', rm1, rm2
	  
	  A01 = (rm1*(VRg-R0*sVg)-VRl*(R0-rm1)) / 
     *       (rm1*rm1*(RRg-2*R0*Rg+R0*R0*Ng)+RRl*(R0-rm1)**2) * 
     *	   (R0*rm1*(R0-rm1)) / 2
	  A02 = (rm2*(VRg-R0*sVg)-VRl*(R0-rm2)) / 
     *       (rm2*rm2*(RRg-2*R0*Rg+R0*R0*Ng)+RRl*(R0-rm2)**2) * 
     *	   (R0*rm2*(R0-rm2)) / 2
	  write(*, *) 'A0 =', A01, A02
	END