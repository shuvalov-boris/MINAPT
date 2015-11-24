$notstrict
c**********************************************************************
c*
c*            S P R E P R
c*
c*       Сплайн-представление контуров 
c*
c* Примечание. Если t-координата первого узла после ее расчета 
c*	программой tprep меньше 0, то перед расчетом сплайна в контур
c*	добавляется точка, лежащая первой и имеющая координаты первого
c*	узла сплайна
c*
c**********************************************************************
C						ОПРЕДЕЛЕНИЕ КОНСТАНТ
      integer			 VMAX
	character	NAMARC*12, RESFILE*12
	parameter (
     *  NMAX   =    200,	! максимальное количество точек контура
     *  NKMAX  =     28,	! максимальное количество узлов сплайна
     *  VMAX   =   1000,	! максимальное количество маркеров скоростей
     *  LENMAX = 100000,	! размер массива idat
     *  NBMAX  =    100,  ! максимальное количество блоков
     *  NAMARC  = 'geodyn.dat ',	! имя файла со входными данными
     *  RESFILE = 'gd_res.dat '	! имя файла с результатами расчетов
     *)
C						СЧИТЫВАЕМЫЕ ДАННЫЕ
	character		  where*38	 ! путь к файлам с геоданными
      real					fi,  ! долгота	  / начало географической
     *					lambda,	 ! широта	  \     системы координат
	*						x0,	 ! "левая" (внутр) граница относительно кон-ра
     *						xr,	 ! ограничение на выбор скоростей 				
	*					   eps,	 ! критерий для PGM сглаживания выбросов
	*					 epsil	 ! точность расчета координат в 10km
	integer*2		 idate0(7),  ! дата начала отсчета времени
     1					 nknot,	 ! количество узлов для сплайна
C										Параметры расчета:
     2					itypcn,	 ! :0 - координаты в 10км, 1 - в градусах
     3					iflgsp,	 ! управляющий параметр для расчета сплайна 
     4					iflgtp,	 ! управляющий параметр для сглаживания кон-ра                   
	5					iflgvr,	 ! управляющий параметр для разложения скор-ей
	6					iflgtd	 ! управляющий параметр ПП TRDIST
C						ЛОКАЛЬНЫЕ ПЕРЕМЕННЫЕ
      real		 dat(LENMAX/2)
	integer*2	  idat(LENMAX),
     *						nb,	 ! текущее количество блоков 
     *					nbspap,  ! номер блока с узлами	и сплайном
     *					nbcont	 ! номер блока с контуром
	integer*4 lentyp(2,NBMAX+1)	 ! начало блока в массиве idat и тип блока
      equivalence (dat,idat)
	real t(NMAX),wk(NMAX*(NKMAX+6)*4),piw(4,2),xy[allocatable](:,:)
	*  , xxyy(2), vg(VMAX), dist(VMAX), vr(VMAX)
	*  , x0sum, x0cur, x0max, x0min
	integer*2 iwk(202)
	integer*4 ibgwk4,nbgknt,nbgcon,nbc,nbk
	character nm
	data maxrec/1000/
C			ЧТЕНИЕ ВХОДНЫХ ДАННЫХ
      open(unit=11, file=NAMARC, action='read', status='old')
      read(11, *, ERR=999) where, fi, lambda, idate0
	
	open(unit=33, file=RESFILE, action='write')

	mr1=maxrec+1
	nv=-1
		
   10 read(11, *, END=999) nknot, itypcn, iflgtp, iflgsp, iflgvr, iflgtd
     1	, x0, xr, eps, epsil
      nb = 0
	nbspap = 0
c
c*		расчет узлов
c
	lnam=LEN_TRIM(where)
	call  inpsfd(fi,lambda,idate0,itypcn,where(:lnam),
     *  nb,lentyp,dat,idat,
     *  NBMAX,LENMAX,piw,MAXREC,iwk,wk,wk(mr1*2+1),wk(mr1*4+1),irc)
	write(6,17) irc
   17 format(3x,'программа inpsfd завершилась с irc=',i4)
	ircis=irc

C			Вывод координат точек контура в файл для визуадизации
      open(unit=14, file='contour.txt', err=999, action='write')
	cibg=4
	do i=1,idat(lentyp(1,1)+3)
	  xxyy(1)=dat(cibg+i)
	  xxyy(2)=dat(cibg+idat(4)+i)
        call cnkmdg(fi,lambda,xxyy,v,1,dx,dy,0,0)   
	  write(14, 15) xxyy   
	enddo
   15 format(2F11.3)
	close(14,err=999)
C			Поиск блока с контуром
	do i = 1, nb
	  if (lentyp(2, i).eq.1) nbcont = i
      enddo

      nbs=nbcont

      call  datwrk(dat,idat,lentyp,lenmax,nb,nbmax,nbs,x,y,t,ibgwk4,
     *  nknot,8,iflg,irc)	
      write(6, *) ' Подпрограмма DATWRK завершилась с irc =', irc
      if(irc.ne.0) call  mescon('  Нет места в архиве',20,nm,0,intg,
     *  0,rl,0)
C			Вывод координат точек узлов в файл для визуадизации
      open(unit=16, file='nodes.txt', err=999, action='write')
	cibg=lentyp(1,nbs)
	nknot=idat(cibg+3)
	cibg=(cibg+6)/2
	do i=1,nknot
	  xxyy(1)=dat(cibg+i)
	  xxyy(2)=dat(cibg+nknot+i)
        call cnkmdg(fi,lambda,xxyy,v,1,dx,dy,0,0)      
 	  write(16, 32) xxyy(1), xxyy(2)
	enddo
	close(16,err=999)
      
	ircdw=irc
	if(irc.ne.0) goto 666
C	write(6, *) ' Удаляем директорию'
C	if(nbspap.eq.0)call  kmon('del '//where(:lnam)//'\*.*'//null_char)
	nbspap=nbs
	nv=0

      write (*, *) 'SPREPR: call MINAPT'
	call minapt(lentyp, idat, dat, 1, nbspap)
c
c	call  prdat(fi,lambda,idate0,nb,lentyp,dat,idat,nbs,-1,0)
c
c*		сплайн аппроксимация контура
c*	  pасчет констант
20	nbgknt=lentyp(1,nbspap)
	nknot=idat(nbgknt+3)
	nbgknt=nbgknt+6
	if(nv.eq.-1) nv=idat(nbgknt-1)
	nbgcon=lentyp(1,nbcont)
	ncon=idat(nbgcon+3)
	nf=idat(nbgcon+4)
	if(nf-nf/2*2.eq.1) nf=nf+1
	nbgcon=nbgcon+6+nf
	nbc=nbgcon/2+1
	nbk=nbgknt/2+1
c*		Расчет t-кооpдиинаты
	call  tprep(dat(nbc),t(2),ncon,dat(nbk),dat(nbk+nknot),dat(nbk
     *  +nknot+nknot),nknot,nbeg,nend,nusl,nusr,ndoubl,eps,wk,iwk,
     *  iflgtp,irc)
	write(6,120) irc,nbeg,nend,nusl,nusr,ndoubl,eps
120	format(3x,' программа tprep. irc=',i3,/,3x,'nbeg,nend,nusl,nusr,
     *ndoubl,eps=',5i5,f9.3)
	print*,'tk=',(dat(nbk+nknot*2+jj),jj=0,nknot-1)
	irctp=irc
c 	write(6,130) (wk(k),k=1,nknot)
c	write(6,135) (iwk(k),k=1,nknot)
c130	format(3x,'wk=',8f9.3)
c135	format(3x,'iwk=',14i5)
c	if(nb.ne.-1) return
	nbeg=nbeg-1
	nconn=nend-nbeg
c
c*		Сплайн-аппpоксимация
c
c*	 Выделяем память и заводим новый массив кооpдинат контуpа
	ib=1
	if(dat(nbk+nknot+nknot).lt.0.) ib=2
	ie=0
	if(dat(nbk+nknot*3-1).gt.t(nend+1)) ie=1
	nconn=nconn+ib-1+ie
	allocate(xy(nconn,2),STAT=irc)
C	if(irc.ne.0) call  mescon('Не хватает опеpативной памяти',29,
C     *  nm,0,intg,0,rl,0)
	if(irc.ne.0) goto 999
c*		Если кpайние узелы вне контуpа - добавляем их в контуp
	if(ib.eq.2) xy(1,1)=dat(nbk)
	if(ib.eq.2) xy(1,2)=dat(nbk+nknot)
	if(ib.eq.2) t(1)=dat(nbk+nknot+nknot)
	if(ie.eq.0) go to 138
c*		Последний узел необходимо внести в контуp
	xy(nconn,1)=dat(nbk+nknot-1)
	xy(nconn,2)=dat(nbk+nknot+nknot-1)
	t(nconn+2-ib)=dat(nbk+nknot*3-1)
c
138	do 140 i=0,nconn-ib-ie
	xy(i+ib,1)=dat(nbc+nbeg+i)
140	xy(i+ib,2)=dat(nbc+nbeg+ncon+i)
C			Вывод координат точек сглаженного контура в файл для визуадизации 
      open(unit=15, file='sp_contour.txt', err=999, action='write')
	cibg=4
	do i=1,nconn
	  xxyy(1)=xy(i,1)
	  xxyy(2)=xy(i,2)
        call cnkmdg(fi,lambda,xxyy,v,1,dx,dy,0,0)      
 	  write(15, 32) xxyy(1), xxyy(2)
	enddo
	close(15, err=999)
c
	print *,'ncon,nconn,nbeg',ncon,nconn,nbeg
	call  spapp(xy,t(3-ib+nbeg),nconn,dat(nbk),dat(nbk+nknot)
     *,dat(nbk+nknot+nknot),nknot,dat(nbk+nknot*3),nknot-1,err,wk,
     *iflgsp,irc)
	write(6,150) irc,err
150	format(3x,'программа spapp. irc,err=',i4,f9.3)
	if(irc.eq.0) nv=nv+1
	if(irc.eq.0) idat(nbgknt-1)=nv
	ircsp=irc
	sperr=err
C		Расчет x,y координат сплайна и их вывод в файл для визуализации
      open(UNIT=17, FILE='spline.txt', ERR=999, ACTION='write')
	do i = 1, nconn
	  call trdist(fi,lambda, xy(i,1), xy(i,2), 
     *    dat(nbk), dat(nbk+nknot), dat(nbk+nknot+nknot), nknot,
     *    dat(nbk+nknot*3),nknot-1, x, y, t,
     *    tga, vrk, distance, epsil, iflgtd)
	  xxyy(1)=x
	  xxyy(2)=y
        call cnkmdg(fi,lambda,xxyy,v,1,dx,dy,0,0)      
 	  write(17, 32) xxyy(1), xxyy(2)
	enddo
   32 format(2F12.4)
      close(17, err=999)
C		Создание файла для визуализации оптимизированных узлов
      open(unit=18, file='sp_nodes.txt', err=999, action='write')
	cibg=lentyp(1,nbs)
	nknot=idat(cibg+3)
	cibg=(cibg+6)/2
	do i=1,nknot
	  xxyy(1)=dat(cibg+i)
	  xxyy(2)=dat(cibg+nknot+i)
        call cnkmdg(fi,lambda,xxyy,v,1,dx,dy,0,0)      
 	  write(18, 32) xxyy(1), xxyy(2)
	enddo
	close(18,err=999)
c*		Освобождаем память
	deallocate(xy,STAT=irc)

C			Расчет радиуса вихря (x0)
      x0sum = 0.0
	x0min = 100000.0 !
	x0max = -x0min
      do i = 0, ncon - 1
	  x0cur = sqrt((dat(nbc + i))**2 + (dat(nbc + ncon + i))**2)
	  x0sum = x0sum + x0cur
	  x0max = max(x0max, x0cur) 
	  x0min = min(x0min, x0cur)
	enddo
	x0 = x0sum / ncon
	xr = x0 * 2.5

C*		РАСЧЕТ ТАНГЕНЦИАЛЬНЫХ КОМПОНЕНТ СКОРОСТИ
      call vprep(fi, lambda, dat, idat, lentyp, 2, nbspap, x0, xr, 
     *  dist,vg,vag,tga,vr,nv,iflgvr,epsil,irc)
	write(6,160) irc
160	format(3x,'программа vprep. irc=',i4)
    
c*		Вывод pезультов pасчета
  666 write(33, 13) ircis, ircdw, irctp, ircsp, irc , 
     *  itypcn, iflgtp, iflgsp, iflgvr, iflgtd, nknot, x0min, x0, x0max
     *  , xr, sperr, nv
   13 format(5I2, ' ', 5I1, I3, 5F8.3, I5)
      goto 10

  999 close(33)
      close(11)
	close(22) 
      end