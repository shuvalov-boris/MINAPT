c*********************************************************************
c*
c*
c*              M O N M I N
c*
c*       Подготовка данных для минимизации и минимизация функции
c*          namfun методом nammet с параметрами из блока nbminf
c*
c*
c*********************************************************************
c	EXTERNAL OUT2,minl5,min03,min05,ralg
	common/gddat/ icninf,istinf,idat,lentyp
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

	integer*4 lentyp(2,NBMAX+1)	 ! начало блока в массиве idat и тип блока
c$include:'geodyn.for'      
	REAL A[allocatable](:),B[allocatable](:),D[allocatable](:),
     *DM[allocatable](:),DS[allocatable](:),X[allocatable](:),
     *FX[allocatable](:),ddmds(3),rl(2),xold(220),prmcv(5),xbst(220)
	INTEGER*2 idate(7),mn(12)
	integer*4 ibg,ibg2
	data lencdt/53/,pi/3.1415926/
	data mn/31,28,31,30,31,30,31,31,30,31,30,31/
c*              Чтение архива данных
	call  rdwrta(namarc,icninf,istinf,lentyp,idat,0,irc)
	if(irc.eq.3) call  mescon('  Ошибки при чтении из архива',29,
     *  nm,0,intg,0,rl,0)
	if(irc.eq.2) call  mescon('  Нет архива в текущей директории',33,
     *  nm,0,intg,0,rl,0)
	if(irc.ne.0) return
c
c*          Пересчет текущей даты в часы
c
	call date(idate(1),idate(2),idate(3),idate(4))
	call time(idate(4),idate(5),idate(6),idate(7))
	id0=0
c*          расчет номера дня в году

	me=idate(2)-1
	if(me.eq.0) go to 23
	  do 22 i=1,me
22      id0=id0+mn(i)
23    id0=id0+idate(1)
	idy0=0
	if(idate(3)-idate(3)/4*4.eq.0) idy0=1
	if(me.ge.2) id0=id0+idy0
c*          пересчет времени дня в часы
	h0=idate0(4)+idate0(5)/60.+idate0(6)/3600.+idate0(7)/360000.
	time0=24.*id0+h0
c
c*              Заведение массивов для минимизации
c
c*          Вытягивание из минимизационного блока информации
	ibg=lentyp(1,nbminf)
	n=idat(ibg+3)
	ibg=ibg+6
	isizea=idat(ibg+30)
	isizeb=idat(ibg+31)
	kodgrd=idat(ibg+32)
	jjxfj=idat(ibg+20)
	isizex=n*jjxfj
	isized=n
	if(kodgrd.eq.2) isized=n+n*(n+1)/2
c
	allocate (a(isizea),b(isizeb),d(isized),dm(n),ds(n),x(isizex),
     *  fx(jjxfj),STAT=irc)
	if(irc.ne.0) call  mescon('  Не хватает оперативной памяти',31,
     *nm,0,intg,0,rl,0)
	if(irc.ne.0) return
c
c*          Запись даты в файл nbminf
	ibg2=ibg/2
	dat(ibg2-2)=time0
c
c*          Присвоение начальных значений массивам
c
c*                D,DM,DS
	ddmds(1)=dat(ibg2+18)
	ddmds(2)=dat(ibg2+19)
	ddmds(3)=dat(ibg2+20)
	do 30 i=1,n
	d(i)=ddmds(1)
	dm(i)=ddmds(2)
30    ds(i)=ddmds(3)
c*                X,FX
	fx(1)=dat(ibg2+lencdt)
	do 40 i=1,n
40    x(i)=dat(ibg2+53+i)
c*          Минимизируем
c       CALL minfun(min05,idat(ibg+20),n,X,FX,D,DM,DS,dat(ibg2+1),
c     *MINL5,OUT2,A,B,IRC,idat(ibg+lencdt*2+2*n),13,1)
c************           Спец. pасчеты: nbminf<0 ************
c     nbminf=-nbminf
c        fxbest=1.e30
       call  fun(x,fx)
c        if(i.ne.375) go to 65
       if(i.ne.375) go to 80
c       fxbest=fx(1)
	   do i=1,n
	   xold(i)=x(i)
	   xbst(i)=x(i)
	   end do
       dxy=0.5
       dfi=pi/18.
       dx0=-2.
       x0=9.
       y0=6.
       do 60 i=0,8
	dx=dx0+dxy*i
	dy0=-2.
       do 60 j=0,8 
	dy=dy0+dxy*j
	do 60 k=0,8
	fip=(k-4)*dfi
	do 60 l=0,8
	akz=1.2**(l-4)     
	    do jj=1,n
	    x(jj)=xold(jj)
	    end do
c*              растяжение
	prmcv(1)=x0
	prmcv(2)=y0
	prmcv(3)=0.
	prmcv(4)=1.
	prmcv(5)=akz
	call convxyt(x,x(n/2+1),t,n/2,prmcv,3)
c
c*              поворот
	prmcv(1)=x0
	prmcv(2)=y0
	prmcv(3)=fip
	call convxyt(x,x(n/2+1),t,n/2,prmcv,2)
c
c*              сдвиг
	prmcv(1)=dx
	prmcv(2)=dy
	call convxyt(x,x(n/2+1),t,n/2,prmcv,1)
c
c        call  mescon('введите приращения по x',23,nm,0,int,0,rl,2)
c        if(rl(1).eq.-111.) go to 60
      call  fun(x,fx)
	print*,'dx,dy,fip,akz,fx,fxbest,dxb,dyb,x=',dx,dy,fip,akz,fx(1),
     *fxbest,dxb,dyb,(x(kk),kk=1,n)
	if(fx(1).ge.fxbest) go to 60
	fxbest=fx(1)
	dxb=dx
	dyb=dy
	fib=fip
	akzb=akz
	    do jj=1,n
	    xbst(jj)=x(jj)
	    end do
60      continue
	do i=1,n
	x(i)=xbst(i)
	end do
	print*,'dx,dy,fip,akz best=',dxb,dyb,fib,akzb
       call fun(x,fx)
c**********************************************************
c*          Сохраняем результат
65      dat(ibg2+lencdt)=fx(1)
	do 70 i=1,n
70    dat(ibg2+lencdt+i)=x(i)
80      call  rdwrta(namarc,icninf,istinf,lentyp,idat,1,irc)
	if(irc.eq.3) call  mescon('  Ошибки при записи в архива',28,
     *  nm,0,intg,0,rl,0)
c*          Освобождаем память
	deallocate (a,b,d,dm,ds,x,fx,STAT=irc)
c
	RETURN
	END
	subroutine fun(x,fx)
c*********************************************************************
c*
c*
c*          Расчет целевой функции
c*
c*
c* Параметры:     x - массив оптимизируемых переменных (вход)
c*            fx - значение целевой функции (выход)
c*
c*
c*********************************************************************
	common/gddat/ icninf,istinf,idat,lentyp
$include:'geodyn.for'
	real optprm(53),ddmds(3),parm(10),x(1),wk[allocatable](:),prm(10)
	real xst(200),yst(200),hm(200),hc(200),ys(200),xold(220)
	integer typesf(2),nblopt(50),kod(14),hddat(6)
	character*8 names(3)
	integer*4ibg,ibg2,ibgo,ibgopt,ibgi,ibgi2,nbg(2,50),ib,ibgc,ibgwk4
	equivalence (optprm,parm),(optprm(11),kod),(optprm(18),
     *ddmds),(optprm(21),names),(optprm(27),typesf),(optprm(28),nblopt),
     *(optprm(53),fxbest),(hddat,time),(hddat(3),iflg),(hddat(4),n),
     *(hddat(5),nbopt),(hddat(6),iblk)
	data lencdt/53/,eps/.1/,ifst/0/,vmin/1./,
     *nbc/2/,np/5/,g/9.7818/,omega/7.292e-5/,iflout/0/,pi/3.1415926/
c!          nbginf - номер блока с гидрологической информацией
	data nbginf/15/
c*          Удаляем последний блок, если вызов не первый
	ibgwk4=0
	nbindl=nb
	if(ifst.eq.1) call datwrk(dat,
     *  idat,lentyp,lenmax,nb,nbmax,nbindl,x,y,t,ibgwk4,0,2,0,irc)
	if(nbminf.lt.0) iflout=1
	if(nbminf.lt.0) nbminf=-nbminf
c           Проверяем: первый ли раз вызывается программа?
	if(ifst.eq.1) go to 100
	ifst=1
c
c*          Инициализация программы расчета целевой функции
c
c*          счиываем справку оптимизационного блока
	ibg=lentyp(1,nbminf)
	do 10 i=1,6
10    hddat(i)=idat(ibg+i-1)
c**      сохраняем начальные значения X
c       do 15 i=1,n
c15      xold(i)=x(i)
c
c*          считываем параметры оптимизации
	ibg2=ibg/2+3
	do 20 i=1,lencdt
20    optprm(i)=dat(ibg2+i)
c*          константы оптимизируемого сплайна
	ibgo=lentyp(1,nbopt)
	ibgopt=ibgo/2+4
	nknot=idat(ibgo+3)
	ibgc=ibgopt+nknot*3
	lenc=(nknot-1)*3
c*          константы блока с гидрологической информацией
	nbgi=0
	do 22 i=1,nb
	if(lentyp(2,i).ne.3.or.i.ne.nbginf) go to 22
	nbgi=i
22    continue
	fdg=2.*omega*sin(lambda*pi/180.)/g*100.
	if(nbgi.eq.0) go to 28
c*          блок с гидрологической информацией существует
	ibgi=lentyp(1,nbgi)
	nst=idat(ibgi+3)
	ibgi2=ibgi/2+3
	do 24 i=1,nst
	xst(i)=dat(ibgi2+i)
	yst(i)=dat(ibgi2+nst+i)
24    hm(i)=dat(ibgi2+nst*3+i)
c
c*          Заполнение массива начал блоков участков оптимизации
c*                и чисел точек в блоке
28    nbo=nblopt(1)
	do 30 i=1,nbo
	nf=0
	if(lentyp(2,nblopt(i+1)).eq.1) nf=idat(lentyp(1,nblopt(i+1))+4)
	if(nf-nf/2*2.eq.1) nf=nf+1
	nbg(1,i)=(lentyp(1,nblopt(i+1))+nf)/2+4
30    nbg(2,i)=idat(nbg(1,i)*2-4-nf)
c*          опpеделение максимального числа вектоpов скоpости 
c*    в одном блоке из списка заказанных
	maxv=0
	do 40 i=1,nbo
	if(lentyp(2,nblopt(i+1)).ne.0) go to 40
	if(maxv.ge.nbg(2,i)) go to 40
	maxv=nbg(2,i)
40    continue
	npntwk=maxv
	if(nknot.gt.maxv) npntwk=nknot
c*          заведение pабочего массива
	allocate(wk(npntwk*6),STAT=irc)
	if(irc.ne.0) call  mescon('  Не хватает оперативной памяти',31,
     *nm,0,intg,0,rl,0)
	if(irc.ne.0) return
c
c*          Расчет целевой функции
c
100   continue
c       print *,'begin'
c*          Обновляем параметры блока nbopt
	do 110 i=1,n
	dat(ibgopt+i-1)=x(i)
110     continue
c*          Замыкание сплайна
c        dat(ibgopt+n/2-1)=x(1)
c        dat(ibgopt+n-1)=x(n/2+1)
c        
c*          Расчет сплайн-коэффициентов
	if((dat(ibgopt).eq.dat(ibgopt+nknot-1)).and.(dat(ibgopt+nknot).
     *   eq.dat(ibgopt+nknot*2-1))) go to 112
c*        сплайн не замкнутый
	call icsccu(dat(ibgopt+nknot+nknot),dat(ibgopt),nknot,
     *  dat(ibgopt+nknot*3),nknot-1,irc)
	if(irc.ne.0) print *,'  Ошибка сплайн-аппроксимации irc=',irc
	call icsccu(dat(ibgopt+nknot+nknot),dat(ibgopt+nknot),nknot,
     *  dat(ibgopt+nknot*3+lenc),nknot-1,irc)
	if(irc.ne.0) print *,'  Ошибка сплайн-аппроксимации irc=',irc
	go to 114
c*         сплайн замкнутый
112     continue
	call icspln(dat(ibgopt+nknot+nknot),dat(ibgopt),nknot,
     *  dat(ibgopt+nknot*3),nknot-1,wk,irc)
	if(irc.ne.0) print *,'  Ошибка сплайн-аппроксимации irc=',irc
	call icspln(dat(ibgopt+nknot+nknot),dat(ibgopt+nknot),nknot,
     *  dat(ibgopt+nknot*3+lenc),nknot-1,wk,irc)
	if(irc.ne.0) print *,'  Ошибка сплайн-аппроксимации irc=',irc
c
114     fx=0.                
	sdf1=0.
c       print *,'spline'
	sdf=0.
c*          пpепаpиpование скоpостей с запоминанием
	nbin=0
	x0=9.
	r=x0+8.
	do 120 i=1,nbo
c      print*,'nbin,ityp=',nblopt(i+1),lentyp(2,nblopt(i+1))
	if(lentyp(2,nblopt(i+1)).ne.0) go to 120
c*        iflg=2 => наклон фронта и кривизну (2) считать; + перебор по 
c*              всем звеньям сплайна (0 в первом бите)
	call  vprep(fi,lambda,dat,idat,lentyp,nblopt(i+1),nbopt,x0,
     *  r,wk(maxv+1),wk(maxv*2+1),wk(maxv*3+1),wk(maxv*4+1),
     *  wk(maxv*5+1),nv,2,eps,irc)
	if(nv.le.0) go to 120
	do 119 j=1,nv
c      wk(maxv*2+j)=wk(maxv*3+j)
c     sdf1=sdf1+wk(maxv*3+j)**2
      wk(maxv*3+j)=0.
119    wk(j)=40.
c119   wk(j)=21.
	ibgwk4=maxv
	call  datwrk(dat,idat,lentyp,lenmax,nb,nbmax,nbin,wk,
     *  wk(maxv*4+1),.123,ibgwk4,nv,6,2,irc)
	nbin=nb
120   continue
c      print *,'препарировали,nv=',nv
c*          Сортировка блока скоростей по возpастанию y-кооpдинаты
c       if(nv.gt.0) return
	ibg=lentyp(1,nb)
	nv=idat(ibg+3)
	ibg2=ibg/2+4
	m=4
	if(idat(ibg+5).eq.2) m=6
	call  sorthm(dat(ibg2),nv,m,2)
c*          Выбрасывание "хвостовых" скоростей. Последняя
c*          выбрасываемая скорость определяет динамическую высоту 
c*          и ширину фронта
c       if(nv.gt.0) return
	ib=ibg2+nv*2
	nvo=nv
	nv=0
	do 130 i=1,nvo
c!!!!!!!!!!!!!!!!!!!!
c       if(dat(ib+i-1).lt.vmin) go to 134
	nv=nv+1
130   continue
134   continue
c        if(i.le.nvo) r=dat(ibg2+nvo+i)
c  Включаем последнюю точку
c     nv=nv+1
c     r=dat(ibg2+nvo+nv)
c     if(nv.eq.nvo) r=x0+5.
      call velprf(lambda,dv2,prm,dat(ibg2+nvo),dat(ibg2+nvo+nvo),nv,
     * 2,eps,irc)
c       prm(1)=0.834
c       prm(2)=1.65+0.8
c       prm(3)=0.
c       prm(4)=11.03
c       prm(5)=0.
c       prm(6)=0.
c       prm(7)=1.80
       print *,'VELPRF:irc,nv.fx,prm=',irc,nv,dv2,(prm(ii),ii=1,7)
	if(nbgi.eq.0) go to 139
c
c*          Расчет динамической топографии между двумя станциями
c*      ys - расстояния от точки измерения до контура
	do 136 i=1,nst
	call  trdist(fi,lambda,xst(i),yst(i),dat(ibgopt),dat(ibgopt+
     *  nknot),dat(ibgopt+nknot*2),nknot,dat(ibgc),nknot-1,xc,yc,tc,
     *  tgac,vrc,dist,eps,1)
136   ys(i)=x0+dist
	if(i.gt.0) go to 1356
c*          подготовка массивов коэффициентов полинома  для расчета
c*                северной границы фронта
c*          ns - степень полинома
1366  ns=np-3+nbc
	if(nbc.eq.0) go to 1345
c*          обнуление nbc коэффициентов
	do 1340 k=1,nbc
1340  wk(ns+2-k)=0.
c*          перекачка ненулевых коэффициентов
1345  do 1350 k=nbc+1,ns+1
1350  wk(ns+2-k)=prm(k-nbc+2)*k
c*********                                    ********
c*          определяем возможную x-координату границы фронта
	wk(ns+1)=0.
	call zrpoly(wk,ns,wk(ns+2),irc)
	if(irc.ne.0) print *,'zrpoly-irc',irc
	xx=r
	dx=1000.
	dxp=1000.
	ns2=ns+ns
	do 1355 k=1,ns2,2
	if(wk(ns+k+2).ne.0) go to 1355
c*          корень не комплексный
	dxs=r-wk(ns+k+1)
	dxn=abs(dxs)
	if(dxn.gt.dx) go to 1355
c*          x-координата ближе к xv, чем предыдущие
	xx=wk(ns+k+1)
	dxp=dxs
	dx=dxn
1355  continue
1356  continue
c*           определение динамических высот
	do 1358 i=1,nst
c     if(ys(i).gt.r) go to 1358
	  call  clcr1(fi,lambda,prm,7,dat(ibgopt),
     *      nknot,ys(i),hc(i),vc,psic,1,eps,irc)
	print*,'x,динам. топогр. ',ys(i),hc(i)
1358  continue
c
c*          Расчет невязок расположения скоростей относительно 
c*    контура. Для скорости v рассчитывается x: 
c*    v=sum(i*a(i)*x**(i-1)). По невязкам между x и положением
c*    маркера скорости v рассчитывается целевая функция
139   nvo1=nvo-1
	nvn=0
c     fx=dv2
c     if(i.ne.333) return
c     x0=prm(1)
c     r=xx
	fx=0.
	vc=0.
c*          Если iflout=1 - готовим инфоpмацию для гpафика
c      if(iflout.eq.1) open(8,FILE='xvvc.dat'//8#0)
	do 190 i=1,nv
	xv=dat(ibg2+nvo1+i)
	v=dat(ibg2+nvo+nvo1+i)
	vv=dat(ibg2+nvo*3-1+i)
	 call  clcr1(fi,lambda,prm,7,dat(ibgopt),
     *      nknot,xv,hightc,vc,psic,2,eps,irc)
c       print *,'xv,v,vc,vv=',xv,v,vc,vv
c      if(iflout.eq.1) write(8,888) xv,v,vc
888   format(4e10.3)
c*          подготовка массивов коэффициентов полинома  ********
c*          ns - степень полинома
c     ns=np-3+nbc
c     if(nbc.eq.0) go to 145
c*          обнуление nbc коэффициентов
c     do 140 k=1,nbc
c140  wk(ns+2-k)=0.
c*          перекачка ненулевых коэффициентов
c145  do 150 k=nbc+1,ns+1
c150  wk(ns+2-k)=prm(k-nbc+2)*k
c*********                                    ********
c*          определяем возможную x-координату маркера скорости
c     wk(ns+1)=v
c     call zrpoly(wk,ns,wk(ns+2),irc)
c     if(irc.ne.0) print *,'zrpoly-irc',irc
c     dx=1000.
c     dxp=1000.
c     ns2=ns+ns
c     do 155 k=1,ns2,2
c     if(wk(ns+k+2).ne.0) go to 155
c*          корень не комплексный
c     dxs=xv-wk(ns+k+1)
c     dxn=abs(dxs)
c     if(dxn.gt.dx) go to 155
c*          x-координата ближе к xv, чем предыдущие
c     xx=wk(ns+k+1)
c     dxp=dxs
c     dx=dxn
c155  continue
c     print *,'dxp,xv=',dxp,xv
c     if(dx.ne.1000.) go to 185
c*          все корни комплексные
c     go to 190
c185  if(xx.lt.x0) xx=x0
c     if(xx.gt.r) xx=r
c     dx=abs(xv-xx)
c     nvn=nvn+1
c     if(typesf(2).eq.0) dx=dx*dx
c     fx=fx+dx
c     fx=fx+vv*vv
	  df1=abs(vv)
c        df=abs(v-vc)
	 df=0.
c        if(v.lt.0..and.xv.le.prm(2)+prm(4)) df=-v*10.
	vabs=sqrt(v*v+vv*vv)
	  if(typesf(2).eq.0) df1=df1*df1
	if(typesf(2).eq.0) df=df*df
c!!!!!            Невязки с весовыми коэффициентами /v/
	  sdf1=sdf1+df1
	  sdf=sdf+df
	fx=fx+df*vabs
	fx=fx+df1*vabs
190   continue
	if(nv.ne.0) fx=fx/nv
      print *,'nv,fx,sdf,sdf1',nv,fx,sdf,sdf1
	if(nv.ne.0) sdf1=sdf1/nv
	if(nv.ne.0) sdf=sdf/nv
c       print *,'fx,df,df1,nv,prm=',fx,sdf,sdf1,nv,(prm(i),i=1,7)
c      if(iflout.eq.1) close(8)
	if(nv.ne.0) return
199   np=0
	do 200 i=1,nbo
	ds=0.
	npnt=0
c     if(lentyp(2,nblopt(i+1)).eq.0) call  dfvel
	if(lentyp(2,nblopt(i+1)).eq.1) call  dfcont(fi,lambda,
     *dat(nbg(1,i)),nbg(2,i),dat(ibgopt),dat(ibgopt+nknot),dat(ibgopt+
     *nknot+nknot),nknot,dat(ibgopt+nknot*3),nknot-1,ds,npnt,typesf(2),
     *eps)
	if(npnt.le.0) go to 200
	fx=fx+ds/npnt
	np=np+1
200   continue
	if(np.gt.0) fx=fx/np
	return
	end
