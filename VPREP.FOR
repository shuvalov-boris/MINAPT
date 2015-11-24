c*********************************************************************
c*
c*
c*
c*        Пересчет положения маркеров скоростей блока nbin относительно
c*      положения контура блока nbcont, а также разложение скорости
c*      на компоненты и расчет наклона и кривизны контура в точках,
c*      ближайших к положениям маркеров, коррекция величин скоростей
c*      с учетом их криволинейного движения по контуру
c*
c*
c* Параметры:   fi,lambda,dat,idat,lentyp - описаны в 
c*                      программе inpsfd
c*              x0 - константа, которая добавляется к расстояниям 
c*                      маркеров относительно контура (x)
c*              xr - точка-ограничитель. Из блока nbin выбираются
c*                      только те v, положения маркеров которых относи-
c*                      тельно контура (x) лежат в диапазоне (x0.-xr)
c*              x - расстояние между маркером скорости и 
c*                      ближайшей точкой на сплайн-кривой в 10км
c*              vg - геострофическая компонента скорости, т.е. каса-
c*                      тельная к контуру.
c*              vag - агеострофическая, т.е. ортоганальная
c*              tga(i),vr(i) - наклон и кривизна контура для маркера i
c*  размерность tga,vr - nvi, где nvi - число скоростей в блоке nbin
c*              nv - число отобраных маркеров
c*              iflag=cab - флаг программы для trdist-двоичное число ab
c*                      b=0 - полный перебор звеньев сплайна при расчете
c*                      ближайшей точке на контуре
c*                      b=1 - расчет для двух ближайших звеньев
c*                      a=0 - tga и vr не считать
c*                      a=1 - считать
c*                      c=0 - коррекцию скоростей не проводить
c*                      с=1 - проводить
c*              eps - точность расчета координат в 10km
c*              irc - код возврата:
c*                      0 - нормальное завершение
c*                      1 - невозможно провести коррекцию скорости
c*                      >1 - ошибки в сплайн-программах - их коды
c*
c* Примечание. Если для блока с узлами nbcont сплайн-аппроксимация не 
c*      проводилась, то расчитываются сплайн-коэффициенты кривой,
c*      проходящей через узлы.
c*              Если точка на контуре, ближайшая к точке положения 
c*      маркера совпадает с концевой, то такой маркер исключается 
c*      из рассмотрения. Исключаются также маркеры с x не принадлежащих
c*      к интервалу (x0-xr)
c*              Коррекция скоростей проводится по круговой модели. 
c*      Рассчитывается радиус кривизны в точке, равный
c*      радиусу кривизны кривой + dist, вводятся поправки на угол
c*      и длину пути, исходя из предположения, что геострофическая и
c*      агеострофическая компоненты скорости были постоянны. Если 
c*      радиус кривизны в конечной точке меняет знак за счет агеостро-
c*      фической компоненты, то такая точка выбрасывается из рассмотрения
c*
c*
c*********************************************************************
	subroutine vprep(fi,lambda,dat,idat,lentyp,
     *nbin,nbcont,x0,xr,x,vg,vag,tga,vr,nv,iflag,eps,irc)
	real fi,lambda,dat(1),x(1),vg(1),vag(1),vr(1),tga(1),tgai(5),vri
     *  , x0, xr, xxyy(2), ffi, llambda, a, b, c, ndst, nvg, gtga(2),
	*  dtg(2), dnr(2), ptan(2), v(1), dvnorm
	integer*2 idat(1), nbin, nbcont
	integer*4 lentyp(2,1),ibg,ix,iy,iv1,iv2,ixk,iyk,itk,ic1,ic2,dd/-2/
	iflg=iflag-iflag/4*4
	ic=iflag/4
C		Вывод входных параметров
      open(unit=23, file='markers.txt', action='write')
	open(unit=24, file='profile.txt', action='write')
	open(unit=35, file='vtan.txt', action='write')
	open(unit=34, file='tga.txt', action='write')
	open(unit=37, file='vnorm.txt', action='write')
	open(unit=36, file='VvsS.txt', action='write')

	xxyy(1)=fi
      xxyy(2)=lambda
	call CONVGD(V,XxYy,WK,1,FI,LAMBDA,TIME,2,0)
	ffi=xxyy(1)
	llambda=xxyy(2)
c
c*              Расчет начал нужных массивов
c
c               nbin (блок скоростей)
	ibg=lentyp(1,nbin)	 !начало блока скоростей в массиве idat
	itypcn=idat(ibg+2)	 !в 10км или градусах?
	nvi=idat(ibg+3)		 !количество маркеров скоростей
	ix=ibg/2+4			 !начало x-координат точек в массиве dat
	time=dat(ix-3)		 !время от начала отсчета в часах
	idt=idat(ibg+4)		 !разница во времени между снимками	в минутах
	iy=ix+nvi			 !начало y-координат точек в массиве dat
	iv1=iy+nvi			 !начало x-координат скоростей в массиве dat
	iv2=iv1+nvi			 !начало y-координат скоростей в массиве dat
c               nknot
	ibg=lentyp(1,nbcont)
	nknot=idat(ibg+3)
	numvar=idat(ibg+5)
	ixk=ibg/2+4
	iyk=ixk+nknot
	itk=iyk+nknot
	ic1=itk+nknot
	ic2=ic1+(nknot-1)*3
	if(numvar.ne.0) go to 10
c*              контур не построен, рассчитываем сплайн-коэффициенты
      write(*, *) 'VPREP FAIL Рассчитываем сплайн коэффициенты'
	call  icsccu(dat(itk),dat(ixk),nknot,dat(ic1),nknot-1,irc)
	if(irc.ne.0) return
	call  icsccu(dat(itk),dat(iyk),nknot,dat(ic2),nknot-1,irc)
	if(irc.ne.0) return
	idat(ibg+5)=1
C	Создание файла для визуализации сплайна	  	
	  call WriteSplineXY(int(nknot+1), dat(cibg), dat(cibg+nknot+1), 
     *  dat(cibg+(nknot+1)*2), dat(cibg+(nknot+1)*3), fi, lambda,
     *  'spline.opt')
c
c*              расчет x,vg,vga,tga и vr c отбраковкой маркеров
c
10      continue
	nv=1
	tkbeg=dat(itk)
	tkend=dat(itk+nknot-1)
	ifl=6
	if(iflg-iflg/2*2.eq.1) ifl=7

	write(*, *) 'vprep: nvi =', nvi

	do 50 i=0,nvi-1
c      print *,'i=',i
c*              считаем x,tga и vr
	call trdist(fi,lambda,dat(ix+i),dat(iy+i),dat(ixk),
     *dat(iyk),dat(itk),nknot,dat(ic1),nknot-1,xc,yc,tc,tgai,
     *vri,dist,eps,ifl)
c      print*,'tkbeg,tkend,tc,dist+x0,xr=',tkbeg,tkend,tc,dist+x0,xr
C	if(tc.eq.tkbeg.or.tc.eq.tkend) then		! зачем ????
	if(tc.lt.tkbeg.or.tc.gt.tkend) then
	  write(*, *) 'line126: GOTO 49'
	  go to 49
	endif
	dist=dist+x0
	if(dist.le.0.or.dist.gt.xr) then
	  write(*, *) 'line134: GOTO 49'
	  goto 49
	endif
C*              точку необходимо вносить в массив
      xxyy(1)=dat(ix+i)
      xxyy(2)=dat(iy+i)
      call cnkmdg(fi,lambda,xxyy,v,1,dx,dy,0,0)

	gtga(1) = xc
	gtga(2) = yc
	call cnkmdg(fi,lambda,gtga,v,1,dx,dy,0,0)

	tga(nv)=tgai(1)
	vr(nv)=vri
	x(nv)=dist
	vi1=dat(iv1+i)
	vi2=dat(iv2+i)
	if(ic.eq.0.or.idt.eq.0.or.vri.eq.0) go to 30
c
c*              корректируем скорости
c        
c*              считаем радиус начальной точки
	radius=1/vri-(dist-x0)
	if(vri.lt.0.) radius=abs(1./vri)+dist-x0
c      print*, 'radius=',radius
	if(radius.le.0.) go to 49
	dt=idt*60.
c*              считаем радиус конечной точки
	x1=dat(ix+i)
	y1=dat(iy+i)
	x2=x1+vi1*dt*1.e-6
	y2=y1+vi2*dt*1.e-6
	dist0=dist-x0
	if(dist0.lt.0.) dist0=-dist0
c      xr=x1+(xc-x1)/dist0*radius*vri/abs(vri)
c      yr=y1+(yc-y1)/dist0*radius*vri/abs(vri)
	tg2=tga(nv)**2
	orty=sqrt(1./(1.+tg2))
	ortx=sqrt(tg2/(1.+tg2))
	if(tga(nv).gt.0..and.vri.gt.0.) orty=-orty
	if(vri.lt.0.) ortx=-ortx
	print *,'ORTX,ORTY=',ORTX,ORTY
	xr0=x1+ortx*radius
	yr0=y1+orty*radius
c       print*,'dx,dy=',xc-x1,yc-y1
	rnew=sqrt((xr0-x2)**2+(yr0-y2)**2)
c*              рассчитываем исправленные компоненты скорости
	vn=(radius-rnew)/dt*vri/abs(vri)*1.e+6
	dfi=((x1-xr0)*(x2-xr0)+(y1-yr0)*(y2-yr0))/radius/rnew
	  print *,'r0,rnew,sin,dist0=',radius,rnew,dfi,dist-x0
c       print*,'vi1,vi1=',vi1,vi2,'x12r=',x1,y1,x2,y2,xr0,yr0,xc,yc
	dfi=acos(dfi)
	if(vn.ne.0.) vp=dfi*vn/alog(rnew/radius)
	if(vn.eq.0.) vp=dfi*radius/dt*1.e+6
c*              рассчитываем vg и vag
30      dg=sqrt(tgai(2)**2+tgai(3)**2)
	xe=1.
	ye=0.
	if (dg.eq.0.) go to 40
	ye=tgai(3)/dg
	xe=tgai(2)/dg
40    vg(nv)=vi1*xe+vi2*ye 
	vag(nv)=-ye*vi1+xe*vi2
	if(vg(nv).lt.0.) vp=-vp
c      print*,'vold,vnew=',vg(nv),vag(nv),vp,vn
	if(ic.eq.1.and.idt.ne.0.and.vri.ne.0) vag(nv)=vn
	if(ic.eq.1.and.idt.ne.0.and.vri.ne.0) vg(nv)=vp
c        print *,'i,xc,yc,tc,dist,v,va=',i,xc,yc,tc,dist,sqrt(
c     * dat(iv1+i)**2+dat(iv2+i)**2),sqrt(vg(nv)**2+vag(nv)**2)
C      goto 666 

C		Нормирование вектора
      a = abs(x(nv) - x0)
	b = sqrt((dat(ix+i)-llambda)**2+(dat(iy+i)-ffi)**2)
	c = sqrt((xc-llambda)**2+(yc-ffi)**2)
	dst = (b**2 + c**2 - a**2) / (2 * c)
      ndst = dst / c
	nvg = vg(nv) / c
C		Расчет тангенциального вектора в геогр. координатах
      svg = abs(vg(nv)) / 10e5 * idt * 60 
      !vangle = 0.
	if (vi1.ne.0) then
	  vangle=atan2(vi2, vi1)
	else 
	  vangle = 3.14/2 * sign(1., vi2)
	endif
      sangle=atan(tga(nv))

	alpha=abs(vangle-sangle)
	alpha1 = alpha
	sgn = sign(1., 3.14/2-alpha)
	sgn1 = sgn
      dtg(1) = dat(ix + i) + sgn * svg * cos(sangle)
	dtg(2) = dat(iy + i) + sgn * svg * sin(sangle)
	call cnkmdg(fi,lambda,dtg,v,1,dx,dy,0,0)
C		Расчет нормального вектора в геогр. координатах
      snd = abs(vag(nv)) / 10e5 * idt * 60
C	vangle = atan2(-vi1, vi2)
      dvnorm = sangle + 3.14/2
	if (dvnorm.gt.3.14/2) dvnorm = dvnorm - 3.14
      alpha=abs(vangle-dvnorm)
	sgn = sign(1., 3.14/2-alpha)
      dnr(1) = dat(ix + i) + sgn * snd * cos(dvnorm)
	dnr(2) = dat(iy + i) + sgn * snd * sin(dvnorm)
	call cnkmdg(fi,lambda,dnr,v,1,dx,dy,0,0)

	ptan(1) = xc
	ptan(2) = yc
	call cnkmdg(fi,lambda,ptan,v,1,dx,dy,0,0)

      write(23, 7) xxyy(1), xxyy(2), !vg(nv), xsgn, ysgn, 
     *  vg(nv)!dist!alpha !atan(tga(nv)) !vangle !sign(1., tga(nv))
	write(36, 7) xxyy, ptan
	write(34, 6) ptan, tga(nv)
	write(35, 7) xxyy, dtg !, vg(nv), svg
	write(37, 7) xxyy, dnr
      write(24, 7) dist, vg(nv)
     *  !ndst,
     *  !c, 
     *  !nvg
    6 format(3F11.3)
    7 format(8F11.3)
    8 format(5F11.3, I7)
    9 format(12F9.3)
666	nv=nv+1
49	write(*, *) nv, i , dd
      xxyy(1)=dat(ix+i)
      xxyy(2)=dat(iy+i)
      call cnkmdg(fi,lambda,xxyy,v,1,dx,dy,0,0)
      if (nv+dd.ne.i) then
        write(*, *) 'vprep: NOT EQUAL!!!', i, xxyy, tc, tkbeg, tkend,
     *    dist, xr, radius      
        dd = dd + 1
	endif
50    continue
	nv=nv-1
	irc=0
	if(ic.eq.1.and.idt.eq.0) irc=1
   99 close(23)
	close(24)
	close(36)
	close(34)
	close(35)
	close(37)
      return
	end
