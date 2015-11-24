$notstrict												
c*********************************************************************
c*
c*
c*     Расчет x,y,t - координат точки сплайн-кривой с минимальным 
c*  расстоянием до точки (xb,yb), кривизны в найденной точке , 
c*  наклона и расстояния между точками.
c*
c* Параметры:   fi,lambda - начало географической системы координат
c*              xb,yb,x,y,t - описаны в заголовке
c*              xk,yx,tk,nknot,c - параметры сплайна (см.,например,
c*                      описание программы spapp)
c*              ic=nknot-1
c*              tga(5) - тангенс угла наклона сплайна в точке (x,y) и
c*                первые и вторые производные x и y по t (dx,dy,ddx,ddy)
c*              vr - кривизна сплайна в точке (x,y).
c*                      vr=1/(радиус кривизны)
c*              dist - расстояние между точками (x,y) и (xb,yb) в 10км
c*                      если dist>0, то точка (xb,yb) лежит слева от
c*                      сплайн-кpивой пpи движении вдоль последней
c*                      по напpавлению возpастания t; <0 - спpава
c*              eps - точность расчета расстояния (вход) в 10 км
c*              iflg=cab (c,a,b - двоичные числа) - флажок:
c*                      b=0 - полный перебор звеньев сплайна при 
c*                      расчете точки (x,y)
c*                      b=1 - расчет для двух звеньев - i и i-1 :
c*                      точка (xk(i),yk(i)) - ближайшая к (xb,yb)
c*                      a=0 - tga(1) и vr не считать
c*                      a=1 - считать
c*                      c=0 - массив tga(2-5) - не заполнять
c*                      c=1 - заполнять (в старой версии программы
c*                        tga была переменная, а не массив)
c*
c* Примечание: xb,yb,x,y - в 10 km в декартовой системе ; fi,lambda - 
c*              в градусах
c* Алгоритм. Для получившейся оптимальной точки в декартовой системе
c*      координат (система описана в программе rdconv) рассчитывается
c*      точность определения dist. Если эта точность превышает
c*      eps, то переходим в географическую систему координат и методом
c*      "золотого сечения" уточняем положение точки (x,y). Направление
c*      поиска по t-координате выбирается в результате численной
c*      оценки градиента функции расстояния между точками с dt=amin(
c*      eps/10.,(tk(nknot)-tk(1))/nknot/1000.). Т.о. считается, что
c*      смещение dt по t-координате мало по сравнению с изменчивостью
c*      функции расстояния. Коэффициенты tga и vr pассчитываются
c*      без попавок на кpиволинейность декаpтовой системы кооpдинат.
c*      Неучет поправок на криволинейность приводит к следующим погреш-
c*      ностям (которые прямо пропорциональны расстоянию по координате x
c*      (долгота)): при x=10000км (около 10 град.) на широте 45 град. 
c*      ошибка в определении угла наклона 1.6 град., ошибка кривизны
c*      будет менее 4.5E-4 (радиус кривизны >2000км)
c*
c*
c*********************************************************************
	subroutine trdist(fi,lambda,xb,yb,xk,yk,tk,nknot,c,ic,x,y,t,tga,
     *  vr,dist,eps,iflg)
	external fndst
	real fi,lambda,xk(nknot),yk(nknot),c(ic,3,2),a(6),xyd(4),
     *  z(2,5),parm(6),wk(10),ts(2),ds(2),tk(nknot),tga(5)
	integer*2 kod(6)
	data kod/2,1,80,100,1000,0/parm/0.,.1,100.,.001,1.,.001/
	data pi/3.1415926/,onedg/11.124/
c*              поиск ближайшего к точке (xb,yb) узла сплайна
	dist=1.e32
	is=0
	do 10 i=1,nknot
	d=(xk(i)-xb)**2+(yk(i)-yb)**2
	if(d.ge.dist) go to 10
c*              найден более близкий узел
	is=i
	dist=d
10      continue
c
	x=xk(is)
	y=yk(is)
	t=tk(is)
c*              поиск x,y,t в декартовой системе координат
	ib=1
	ie=nknot
	if(iflg-iflg/2*2.eq.0) go to 20
c*              поиск не по всем звеньям
	ib=is-1
	ie=is
	if(ib.lt.1) ib=1
	if(ie.ge.nknot) ie=nknot-1
c
c*              цикл по звеньям
c
20      do 70 i=ib,ie
c*              расчет коэффициентов полинома - производной от
c*              функции расстояния между точками
		do 30 j=1,6
30              a(j)=0.
	xyd(1)=xk(i)-xb
	xyd(2)=yk(i)-yb
		do 40 j=1,2
		a(1)=a(1)+3.*c(i,3,j)**2
		a(2)=a(2)+5.*c(i,2,j)*c(i,3,j)
		a(3)=a(3)+4.*c(i,3,j)*c(i,1,j)+2.*c(i,2,j)**2
		a(4)=a(4)+3.*c(i,3,j)*xyd(j)+3.*c(i,2,j)*c(i,1,j)
		a(5)=a(5)+2.*c(i,2,j)*xyd(j)+c(i,1,j)**2
		a(6)=a(6)+c(i,1,j)*xyd(j)
40              continue
c*              расчет корней полинома
	ns=5
50    call  zpolr(a(6-ns),ns,z,ier)
	ns=ns-1
	if(ier.eq.130.and.ns.eq.0) go to 70
	if(ier.eq.130) go to 50
c*              проверяем получившиеся точки на оптимальность dist
	ns=ns+1
		do 60 j=1,ns
		if(z(2,j).ne.0.) go to 60
c*              корень полинома не комплексный
		if(z(1,j).lt.0..or.z(1,j).gt.tk(i+1)-tk(i)) go to 60
c*              t-координата (z(1,j)) - в нужном диапозоне
		call  icsevu(tk,xk,nknot,c,ic,z(1,j)+tk(i),xx,1,ier)
		call  icsevu(tk,yk,nknot,c(1,1,2),ic,z(1,j)+tk(i),yy,
     *          1,ier)
		d=(xx-xb)**2+(yy-yb)**2
		if(d.gt.dist) go to 60
		t=tk(i)+z(1,j)
		x=xx
		y=yy
		dist=d
60              continue
70      continue
c
	dist=sqrt(dist)
	if(6.284*dist.le.eps) go to 120
c*              оценка точности расчета расстояния
	xyd(1)=xb
	xyd(3)=yb
	xyd(2)=x
	xyd(4)=y
	call  ckd(fi,lambda,xyd,dx)
	epsx=abs(x-xb)*dx/dist
	xbd=xyd(1)
	ybd=xyd(3)
	if(epsx.le.eps) go to 120
c
c*              необходимо уточнить положение точки (x,y)
c
	g=1.
	dt=amin1(eps/2.,(tk(nknot)-tk(1))/nknot/1000.)
	if(dt.le.0.) dt=amax1(eps/2.,(tk(nknot)-tk(1))/nknot/1000.)
	t1=t+dt
	if(t1.gt.tk(nknot)) g=-1.
	t1=t+dt*g
	call  fndst(t,d,xbd,ybd,fi,lambda,xk,yk,tk,nknot,c,ic)
	call  fndst(t1,d1,xbd,ybd,fi,lambda,xk,yk,tk,nknot,c,ic)
	if(d1.le.d) go to 95
c*              в выбранном направлении функция расстояния растет.
c*              меняем направление
	g=-g
	t1=t+dt*g
	if(t1.lt.tk(1).or.t1.gt.tk(nknot)) go to 120
	call  fndst(t1,d1,xbd,ybd,fi,lambda,xk,yk,tk,nknot,c,ic)
	if(d1.ge.d) go to 120
c*              уточняем положение точки (x,y)
95      stp=100.*dt
c       print *,'dt,g,t,t1',dt,g,t,t1
	if(t+stp*g.gt.tk(nknot)) stp=tk(nknot)-t
	if(t+stp*g.lt.tk(1)) stp=t-tk(1)
	ts(1)=t
	parm(2)=stp
	parm(3)=tk(nknot)-tk(1)
	parm(4)=dt
	parm(6)=dt
	kod(4)=100
	kod(2)=1
	ds(1)=d
	call ml5(kod,1,ts,ds,g,fndst,parm,wk,ier,xbd,ybd,fi,lambda,
     *  xk,yk,tk,nknot,c,ic)
	dist=ds(1)
	t=ts(1)
	if(t.gt.tk(nknot)) t=tk(nknot)
	if(t.lt.tk(1)) t=tk(1)
	call  icsevu(tk,xk,nknot,c,ic,t,x,1,ier)
	call  icsevu(tk,yk,nknot,c(1,1,2),ic,t,y,1,ier)
c
120     continue
c*              рассчитываем кривизну и наклон
	call  dcsevu(tk,xk,nknot,c,ic,t,dx,1,ddx,1,ier)
	call  dcsevu(tk,yk,nknot,c(1,1,2),ic,t,dy,1,ddy,1,ier)
c*              опpеделение знака dist
	sgn=dy*(x-xb)-dx*(y-yb)
	if(sgn.lt.0.) dist=-dist
	if(iflg/4-iflg/8*2.eq.0) go to 125
c*          Заполняем массив tga
	tga(2)=dx
	tga(3)=dy
	tga(4)=ddx
	tga(5)=ddy
125   if(iflg/2-iflg/4*2.eq.0) return
c*          Считаем тангенс и кривизну
	tga(1)=1.e+32
	if(abs(dx).lt.1.e-32) go to 130
	tga(1)=dy/dx
130     sum=sqrt(dx*dx+dy*dy)**3
	if(sum.lt.1.e-30) vr=1.e30
	if(sum.lt.1.e-30) return
	vr=(dx*ddy-dy*ddx)/sum
	return
90    write(*, *) 'TRDIST: Opening file error'
      RETURN
92    write(*, *) 'TRDIST: Closing file error'
      RETURN
	end
C*********************************************************************
C*
C*
C*         программа, обратная п/п convgd в программе rdconv.
C*        осуществляет перевод информации в
C*        географическую систему координат
C*
C* параметры: fi,lambda - начало координат (географ.)
C*            xyd(4) - координаты двух точек  в 10км.
C*              на выходе - в географической системе координат
C*            dx - точность расчета растояния между точками по х- 
C*              координате в 10км (выход)
C*
C*
C*********************************************************************
	subroutine ckd(fi,lambda,xyd,dx)
	real*4 xyd(2,2),lami,lambda,ri(2),fi2(2)
	data onedg/111.24/,r0/6378.245/,pi/3.1415926/
	do 50 i=1,2
c*              преобразование данных из декартовой системы координат в
c*              географическую
10      call  dgkm(x0,y0,fi,lambda)
	y=xyd(i,2)*10.+y0
	lami=y/onedg
	rearth=r0*(1.-0.00335*sin(lami*pi/180.)**2)
	rlitl=rearth*cos(lami*pi/180.)
	ri(i)=rlitl
	fii=fi+xyd(i,1)*1800./pi/rlitl
	fi2(i)=fii
	if(fii.gt.180.) fii=fii-360.
	xyd(i,1)=fii
	xyd(i,2)=lami
50      continue
c*              расчет dx
	dfi=abs(fi2(2)-fi2(1))
	if(dfi.gt.180.) dfi=360.-dfi
	dx=pi*abs(ri(2)-ri(1))*dfi/1800.
	return
	end     
C#############################################################################
	SUBROUTINE ML5(KOD,N,X,FX,S,F,EST,Y,IER,xbd,ybd,fi,lambda,xk,yk,
     *tk,nknot,cc,ic)
	real X(N),Y(N),S(N),cc(ic,3,2),xk(nknot),yk(nknot),tk(nknot)
	real EST(1),lambda
	integer*2 kod(1)
	DATA C/1.618034/
	DATA C1/0.381966/
	DATA C2/0.618034/
C      BЫБOP ДЛИHЫ OTPEЗKA
	LIM=MIN0(KOD(3),KOD(4))
	L1=LIM
	H=0.
	DO 66 I=1,N
   66 H=H+S(I)*S(I)
	GN=SQRT(H)
	H=EST(2)*EST(5)/GN
	DL=EST(4)/GN
	A=0.
	B=H
	H1=0.
	K=LIM-1
    5 DO 2 I=1,N
    2 Y(I)=X(I)+B*S(I)
	CALL F(Y,FY,xbd,ybd,fi,lambda,xk,yk,tk,nknot,cc,ic)
	LIM=LIM-1
	IF(LIM)50,50,55
   55 IF(FY-FX)3,4,4
    3 A=A+H1
	ALFA1=A+H
	H1=H
	H=C*H
	FX=FY
	B=ALFA1+H
	GO TO 5
    4 RAZN=B-A
	IF(K-LIM)101,6,21
    6 FZ=FX
   19 ALFA1=A+C1*RAZN
	DO 8 I=1,N
    8 Y(I)=X(I)+ALFA1*S(I)
	CALL F(Y,FY,xbd,ybd,fi,lambda,xk,yk,tk,nknot,cc,ic)
	LIM=LIM-1
	IF(LIM.LE.0.AND.FX.LT.FZ)GO TO 53
	IF(LIM)17,17,18
   18 IF(FX-FZ)21,16,16
   16 B=ALFA1
	RAZN=B-A
	IF(RAZN-DL)17,17,19
   17 FX=FZ
	GO TO 44
   21 ALFA2=A+C2*RAZN
	DO 10I=1,N
   10 Y(I)=X(I)+ALFA2*S(I)
	CALL F(Y,FY,xbd,ybd,fi,lambda,xk,yk,tk,nknot,cc,ic)
	LIM=LIM-1
	IF(LIM)40,40,7
C      ПOИCK MИHИMУMA ФУHKЦИИ METOДOM ЗOЛOTOГO CEЧEHИЯ
    7 IF(FX-FY)12,11,11
   11 IF(RAZN-DL)40,40,20
   20 RAZN=C2*RAZN
	A=ALFA1
	ALFA1=ALFA2
	FX=FY
	GO TO 21
   12 IF(RAZN-DL)53,53,15
   15 RAZN=C2*RAZN
	B=ALFA2
	ALFA2=ALFA1
	FY=FX
	ALFA1=A+C1*RAZN
	DO 13 I=1,N
   13 Y(I)=X(I)+ALFA1*S(I)
	CALL F(Y,FY,xbd,ybd,fi,lambda,xk,yk,tk,nknot,cc,ic)
	LIM=LIM-1
	IF(LIM)53,53,7
   40 D0=ALFA2
  100 FX=FY
   58 DO 41 I=1,N
   41 X(I)=X(I)+D0*S(I)
  103 EST(2)=D0*GN
  101 KOD(4)=KOD(4)-L1+LIM
	IER=0
	RETURN
C      BЫXOД ПO ЛИMИTИPOBAHHOMЫ KOЛИЧECTBУ BЫЧИCЛEHИЙ ФУHKЦИИ
   50 IF(K-LIM)101,51,45
   45 IF(FX-FY)53,57,57
   57 D0=B
	GO TO 100
   53 D0=ALFA1
	GO TO 58
   51 IF(FX-FY)44,42,42
   42 D0=H
	GO TO 100
   44 EST(2)=0.
	GO TO 101
	END
c********************************************************************
c*
c*
c*      целевая функция - расчет расстояния от точки (xbd,ybd) до
c*  точки на сплайн-кривой с координатой t
c*
c* параметры: d - искомое расстояние в 10км
c*            fi,lambda - начало геогафической системы координат
c*            xk,yk,tk,c,ic - параметры сплайна
c*
c* примечание. если t лежит вне диапозона (tk(1)-tk(nknot)), то
c*      к d прибавляется штраф.
c*
c*
c********************************************************************
	subroutine fndst(t,d,xbd,ybd,fi,lambda,xk,yk,tk,nknot,c,ic)
	real lambda,xk(nknot),yk(nknot),tk(nknot),c(ic,3,2)
	real xy(2)
	ts=t
	fine=0.
	if(t.lt.tk(1)) ts=tk(1)
	if(t.lt.tk(1)) fine=tk(1)-t
	if(t.gt.tk(nknot)) ts=tk(nknot)
	if(t.gt.tk(nknot)) fine=t-tk(nknot)
	call  icsevu(tk,xk,nknot,c,ic,ts,xy,1,ier)
	call  icsevu(tk,yk,nknot,c(1,1,2),ic,ts,xy(2),1,ier)
	call  cnkmdg(fi,lambda,xy,v,1,dx,dy,0,0)
	call  distdg(xy,xy(2),xbd,ybd,d)
	d=d+fine
	return
	end