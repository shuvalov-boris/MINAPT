c*********************************************************************
c*
c*
c*          работа с данными файла dat(idat)
c*
c*
c* Параметры:   dat,idat,lentyp,lenmax,nb,nbmax, - описаны в inpsfd
c*          nbin - номер текущего блока
c*          x,y,t - массивы данных (свои для каждого вида работ)
c*          ibgwd - порядковый номер начала вводимых данных в 
c*             массиве данных dat блока nbin
c*          n - число вводимых элементов данных (при iwk=4 - число
c*             оптимизируемых переменных) 
c*          iwk - тип работы:
c*             1 - сменить тип блока с 1 на 2
c*             2 - удалить блок
c*             3 - добавить (n>0) или удалить (n<0) узлы или скорости
c*             4 - создание или модернизация оптимизационного блока
c*             5 - создание копии файла
c*             6 - создание блока 'препарированных' скоростей или
c*                добавление к блоку
c*             7 - создание блока фрагментов контура
c*			 8 - создание блока узлов по блоку контура
c*             9 - изменение блока узлов сплайна
c*          iflg - флажек управления
c*          irc - код возврата:
c*             0 - нормальное завершение
c*             1 - нет места в массиве dat(idat)
c*             2 - нет свободных блоков
c*
c* Примечание. Типы входных параметров для разных видов работ разные,
c*    поэтому они иногда описываются перед блоком, выполняющим 
c*    указанную работу.
c*
c*
c*********************************************************************
	subroutine datwrk(dat,idat,lentyp,lenmax,nb,nbmax,nbin,x,y,t,
     *  ibgwd,n,iwk,iflg,irc)
	real*4 dat(1),x(1),y(1),t(1),rl
	integer*2 idat(1),intg(2), nb, nbmax, nbin, nknot, nc
	integer*4 lentyp(2,1),lenmax,ibgwk4,ibeg,ibg,ibg2,len,
     *  ibgfr,ib
	equivalence (rl,intg)
	data lencdt/53/
	go to (10,20,30,60,90,100,140, 200, 210),iwk
c
c*          смена типа блока (точки контура -> узлы сплайна)
c
10    irc=1
	ibeg=lentyp(1,nbin)+6
	nknot=idat(ibeg-3)
	nwd=idat(ibeg-2)
	if(nwd-nwd/2*2.eq.1) nwd=nwd+1
	dl=(nknot-1)*14+2-nwd
C	write(6,1) lenmax, nb, dl
C    1 format(' lenmax = ', I9, ', nb = ', I4, ', dl = ', I7)
	if(lenmax.lt.lentyp(1,nb+1)+dl-1) return
c*          место в массиве idat(dat) есть
	nwd=-nwd
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibeg,nwd,irc)
	ibeg=ibeg+nknot*4
	nwd=dl-nwd
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibeg,nwd,irc)
c*          обнуляем информацию в справке блока
	ibeg=lentyp(1,nbin)
	idat(ibeg)=0
	idat(ibeg+1)=0
	idat(ibeg+4)=0
	idat(ibeg+5)=0
	lentyp(2,nbin)=2
	return
c
c*          удаление блока
c
20    ibeg=lentyp(1,nbin)
	nwd=ibeg-lentyp(1,nbin+1)
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibeg,nwd,irc)
	return
c
c*          добавление/удаление узлов или скоростей к блоку
c*     в случае добавления узлов переменные x,y,t - одномерные массивы
c*     в случае со скоростями - x(n,2) - координаты, y(n,2) - 
c*     скорости, t(n,2) - наклон фронта и кривизна (только для 
c*     "препарированных" скоростей)
c
30    nknot=idat(lentyp(1,nbin)+3)
	itpblk=lentyp(2,nbin)
	if(n.le.0) go to 35
c*          узлы добавляем
	irc=1
	dl=n*8
	if(idat(lentyp(1,nbin)+5).eq.2) dl=dl+n*4
	if(itpblk.eq.2) dl=n*18
	if(lentyp(1,nb+1)+dl-1.gt.lenmax) return
c*          место в массиве dat(idat) еще есть
35    ibeg=lentyp(1,nbin)+6+(ibgwd-1)*2
	n2=n*2
	ntimes=3
	if(itpblk.eq.0) ntimes=4
	if(itpblk.eq.0.and.dl.eq.n*12) ntimes=6
	do 50 i=1,ntimes
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibeg,n2,irc)
	if(n.lt.0) go to 50
	if(itpblk.eq.0) go to 45
c*          добавляем узлы
	  do 40 j=1,n
	  if(i.eq.1) dat(ibeg/2+j)=x(j)
	  if(i.eq.2) dat(ibeg/2+j)=y(j)
	  if(i.eq.3) dat(ibeg/2+j)=t(j)
40      continue
	  go to 50
c*          добавляем скорости
45      continue
	  do 47 j=1,n
	  if(i.eq.1) dat(ibeg/2+j)=x(j)
	  if(i.eq.2) dat(ibeg/2+j)=x(j+n)
	  if(i.eq.3) dat(ibeg/2+j)=y(j)
	  if(i.eq.4) dat(ibeg/2+j)=y(j+n)
	  if(i.eq.5) dat(ibeg/2+j)=t(j)
	  if(i.eq.6) dat(ibeg/2+j)=t(j+n)
47      continue
50    ibeg=ibeg+n2+nknot*2
	if(itpblk.eq.0) go to 55
c
c*          меняем размер массива под коэффициенты сплайна
	ibeg=ibeg-(ibgwd-1)*2
	nwd=n*12
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibeg,nwd,irc)
c*          корректируем блок справки
55    idat(lentyp(1,nbin)+3)=nknot+n
	if(itpblk.eq.2) idat(lentyp(1,nbin)+5)=0
	return
c
c*          Создаем или модернизируем оптимизационный блок
c
c*  На входе: x(n) - оптимизируемые переменные,
c*          y(n) - их статусы (целые, integer*2)
c*          t(52) - параметры оптимизации
c*          iflg - флаг имени целевой функции
c*          ibgwd - номер блока, чьи параметры оптимизируются
c
c*          Проверка возможности записи в архив
60    irc=2
	if(nbin.eq.0.and.nb.eq.nbmax) return
	irc=1
	leninf=lencdt*2+3*n+6
	lenold=0
	if(nbin.ne.0) lenold=lencdt*2+idat(lentyp(1,nbin)+3)*3+6
	idlen=leninf-lenold
	if(lenmax.lt.lentyp(1,nb+1)+idlen-1) return
	irc=0
c
	if(lenold.eq.0) nb=nb+1
	if(lenold.eq.0) nbin=nb
	if(lenold.eq.0) lentyp(1,nb+1)=lentyp(1,nb)+leninf
	if(lenold.eq.0) go to 65
c*          Необходимо сдвинуть данные в архиве
	ibgwk4=112+lentyp(1,nbin)
	call  movdat(idat,lentyp,lenmax,nb,nbin,ibgwk4,idlen
     *  ,irc)
c
c*          Заносим данные в блок
c
65    lentyp(2,nbin)=4
c*          Корректируем справочник блока
	ibg=lentyp(1,nbin)+6
	ibg2=ibg/2
	dat(ibg2-2)=0.123
	idat(ibg-3)=n
	idat(ibg-4)=iflg
	idat(ibg-2)=ibgwd
	idat(ibg-1)=0
c*          Ввод параметров оптимизации
	do 70 i=1,lencdt
70    dat(ibg2+i)=t(i)
c*          Ввод переменных и задание fxbest-начального
	ibg2=ibg2+lencdt
	dat(ibg2)=1.e30
	do 75 i=1,n
75    dat(ibg2+i)=x(i)
c*          Ввод статусов переменных
	n2=n/2
	ibg2=ibg2+n
	if(n2.eq.0) go to 85
	do 80 i=1,n2
80    dat(ibg2+i)=y(i)
85    if(n-n2*2.eq.0) return
c*          n - нечетное число
	rl=y(n2+1)
	idat(ibg+lencdt*2-1+n*3)=intg(2)
	return
c
c*          Создание копии файла
c
90    irc=2
	if(nb.ge.nbmax) return
	ibg=lentyp(1,nbin)
	len=lentyp(1,nbin+1)-lentyp(1,nbin)
	irc=1
	ibgfr=lentyp(1,nb+1)
	if(ibgfr+len-1.gt.lenmax) return
c*          В архиве есть место для копирования
	irc=0
	do 95 i=0,len-1
95    idat(ibgfr+i)=idat(ibg+i)
	nb=nb+1
	lentyp(1,nb+1)=lentyp(1,nb)+len
	lentyp(2,nb)=lentyp(2,nbin)
	return
c
c*          Заведение блока скоростей или добавление к блоку
c
c*
c* Параметры: структура массивов x(ibgwd,4),y(ibgwd,2):
c*    в x - x и y координаты маркеров и две компоненты скорости
c*    в y - наклон и кривизна контура в точках, ближайших к точкам
c*    положений маркеров
c*    n - число маркеров, n<=ibgwd
c*    t - время, вносимое в справку блока
c*    iflg - значение 6-го элемента справки
c*    nbin - <=0:
c*          тип координат (itypcn=iabs(nbin))
c*          >0:
c*          номер блока к которому добавляются данные
c*
c*
100   continue
c*          проверка наличия свободных блоков и места
	irc=2
	if(nb.ge.nbmax.and.nbin.le.0) return
	ibg=lentyp(1,nb+1)
	if(nbin.gt.0) ibg=lentyp(1,nbin)
	len=n*8
	if(nbin.le.0) len=len+6
	if(nbin.gt.0) iflg=idat(ibg+5)
	if(iflg.eq.2) len=len+n*4
	irc=1
	if(lentyp(1,nb+1)+len-1.gt.lenmax) return
	irc=0
	no=0
	if(nbin.gt.0) no=idat(ibg+3)
	nn=no+n
	if(nbin.gt.0) go to 105
c*          формируем справку блока
	dat(ibg/2+1)=t(1)
	idat(ibg+2)=iabs(nbin)
	idat(ibg+3)=n
	idat(ibg+4)=0
	idat(ibg+5)=iflg
	go to 109
c*          Раздвигаем блок nbin
105   ntime=4
	if(iflg.eq.2) ntime=6
	ib=no+no+ibg+6
	do 107 i=1,ntime
	call  movdat(idat,lentyp,lenmax,nb,nbin,ib,n+n,irc)
107   ib=ib+nn+nn
c*          формируем блоки
109   iy=ibgwd
	iv1=iy+ibgwd
	iv2=iv1+ibgwd
	ibg2=ibg/2+3+no
	n2=nn+nn
	n3=n2+nn
	n4=n3+nn
	n5=n4+nn
	do 110 i=1,n
	dat(ibg2+i)=x(i)
	dat(ibg2+nn+i)=x(iy+i)
	dat(ibg2+n2+i)=x(iv1+i)
110   dat(ibg2+n3+i)=x(iv2+i)
c
	if(iflg.ne.2) go to 130
c*          Заносим наклон и кривизну
	do 120 i=1,n
	dat(ibg2+n4+i)=y(i)
120   dat(ibg2+n5+i)=y(ibgwd+i)
c
130   if(nbin.gt.0) idat(ibg+3)=nn
	if(nbin.gt.0) return
	nb=nb+1
	lentyp(2,nb)=0
	lentyp(1,nb+1)=lentyp(1,nb)+len
	return
c*
c*          заведение блока фрагментов контура
c*
c* параметры: x(n),y(n) - координаты точек контура
c*          t(1) - время, вносимое в справку блока
c*          iflg - значение 6-го элемента справки
c*          nbin - тип координат (itypcn)
c*          
c* примечание. совпаданеие координат двух последовательно идущих
c*    в списке точек означает окончание фрагмента. одна из таких 
c*    точек выбрасывается
c*
140   continue
c*          проверка наличия свободных блоков и места в архиве
	irc=2
	if(nb.ge.nbmax) return
	irc=1
	ibg=lentyp(1,nb+1)
	len=n*4+6+2
	if(ibg+len-1.gt.lenmax) return
c*          в архив писать можно
	irc=0
	np=1
	nf=1
	idat(ibg+6)=1
	if(n.le.1) go to 170
c
c*          разбиваем точки контура на фрагменты
c
	nplook=0
150   nplook=nplook+1
	if(nplook.eq.n) go to 170
	np=np+1
	if(x(np-1).ne.x(np).or.y(np-1).ne.y(np)) go to 150
c*          конец фрагмента
	if(nplook.eq.n-1) np=np-1
	if(nplook.eq.n-1) go to 170
c*          запоминаем положение начала нового фрагмента
	nf=nf+1
	idat(ibg+5+nf)=np
c*          сдвиг массивов координат
	ie=n-nplook-nf+1
	do 160 i=1,ie
	x(np+i-1)=x(np+i)
160   y(np+i-1)=y(np+i)
c
c*          коррекция счетчиков числа точек и числа просмотренных 
c*                      точек
	nplook=nplook+1
	go to 150
c
170   continue
c
c*          формируем справку блока
c
	dat(ibg/2+1)=t(1)
	idat(ibg+2)=nbin
	idat(ibg+3)=np
	idat(ibg+4)=nf
	idat(ibg+5)=iflg
c
c*          внесение координат в блок
	ibg=(ibg+nf+6)/2
	do 180 i=1,np
	dat(ibg+i)=x(i)
180   dat(ibg+i+np)=y(i)
c*          коррекция таблицы блоков
	nb=nb+1
	lentyp(1,nb+1)=lentyp(1,nb)+np*4+6+nf+nf-nf/2*2
	lentyp(2,nb)=1
	return
C
C       Создание блока узлов размером n, распределенных
C	  равномерно (по точкам) по всему контуру (блок nbin)
C
200   nknot = n + 1
      irc=3
	if (n.le.1) goto 222
      irc=2
	if(nb+1.ge.nbmax) goto 222
	ibg=lentyp(1,nbin)
	irc=1
	ibgfr=lentyp(1,nb+1)
	len=6+nknot*6+(nknot-1)*12
	if(ibgfr+len-1.gt.lenmax) goto 222
c*          В архиве есть место для копирования
	irc=0
	NC=idat(ibg+3)
	idat(ibgfr)  =0
	idat(ibgfr+1)=0
	idat(ibgfr+2)=idat(ibg+2)
	idat(ibgfr+3)=nknot
	idat(ibgfr+4)=0
	idat(ibgfr+5)=0

       nf=idat(ibg+4)
	 nf1=nf
	 if(nf-nf/2*2.eq.1) nf1=nf1+1
	ibg=ibg+6+nf1
	ibgfr=ibgfr+6

	do i = 1, nknot - 1
	  idx = ((NC - 1) * (i - 1) / (nknot-1)) + 1
	  dat(ibgfr/2+i) = dat(ibg/2+idx)
	  dat(ibgfr/2+i+nknot) = dat(ibg/2+idx+NC)
	enddo
	dat(ibgfr/2+nknot) = dat(ibgfr/2 + 1)
	dat(ibgfr/2+nknot*2) = dat(ibgfr/2 + nknot + 1)	
	call Parametrize(nknot, dat(ibgfr/2+1), dat(ibgfr/2+nknot+1),
     *  dat(ibgfr/2+nknot*2+1))
	write(*, *) 'parametrize', dat(ibgfr/2+1), dat(ibgfr/2+nknot+1),
     *  dat(ibgfr/2+nknot*2+1) 
      nb=nb+1
      lentyp(1,nb+1)=lentyp(1,nb)+6+nknot*6+(nknot-1)*12
      lentyp(2,nb)=2
	nbin=nb
      return
C
C	изменение блока узлов и его параметризация
C*		x(n),y(n) - координаты узлов
c*          nbin - тип координат (itypcn)
C
  210	  nknot = n
  	  write(*, *) 'DATWRK: nb =', nb, ', n =', n, ', nbin =', nbin,
     *    ', len =', len, ', LT(1) =', lentyp(1, nb),
     *    ', LT(2) =', lentyp(2, nb)
        irc = 3
	  if (n.le.1.or.nbin.le.0) goto 222
c*          В архиве есть место для копирования
        irc=0
	  ibg=lentyp(1,nbin)/2+3
        write(*, *) 'datwrk: изменение данных'
	  do i = 1, nknot
	    dat(ibg+i) = x(i)
	    dat(ibg+i+nknot+1) = y(i)
	  enddo
	  dat(ibg+nknot+1) = x(1)
	  dat(ibg+(nknot+1)*2) = y(1)
	  nknot=nknot+1
	  call Parametrize(nknot, dat(ibg+1), 
     *    dat(ibg+nknot+1), dat(ibg+nknot*2+1))
        write(*, *) 'DW DAT: ' 
	  write(*, '(3F12.6)') (dat(ibg+i), dat(ibg+nknot+i)
     *    , dat(ibg+nknot*2+i), i = 1, nknot)
        return

  222 end


C**********************************************************************
C*
C*						Parametrize
C*
C*		Параметризация массива точек в порядке записи в массиве
C*
C*	Параметры:
C*		N - количество точек
C*		X(N) - x-координаты точек
C*		Y(N) - y-координаты точек
C*		T(N) - t-координаты точек (выход)
C*
C**********************************************************************
      SUBROUTINE Parametrize(N, X, Y, T)
	  IMPLICIT NONE
	  INTEGER*2, INTENT(IN) :: N
	  REAL, INTENT(IN)    :: X(N), Y(N)
	  REAL, INTENT(OUT)   :: T(N)

	  INTEGER i

	  T(1) = 0.
	  DO i = 2, N
	    T(i) = T(i-1) + sqrt(
     *	  (x(i) - x(i - 1))**2 + (y(i) - y(i - 1))**2)
	  END DO

	END