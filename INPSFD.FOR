C**********************************************************************
C*
C*
C*                   ПРОГРАММА ввода контурно-векторной
C*        и гидрологической информации 
C*
C*
C* Параметры: fi,lambda - начало географической системы координат (дол-
C*                       гота, шиpота)
C*            idate0(7) - дата начала отсчета времени: год, месяц, день,
C*                       час, минута, секунда, тик
C*            itypcn - :0 - координаты в 10км, 1 - в градусах
C*            where - директория, где находятся считываемые файлы
C*            nb - число блоков данных уже имеющихся (на входе), на 
C*                       выходе - общее число блоков
C*            lentyp(2,nbmax+1) - начала блоков в массиве idat(dat) и
C*                       типы блоков
C*            nbmax,lenmax - ограничения на массивы lentyp и idat
C*            paspim - паспоpт изобpажения (см. коммент. пpог. rdconv)
C*            maxrec - максимальное число записей в считываемом блокe
C*            metpos(maxrec+2),xy(2*maxrec+2),v(2*maxrec+2),wk(2*maxrec+
C*                      2) - рабочие массивы
C*            irc - код возврата:
C*                      0 - нормальное завершение
C*                      1 - были ошибки при вводе некоторых блоков (не
C*                  все блоки ввелись)
C*                      2 - превышены пределы nbmax или lenmax
C*                      3 - ничего не было введено
C*
C*
C* Алгоритм. Файлы с информацией должны быть помещены в директорию 
C*   where, не совпадающей с текущей. Считанные блоки данных 
C*   имеют в массиве dat (idat) следующую структуру (обозначения: i -
C*   integer*2, r - real*4):
C*           Тип данных 0 - скорости
C*     Справка блока:
C*   Время time (r) в час.; itypcn (i); число точек (i); интервал между сним-
C*   ками в мин; тип скоpости (i): 0-исходные данные,>0-обpаботанные:1-данных о 
C*   наклоне и кpивизне нет, 2-есть
C*     Данные(r): x-координаты точки; y-коор.; x-коор. скоростей; y-коор
C*   скоpостей; тангенсы наклонов контуpа и кpивизна в точках контуpа,
C*   ближайших к точкам (x,y)
C*           Тип данных 1 - точки фрагментов контура
C*     Справка блока:
C*   Время time (r) в час.; itypcn (i); число точек (i); число
C*   фрагментов контура (i); резерв (i)
C*     Данные(r): начала фpагментов (i); x-координаты точки; y-коор.
C*           Тип данных 2 - узлы сплайна.
C*     Справка блока:
C*   Резерв (2i); itypcn (i); число узлов (i); резерв (i); номер вариан
C*   та(i):0-сплайн-аппроксимация не проводилась
C*     Данные(r): x-координаты узлов; y-коор.; t-коор.; коэффициенты
C*   сплайна
C*           Тип данных 3 - гидрологическая информация
C*     Справка блока:
C*   Номер разреза (i); 1-кривизна считалась,0-нет (i); itypcn (i); 
C*   число точек (i); резерв (2i)
C*     Данные(r): x-координаты точки; y-коор.; время пост. станций в ч.;
C*   топография в м; радиусы кривизны в 10км; номера станций (i)
C*           Тип данных 4 - оптимизационная информация
C*     Справка блока:
C*   Время начала оптимизации (r);0-имя целевой функции fun,1-нет (i);
C*    n-число минимизируемых переменных (i); номер блока,параметры кото-
C*   рого оптимизируются (i); резерв (i)
C*     Данные: (r): parm(8),pr2,pr7; (i): kod(8),kod4,kod5,sizea,sizeb,
C*   kodgrd (порядок метода оптимизации),nicp (интервал контрольной
C*   точки); (r):d.dm,ds; (ch*1): nammet(8),namout(8),nammnl(8); 
C*   (i): тип целевой функции,тип оптимизационного критерия;nbopt(50)-
C*   номера блоков с данными для оптимизации;
C*   (r): fxbest,x(n), (i):statx(n)
C*
C* Примечания: 
C*     ЕСЛИ ВРЕМЯ В ФАЙЛЕ РАВНО .123, ТО ЭТО ОЗНАЧАЕТ, ЧТО ОНО НЕ 
C*  ВВОДИЛОСЬ
C*
C*
C**********************************************************************
	subroutine inpsfd(fi,lambda,idate0,itypcn,where,nb,lentyp,dat,
     *idat,nbmax,lenmax,paspim,maxrec,metpos,xy,v,wk,irc)
	integer*2 idate0(7),idat(1),metpos(1), nb
	integer*4 lentyp(2,1),lenmax,nbg
	real dat(1),fi,lambda,paspim(4,2),xy(1),wk(1),v(1),time(2)
	CHARACTER LINE1(80),LINE80*80,FL*80,FLI(80),where*(*),null_char*1
     *  , LTMP(180), LTMPI*180
	EQUIVALENCE (LINE1,LINE80),(FLI,FL), (LTMPI, LTMP)
	if(nb.eq.0) lentyp(1,1)=1
	nbs=nb
C*          ПОСТРОЕНИЕ ТАБЛИЦЫ ИМЕН
	lnw=LEN_TRIM(where)
	null_char=8#0
	LINE80='DIR '//where(:lnw)//'\*.* >tabl.nam'//null_char
	CALL  KMON(LINE80)
	irc=0
	iflg=1
	LINE80='TABL.NAM'
	LINE1(9)=8#0 
	OPEN(9,ERR=100,FILE=LINE80,STATUS='OLD')
	IFLG=0
C*      ПРОПУСК семи СТРОК
	 DO 5 I=1,7
	 READ(9,4)
4      FORMAT(80A1)
 5     CONTINUE
	 iflrc=1
c
c*          Цикл по именам массивов
c
c*    Удаляем пpедыдущий считанных файл, если он считан без ошибок
9     if(iflrc.eq.0) call  kmon('del '//fl(:ind+1))
	READ (9,6) (LTMP(I), I=1,180)
    6 format(180A1)
C	 IF(Ltmp(1).EQ.' ') CLOSE(9)
	 IF(Ltmp(1).EQ.' ') then
	   CLOSE(9)
	   CALL KMON('DEL TABL.NAM'//null_char)
	   go to 100
	 endif
C	 IF(Ltmp(1).EQ.' ') go to 100
      do i=1,80
	  line1(i)=ltmp(i+36)
	enddo
	lnam=LEN_TRIM(line80)

C*           СОЗДАЕМ ИМЯ ФАЙЛА
C*Было  FL=where(:lnw)//'\'
	 FL=where(:lnw)//'/'
	 ind=lnw+1
C	 DO 10 I=1,8
C	 IF(LINE1(I).EQ.' ') GO TO 11
C	 FLI(I+ind)=LINE1(I)
C10     CONTINUE
C	 i=9
C11     ind=ind+i
      do i=1,lnam
	  fli(i+ind)=line1(i)
	enddo
      ind=ind+lnam

C	 if(line1(10).eq.' ') go to 25
C	 FLI(ind)='.'
C	 DO 20 J=1,3
C	 FLI(ind+J)=LINE1(J+9)
C20     CONTINUE
C	 ind=ind+4
C
c*              чтение блока данных
25     fli(ind+1)=8#0
	 iflrc=0
	 iflag=0
	 nbf=0
30     if(nb.ge.nbmax) go to 105
	 call  rdconv(fl,paspim,fi,lambda,idate0,xy,v,n,metpos,time,wk,
     * iflag,itypcn,maxrec,ir)
	 nbf=nbf+1
	 if(ir.le.3) go to 50
c*              считывание завершилось аварийно
	 if(ir.ne.5) irc=1
	 if(ir.eq.4) print *,'Файл ',(fli(i),i=1,ind),' не найден'
	 if(ir.eq.4) iflrc=1
	 if(ir.eq.4) go to 9
	 if(ir.eq.6) print *,'Блок ',nbf,' файла ',(fli(i),i=1,ind),
     * ' большеват'
	 if(ir.eq.7) print *,'В блоке ',nbf,' файла ',(fli(i),i=1,ind),
     * ' отсутствуют даты'
	 if(ir.ge.6) iflrc=1
	 if(ir.ge.6) go to 9
	 go to 80
c*              блок данных успешно считан
50     nf=metpos(1)
	 nf1=nf
	 if(nf-nf/2*2.eq.1) nf1=nf1+1
	 lentyp(1,nb+2)=4*n+lentyp(1,nb+1)+6
	 if(ir-ir/2*2.eq.0) lentyp(1,nb+2)=lentyp(1,nb+2)+nf1
	 if(ir-ir/2*2.eq.1) lentyp(1,nb+2)=lentyp(1,nb+2)+4*n
	 if(lenmax.lt.lentyp(1,nb+2)-1) then
	   print*, ' Не хватает места в массиве dat (line196)'
	   go to 107
	 endif
c
c*              занасение данных в массив dat
c
	nb=nb+1
	lentyp(2,nb)=1
	if(ir-ir/2*2.eq.1) lentyp(2,nb)=0
	if(ir-ir/2*2.eq.0) write(6,51) (fli(i),i=1,ind)
51    format(' считан блок контурной информации из файла ',80a1)
	if(ir-ir/2*2.eq.1) write(6,52) (fli(i),i=1,ind)
52    format(' считан блок векторов из файла ',80a1)
	write(6,522) nb
522   format(' блоку пpисвоен номеp ',i4)
	nbg=lentyp(1,nb)
c*              занесение информации в справку файла
	dat(nbg/2+1)=time(1)
	idat(nbg+2)=itypcn
	idat(nbg+3)=n
	idat(nbg+4)=0
	if(ir-ir/2*2.eq.1) idat(nbg+4)=time(2)
	inbg=nbg
	idat(nbg+5)=0
	nbg=nbg+6
	if(lentyp(2,nb).eq.0) go to 55
c*              запись числа фрагментов контура и их начала
	idat(nbg-2)=nf
	do 53 i=1,nf
53    idat(nbg+i-1)=metpos(i+1)
	nbg=nbg+nf1
c*              заносим координаты точек
55    do 60 i=1,n
       dat(nbg/2+i)=xy(i)
60     dat(nbg/2+i+n)=xy(i+n)
	nbg=nbg+n*4
	if(lentyp(2,nb).eq.1) go to 70
c*              заносим компоненты векторов
	do 65 i=1,n
      dat(nbg/2+i)=v(i)
65    dat(nbg/2+i+n)=v(i+n)
c*              проверка окончания файла
70    iflag=1
	if(ir.gt.1) go to 9
	go to 30
c
c*              считаваем гидрологическую информацию
c
80     iflag=0
	 nbf=0
88     if(nb.ge.nbmax) go to 105
	 call  rdginf(fl,fi,lambda,idate0,xy,v,n,metpos,wk,iflag,
     * itypcn,maxrec,ir)
	 nbf=nbf+1
	 if(ir.le.3) go to 81
c*              считывание завершилось аварийно
	 irc=1
	 if(ir.eq.4) print *,'Файл ',(fli(i),i=1,ind),' не найден'
	 if(ir.eq.5) print *,'Ошибка при чтении файла ',(fli(i),i=1,ind)
	 if(ir.eq.6) print *,'Блок ',nbf,' файла ',(fli(i),i=1,ind),
     * 'большеват'
	 iflrc=1
	 go to 9
c*              блок данных успешно считан
81     lentyp(1,nb+2)=11*n+lentyp(1,nb+1)+6+n-n/2*2
	 if(lenmax.lt.lentyp(1,nb+2)-1) go to 107
c
c*              занасение данных в массив dat
c
	nb=nb+1
	lentyp(2,nb)=3
	write(6,82) (fli(i),i=1,ind)
82    format(' считан блок гидрологической информации из файла ',36a1)
	write(6,822) nb
822   format(' блоку пpисвоен номеp ',i4)
	nbg=lentyp(1,nb)
      write(6, 8) ' INPSFD: nbg = ', nbg
    8 format(A, I7)
c*              занесение информации в справку файла
	idat(nbg)=metpos(1)
	idat(nbg+1)=0
	idat(nbg+2)=itypcn
	idat(nbg+3)=n
	idat(nbg+4)=0
	idat(nbg+5)=0
	nbg=nbg+6
c*              заносим координаты точек
	do 83 i=1,n
	dat(nbg/2+i)=xy(i)
83    dat(nbg/2+i+n)=xy(i+n)
	nbg=nbg+n*4
c*              заносим время и динамическую топографию
	do 84 i=1,n
	dat(nbg/2+i)=v(i)
84    dat(nbg/2+i+n)=v(i+n)
	nbg=nbg+n*6
c*              заносим номера станций
	do 85 i=1,n
85    idat(nbg+i-1)=metpos(i+1)
c*              проверка окончания файла
	iflag=1
	if(ir.gt.1) go to 9
	go to 88
c
C*         ОБРАБОТКА ОШИБОК
100    IF(IFLG.EQ.1) STOP 'ФАЙЛ TABL.NAM НЕ ОТКРЫТ'
	 if(nbs.eq.nb) irc=3
	 if(irc.eq.3) print *,'Ничего не считано'
	 go to 110
105    print *,'Число блоков больше максимального'
107    print *,'Мал массив под данные'
	 irc=2
110   return
	 END
