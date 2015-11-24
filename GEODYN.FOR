c*              G E O D Y N . F O R
c*              стуктуpа аpхива
      character	NAMARC*12, RESFILE*12
      integer VMAX
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
      integer*2	  idat(LENMAX), nb, nbvel
	real		 dat(LENMAX/2)

	character where*38,img*6,from*30,namfun*8
	real fi, lambda, x0, xr, eps, epsil, paspim(4,2)
	integer*2 idate0(7), nknot, itypcn, iflgsp, iflgtp, iflgvr, 
     *  iflgtd, icninf(130), istinf(40) 
	
	equivalence (icninf(1),nbgst),(icninf(2),nbgvi),
     *  (icninf(5),where),(icninf(24),img),(istinf,nbgdat),
     *  (istinf(3),fi),(istinf(5),lambda),
     *  (istinf(8),idate0),(istinf(15),paspim),(istinf(31),itypcn),
     *  (icninf(81),iretcd),(icninf(27),from),(icninf(43),eps),
     *  (istinf(32),namfun),(istinf(36),nk),
     *  (istinf(39),iflgsp)

      equivalence (dat,idat)

c
C     dat real (50000) : idat int (100000)
C     
C     icninf: (130) integer
C        (1)     nbgst    -
C        (2)     nbgvi    -
C        (3)     lenmax int*4     ограничение на массив idat
C        (5)     where  char*38   for INPSFD директория, где находятся считываемые файлы
C       (24)     img    char*6    -
C       (27)     from   char*30   -
C       (43)     eps              -
C       (81)     iretcd           -
C     
C     istinf:  (40) integer
C        (2)     nbmax               ограничение на массив lentyp
C        (3)     fi     real         широта
C        (5)     lambda real         долгота
C        (7)     nb                  число блоков данных уже имеющихся
C        (8)     idate0 (7) int      дата начала отсчета времени
C       (15)     paspim (4, 2) real  паспорт изображения
C       (31)     itypcn              0: координаты в 10км; 1 - в градусах
C       (32)     namfun char*8       -
C       (36)     nbminf	int*2	     количество узлов (nknot <= 28)
C       (37)     nbspap              0: узлы необходимо ввести
C       (38)     nbcont              =1? #TODO auto
C       (39)     iflgsp              =0