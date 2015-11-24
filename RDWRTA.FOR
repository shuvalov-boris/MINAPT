$NOTSTRICT
c*********************************************************************
c*
c*
c*          Чтение и запись данных в архив. Архив считывается или 
c*    записывается полностью.
c*
c* Параметры:   name - имя архива (вместе с директорией)
c*          icninf - массив информации по управлению счетом
c*          istinf - массив установочной информации
c*          lentyp - типы и начала блоков в массиве idat
c*          idat - блоки данных
c*          iflg - флажек:
c*                0 - читать
c*                1 - писать
c*                2 - писать с заведением архива
c*          irc - код возврата:
c*                0 - нормальное завершение
c*                1 - при iflg=2 архив с именем name уже
c*                 существует. запись не производится
c*                2 - при iflg<2 архива нет или ошибка при его 
c*                открытии
c*                3 - ошибка чтения или записи. работа 
c*                с архивом прекращается
c*                4 - нет места под запись на внешнем носителе
c*
c*
c*********************************************************************
	subroutine rdwrta(name,icnin,istin,lentyp,idat,iflg,irc)

	character where*38
	real fi, lambda
      integer*2 idat(100000),idate0(7),icninf(130), istinf(40), nknot
     *  ,icnin(130), istin(40), nb, itypcn,iflgsp
	integer*4 lentyp(2,101)

	integer*2 bcio(2)
	integer*4 len,nbrw
	character*(*) name
	data maxlen/16384/
	equivalence 
     *  (icninf(5),where),(istinf(3),fi),(istinf(5),lambda),
     *  (istinf(8),idate0),(istinf(31),itypcn),
     *  (istinf(36),nknot),(istinf(39),iflgsp)
c
c*          открытие архива
c
C	call  bopen(bcio,name,0,ir)
    8 format(A, A)
      open(unit=12, file=name, err=200, action='read', status='old')
C	if(iflg.ne.2) go to 200
	if(iflg.eq.0) go to 10
	goto 190
C	if(iflg.eq.2) go to 210
C	open(unit=12, file=name, err=200, action='write', status='old')
C	go to 100
c
c*          чтение данных
10    read(12, *, ERR=91) where
2     format(A38)
      read(12, *, ERR=91) fi, lambda
      read(12, *, ERR=91) idate0
      read(12, *, ERR=91) itypcn, iflgsp
      read(12, *, ERR=91) nknot

	nb=istinf(7)
	len=(lentyp(1,nb+1)-1)*2
	irc=0
	if(len.eq.0) go to 190
c*        пересчет длины данных для архива с длиной более maxlen байт
	nbrw=len/maxlen
	lenl=len-nbrw*maxlen
	if(nbrw.le.0) go to 190
C	do 20 i=1,nbrw
C	call  bread(bcio,idat,maxlen,irc)
C	if(irc.ne.0) go to 220
C20    continue        
C50    if(lenl.le.0) go to 190
C	call  bread(bcio,idat,lenl,irc)
C	if(irc.ne.0) go to 220
C	go to 190
c
c*          запись в архив
c
100   call  bwrite(bcio,icninf,icninf(1)*2,ir)
	if(ir.eq.-1) go to 230
	if(ir.ne.0) go to 220
	call  bwrite(bcio,istinf,istinf(1)*2,ir)
	if(ir.eq.-1) go to 230
	if(ir.ne.0) go to 220
C	call  bwrite(bcio,lentyp,(istinf(2)+1)*8,ir)
	if(ir.eq.-1) go to 230
	if(ir.ne.0) go to 220
	nb=istinf(7)
	len=(lentyp(1,nb+1)-1)*2
	irc=0
	if(len.eq.0) go to 190
c*        пересчет длины данных для архива с длиной более maxlen байт
	nbrw=len/maxlen
	lenl=len-nbrw*maxlen
	if(nbrw.le.0) go to 150
	do 120 i=1,nbrw
	call  bwrite(bcio,idat,maxlen,irc)
	if(irc.eq.-1) go to 230
	if(irc.ne.0) go to 220
120   continue        
150   if(lenl.le.0) go to 190
	call  bwrite(bcio,idat,lenl,irc)
	if(irc.eq.-1) go to 230
	if(irc.ne.0) go to 220
C   call  bclose(bcio,ir)
190     close(12, err=92)
      istin=istinf
	icnin=icninf
	return
200   irc=2
	return
210   irc=1
	return
220   irc=3
	return
230   irc=4
	return
   91 write(*, *) 'RDWRTA: Ошибка чтения данных'
      return
   92 write(*, *)	'RDWRTA: Ошибка закрытия файла'
	end