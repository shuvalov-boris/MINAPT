c*********************************************************************
c*
c*
c*	Сохранение оптимальных значений целевой функции на 
c*		внешнем носителе в символьной форме
c*
c* Паpаметpы:   irc-код возвpата: 0-ноpмальное завеpшение,
c*			1-ошибка пpи откpытии, или пpи записи
c*
c* Примечание. Оптимальные значения выбираются из common-областей 
c*		оптимизационного комплекса
c*
c*
c*********************************************************************
	subroutine savmui(x, fx)
	common/mfcmin/ kod1(8),iw,nnow,n,m,is,maxit,maxfx,ngrad,fxopt,
     *  ngrd,ityp1,istps,kd4,kd5
	common /mfcom/ xopt(500)
	irc=0
	open (19,FILE='optxfx.sav',ERR=30)
	write(19,10,ERR=30) n,fxopt
10	format(i6,e13.6)
	write(19,20,ERR=30) (xopt(i),i=1,n)
20	format(6e13.6)
	close(19)
	return
30	irc=1
	return
	end
