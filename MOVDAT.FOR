c*********************************************************************
c*
c*
c*		сдвиг части массива данных idat
c*
c* параметры:   idat,lentyp,lenmax,nb - описаны в программе inpsfd
c*		nbin - блок, в котором производятся изменения
c*		ibeg - начало изменяемого фрагмента
c*		nwd - число слов (<0 - удаляемых,>0 - добывляемых)
c*		irc : 0 - нормальное завершение
c*		      1 - нет места для добавляемых данных
c*
c* Примечание. Если число слов, на которое происходит уменьшение
c*	массива, равно его длине, то блок удаляется полностью и
c*	число блоков становится на 1 меньше
c*
c*
c*********************************************************************
	subroutine movdat(idat,lentyp,lenmax,nb,nbin,ibeg,nwd,irc)
	integer*2 idat(1), nb, nbin, irc
	integer*4 lentyp(2,1),ibeg,lastwd,ibg,lenmax, nwd
	irc=0
	lastwd=lentyp(1,nb+1)
	if(nwd.gt.0) go to 50
c
c*		удаляем данные из массива
c
	nbo=nb
	if(lentyp(1,nbin)-lentyp(1,nbin+1).ne.nwd) go to 20
c*		удаляем весь блок
	nb=nb-1
	if(nb.lt.nbin) return
	do 10 i=nbin,nb+1
	lentyp(2,i)=lentyp(2,i+1)
10	lentyp(1,i)=lentyp(1,i+1)+nwd
c*		сдвигаем массив данных
20	ibg=ibeg-nwd
	do 30 i=ibg,lastwd-1
30	idat(i+nwd)=idat(i)
	if(nbo.ne.nb) return
	go to 70
c
c*		раздвигаем массив
c
50	if(lastwd+nwd-1.gt.lenmax) irc=1
	if(lastwd+nwd-1.gt.lenmax) return
	if(ibeg.eq.lastwd) go to 70 
	do 60 i=1,lastwd-ibeg
60	idat(lastwd+nwd-i)=idat(lastwd-i)
c*		коррекщия длин блоков
70	do 80 i=nbin+1,nb+1
80	lentyp(1,i)=lentyp(1,i)+nwd
	return
	end