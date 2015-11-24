$NOTSTRICT
$STORAGE:2
      subroutine mescon(ms,ln,name,lnm,intg,ni,rl,nr)
C********************************************************************** 
C*
C*
C*      вывод сообщения на экран и ввод данных с экрана
C*
C*
C* параметры: ms(ln) - текстовое сообщение
C*               вводимые данные
C*            name(ln) - символьная информация
C*            intg(ni) - целые числа
C*            rl(nr) - реальные числа
C*
C*            irc - внутренний код возврата:
C*            irc=0 - нормальное завершение
C*            irc=10000 - конец экранного файла или ошибка ввода
C*            0<irc<10000 - беличина парамета iostat
C*
C* примечания: 1 при нулевой длине массивов для ввода, считывания нет
C*             2 последовательность ввода с экрана: символьная информа-
C*             ция, целые числа, реальные числа
C*
C*
C********************************************************************* 
      character*1 ms(ln),name*(*),line80*80,line(80)
      integer intg(ni)
      real*8 rl(nr)
      equivalence(line80,line)
      irc=0
c*         определение вида оператора read
      nrd=0
C*       надо ли читать в name ?
      if(lnm.gt.0) nrd=1
C*       надо ли читать в intg ?
      if(ni.gt.0) nrd=nrd+2
C*       надо ли читать в rl ?
      if(nr.gt.0) nrd=nrd+4
      nrd=nrd+1
c
3     write(*,5) ms(:ln)
5     format(1x, 77a)
      return
C*       считывание информации с экрана
      go to (10,20,30,40,50,60,70,80),nrd
C*    
10    return
C*    
C*    Было: read(*,25,err=100,iostat=irc,end=100) 
20    read(*,25,err=100,end=100) (line(i),i=1,lnm)
25    format(80a1)
      name(:lnm)=line80(:lnm)
      go to 90
C*    
30    read(*,*,err=100,end=100) intg
      go to 90
C*    
40    read(*,*,err=100,end=100) (line(i),i=1,lnm),intg
      name(:lnm)=line80(:lnm)
      go to 90
C*    
50    read(*,*,err=100,end=100) rl
      go to 90
C*    
60    read(*,*,err=100,end=100) (line(i),i=1,lnm),rl
      name(:lnm)=line80(:lnm)
      go to 90
C*    
70    read(*,*,err=100,end=100) intg,rl
      go to 90
C*    
80    read(*,*,err=100,end=100) (line(i),i=1,lnm),intg,rl
      go to 90
C*       анализ ошибок и возврат
90    if(irc.ne.0) go to 102
      return
C*      конец файла или ошибка ввода     
100   irc=10000
102   write(*,105) irc
105   format(3x,'некорректный ввод. irc=',i6)
      go to 3
      end