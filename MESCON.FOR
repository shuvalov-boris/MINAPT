$NOTSTRICT
$STORAGE:2
      subroutine mescon(ms,ln,name,lnm,intg,ni,rl,nr)
C********************************************************************** 
C*
C*
C*      ����� ��������� �� ����� � ���� ������ � ������
C*
C*
C* ���������: ms(ln) - ��������� ���������
C*               �������� ������
C*            name(ln) - ���������� ����������
C*            intg(ni) - ����� �����
C*            rl(nr) - �������� �����
C*
C*            irc - ���������� ��� ��������:
C*            irc=0 - ���������� ����������
C*            irc=10000 - ����� ��������� ����� ��� ������ �����
C*            0<irc<10000 - �������� �������� iostat
C*
C* ����������: 1 ��� ������� ����� �������� ��� �����, ���������� ���
C*             2 ������������������ ����� � ������: ���������� �������-
C*             ���, ����� �����, �������� �����
C*
C*
C********************************************************************* 
      character*1 ms(ln),name*(*),line80*80,line(80)
      integer intg(ni)
      real*8 rl(nr)
      equivalence(line80,line)
      irc=0
c*         ����������� ���� ��������� read
      nrd=0
C*       ���� �� ������ � name ?
      if(lnm.gt.0) nrd=1
C*       ���� �� ������ � intg ?
      if(ni.gt.0) nrd=nrd+2
C*       ���� �� ������ � rl ?
      if(nr.gt.0) nrd=nrd+4
      nrd=nrd+1
c
3     write(*,5) ms(:ln)
5     format(1x, 77a)
      return
C*       ���������� ���������� � ������
      go to (10,20,30,40,50,60,70,80),nrd
C*    
10    return
C*    
C*    ����: read(*,25,err=100,iostat=irc,end=100) 
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
C*       ������ ������ � �������
90    if(irc.ne.0) go to 102
      return
C*      ����� ����� ��� ������ �����     
100   irc=10000
102   write(*,105) irc
105   format(3x,'������������ ����. irc=',i6)
      go to 3
      end