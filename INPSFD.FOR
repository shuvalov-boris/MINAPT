C**********************************************************************
C*
C*
C*                   ��������� ����� ��������-���������
C*        � ��������������� ���������� 
C*
C*
C* ���������: fi,lambda - ������ �������������� ������� ��������� (���-
C*                       ����, ��p���)
C*            idate0(7) - ���� ������ ������� �������: ���, �����, ����,
C*                       ���, ������, �������, ���
C*            itypcn - :0 - ���������� � 10��, 1 - � ��������
C*            where - ����������, ��� ��������� ����������� �����
C*            nb - ����� ������ ������ ��� ��������� (�� �����), �� 
C*                       ������ - ����� ����� ������
C*            lentyp(2,nbmax+1) - ������ ������ � ������� idat(dat) �
C*                       ���� ������
C*            nbmax,lenmax - ����������� �� ������� lentyp � idat
C*            paspim - �����p� ����p������ (��. �������. �p��. rdconv)
C*            maxrec - ������������ ����� ������� � ����������� ����e
C*            metpos(maxrec+2),xy(2*maxrec+2),v(2*maxrec+2),wk(2*maxrec+
C*                      2) - ������� �������
C*            irc - ��� ��������:
C*                      0 - ���������� ����������
C*                      1 - ���� ������ ��� ����� ��������� ������ (��
C*                  ��� ����� �������)
C*                      2 - ��������� ������� nbmax ��� lenmax
C*                      3 - ������ �� ���� �������
C*
C*
C* ��������. ����� � ����������� ������ ���� �������� � ���������� 
C*   where, �� ����������� � �������. ��������� ����� ������ 
C*   ����� � ������� dat (idat) ��������� ��������� (�����������: i -
C*   integer*2, r - real*4):
C*           ��� ������ 0 - ��������
C*     ������� �����:
C*   ����� time (r) � ���.; itypcn (i); ����� ����� (i); �������� ����� ����-
C*   ���� � ���; ��� ���p���� (i): 0-�������� ������,>0-��p���������:1-������ � 
C*   ������� � �p������ ���, 2-����
C*     ������(r): x-���������� �����; y-����.; x-����. ���������; y-����
C*   ���p�����; �������� �������� �����p� � �p������ � ������ �����p�,
C*   ��������� � ������ (x,y)
C*           ��� ������ 1 - ����� ���������� �������
C*     ������� �����:
C*   ����� time (r) � ���.; itypcn (i); ����� ����� (i); �����
C*   ���������� ������� (i); ������ (i)
C*     ������(r): ������ �p�������� (i); x-���������� �����; y-����.
C*           ��� ������ 2 - ���� �������.
C*     ������� �����:
C*   ������ (2i); itypcn (i); ����� ����� (i); ������ (i); ����� ������
C*   ��(i):0-������-������������� �� �����������
C*     ������(r): x-���������� �����; y-����.; t-����.; ������������
C*   �������
C*           ��� ������ 3 - ��������������� ����������
C*     ������� �����:
C*   ����� ������� (i); 1-�������� ���������,0-��� (i); itypcn (i); 
C*   ����� ����� (i); ������ (2i)
C*     ������(r): x-���������� �����; y-����.; ����� ����. ������� � �.;
C*   ���������� � �; ������� �������� � 10��; ������ ������� (i)
C*           ��� ������ 4 - ��������������� ����������
C*     ������� �����:
C*   ����� ������ ����������� (r);0-��� ������� ������� fun,1-��� (i);
C*    n-����� �������������� ���������� (i); ����� �����,��������� ����-
C*   ���� �������������� (i); ������ (i)
C*     ������: (r): parm(8),pr2,pr7; (i): kod(8),kod4,kod5,sizea,sizeb,
C*   kodgrd (������� ������ �����������),nicp (�������� �����������
C*   �����); (r):d.dm,ds; (ch*1): nammet(8),namout(8),nammnl(8); 
C*   (i): ��� ������� �������,��� ���������������� ��������;nbopt(50)-
C*   ������ ������ � ������� ��� �����������;
C*   (r): fxbest,x(n), (i):statx(n)
C*
C* ����������: 
C*     ���� ����� � ����� ����� .123, �� ��� ��������, ��� ��� �� 
C*  ���������
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
C*          ���������� ������� ����
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
C*      ������� ���� �����
	 DO 5 I=1,7
	 READ(9,4)
4      FORMAT(80A1)
 5     CONTINUE
	 iflrc=1
c
c*          ���� �� ������ ��������
c
c*    ������� �p�������� ��������� ����, ���� �� ������ ��� ������
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

C*           ������� ��� �����
C*����  FL=where(:lnw)//'\'
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
c*              ������ ����� ������
25     fli(ind+1)=8#0
	 iflrc=0
	 iflag=0
	 nbf=0
30     if(nb.ge.nbmax) go to 105
	 call  rdconv(fl,paspim,fi,lambda,idate0,xy,v,n,metpos,time,wk,
     * iflag,itypcn,maxrec,ir)
	 nbf=nbf+1
	 if(ir.le.3) go to 50
c*              ���������� ����������� ��������
	 if(ir.ne.5) irc=1
	 if(ir.eq.4) print *,'���� ',(fli(i),i=1,ind),' �� ������'
	 if(ir.eq.4) iflrc=1
	 if(ir.eq.4) go to 9
	 if(ir.eq.6) print *,'���� ',nbf,' ����� ',(fli(i),i=1,ind),
     * ' ���������'
	 if(ir.eq.7) print *,'� ����� ',nbf,' ����� ',(fli(i),i=1,ind),
     * ' ����������� ����'
	 if(ir.ge.6) iflrc=1
	 if(ir.ge.6) go to 9
	 go to 80
c*              ���� ������ ������� ������
50     nf=metpos(1)
	 nf1=nf
	 if(nf-nf/2*2.eq.1) nf1=nf1+1
	 lentyp(1,nb+2)=4*n+lentyp(1,nb+1)+6
	 if(ir-ir/2*2.eq.0) lentyp(1,nb+2)=lentyp(1,nb+2)+nf1
	 if(ir-ir/2*2.eq.1) lentyp(1,nb+2)=lentyp(1,nb+2)+4*n
	 if(lenmax.lt.lentyp(1,nb+2)-1) then
	   print*, ' �� ������� ����� � ������� dat (line196)'
	   go to 107
	 endif
c
c*              ��������� ������ � ������ dat
c
	nb=nb+1
	lentyp(2,nb)=1
	if(ir-ir/2*2.eq.1) lentyp(2,nb)=0
	if(ir-ir/2*2.eq.0) write(6,51) (fli(i),i=1,ind)
51    format(' ������ ���� ��������� ���������� �� ����� ',80a1)
	if(ir-ir/2*2.eq.1) write(6,52) (fli(i),i=1,ind)
52    format(' ������ ���� �������� �� ����� ',80a1)
	write(6,522) nb
522   format(' ����� �p������ ����p ',i4)
	nbg=lentyp(1,nb)
c*              ��������� ���������� � ������� �����
	dat(nbg/2+1)=time(1)
	idat(nbg+2)=itypcn
	idat(nbg+3)=n
	idat(nbg+4)=0
	if(ir-ir/2*2.eq.1) idat(nbg+4)=time(2)
	inbg=nbg
	idat(nbg+5)=0
	nbg=nbg+6
	if(lentyp(2,nb).eq.0) go to 55
c*              ������ ����� ���������� ������� � �� ������
	idat(nbg-2)=nf
	do 53 i=1,nf
53    idat(nbg+i-1)=metpos(i+1)
	nbg=nbg+nf1
c*              ������� ���������� �����
55    do 60 i=1,n
       dat(nbg/2+i)=xy(i)
60     dat(nbg/2+i+n)=xy(i+n)
	nbg=nbg+n*4
	if(lentyp(2,nb).eq.1) go to 70
c*              ������� ���������� ��������
	do 65 i=1,n
      dat(nbg/2+i)=v(i)
65    dat(nbg/2+i+n)=v(i+n)
c*              �������� ��������� �����
70    iflag=1
	if(ir.gt.1) go to 9
	go to 30
c
c*              ��������� ��������������� ����������
c
80     iflag=0
	 nbf=0
88     if(nb.ge.nbmax) go to 105
	 call  rdginf(fl,fi,lambda,idate0,xy,v,n,metpos,wk,iflag,
     * itypcn,maxrec,ir)
	 nbf=nbf+1
	 if(ir.le.3) go to 81
c*              ���������� ����������� ��������
	 irc=1
	 if(ir.eq.4) print *,'���� ',(fli(i),i=1,ind),' �� ������'
	 if(ir.eq.5) print *,'������ ��� ������ ����� ',(fli(i),i=1,ind)
	 if(ir.eq.6) print *,'���� ',nbf,' ����� ',(fli(i),i=1,ind),
     * '���������'
	 iflrc=1
	 go to 9
c*              ���� ������ ������� ������
81     lentyp(1,nb+2)=11*n+lentyp(1,nb+1)+6+n-n/2*2
	 if(lenmax.lt.lentyp(1,nb+2)-1) go to 107
c
c*              ��������� ������ � ������ dat
c
	nb=nb+1
	lentyp(2,nb)=3
	write(6,82) (fli(i),i=1,ind)
82    format(' ������ ���� ��������������� ���������� �� ����� ',36a1)
	write(6,822) nb
822   format(' ����� �p������ ����p ',i4)
	nbg=lentyp(1,nb)
      write(6, 8) ' INPSFD: nbg = ', nbg
    8 format(A, I7)
c*              ��������� ���������� � ������� �����
	idat(nbg)=metpos(1)
	idat(nbg+1)=0
	idat(nbg+2)=itypcn
	idat(nbg+3)=n
	idat(nbg+4)=0
	idat(nbg+5)=0
	nbg=nbg+6
c*              ������� ���������� �����
	do 83 i=1,n
	dat(nbg/2+i)=xy(i)
83    dat(nbg/2+i+n)=xy(i+n)
	nbg=nbg+n*4
c*              ������� ����� � ������������ ����������
	do 84 i=1,n
	dat(nbg/2+i)=v(i)
84    dat(nbg/2+i+n)=v(i+n)
	nbg=nbg+n*6
c*              ������� ������ �������
	do 85 i=1,n
85    idat(nbg+i-1)=metpos(i+1)
c*              �������� ��������� �����
	iflag=1
	if(ir.gt.1) go to 9
	go to 88
c
C*         ��������� ������
100    IF(IFLG.EQ.1) STOP '���� TABL.NAM �� ������'
	 if(nbs.eq.nb) irc=3
	 if(irc.eq.3) print *,'������ �� �������'
	 go to 110
105    print *,'����� ������ ������ �������������'
107    print *,'��� ������ ��� ������'
	 irc=2
110   return
	 END
