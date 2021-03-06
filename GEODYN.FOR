c*              G E O D Y N . F O R
c*              �����p� �p娢�
      character	NAMARC*12, RESFILE*12
      integer VMAX
      parameter (
	*  NMAX   =    200,	! ���ᨬ��쭮� ������⢮ �祪 ������
     *  NKMAX  =     28,	! ���ᨬ��쭮� ������⢮ 㧫�� ᯫ����
     *  VMAX   =   1000,	! ���ᨬ��쭮� ������⢮ ��થ஢ ᪮��⥩
     *  LENMAX = 100000,	! ࠧ��� ���ᨢ� idat
     *  NBMAX  =    100,  ! ���ᨬ��쭮� ������⢮ ������
     *  NAMARC  = 'geodyn.dat ',	! ��� 䠩�� � �室�묨 ����묨
     *  RESFILE = 'gd_res.dat '	! ��� 䠩�� � १���⠬� ���⮢
	*)

	integer*4 lentyp(2,NBMAX+1)	 ! ��砫� ����� � ���ᨢ� idat � ⨯ �����
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
C        (3)     lenmax int*4     ��࠭�祭�� �� ���ᨢ idat
C        (5)     where  char*38   for INPSFD ��४���, ��� ��室���� ���뢠��� 䠩��
C       (24)     img    char*6    -
C       (27)     from   char*30   -
C       (43)     eps              -
C       (81)     iretcd           -
C     
C     istinf:  (40) integer
C        (2)     nbmax               ��࠭�祭�� �� ���ᨢ lentyp
C        (3)     fi     real         ���
C        (5)     lambda real         ������
C        (7)     nb                  �᫮ ������ ������ 㦥 ��������
C        (8)     idate0 (7) int      ��� ��砫� ����� �६���
C       (15)     paspim (4, 2) real  ��ᯮ�� ����ࠦ����
C       (31)     itypcn              0: ���न���� � 10��; 1 - � �ࠤ���
C       (32)     namfun char*8       -
C       (36)     nbminf	int*2	     ������⢮ 㧫�� (nknot <= 28)
C       (37)     nbspap              0: 㧫� ����室��� �����
C       (38)     nbcont              =1? #TODO auto
C       (39)     iflgsp              =0