c*********************************************************************
c*
c*
c*  ������ ������� ����, ������������ ������, �������� � ������������� 
c*      (��� �������, ������� � ���������� ������� ������������ �������)
c*      ��� ������� � ���������� �������� ����:
c*                    { 1-alpha*d**2       0<d<R
c*  hight=F(d)=A*f(d)={
c*                    {  beta*(R+P-d)**2   R+P>d>R
c*                    {  0                   d>R+P
c*       �������: �����, ����� �������� �������� ��������, ��� �����
c*              ������������� �����
c*              
c*
c* ���������:   fi,lambda - ������ �������������� ������� ���������
c*              prm(np) - ��������� ������������������ ��������
c*              (1) - A0-������������ ������� ������� ���� (���������)
c*              (2,3) - x0,y0 - ���������� ��������� �������:
c*    ��� ������� - ��� �����, ��� ������� y0=0, � x0 - ����������
c*    �� ����� ������ �� ��� ������ �������
c*              (4,5) - ��������� a � b: ����������� ������ 
c*      R=sqrt(a**2+b**2). ��� ����� (a/r,b/r) - ������ �����������
c*      ����� ������� �������. ���� ������ �� �����, �� b=0
c*              (6) - �������� ��� ������� �������������� e
c*              (7) - ���������� �� ����� ������������ �������� ������
c*                    �� ��� "��������" �������
c*              dop,idop - �������� � ����� ���������� �����������������
c*                              �������
c*              idop(1)>0 - ����� ����� ������� (������ - �����), dop - 
c*                      ���� ������� � ��� ������������
c*              idop(1)=0 - ������ - ����� ������������� �����
c*              xy - �����, � ������� ���������� ���������� ����������
c*                      (�������� � ��������� ����������� t)
c*              dist - ���������� �� ����� xy �� ������ ������� (�����)
c*                ��� ����������� �� ����� xy �� ����� ������ +x0 (�����)
c*              hight - ������� ���� ��� ������������ ������ � ����� xy � � 
c*                 (�������,  ����������)
c*              v - �������� � ����� xy � ��/��� (������)
c*              psi - ������������� � 1/��� (������ �����������)
c*              ifl=cab, ab - ����������� ����� - ��������� �������������
c*                          �������� � ���� ������� ����������:
c*                 b    b=1,3,5,7 - ������� hight
c*                      b=2,3,6,7 - ������� v
c*                      b=4,5,6,7 - ������� psi
c*                      b=8 - hight - ���������� �� ����� �� "������"
c*                      ������� � ������ ������ � "����������������" d 
c*                      � ������ �������������� �����
c*                 a    a=0 - xy - ��� dist - ���� ����������
c*                      a=1,3 - ������ (x,y)
c*                      a=2,4 - ������ (x,y,t)
c*                      a>2 - ������� ������������ ����������, 
c*                        <3 - ������� ����
c*                 c=0 - ������� �������� � ������������� �� �������
c*              ����; 1 - ������� �������������� �������
c*              eps - �������� �������� ������ (� 10��)
c*              irc �� ������������ ����
c*
c* ����������. ��������� prm �� ������� ������� ����������� ���
c*      �������������� ������������ ������� �������� �����������
c*             ��� ������� ������������ ������ � �������� ������� ����
c*      ����������� ��������: 1-� � 2-� - ��������� ��������������
c*      ����� U**2*vr/g (��� ������� ��������������, ��� ����� �������
c*      ����� �� �����, ���������� ����� ������ (a,b) � ����� �����), 
c*      3-� - ������ ��� �������, ����� ����� �� �� �������������� �����
c*      (�������������� �� ������ ������ ���������� ������� ��������
c*      ����� ���������� �������� ����)
c*
c*
c*********************************************************************
	subroutine clcr1(fi,lambda,prm,np,dop,idop,xy,dist,hight,v,psi,
     *ifl,eps,irc)
	integer idop(1)
	real prm(np),dop(1),xy(1),v(2),lambda,tgdxdy(5)
	data W/7.29e-5/,pi/3.141593/,g/9.81/
	iflg=ifl-ifl/10*10
	inpf=ifl/10-ifl/100*10
	iflcon=ifl/100-ifl/1000*100
c
c*              ������ ���������� F(r)
c
	f=2.*W*sin(lambda*pi/180.)
	constv=g/f/100.
	nknot=idop(1)
	a0=prm(1)
	x0=prm(2)
	y0=prm(3)
	a=prm(4)
	b=prm(5)
	e=1./(1.+prm(6)**2)
	r2=a*a+b*b
	r=sqrt(r2)
	if(inpf.ne.2) go to 10
c*              ������������ t-����������
	t=xy(3)
10      p=prm(7)**2
	alpha=1./r/(r+p)
	beta=1./p/(r+p)
c        print *,'a0,r,p,alpha,beta,e=',a0,r,p,alpha,beta,e
c*              ������ dist (�.�. d)
	dx=xy(1)-x0
	dy=xy(2)-y0
	dist=dx
	if(inpf.ne.0.and.nknot.eq.0) dist=sqrt((dx*a+dy*b)**2+e*(-dx*b
     *  +dy*a)**2)/r
	if(nknot.le.0.or.inpf.eq.0) go to 20
c*              ������� ���������� �� �������
	iflsp=0
	if(iflg.le.8.and.iflg.gt.1) iflsp=6
	call  trdist(fi,lambda,xy(1),xy(2),dop,dop(nknot+1),
     *  dop(nknot*2+1),nknot,dop(nknot*3+1),nknot-1,xp,yp,tp,tgdxdy,vr,
     *  dist,eps,iflsp)
c*              ���������� ������ �� ����� ��� �������� ����� ������-
c*        ������ (t ������) ������������� (�������� ��������� trdist)
	dist=dist+x0
20      if(iflg.eq.8) hight=dist
	if(iflg.eq.8) return
c
	if(iflg/2-iflg/4*2.eq.0.and.inpf.eq.0) go to 40
	if(iflg.eq.1.and.inpf.le.2) go to 40
c*              ������� v
	v1=0.
	if(iflg/2-iflg/4*2.eq.1) v(1)=0.
	v2=0.
	if(iflg/2-iflg/4*2.eq.1.and.inpf.ne.0) v(2)=0.
	vabs=0.
c        print *,'dist=',dist
	if(dist.le.0..or.dist.gt.p+r) go to 40
c*          vabs=-A*f'(d) �� d (-Df/Dd)
	if(dist.le.r) vabs=2.*a0*alpha*dist
	if(dist.gt.r) vabs=2.*a0*beta*(r+p-dist)
	if(iflcon.eq.0) vabs=vabs*constv
	v1=vabs
	if(iflg/2-iflg/4*2.eq.1) v(1)=v1
	if(inpf.eq.0) go to 40
	if((iflg.eq.1.or.iflg.eq.5).and.inpf.le.2) go to 40
c*              ���������� ������� ��� ���������� ��������
	if(nknot.gt.0) go to 35
c*              ������� d' �� x � y (�������� ������� ���������)    
	ycomp=dx*a*b*(1.-e)+dy*(b*b+e*a*a)
	xcomp=dx*(a*a+e*b*b)+dy*a*b*(1.-e)
	ry=ycomp/(dist*r*r)
	v1=vabs*ry
	if(iflg/2-iflg/4*2.eq.1) v(1)=vabs*ry
	rx=xcomp/(dist*r*r)
	v2=-vabs*rx
	if(iflg/2-iflg/4*2.eq.1) v(2)=-vabs*rx
	go to 40       
c*          ������ �������� ��� ������ (���������� ������� ���������)
35    v2=vabs*vabs
	tga=tgdxdy(1)
	tga2=tga*tga
	v1=sqrt(v2/(1.+tga2))
	if(tgdxdy(2).lt.0.) v1=-v1
	if(iflg/2-iflg/4*2.eq.1) v(1)=v1
	v2=sqrt(v2*tga2/(1.+tga2))
	if(tgdxdy(3).lt.0.) v2=-v2
	if(iflg/2-iflg/4*2.eq.1) v(2)=v2
c
40      if(iflg-iflg/2*2.eq.0) go to 45
c*              ������� hight
	hight=a0
	if(dist.le.0) go to 41
	hight=0.
	if(dist.gt.p+r) go to 45
	if(dist.le.r) hight=a0*(1.-alpha*dist*dist)
	if(dist.gt.r) hight=a0*beta*(r+p-dist)**2
c
41      if(inpf.le.2) go to 45
c*         ��������� ������� ������������ ����������, � �� ������� ����
c*              R+P >=dist>=R. ������� �������� � ������
	dr=dist
	if(dist.le.r) dr=r
	cv2=constv*constv
	sa0=1.
c        print*,'dist,r,hight,vabs=',dist,r,hight,vabs
	if(a0.lt.0.) sa0=-1.
	if(nknot.eq.0) hight=hight+4.*a0*a0*beta*beta*sa0/g*
     *((r+p)**2*alog(dr/(r+p))+2.*(r+p)*(r+p-dr)+dr*dr/2.-(r+p)**2/2.)*
     * cv2*1.e-4
	if(nknot.gt.0) hight=hight+4./3.*vr*a0*a0*beta*beta/g*
     *(r+p-dr)**3*cv2
c        print *,hight
	if(dist.ge.r) go to 43
c*              dist <r
	dr=dist
	if(nknot.eq.0) hight=hight-2.*a0*a0*alpha*alpha*sa0/g*
     * (r*r-dr*dr)*cv2*1.e-4
	if(nknot.gt.0) hight=hight+4./3.*vr*a0*a0*alpha*alpha/g*(r**3-
     *dr**3)*cv2
c        print *,hight
c*              3-� �������� ��� �������������� �����
43      if(nknot.eq.0) hight=hight+vabs*vabs*(1.-rx*rx-ry*ry)/g*0.5e-4
c        print *,hight
c
45      if(iflg/4-iflg/8*4.eq.0) return
c*              ������� psi
	psi=0.
	if(dist.le.0..or.dist.gt.p+r) return
c*        ������ ������������� ������ (���������� ������ � �����)
	psi=-2.*alpha
	if(dist.gt.r) psi=2.*beta
	if(iflcon.eq.0) psi=psi*constv*1.e-4
	if(inpf.eq.0) return
	if(nknot.eq.0) go to 50
c*              ������ ������������� �������� ����� (��� ������)
	psir=vabs*vr
	if(iflcon.eq.0) psir=psir*1.e-4
	psi=psi+psir
	return
c*              ������ ������������� ��� �������������� �����
50    rxx=(a*a+e*b*b)/(dist*r*r)-rx*rx*r
	ryy=(b*b+e*a*a)/(dist*r*r)-ry*ry*r
	if(iflcon.eq.0) vabs=vabs*1.e-4
	psi=psi*rx*rx-vabs*rxx+psi*ry*ry-vabs*ryy
c
	return
	end
