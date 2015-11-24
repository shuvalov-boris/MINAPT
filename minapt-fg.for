	if (nknot.le.1) goto 999
	if(nb+1.ge.nbmax) goto 999
	ibgfr=lentyp(1,nb+1)
	len=6+nknot*6+(nknot)*12
	if(ibgfr+len-1.gt.lenmax) goto 999
c*          A a?oeaa anou ianoi aey eiie?iaaiey
	NC=idat(ibg+3)
	idat(ibgfr)  =0
	idat(ibgfr+1)=0
	idat(ibgfr+2)=itypcn
	idat(ibgfr+3)=nknot + 1
	idat(ibgfr+4)=0
	idat(ibgfr+5)=0

	ibgfr=ibgfr+6

	dat(ibgfr/2+(nknot+1)*2) = 0.

	do i = 1, nknot
	  dat(ibgfr/2+i) = X(i + 2)
	  dat(ibgfr/2+i+nknot + 1) = X(i + nknot + 2)
	  if (i.gt.1) dat(ibgfr/2+i+(nknot+1)*2) = 
     *	  dat(ibgfr/2 + i - 1 + (nknot+1)*2) + sqrt(
	*	  (dat(ibgfr/2 + i) - dat(ibgfr/2 + i - 1))**2 +
	*	  (dat(ibgfr/2+i+nknot+1) - dat(ibgfr/2+i+nknot))**2)
	enddo

	dat(ibgfr/2 + nknot + 1) = X(3)
	dat(ibgfr/2 + 2*nknot + 2) = X(3 + nknot)
	dat(ibgfr/2+(nknot+1)*3) = 
     *	dat(ibgfr/2 + (nknot+1)*3 - 1) + sqrt(
     *	(dat(ibgfr/2 + nknot + 1) - dat(ibgfr/2 + nknot))**2 +
	*	(dat(ibgfr/2+nknot*2+2) - dat(ibgfr/2+nknot*2+1))**2)	  



      nb=nb+1
      lentyp(1,nb+1)=lentyp(1,nb)+6+(nknot+1)*6+(nknot)*12
      lentyp(2,nb)=2
	nbs=nb
	nbspap=nbs

 





















c*		nieaei aii?ieneiaoey eiioo?a
c*	  pan?ao eiinoaio
      print*, 'MINAPT: nbspap =', nbspap
20	nbgknt=lentyp(1,nbspap)
	nknot=idat(nbgknt+3)
	nbgknt=nbgknt+6
	nv=idat(nbgknt-1)
	nbgcon=lentyp(1,nbcont)
	ncon=idat(nbgcon+3)
	nf=idat(nbgcon+4)
	if(nf-nf/2*2.eq.1) nf=nf+1
	nbgcon=nbgcon+6+nf
	nbc=nbgcon/2+1
	nbk=nbgknt/2+1
	write(*,*) 'MINAPT: LENTYP:'
	write(*, *) ((lentyp(j, i), j = 1, 2), i = 1, nb)
	print*, 'MINAPT: nbc=', nbc, ',ncon=', ncon, ',nbk=', nbk,
     *  ',nknot=', nknot, ',eps=', eps, ',iflgtp=', iflgtp
      print*, 'dat(nbc) =', dat(nbc), ', t(2) =', t(2), ', dat(nbk) ='
     *  , dat(nbk), ', dat(nbk+nknot) =', dat(nbk+nknot)
	*  , ', dat(nbk+nknot*2) =', dat(nbk+nknot*2) 
c*		?an?ao t-eiipaeeiaou
	call  tprep(dat(nbc),t(2),ncon,dat(nbk),dat(nbk+nknot),dat(nbk
     *  +nknot+nknot),nknot,nbeg,nend,nusl,nusr,ndoubl,eps,wk,iwk,
     *  iflgtp,irc)
c	goto 666
	write(*,120) irc,nbeg,nend,nusl,nusr,ndoubl,eps
120	format(3x,' i?ia?aiia tprep. irc=',i3,/,3x,'nbeg,nend,nusl,nusr,
     *  ndoubl,eps=',5i5,f9.3)
	print*,'tk=',(dat(nbk+nknot*2+jj),jj=0,nknot-1)
	irctp=irc
	
c 	write(6,130) (wk(k),k=1,nknot)
c	write(6,135) (iwk(k),k=1,nknot)
c130	format(3x,'wk=',8f9.3)
c135	format(3x,'iwk=',14i5)
c	if(nb.ne.-1) return
	nbeg=nbeg-1
	nconn=nend-nbeg
c
c*		Nieaei-aiipieneiaoey
c
c*	 Auaaeyai iaiyou e caaiaei iiaue iannea eiipaeiao eiioopa
	ib=1
	if(dat(nbk+nknot+nknot).lt.0.) ib=2
	ie=0
	if(dat(nbk+nknot*3-1).gt.t(nend+1)) ie=1
	nconn=nconn+ib-1+ie
	allocate(xy(nconn,2),STAT=irc)
C	if(irc.ne.0) call  mescon('Ia oaaoaao iiapaoeaiie iaiyoe',29,
C     *  nm,0,intg,0,rl,0)
	if(irc.ne.0) goto 999
c*		Anee epaeiea ocaeu aia eiioopa - aiaaaeyai eo a eiioop
	if(ib.eq.2) xy(1,1)=dat(nbk)
	if(ib.eq.2) xy(1,2)=dat(nbk+nknot)
	if(ib.eq.2) t(1)=dat(nbk+nknot+nknot)
	if(ie.eq.0) go to 138
c*		Iineaaiee ocae iaiaoiaeii aianoe a eiioop
	xy(nconn,1)=dat(nbk+nknot-1)
	xy(nconn,2)=dat(nbk+nknot+nknot-1)
	t(nconn+2-ib)=dat(nbk+nknot*3-1)
c
138	do 140 i=0,nconn-ib-ie
	xy(i+ib,1)=dat(nbc+nbeg+i)
140	xy(i+ib,2)=dat(nbc+nbeg+ncon+i)
C			Auaia eii?aeiao oi?ae naea?aiiiai eiioo?a a oaee aey aecoaaecaoee 
      open(unit=15, file='sp_contour.txt', err=999, action='write')
	cibg=4
	do i=1,nconn
	  xxyy(1)=xy(i,1)
	  xxyy(2)=xy(i,2)
        call cnkmdg(fi,lambda,xxyy,v,1,dx,dy,0,0)      
 	  write(15, 32) xxyy(1), xxyy(2)
	enddo
	close(15, err=999)
c
	print *,'ncon,nconn,nbeg',ncon,nconn,nbeg
	call  spapp(xy,t(3-ib+nbeg),nconn,dat(nbk),dat(nbk+nknot)
     *,dat(nbk+nknot+nknot),nknot,dat(nbk+nknot*3),nknot-1,err,wk,
     *iflgsp,irc)
	write(6,150) irc,err
150	format(3x,'i?ia?aiia spapp. irc,err=',i4,f9.3)
	if(irc.eq.0) nv=nv+1
	if(irc.eq.0) idat(nbgknt-1)=nv
	ircsp=irc
	sperr=err



C		?an?ao x,y eii?aeiao nieaeia e eo auaia a oaee aey aecoaeecaoee
      write(*, *) 'nconn =', nconn
      open(UNIT=17, FILE='spline.txt', ERR=999, ACTION='write')
	do i = 1, nconn
	  call trdist(fi,lambda, xy(i,1), xy(i,2), 
     *    dat(nbk), dat(nbk+nknot+1), dat(nbk+(nknot+1)*2), nknot+1,
     *    dat(nbk+(nknot+1)*3),nknot, xx, yy, t,
     *    tga, vrk, distance, epsil, iflgtd)
	  xxyy(1)=xx
	  xxyy(2)=yy
        call cnkmdg(fi,lambda,xxyy,v,1,dx,dy,0,0)      
 	  write(17, 32) xxyy(1), xxyy(2)
	enddo
   32 format(3F12.4)
      close(17, err=999)

c*		Inaiai?aaai iaiyou
      write(*, *) 'Inaiai?aaai iaiyou xy, irc =', irc
	deallocate(xy,STAT=irc)
	write(*, *) 'Inaiaiaeee iaiyou xy'





















CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC^CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC/^\CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//^\\CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	SUBROUTINE SeidelMethod(N, Xcur, GetNewX, EPS)
	INTEGER N
      REAL EPS, Xprev(N), Xcur(N), Xnew, dx

	  DO i = 1, 50
	    Xprev = XCur

	    write(*, *) 'Caeaaeu ', i
	    
		DO j = 1, N
	      CALL GetNewX(Xprev, Xnew, j, irc)
	      Xcur(j) = Xnew
	    ENDDO

	    write(*, *) 'Xcur is', Xcur(1), Xcur(2)

	    dx = 0.
	    DO j = 1, N
	      dx = max(dx, abs(Xcur(j) - Xprev(j)))
	    ENDDO
	    write(*, *) 'iaen. iaaycea ', dx
	    IF (dx.le.EPS) RETURN

	  ENDDO

	END



	SUBROUTINE MinA0R0(X, Xnew, index, irc)
	COMMON /SumData/ SumVl, SumVg, SumRl, SumRg, SumVRl,
     *  SumVRg, SumRRl, SumRRg, NumLess, NumGreater, Rm		 
      INTEGER index, irc
	REAL X(1), Xnew, A0, R0
	  	write(*, *) 'Ex SumVl =', SumVl
		write(*, *) 'Ex SumRl =', SumRl
		write(*, *) 'Ex SumVRl =', SumVRl
		write(*, *) 'Ex SumRRl =', SumRRl
		write(*, *) 'Ex SumVg =', SumVg
		write(*, *) 'Ex SumRg =', SumRg
		write(*, *) 'Ex SumVRg =', SumVRg
		write(*, *) 'Ex SumVRg =', SumRRg
		write(*, *) 'Ex NumLess =', NumLess
		write(*, *) 'Ex NumGreater =', NumGreater
		write(*, *) 'Ex Rm =', Rm

	  irc = 1
	  
	  A0 = X(1)
	  R0 = X(2)
		write(*, *) 'Ex A0 =', A0
		write(*, *) 'Ex R0 =', R0

	  SELECT CASE (index)
	  CASE (1)	!A0
	    Xnew = ((SumVRg - Rm * SumVg) / (R0 - Rm) - 
     *				SumVg - SumVRl / Rm) / 
     *			(2 / (R0 * (R0 - Rm)) * 
     *				((SumRRg - 2 * Rm * SumRg + Rm**2) / (R0 - Rm) -
     *				 2 * (SumRg - NumGreater * Rm)) +
     *			2 * SumRRl / (R0 * Rm**2) + 2 * NumGreater / R0)
	    irc = 0
        CASE (2)  !R0               
          Xnew = (4 * A0 * (SumRRl / Rm**2 + NumGreater) - 
     *			4 * (Rm - 2 * R0) / (R0 - Rm)**3 * 
     *				(SumRRg - 2 * SumRg * Rm + Rm**2) + 
     *			4 * (Rm - 2 * R0) / (R0 - Rm)**2 * 
     *				(SumRg - NumGreater * Rm) +
     *			2 * R0 * (Rm - 2 * R0) / A0 * (R0 - Rm) * 
     *				(SumVRg - Rm * SumVg) - 
     *			4 * A0 / (R0 - Rm) * (SumRg - NumGreater * Rm)) / 
     *			(-2 * (SumVRl / Rm + SumVg))
	    irc = 0
        CASE DEFAULT 
	    irc = 1
	  END SELECT

	END


	SUBROUTINE Example(X, Xnew, index, irc)
	INTEGER index, irc
	REAL X(1), Xnew

	  irc = 1

	  SELECT CASE (index)
	  CASE (1)
          Xnew = sqrt((X(1) * (X(2) + 5) - 1) / 2)
	    irc = 0
        CASE (2)                  
          Xnew = sqrt(X(1) + 3 * log10(X(1)))
	    write(*, *) 'log(', X(1), ') =', log(X(1)) 
	    irc = 0
        CASE DEFAULT 
	    irc = 1
	  END SELECT

	END












      SUBROUTINE InitArraysSortedByDr(N, X, Y, V, CX, CY, VX, VY, VG, D)
      REAL X(N), Y(N), V(N), CX, CY, VX(N), VY(N), VG(N), D(N)
	INTEGER N
      REAL dist(N)

        DO i = 1, N
	    dist(i) = sqrt((CX - X(i))**2 + (CY - Y(i))**2) !i?iaa?eou ia?aea ianneaia X e Y
	  END DO

	  DO i = 1, N
	    DO j = i, N	
	      k=0
	    END DO
	  END DO

	END
