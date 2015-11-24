      SUBROUTINE GESSE(X,FX,GR,NM,STP,DGR,MULT)
C*        ÏPOÃPAMMA PAC×ETA MATPÈÖÛ ÃECCE
      COMMON/INDXDR/ INDXD
     */grcom/ig,ing,ifst,ngrad,mlt,ibinar(100)
      integer*4 ibinar 
      REAL*4 X(1),FX(1),GR(1),STP(1),DGR(1),MULT(1),mlt
      INTEGER NM(2),NNMM(2)
c	write (6,66) ig,indxd
c66	format(3x,'ig,ind=',2i5)
      if(ig.eq.1) return
      NNMM(1)=NM(1)
      NNMM(2)=4
      INDXD=0
      CALL  GRAD(X,FX,GR,NNMM,STP,DGR,MULT)
      INDXD=1
      RETURN
      END
      BLOCK DATA
      COMMON /INDXDR/ IND
      data IND/0/
      END
