C**********************************************************************
      SUBROUTINE PROJECT (SURF,X,NMOL,XM)
C**********************************************************************
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
C
      
      COMPLEX*8 SURF(NP,NP,NP)
      REAL XOUT
      DIMENSION X(MAXATM,5),y(3)
      CHARACTER*3 CRES

C
      TR=1.0
      TI=0.0
! create a neitive grid for molecule 1 and 0 grid for molecule 2     
      XOUT=RHO
      if (NMOL.EQ.2) XOUT=0.0
c
      ADD=XM*STEPSZ
c      print *,' Surface of A extended by ',xm,' x step( ',add,')'
c      IF (XM.EQ.10.) ADD=3.3862
c      if (xm.eq.0) add=
      DO 20 I=1,NP
      DO 20 J=1,NP
      DO 20 K=1,NP
         SURF(I,J,K)=CMPLX(XOUT,0.)

   20 CONTINUE
      xnp2=stepsz*np2
      ISTOP=NA
      IF (NMOL.EQ.2) ISTOP=NB
C SITE OPTION
      IF (NSITA.NE.0.AND.NSITB.EQ.0.AND.NMOL.EQ.2) TI=1.
      IF (NSITA.EQ.0.AND.NSITB.NE.0.AND.NMOL.EQ.1.AND.XM.EQ.10.)
     1     TI= 1.0
c     print *,'isita,isitb,nmol,xm',isita,isitb,nmol,xm
c     print *,'tr=',tr
c     print *,'ti=',ti
C
      DO K=1,ISTOP
         IF (NMOL.EQ.1.AND.NSITA.NE.0.AND.XM.EQ.10.) THEN
            TI=0.
            DO L=1,NATSITA
               IF (K.EQ.IATSITA(L,1)) THEN
                  IDXSITA = IATSITA(L,2)
                  TI=XSITA(IDXSITA)
               END IF
            END DO
c        print *,'k,ti=',k,ti
         END IF
         IF (NMOL.EQ.2.AND.NSITB.NE.0) THEN
            TI=0.
            DO L=1,NATSITB
               IF (K.EQ.IATSITB(L,1)) THEN 
                  IDXSITB = IATSITB(L,2)
                  TI=XSITB(IDXSITB)
               END IF
            END DO
c     print *,'k,ti=',k,ti
         END IF
C
         RA=(X(K,4)+ADD)/STEPSZ
         RA2=RA**2.
         do ii=1,3
            y(ii)=(x(k,ii)+xnp2)/stepsz
         end do
         dy1=y(1)-ra
         MINX=INT(dy1)
         if (float(minx).ne.dy1) minx=minx+1
         MAXX=INT((y(1)+RA))
         DO NX=MINX,MAXX
            XX=(y(1)-NX)**2
            TMP=RA2-XX
            IF (TMP.GE.0.) THEN
               DD=SQRT(RA2-XX)
               dy2=y(2)-dd
               MINY=INT(dy2)
               if (float(miny).ne.dy2) miny=miny+1
               MAXY=INT((y(2)+DD))
               if(miny.le.maxy) then     
                  DO NY=MINY,MAXY
                     YY=(y(2)-NY)**2
                     TMP=RA2-XX-YY
                     IF (TMP.GE.0.)  THEN
                        DD=SQRT(RA2-XX-YY)        
                        dy3=y(3)-dd
                        MINZ=INT(dy3)
                        if (float(minz).ne.dy3) minz=minz+1
                        MAXZ=INT((y(3)+DD))
                        if (minz.le.maxz) then
                           DO NZ=MINZ,MAXZ
c     TR=1.
c                 TI=0.
                              SURF(NX,NY,NZ)=CMPLX(TR,TI)
                           END DO
                        end if
                     END IF
                  END DO
               end if
            END IF
         END DO
      END DO
      RETURN
      END
