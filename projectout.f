C**********************************************************************
      SUBROUTINE PROJECTOUT (SURF,X,NMOL,XM)
C**********************************************************************
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
C
      COMPLEX*8 SURF(NP,NP,NP)
      CHARACTER*3 CRES
      DIMENSION X(MAXATM,5),y(3)
C	
c     icout=0
      ADD=XM*STEPSZ
      IF (XM.EQ.10.) ADD=1.4
      xnp2=stepsz*np2
      ISTOP=NA
      IF (NMOL.EQ.2) ISTOP=NB
      DO K=1,ISTOP
         if (x(k,5).GT.0.) then
c     icout=icout+1
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
                           IF(NMOL.EQ.1)THEN
C     OPTION OMTA
                              IF (X(K,5).EQ.1.) then
                     TR=OMITRHO
                     TI=0.
                  ELSE
C     OPTION TRIM 
                     IF (X(K,5).EQ.2.) then
                       TR=REAL(SURF(NX,NY,NZ))
                       TI=AIMAG(SURF(NX,NY,NZ))
                       IF (TR.NE.RHO) THEN
                          TR=XCUT
c                         TI=0.
                       end if
                     END IF
                   END IF
                ELSE
C     OPTION OMTB
                   IF (X(K,5).EQ.1.) then
                      TR=0.01
                      TI=0.
C     THE ELSE OPTION HERE RELATES TO X(K,5)=2., FOR OPTION CUT for B
                   ELSE
                      TR=XCUT
                      TI=AIMAG(SURF(NX,NY,NZ))
		   ENDIF
                 END IF
		 SURF(NX,NY,NZ)=CMPLX(TR,TI)
                END DO
             end if
             END IF
          END DO
       end if
      END IF
      END DO
      end if
      END DO
c     print *,'In projectout: icout=',icout
      RETURN
      END
      
