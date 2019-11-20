C**********************************************************************
      SUBROUTINE PROJPOT (SURF,XYZPOT,ISURF,NMOL)
C**********************************************************************
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
      INCLUDE 'centr.h'
C     
      COMPLEX*8 SURF(NP,NP,NP)
      CHARACTER*3 CRES
      DIMENSION ISURF(np,np,np)
      DIMENSION XYZPOT(MAXPOT,4),y(3)
c
c Counter for the points that are included into the potential calculations.
      iproj=0
c Counter for the points that are excluded from the potential calculations.
      iprho=0
      ipall=0
C     DO_1
      DO 10 I=1,NP
	 DO 10 J=1,NP
            DO 10 K=1,NP
               ISURF(I,J,K)=0
C     !
               TI=AIMAG(SURF(I,J,K))
               TR=REAL(SURF(I,J,K))
               TI=0.0
               SURF(I,J,K)=CMPLX(TR,TI)
c!
 10         CONTINUE  
C     DO_1
        xnp2=stepsz*np2
C     IF_1
        if (nmol.eq.1) then
           ISTOP=NPOTA
           RA=RPOTA
           WVL=RHO
      end if
C     IF_1
C
C IF_1
      IF (NMOL.EQ.2)then
         ISTOP=NPOTB
         RA=RPOTB
         WVL=0.01
      ENDIF
C IF_1 
      RA2=RA**2
C     DO_1
      DO K=1,ISTOP
C DO_2
         do ii=1,3
            y(ii)=(xyzpot(k,ii)+xnp2)/stepsz
         end do
C DO_2
         dy1=y(1)-ra
         MINX=INT(dy1)
         if (float(minx).ne.dy1) minx=minx+1
         MAXX=INT((y(1)+RA))
C     IF_1
         if(minx.le.maxx) then
C DO_2
            DO NX=MINX,MAXX
               XX=(y(1)-NX)**2
               TMP=RA2-XX
C     IF_2
               IF (TMP.GE.0.) THEN
                  DD=SQRT(RA2-XX)
                  dy2=y(2)-dd
                  MINY=INT(dy2)
                  if (float(miny).ne.dy2) miny=miny+1
                  MAXY=INT((y(2)+DD))
C IF_3
                  if(miny.le.maxy) then     
C     DO_3
                     DO NY=MINY,MAXY
                        YY=(y(2)-NY)**2
                        TMP=RA2-XX-YY
C     IF_4
                        IF (TMP.GE.0.)  THEN
                           DD=SQRT(RA2-XX-YY)        
                           dy3=y(3)-dd
                           MINZ=INT(dy3)
                           if (float(minz).ne.dy3) minz=minz+1
                           MAXZ=INT((y(3)+DD))
C     IF_5
                           if (minz.le.maxz) then
C     DO_4
                              DO NZ=MINZ,MAXZ
                                 TR=REAL(SURF(NX,NY,NZ))
                                 TI=AIMAG(SURF(NX,NY,NZ))
                                 IF (TR.NE.WVL) THEN
                                    TI=TI+XYZPOT(K,4)
                                    ISURF(NX,NY,NZ)=ISURF(NX,NY,NZ)+1
c     ELSE
c     iprho=iprho+1 
                                 END IF
                                 SURF(NX,NY,NZ)=CMPLX(TR,TI)
C     DO_4
                              END DO
C     IF_5
                           end if
C IF_4
                        END IF
C     DO_3
                     END DO
C     IF_3
                  end if
C     IF_2
               END IF
C     DO_2
            END DO
C     IF_1
         end if
C     DO_1
      END DO
      icmax=0
c     ncc=0
	DO NX=1,NP
           DO NY=1,NP
              DO NZ=1,NP
                 TI=AIMAG(SURF(NX,NY,NZ))
                 TR=REAL(SURF(NX,NY,NZ))
                 IF(TI.NE.0.0) THEN
c               NCC=NCC+1
                    IIC=ISURF(NX,NY,NZ)
	        TI=TI/IIC
             END IF
	    SURF(NX,NY,NZ)=CMPLX(TR,TI)
         ENDDO
      ENDDO
      ENDDO
300   format('HETATM',i5,2x,'A1 ',1x,'POS',i6,4x,f8.3,f8.3,f8.
     1     3,2x,f10.5)
 301  format(f7.3)
	
      RETURN
      END
