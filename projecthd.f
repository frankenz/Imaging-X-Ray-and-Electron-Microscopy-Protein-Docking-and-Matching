C**********************************************************************
      SUBROUTINE PROJECTHD (SURF,X,IHD,NMOL,XM)
C**********************************************************************
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
C
      COMPLEX*8 SURF(NP,NP,NP)
      CHARACTER*3 CRES
      DIMENSION X(MAXATM,5),y(3)
      DIMENSION IHD(MAXATM)
C
C
      hdd=hd
      if (nmol.eq.2) hdd=1.

	  
      nhyd=0
      ADD=XM*STEPSZ
      IF (XM.EQ.10.) ADD=1.4
C
      xnp2=stepsz*np2
      ISTOP=NA
      if (nmol.eq.2) istop=nb


      DO K=1,ISTOP
C     
         IF (IHD(K).EQ.1) THEN
            nhyd=nhyd+1
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
                                 tr=real(surf(nx,ny,nz))
                                 if (tr.ge.0.) then
                                    SURF(NX,NY,NZ)=CMPLX(tr,HDD)
                                 end if
                              END DO
                           end if
                        END IF
                     END DO
                  end if
               END IF
            END DO
         END IF
c         if (nmol.eq.2.and.ihd(k).eq.1) then
c         print *,k,nhyd
c         call pixels(np,surf,1.,0.)
c         end if
      END DO
CC      print *,'in projechd: nmol=',nmol,'  nhyd=',nhyd
CC      print *,'in projechd:call pixels for (np,surf,hdd,0.0)'
CC      call pixels (np,surf,hdd,0.0)
CC      print *,'in projechd:call pixels for (np,surf,0.0,0.0)'
CC      call pixels (np,surf,0.0,0.0)
      RETURN
      END
