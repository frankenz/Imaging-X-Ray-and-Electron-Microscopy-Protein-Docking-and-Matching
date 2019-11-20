C Algorithm developed by: Ziv Frankenstein, PhD

C
C
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
      INCLUDE 'main.h'
      INCLUDE 'title.h'
      INCLUDE 'cosin.h'
      INCLUDE 'centr.h'
C
      COMPLEX*8 SURFA,SURFB
      CHARACTER*3 CRES
C      INTEGER iclock1, iclock2, i_Time

C
      DIMENSION SURFA(MAXNP,MAXNP,MAXNP),SURFB(MAXNP,MAXNP,MAXNP)
      DIMENSION ISURF(MAXNP,MAXNP,MAXNP) 
      DIMENSION XA(MAXATM,5),XB(MAXATM,5),XBB(MAXATM,3)
c     DIMENSION ORIGA(3),RT(3),TRNS(3)
      DIMENSION RT(3),TRNS(3)
      DIMENSION XYZPOT(MAXPOT,4),XYZPOTT(MAXPOT,3)
c     DIMENSION TRNS(3)
C


      PRINT*,'THIS VERSION OF MOLFIT ALLOWS CUTTING SIDE CHAINS '
      IOUT=15
C
c
C  READ INSTRUCTIONS IN FREE FORMAT
C
      iclock1 = mclock()
      CALL XDRED
      iclock2 = mclock()
      iTime = iclock2-iclock1
      write (26, 800) iTime
 800  format (' timing routine xdred CPUtime= ', I6)


C
c     IF (ISITA.EQ.0.AND.ISITB.EQ.0) THEN
c       CALL RENUM_PDB (1)
c       CALL RENUM_PDB (2)
c     END IF
C! Reads how many side chaine atoms to keep
      IF (KCUT.GT.0) THEN
        PRINT 2
2	FORMAT(//' OMITTING ATOMS AS FOLLOWS: RES     TRIMMING AFTER ATOM')
        DO IJ=1,KCUT
          PRINT 3, CRES(IJ),NCRES(IJ)
        END DO
      END IF
3     FORMAT (31X,A3,13X,I3)
C
C  READ COORDINATES OF MOLECULE A
C
      iclock1 = mclock()
      CALL READATOM (10,1,NA,XA)
      iclock2 = mclock()
      iTime = iclock2-iclock1
      write (26, 801) iTime
 801  format (' timing routine readatom CPUtime= ', I6)
C	PRINT*,'XCUT IN MAIN AFTER READ ATOM',XCUT
C
C  READ COORDINATES OF MOLECULE B
C
      CALL READATOM (11,2,NB,XB)
c	PRINT*,'XCUT IN MAIN AFTER READ ATOM',XCUT
C DETERMINE MOLECULAR SIZES
C
      IF (IORIGN.EQ.0) THEN
         iclock1 = mclock()
         CALL MOLSIZE (XA,CENA,RA,NA,1)
         iclock2 = mclock()
         iTime = iclock2-iclock1
         write (26, 802) iTime
 802     format (' timing routine molsize for A CPUtime= ', I6)
         CALL MOLSIZE (XB,CENB,RB,NB,2)
      END IF
C
C DETERMINE ORIGIN OF MOLECULE A AND CHECK IF AN INITIAL ROTATION 
C IS REQUESTED
C
      iclock1 = mclock()
      CALL ORIGIN (XA,1)
      iclock2 = mclock()
      iTime = iclock2-iclock1
      write (26, 803) iTime
 803  format (' timing routine origin for A CPUtime= ', I6)

c     DO I=1,3
c       ORIGA(I) = XYZ0(I)
c     END DO
C
      IF (IROTA.EQ.1) THEN
        PRINT 211,ANGSA
  211 FORMAT (/,1X,'INITIAL ROTATION FOR MOLECULE A BY EULER ANGELS
     1:',3(F8.3,2X))
        DO IA=1,NA
          DO IK=1,3
            XBB(IA,IK)=XA(IA,IK)
          END DO
        END DO

        iclock1 = mclock()
        CALL COSINS (ANGSA,RAD)
        iclock2 = mclock()
        iTime = iclock2-iclock1
        write (26, 804) iTime
 804    format (' timing routine cosins CPUtime= ', I6)

        iclock1 = mclock()
        CALL ROTATE (XBB,XA,NA)
        iclock2 = mclock()
        iTime = iclock2-iclock1        
        write (26, 805) iTime
 805    format (' timing routine rotate CPUtime= ', I6)
      END IF
C
C DETERMINE ORIGIN FOR MOLECULE B
C
      iclock1 = mclock()
      CALL ORIGIN (XB,2)
      iclock2 = mclock()
      iTime = iclock2-iclock1
      write (26, 806) iTime
 806  format (' timing routine origin for B CPUtime= ', I6)

      DO 210 IA=1,NB
      DO 210 IK=1,3
        XBB(IA,IK)=XB(IA,IK)
  210 CONTINUE
C
      IF (NSURF.EQ.0) GO TO 206
C
C PREPARE SURFACE LAYER FOR MOLECULE A
C     print *,'in main, xresl=',xresl

      iclock1 = mclock()
      IF(IOUTC.EQ.0)CALL SURFACE(SURFA,SURFB,XYZPOT,ISURF,XA,1)
      iclock2 = mclock()
      iTime = iclock2-iclock1
      write (26, 807) iTime
 807  format (' timing routine surface 1st time CPUtime= ', I6)


c     istp=1
c     if (istp.eq.1) stop
      IF (IPOT.EQ.1) THEN
        DO IB=1,NPOTB
          DO JB=1,3
            XYZPOTT(IB,JB)=XYZPOT(IB,JB)
          END DO
        END DO
      END IF
c     istp=1
c     if (istp.eq.1) stop
  206 CONTINUE
      IF (IOUTC.NE.0) GO TO 662
      IF (NSURF.EQ.0.AND.IOUTC.EQ.0) STOP
C
      IF (ISYMM.GT.1) THEN
        ASYMM=RAD*360./ISYMM
        DSYMM=RAD*DSYMM
c       print *,'found ISYMM>1 in main'
      END IF
C
      IF (IROT.EQ.0) THEN
        NTHE=1
        NPHI=1
        NPSI=1
        DANG=0.0
      ELSE
        PRINT 220,DANG,ANG1,ANG2,ANG3
  220   FORMAT(/,1X,'MOLECULE B WILL BE ROTATED AT INTERVALS OF ',
     1  F6.2,' DEGREES IN THETA AND PSI.',/,1X
     2  ,'THE ROTATION IN PHI DEPENDS ON THETA.',/,1X
     3  ,'THE INITIAL THETA, PHI AND PSI ANGLES ARE: '
     4  ,3(F7.2,1X),'DEGREES.')
        DANG=RAD*DANG
        ANG1=RAD*ANG1
        ANG2=RAD*ANG2
        ANG3=RAD*ANG3
        IF (IIROT.GT.0) THEN
         NTHE=JROT(1)
         NPHI=JROT(2)
         NPSI=JROT(3)
	 PRINT 231,NTHE,NPHI,NPSI
  231    FORMAT(1X,'THE NUMBER OF STEPS IN THETA,PHI AND PSI IS: ',3I4)
        ELSE
         NTHE=NINT((PI-ANG1)/DANG)+1
         NPSI=NINT((2.*PI-ANG3)/DANG)
c        IF (IIROT.EQ.2) THEN
c          NTHE=JROT(1)
c          IF(JROT(2).NE.0) NPHI=JROT(2)
c          NPSI=JROT(3)
c        END IF
         PRINT 230,NTHE,NPSI
  230    FORMAT(1X,'THE NUMBER OF STEPS IN THETA AND PSI IS: ',2I4,'.',
     1         /1X,'THE NUMBER OF STEPS IN PHI DEPENDS ON THETA.')
	END IF
      END IF
      SIANG2=SIN(DANG/2.)
      PRINT 240, NPEAK,TPEAK
  240 FORMAT(/,1X,I3,' HIGHEST PEAKS (MINIMUM HEIGHT =',F7.1,
     1') KEPT FOR EACH ORIENTATION.')
      HPEAK=0.
c
      IF (IREST.EQ.1) THEN
        TSTART=TREST*RAD
        PHSTART=PHREST*RAD
        PSSTART=PSREST*RAD
        NTSTART=NINT(TSTART/DANG)+1
        NPSSTART=NINT(PSSTART/DANG)+1
      END IF
      NUMROT=0
      DO ITHE=1,NTHE
        THEANG=ANG1+(ITHE-1)*DANG
	TT=THEANG/RAD
        COTHE=COS(THEANG)
	SITHE=SIN(THEANG)
        IF (IIROT.GT.0) THEN
          DPHI=DANG
          GO TO 250
        END IF
	IF (ABS(SITHE).LT.0.00001) THEN
	  NPHI=1
	  NPHSTART=1
        ELSE
	  DPHI=2.*ASIN(SIANG2/SITHE)
	  NPHI=INT(2.*PI/DPHI)+1
          DPHI=2.*PI/NPHI
C!!!! eliminate after capri1
c       dphi=dang
c       nphi=NINT((2.*PI-ANG2)/DANG)
c!!!!
	  IF (IREST.EQ.1) NPHSTART=NINT(PHSTART/DPHI)+1
        END IF
  250   CONTINUE
        PRINT 260,NPHI,TT,DPHI/RAD
  260   FORMAT (1X,I5,' PHI STEPS FOR THETA = ',F7.3,'; DELTA=',F7.3)
        DO IPHI=1,NPHI
          PHIANG=ANG2+(IPHI-1)*DPHI
          COPHI=COS(PHIANG)
          SIPHI=SIN(PHIANG)
          DO IPSI=1,NPSI
            PSIANG=ANG3+(IPSI-1)*DANG
	    IF (IREST.GT.0) THEN
	      IF (ITHE.LT.NTSTART) THEN
                DO IPK=1,NPEAK
                  READ (12,9999) ITEMP
                END DO
 9999           FORMAT (I10)
		GO TO 500
              END IF
	      IF (ITHE.EQ.NTSTART.AND.IPHI.LT.NPHSTART) THEN
       	        DO IPK=1,NPEAK
         	  READ (12,9999) ITEMP
                END DO
	        GO TO 500
              END IF
	      IF (ITHE.EQ.NTSTART.AND.IPHI.EQ.NPHSTART.AND
     1           .IPSI.LT.NPSSTART) THEN
                DO IPK=1,NPEAK
         	  READ (12,9999) ITEMP
                END DO
		GO TO 500
              END IF
	    END IF
            COPSI=COS(PSIANG)
            SIPSI=SIN(PSIANG)
            IF (IIROT.LT.2)THEN
              CALL ROTATE (XBB,XB,NB,ITEST)
              IF (ITEST.EQ.1) GO TO 500
              IF (IPOT.EQ.1) CALL ROTATEPOT (XYZPOTT,XYZPOT,NPOTB)
            END IF
            IF (IIROT.EQ.2) THEN
              CALL ROTXYZ (XBB,XB,NB)
              IF (IPOT.EQ.1) CALL ROTXYZPOT (XYZPOTT,XYZPOT,NPOTB)
            END IF
            NUMROT=NUMROT+1


            iclock1 = mclock()
            CALL PROJECT (SURFB,XB,2,1.)
            iclock2 = mclock()
            iTime = iclock2-iclock1
            write (26, 808) iTime
 808        format (' timing routine project time CPUtime= ', I6)


            ifft=0
            if (xresl.gt.0.0) then
              CALL mFFT (1,SURFB,NP,NP,NP,IER)
              ifft=1
              if (iomtb.eq.1.or.kcut.gt.0.or.ipot.eq.1) ifft=0
              call resl (surfb,np,np2,stepsz,xresl,ifft,2)
            end if
            IF (IOMTB.EQ.1.OR.KCUT.GT.0) THEN
              CALL PROJECTOUT (SURFB,XB,2,1.)
            END IF
C	CALL PRTMOL (SURFB,0.1,STEPSZ,NP,NP2,2)
c           call pixels (np,surfb,1.0,0.0)
c           call pixels (np,surfb,0.0,0.0)
C read potntial for molecule B and put in imaginary or SURFB
	    IF(IPOT.EQ.1) then
              iclock1 = mclock()
	      CALL PROJPOT (SURFB,XYZPOT,ISURF,2)
              iclock2 = mclock()
              iTime = iclock2-iclock1
              write (26, 813) iTime
 813          format (' timing routine projpot time CPUtime= ', I6)

C	CALL PRTPOT (SURFB,STEPSZ,NP,NP2,2)
c              call pixelreal (np,surfb,1.0)
c              call pixelreal (np,surfb,0.0)
c              call pixelimag(np,2,surfb)
            end if
C
c           istp=1
c           if (istp.eq.1) stop

            iclock1 = mclock()
	    CALL PRODUCT (SURFA,SURFB,stepsz,ifft,NP)
            iclock2 = mclock()
            iTime = iclock2-iclock1
            write (26, 809) iTime
 809        format (' timing routine product time CPUtime= ', I6)


C    CALL PIXELS (NP,SURFB,1.0,0.0)
C    CALL PIXELS (NP,SURFB,RHO,0.0)
C    CALL PIXELS (NP,SURFB,0.0,0.0)

            iclock1 = mclock()
            CALL PEAKS(12,THEANG,PHIANG,PSIANG,HPEAK,SURFB,ISURF)
            iclock2 = mclock()
            iTime = iclock2-iclock1
            write (26, 810) iTime
 810        format (' timing routine peaks time CPUtime= ', I6)


c           istp=1
c           if (istp.eq.1) stop

  500       CONTINUE
          END DO
        END DO
      END DO
      PRINT 693,NUMROT
 693  FORMAT (' SCAN DONE.',I5,' ORIENTATIONS TESTED.')
C
C PRINT LIST OF HIGHEST PEAKS, THEIR POSITIONS AND THE ORIENTATION
C
      iclock1 = mclock()
      CALL PEAKLST (12,1,HPEAK,RT1,RT2,RT3,ITX,ITY,ITZ,IVECS)
      iclock2 = mclock()
      iTime = iclock2-iclock1
      write (26, 811) iTime
 811  format (' timing routine peaklist time CPUtime= ', I6)

  662 CONTINUE 
      IF (IOUTC.EQ.0) THEN
        ROTT(1)=RT1
        ROTT(2)=RT2
        ROTT(3)=RT3
        ITRN(1)=ITX
        ITRN(2)=ITY
        ITRN(3)=ITZ
      ENDIF
      DO I=1,3
        TRNS(I)=ITRN(I)*STEPSZ
      END DO
      PRINT 726, (ROTT(I),I=1,3),(ITRN(II),II=1,3)
  726 FORMAT (////1X,'OUTPUT COORDINATES FOR MOLECULE B ROTATED BY '/,
     1  1X,'THETA, PHI, PSI =',3F10.2,/,1X,
     2'AND TRANSLATED BY ',3I5,' STEPS.')
      CALL COSINS (ROTT,RAD)
      IF (IIROT.LT.2) CALL ROTATE (XBB,XB,NB,ITEST)
      IF (IIROT.EQ.2) CALL ROTXYZ (XBB,XB,NB)
      IF (IPOT.EQ.1) THEN
        IF (IIROT.LT.2) CALL ROTATEPOT (XYZPOTT,XYZPOT,NPOTB)
        IF (IIROT.EQ.2) CALL ROTXYZPOT (XYZPOTT,XYZPOT,NPOTB)
      END IF
      iclock1 = mclock()
      CALL WRTC (XB,TRNS,XYZPOT,11)
      iclock2 = mclock()
      iTime = iclock2-iclock1
      write (26, 812) iTime
 812  format (' timing routine wrtc time CPUtime= ', I6)
C
C END OF SCAN
C
      STOP
      END 































