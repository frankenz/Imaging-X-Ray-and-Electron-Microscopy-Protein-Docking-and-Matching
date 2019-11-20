C Algorithm developed by: Ziv Frankenstein, PhD
C
C
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
      INCLUDE 'main.h'
      INCLUDE 'title.h'
      INCLUDE 'cosin.h'
      INCLUDE 'centr.h'
      INCLUDE 'table.h'
C
      COMPLEX*8 SURFA,SURFB
      CHARACTER*3 CRES
      INTEGER iNumOri, iMaxNumOri
C
      DIMENSION SURFA(MAXNP,MAXNP,MAXNP),SURFB(MAXNP,MAXNP,MAXNP)
      DIMENSION ISURF(MAXNP,MAXNP,MAXNP) 
      DIMENSION XA(MAXATM,5),XB(MAXATM,5),XBB(MAXATM,3)
c     DIMENSION ORIGA(3),RT(3),TRNS(3)
      DIMENSION RT(3),TRNS(3)
      DIMENSION XYZPOT(MAXPOT,4),XYZPOTT(MAXPOT,3)
      DIMENSION IHDA(MAXATM),IHDB(MAXATM)
c     DIMENSION TRNS(3)
C
c
c
      iMaxNumOri= 20000
      PRINT*,'FitEM2EMIN VERSION 1.0 '
      PRINT*,'THIS VERSION ALLOWS SURFACE GEOMETRIC
     1 FITTING'

      PRINT*
      IOUT=15

C
C Initialize the arrays that will hold residues details if the OMTA
C or OMTB options are on
      DO J=1, MAXOMIT
         CID_omitA(J)=''
         RID_omitA(J)=''
         rNUM_omitA(J)=0
         CID_omitB(J)=''
         RID_omitB(J)=''
         rNUM_omitB(J)=0
      ENDDO
C
c
C  READ INSTRUCTIONS IN FREE FORMAT
      CALL ReadData

C
C Print a messege if this is a refinement run and not docking scan
      if (IREFN .GT. 0) then
        print 2222,irefn
      end if
 2222 format (/,'REFINEMENT FOR',i3,' SOLUTIONS.',/,
     1       'IT IS ASSUMED THAT A SCAN WAS PERFORMED BEFORE.',/)
C
C Check that the user do not tries to perform both geo-elec and
C geo-hyd scan at the same docking run. If he does stop the scan.
C The same for operating geo-elec and site together or site with
C geo-hyd
      IF ((IHDRF .EQ. 1) .AND. (IPOT .EQ.1)) THEN 
         PRINT*,'ERROR - YOU CAN DO EITHER GEOMETRIC-ELECTROSTATIC SCAN
     1OR GEOMETRIC-HYDROPHOBIC SCAN. YOU CAN NOT DO BOTH AT THE SAME
     2SCAN.'
         STOP
      ELSEIF ((IHDRF .EQ. 1) .AND. ((NSITA.GT.0) .OR. 
     1        (NSITB.GT.0))) THEN 
         PRINT*,'ERROR - YOU CAN DO EITHER GEOMETRIC-ELECTROSTATIC SCAN
     1OR WEIGHTED GEOMETRIC SCAN. YOU CAN NOT DO BOTH AT THE SAME SCAN.'
         STOP
      ELSEIF (((NSITA.GT.0) .OR. (NSITB .GT. 0)) .AND. 
     1        (IPOT.EQ.1)) THEN 
         PRINT*,'ERROR - YOU CAN DO EITHER GEOMETRIC-ELECTROSTATIC SCAN
     1OR WEIGHTED GEOMETRIC SCAN. YOU CAN NOT DO BOTH AT THE SAME
     2SCAN.'
         STOP
      ELSE
         CONTINUE
      ENDIF
C
C Check that the user did not inserted more than one site option for
C both molecule A and B. 
      IF (((NSITA .GT. 1) .AND. (NSITB .GT. 0)) .OR.
     1     ((NSITA .GT. 0) .AND. (NSITB .GT. 1))) PRINT 10 
 10   FORMAT ('WARNING: NO MORE THAN ONE WIEGHT IS ALLOWED WHEN OPERATIN
     1G BOTH SITA AND SITB TOGETHER ',/)

C If only option OMTA is on (OMTB is off), change the value of
C OMITRHO (the value that is inserted into the omitted voxels)
C into -15. (When both molecules are omitted than there is a 
C need to give the omitted voxels in A a very negative value).
      IF ((IOMTA .EQ. 1) .AND. (IOMTB .EQ. 0)) OMITRHO = RHO
c
C
c     IF (ISITA.EQ.0.AND.ISITB.EQ.0) THEN
c       CALL RENUM_PDB (1)
c       CALL RENUM_PDB (2)
c     END IF
C! Reads how many side chaine atoms to keep
      IF (KCUT.GT.0) THEN
        PRINT 2
2	FORMAT(//' OMITTING ATOMS AS FOLLOWS: RES TRIMMING AFTER ATOM')
        DO IJ=1,KCUT
          PRINT 3, CRES(IJ),NCRES(IJ)
        END DO
        PRINT*
      END IF
3     FORMAT (31X,A3,13X,I3)
C

C  READ COORDINATES OF MOLECULE A
      CALL READEM (10,1,NA,XA,IHDA)

C      CALL READATOM (10,1,NA,XA,IHDA)

C
C  READ COORDINATES OF MOLECULE B
C
      CALL READEM (11,2,NB,XB,IHDB)

C      CALL READATOM (11,2,NB,XB,IHDB)


C
C DETERMINE MOLECULAR SIZES
C
      IF (IORIGN.EQ.0) THEN
         CALL MOLSIZE (XA,CENA,RA,NA,1)
         CALL MOLSIZE (XB,CENB,RB,NB,2)
      END IF
C
C DETERMINE ORIGIN OF MOLECULE A AND CHECK IF AN INITIAL ROTATION 
C IS REQUESTED
      CALL ORIGIN (XA,1)


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

        CALL COSINS (ANGSA,RAD)
        CALL ROTATE (XBB,XA,NA,ITEST)
      END IF
C
C DETERMINE ORIGIN FOR MOLECULE B
      CALL ORIGIN (XB,2)

      DO 210 IA=1,NB
      DO 210 IK=1,3
        XBB(IA,IK)=XB(IA,IK)
  210 CONTINUE
C
      IF (NSURF.EQ.0) GO TO 206
C
C If RESL option is on print a message to the output file
 215  FORMAT(/,' APPLYING RESOLUTION MODIFICATION (RESL=',f4.0,
     1     ', DT=', F5.2,')')
      IF (XRESL .GT. 0.0) PRINT 215, xresl,dt
C
C PREPARE SURFACE LAYER FOR MOLECULE A

      IF(IOUTC.EQ.0)CALL SURFACE(SURFA,SURFB,XYZPOT,ISURF,XA,IHDA,1)
C      iflag = 1
C      if (iflag .EQ. 1) stop

C     istp=1
C     if (istp.eq.1) stop
C Copy the xyzpot array to another array - working array
C (simillar to what we do with the B molecule geo arrays).
      IF (IPOT.EQ.1) THEN
        DO IB=1,NPOTB
          DO JB=1,3
            XYZPOTT(IB,JB)=XYZPOT(IB,JB)
          END DO
        END DO
      END IF
C     istp=1
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
c REFN option
      if (irefn.gt.0) go to 1000
C
      IF (IROT.EQ.0) THEN
        NTHE=1
        NPHI=1
        NPSI=1
        DANG=0.0
      ELSE
        PRINT 220,DANG,ANG1,ANG2,ANG3
  220   FORMAT(//,1X,'MOLECULE B WILL BE ROTATED AT INTERVALS OF ',
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
  240 FORMAT(/,I4,' HIGHEST PEAKS (MINIMUM HEIGHT =',F11.1,
     1') ARE KEPT FOR EACH ORIENTATION.'///)
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

            CALL PROJECT (SURFB,XB,2,1.)
c              print*, 'project of molecule B'
c             call pixels (np,surfb,1.0,0.0)
c             call pixels (np,surfb,0.0,0.0)
c             call pixels (np,surfb,rho,0.0)
C
C     Resolution option. Enter this option only if the user has used
C     the RESL option (which means that xresl > 0.0).
            IF (xresl. GT. 0.0) THEN
	RESOL=xresl
        xresl= 0.577*(NP-1)*STEPSZ/RESOL
c        PRINT*,'APPLYING RESOLUTION MODIFICATION ON B (RESL=',RESOL
c     1 ,xresl,'DT=',DT,')'
c     do forward fft (for the resolution)
               CALL MFFT(1,SURFB,np,np,np,ier)
C     Do the resolution and dichotomization
               CALL RESOLUTION (np,surfb,xresl)
               CALL DICHO (np,DT,rangeu,surfb)
            END IF
c            print*, 'after resolution'
c           call pixels (np,surfb,1.0,0.0)
c           call pixels (np,surfb,0.0,0.0)

c 	CALL PRTMOL (SURFB,0.1,STEPSZ,NP,NP2,2)

C
C Projection of OMIT and/or TRIM options
            IF (IOMTB.EQ.1.OR.KCUT.GT.0) THEN
              CALL PROJECTOUT (SURFB,XB,2,1.)
            END IF
c	CALL PRTMOL (SURFB,0.1,STEPSZ,NP,NP2,2)
c           call pixels (np,surfb,1.0,0.0)
c           call pixels (np,surfb,0.0,0.0)
c           call pixels (np,surfb,0.01,0.0)
C
C Projects hydrophobicity for molecule B and put in imaginary of SURFB
            IF (IHDRF .EQ. 1) THEN
               CALL PROJECTHD (SURFB,XB,IHDB,2,1.)
            END IF
C
C read potntial for molecule B and put in imaginary of SURFB
	    IF(IPOT.EQ.1) then
	      CALL PROJPOT (SURFB,XYZPOT,ISURF,2)

C	CALL PRTPOT (SURFB,STEPSZ,NP,NP2,2)
c              call pixelreal (np,surfb,1.0)
c              call pixelreal (np,surfb,0.0)
c              call pixelimag(np,2,surfb)
            end if
C
C Do fft for molecule B, conjugation and back fft
            ifft = 0
	    CALL PRODUCTNEW (SURFA,SURFB,stepsz,ifft,NP)
c
            CALL PEAKS(12,THEANG,PHIANG,PSIANG,HPEAK,SURFB,ISURF)
c
c           istp=1
c           if (istp.eq.1) stop

  500       CONTINUE
          END DO
        END DO
      END DO
      PRINT 693,NUMROT
 693  FORMAT (' SCAN DONE.',I5,' ORIENTATIONS TESTED.')
C
C Check if the total number of peaks is higher than 20000. 
C If so let the user know he sould use another sorting function.
      iNumOri = NUMROT * npeak
      IF (iNumOri .GT. iMaxNumOri) THEN
         print* 
         print*,'*****************************************'
         print*, 'USE "SORTALL" ROUTINE TO SORT THE PEAKS'
         goto 662
      END IF
C
C PRINT LIST OF HIGHEST PEAKS, THEIR POSITIONS AND THE ORIENTATION
C
      CALL PEAKLST (12,1,HPEAK,RT1,RT2,RT3,ITX,ITY,ITZ,IVECS)
c
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

      CALL WRTC (XB,TRNS,XYZPOT,45)
c
C
C END OF SCAN
C
 1000 continue
c REFN option
      npeak=1
      iirot=2
      do j=1,irefn
        read (iout,660) ij,p1,(rt(i),i=1,3),kkx,kky,kkz
  660   format(i7,12x,f8.3,2x,3f8.2,2x,3i6)
        print 661,ij,(rt(i),i=1,3),kkx,kky,kkz,p1
  661   format(//,'Refining solution',i3,'.',/,' Starting from:',
     1         3f8.2,2x,3i6,'; Score=',f5.0)
        CALL COSINS(RT,rad)
        CALL ROTATE (XBB,XB,NB,ITEST)
c       print *,'j=',j
        CALL REFINE (SURFA,SURFB,ISURF,
     1               IJ,KKX,KKY,KKZ,RT,XB,P1)
c       print *,'back in main: refine peak',ij,':',p1,rt,kkx,kky,kkz
        write (25,660) ij,p1,rt,kkx,kky,kkz
      end do
C
C
      STOP
      END 
































