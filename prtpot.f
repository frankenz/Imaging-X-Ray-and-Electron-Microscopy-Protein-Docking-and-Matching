C***********************************************************************
      SUBROUTINE PRTPOT (SURF,STEPSZ,NP,np2,NMOL)
C***********************************************************************
C print-out a pdb-like file of the non-zero pixels
C
	COMPLEX*8 SURF(NP,NP,NP)
	INCLUDE 'param.h'
c       INCLUDE 'blk2_3.h'
	INCLUDE 'centr.h'
	DIMENSION XYZ00(3)
	CHARACTER*3 CRES

	INTEGER OUT
c
	IF(NMOL.EQ.1) THEN
	DO I=1,3
	XYZ00(I)=CENA(I)
	END DO

c        CRHO=rho

	OUT=36
	END IF
c
	IF(NMOL.EQ.2) THEN
	DO I=1,3
	XYZ00(I)=CENB(I)
	END DO

c       CRHO=0.1

	OUT=37
	END IF
C
       print *,'inside prtmol'
c        print *,(xyz00(i),i=1,3)
      print *,'PRINT PRTMOL NP IS ',NP
      IC1=0
      NAT=0
      NRES=1
      DO I=1,NP
        DO J=1,NP
          DO K=1,NP
	    TI=AIMAG(SURF(I,J,K))
            TR=REAL(SURF(I,J,K))
c	IF (NMOL.EQ.2)TR=TR*(-1)
            IF (TR.EQ.1.AND.TI.GT.0.) THEN 
              X=(I-NP2)*STEPSZ+XYZ00(1)
              Y=(J-NP2)*STEPSZ+XYZ00(2)
              Z=(K-NP2)*STEPSZ+XYZ00(3)
              NAT=NAT+1
              IC1=IC1+1
              IF (IC1.EQ.31) THEN
                NRES=NRES+1
                IC1=0
               END IF
              WRITE (OUT,999) NAT,NRES,X,Y,Z
            END IF
          END DO
        END DO
      END DO
      print *,' nres after edge pixels=',nres
      ic1=0
      nres=nres+1
      DO I=1,NP
        DO J=1,NP
          DO K=1,NP
	    TI=AIMAG(SURF(I,J,K))
            TR=REAL(SURF(I,J,K))
c	IF (NMOL.EQ.2)TR=TR*(-1)
            IF (TR.EQ.1.AND.TI.LT.0.) THEN
              X=(I-NP2)*STEPSZ+XYZ00(1)
              Y=(J-NP2)*STEPSZ+XYZ00(2)
              Z=(K-NP2)*STEPSZ+XYZ00(3)
              NAT=NAT+1
              IC1=IC1+1
              IF (IC1.EQ.31) THEN
                NRES=NRES+1
                IC1=0
              END IF
              WRITE (OUT,998) NAT,NRES,X,Y,Z
            END IF
          END DO
        END DO
      END DO
 999  FORMAT('HETATM',I5,2X,'H1',2X,'HT1',I6,4X,3F8.3)
 998  FORMAT('HETATM',I5,2X,'H2',2X,'HT2',I6,4X,3F8.3)
      RETURN
      END 
