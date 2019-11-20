C***********************************************************************
      SUBROUTINE PRTMOL (SURF,RHO,STEPSZ,NP,np2,NMOL)
C***********************************************************************
C print-out a pdb-like file of the non-zero pixels
C
      DIMENSION XYZ00(3)
      INCLUDE 'param.h'
c     INCLUDE 'blk2_3.h'
      INCLUDE 'centr.h'
      COMPLEX*8 SURF(NP,NP,NP),CRHO
      CHARACTER*3 CRES
      CHARACTER*1 CI
      INTEGER OUT,NACI

      NACI=1
      CI='A'
c
c
	IF(NMOL.EQ.1) THEN
	DO I=1,3
	XYZ00(I)=CENA(I)
	END DO
	OUT=34
	END IF
c
	IF(NMOL.EQ.2) THEN
	DO I=1,3
	XYZ00(I)=CENB(I)
	END DO
	OUT=35
	END IF
C
      IC1=0
      NAT=0
      NRES=1
      CRHO=CMPLX(RHO,0.)
      DO I=1,NP
        DO J=1,NP
          DO K=1,NP
C	     TR=REAL(SURF(I,J,K))
             IF(SURF(I,J,K).EQ.(1.,0.)) THEN
c            IF(SURF(I,J,K).EQ.(0.,0.)) THEN
              X=(I-NP2)*STEPSZ+XYZ00(1)
              Y=(J-NP2)*STEPSZ+XYZ00(2)
              Z=(K-NP2)*STEPSZ+XYZ00(3)
              NAT=NAT+1
              IC1=IC1+1
              IF (IC1.EQ.51) THEN
                NRES=NRES+1
                IC1=0
              END IF
              WRITE (OUT,999) NAT,CI,NRES,X,Y,Z

                IF (NAT.EQ.99999) THEN
		  NAT=0
                  NRES=1
                  IC1=0
                  NACI=NACI+1
		  IF (NACI.EQ.2) CI='B'
		  IF (NACI.EQ.3) CI='C'
		  IF (NACI.EQ.4) CI='D'
		  IF (NACI.EQ.5) CI='E'
                END IF

            END IF
          END DO
        END DO
      END DO
      print *,' nres after edge pixels=',nres
      ic1=0
      nres=nres+1
c      DO I=1,NP
c        DO J=1,NP
c          DO K=1,NP
c	    IF (SURF(I,J,K).EQ.CRHO)THEN
c              X=(I-NP2)*STEPSZ+XYZ00(1)
c              Y=(J-NP2)*STEPSZ+XYZ00(2)
c              Z=(K-NP2)*STEPSZ+XYZ00(3)
c              NAT=NAT+1
c              IC1=IC1+1
c              IF (IC1.EQ.31) THEN
c                NRES=NRES+1
c                IC1=0
c              END IF
c              WRITE (OUT,998) NAT,NRES,X,Y,Z
c            END IF
c          END DO
c        END DO
c      END DO
 999  FORMAT('HETATM',I5,2X,'H1',2X,'HT1',1X,A1,I4,4X,3F8.3)
 998  FORMAT('HETATM',I5,2X,'H2',2X,'HT2',1X,A1,I4,4X,3F8.3)
      RETURN
      END 
