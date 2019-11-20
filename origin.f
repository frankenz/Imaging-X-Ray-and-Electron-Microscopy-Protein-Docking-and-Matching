C***********************************************************************
      SUBROUTINE ORIGIN (X,NMOL)
C**********************************************************************
C FINDS THE ORIGIN AND STEP SIZE OF THE GRID
C
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
      INCLUDE 'centr.h'
C
      DIMENSION X(MAXATM,5)
      CHARACTER*3 CRES
      INTEGER NP_tmp
C
      IF (IORIGN.NE.0) THEN
        DO I=1,3
          XYZ0(I)=ORIG(I)
        END DO
      ELSE
C  FIND STEP SIZE AND ORIGIN OF GRID BASED ON THE SIZES OF THE TWO
C  MOLECULES
        IF (ISTEP.EQ.0) THEN
          IF (NP.EQ.0) STOP
          STEPSZ=2.*(RA+RB)/(NP-5)
        ELSE
          IF (NP.EQ.0) THEN
            NP=NINT(2.*(RA+RB)/STEPSZ)+5
            NSURF=1
c           if (ipot.eq.1) then
c           knp=int(np/2)
c           if (knp*2.eq.np) then
c           print *,' NP IS CHANGED TO NP-1 TO MAKE IT ODD'
c           np=np-1
c           end if
c           end if


            CALL findNP(NP)
            IF (NP .EQ. -1) THEN
               print*,'******************************************'
               print*, 'ERROR, THE STEP YOU INSERTED IS TOO SMALL.'
               PRINT*,'THE MINIMAL STEP SIZE IS: ',2.*(RA+RB)/(432-5)
               STOP
            END IF
          ELSE
            NP_tmp=NINT(2.*(RA+RB)/STEPSZ)+5          
            IF (NP .LT. NP_tmp) THEN
               print*,'******************************************'
               print*, 'ERROR, THE STEP YOU INSERTED IS TOO SMALL.'
               print 50, NP,2.*(RA+RB)/(NP-5)
 50            FORMAT(' THE MINIMAL STEP SIZE FOR SURF= '
     1                      ,I4,' IS: ', F6.4)
               STOP
            END IF
          END IF
c      IF (NMOL .EQ. 1) PRINT 90,NP,NP,NP
c          PRINT 90,NP,NP,NP
c   90     FORMAT (/' THE MOLECULES ARE PROJECTED INTO AN ARRAY'
c     *           ,' OF ',I3,'x',I3,'x',I3,'  POINTS')
        END IF
        IF (NMOL.EQ.1) THEN
          DO I=1,3
            XYZ0(I)=CENA(I)
          END DO
        ELSE
          DO I=1,3
            XYZ0(I)=CENB(I)
          END DO
        ENDIF
      ENDIF
c
      IF (NMOL .EQ. 1) PRINT 90,NP,NP,NP
c          PRINT 90,NP,NP,NP
   90     FORMAT (/' THE MOLECULES ARE PROJECTED INTO AN ARRAY'
     *           ,' OF ',I3,'x',I3,'x',I3,'  POINTS')

c PRINT THE STEP SIZE ONLY ONCE TO THE OUTPUT FILE
      IF (NMOL .EQ. 1) PRINT 120,STEPSZ
 120  FORMAT (/' GRID STEP SIZE = ',F8.5,' ANGSTROMS.')
C      PRINT 120,NMOL,(XYZ0(I),I=1,3),STEPSZ
C 120  FORMAT (//1X,'ORIGIN OF PROJECTION FOR MOLECULE ',I1,
C     1       /1X,'AT X= ',F7.3,', Y= ',
C     2  F7.3,', Z= ',F7.3,'.'/1X,'GRID STEP SIZE = ',F8.5,
C     3  ' ANGSTROMS.')
      ISTOP=NA
      IF (NMOL.EQ.2) ISTOP=NB
      DO IA=1,ISTOP
        DO J=1,3
          X(IA,J)=X(IA,J)-XYZ0(J)
        END DO
      END DO
      ISTEP=1
      RETURN
      END
