C**********************************************************************
      SUBROUTINE MOLSIZE (X,XYZSUM,RMAX,ISTOP,NMOL)
C**********************************************************************
      INCLUDE 'param.h'
C
      DIMENSION X(MAXATM,5),XYZSUM(3)
C
C DETERMINE MOLECULAR CENTER OF MASS; ALL ATOMS HAVE IDENTICAL WEIGHTS.
      DO J=1,3
        XYZSUM(J)=0.
        DO I=1,ISTOP
          XYZSUM(J)=XYZSUM(J)+X(I,J)
	END DO  
	XYZSUM(J)=XYZSUM(J)/ISTOP
      END DO  
C DETERMINE MOLECULAR RADIUS
      RMAX=0.
      DO I=1,ISTOP
        R=0.
        DO J=1,3
          R=R+(X(I,J)-XYZSUM(J))**2
        END DO   
        IF (R.GT.RMAX) RMAX=R
      END DO  
      RMAX=SQRT(RMAX)
      PRINT 200,NMOL
  200 FORMAT (/1X,'MOLECULE ',I1,':')
      PRINT 210,(XYZSUM(I),I=1,3),RMAX
  210 FORMAT (1X,'POSITION OF MOLECULAR CENTER: ',3(F8.3,',',2X)
     1   /1X,'THE MAXIMUM RADIUS OF THE MOLECULE IS ',F8.3,' ANGSTROM.')
      RETURN
      END
