C***********************************************************************
      SUBROUTINE EDGE (NP,RHO,X1)
C***********************************************************************
      INCLUDE 'param.h'
      COMPLEX*8 X1(NP,NP,NP)
C
C     EDGING   (SINGLE EDGE, RHO INSIDE)
      DO 500 I=2,NP-1
      DO 500 J=2,NP-1
      DO 500 K=2,NP-1
C
C IF REAL(X1(I,J,K))=1. AND NONE OF ITS NEIGHBORS IS 1.,
C THEN X3(I,J,K) IS SET TO (RHO,0.)
C
        IF(X1(I,J,K).EQ.(1.,0.)) THEN
          IF(X1(I-1,J,K).EQ.CMPLX(RHO,0.))  GO TO 500
          IF(X1(I,J-1,K).EQ.CMPLX(RHO,0.))  GO TO 500
          IF(X1(I,J,K-1).EQ.CMPLX(RHO,0.))  GO TO 500
          IF(X1(I+1,J,K).EQ.CMPLX(RHO,0.))  GO TO 500
          IF(X1(I,J+1,K).EQ.CMPLX(RHO,0.))  GO TO 500
          IF(X1(I,J,K+1).EQ.CMPLX(RHO,0.))  GO TO 500
C
C         INTERIOR
C
          X1(I,J,K)=CMPLX(0.,0.)
          GO TO 500
        END IF
C
  500 CONTINUE
C 
       RETURN
       END
