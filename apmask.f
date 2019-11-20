C***********************************************************************
      SUBROUTINE APMASK (NP,X1,X2)
C***********************************************************************
      INCLUDE 'param.h' 
C
      COMPLEX*8 X1(NP,NP,NP),CRHO
      COMPLEX*8 X2(NP,NP,NP)
C
      CRHO=CMPLX(0.,0.)      
      DO 500 I=1,NP
      DO 500 J=1,NP
      DO 500 K=1,NP
        IF (X1(I,J,K).EQ.CRHO) X2(I,J,K)=CRHO
  500 CONTINUE
      DO I=1,NP
        DO J=1,NP
          DO K=1,NP
            X1(I,J,K)=X2(I,J,K)
          END DO
        END DO
      END DO
      RETURN
      END
