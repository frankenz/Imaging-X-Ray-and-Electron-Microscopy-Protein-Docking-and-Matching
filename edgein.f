C***********************************************************************
      SUBROUTINE EDGEIN(NP,NN,RHO1,X1)
C***********************************************************************
      INCLUDE 'param.h'
C
      COMPLEX*8 X1(NP,NP,NP),CRHO,CRHO1
c
      CRHO=CMPLX(0.,0.)
      CRHO1=CMPLX(RHO1,0.)
      NLOOP=0
   10 NLOOP=NLOOP+1
      IF (NLOOP.GT.NN) RETURN   
      DO 500 I=2,NP-1
      DO 500 J=2,NP-1
      DO 500 K=2,NP-1
        IF(X1(I,J,K).EQ.CRHO)  THEN
          IF(X1(I-1,J,K).EQ.(1.,0.))  GO TO 602
          IF(X1(I,J-1,K).EQ.(1.,0.))  GO TO 602
          IF(X1(I,J,K-1).EQ.(1.,0.))  GO TO 602
          IF(X1(I+1,J,K).EQ.(1.,0.))  GO TO 602
          IF(X1(I,J+1,K).EQ.(1.,0.))  GO TO 602
          IF(X1(I,J,K+1).EQ.(1.,0.))  GO TO 602
          GO TO 500
  602     X1(I,J,K)=(9.,0.)
        ELSE
          CONTINUE
        ENDIF  
  500 CONTINUE
      DO 600 I=2,NP-1
      DO 600 J=2,NP-1
      DO 600 K=2,NP-1
        IF(X1(I,J,K).EQ.(9.,0.)) THEN
           X1(I,J,K)=CRHO1
        ELSE
          CONTINUE
        ENDIF  
  600 CONTINUE
      GO TO 10
      END
