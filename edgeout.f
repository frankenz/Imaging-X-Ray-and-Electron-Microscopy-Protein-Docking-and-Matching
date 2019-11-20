C***********************************************************************
      SUBROUTINE EDGEOUT (NP,NN,RHO,RHO1,X1)
C***********************************************************************
      INCLUDE 'param.h'
C
      REAL NN
      COMPLEX*8 X1(NP,NP,NP),CRHO,CRHO1
c
      print*,'Surface of A extended by',NN,' grid steps outwards'
      PRINT*,'IN EDGEOUT CHANGE BACK NN (XXM) AND ILOOP TO INTEGERS!!!'


      CRHO=CMPLX(RHO,0.)
      CRHO1=CMPLX(RHO1,0.)
      nnn= nint(nn)
      ILOOP=0
   10 ILOOP=ILOOP+1
      IF (ILOOP.GT.NNn) RETURN
   
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
      print *,'inside edgeout 1'
      call pixels (np,x1,1.0,0.0)
      call pixels (np,x1,-15.,0.0)
      call pixels (np,x1,9.0,0.0)
      call pixels (np,x1,0.0,0.0)
      DO 600 I=2,NP-1
      DO 600 J=2,NP-1
      DO 600 K=2,NP-1
        IF(X1(I,J,K).EQ.(9.,0.)) THEN
           X1(I,J,K)=CRHO1
        ELSE
          CONTINUE
        ENDIF  
  600 CONTINUE
       print *,'inside edgeout 2'
      call pixels (np,x1,1.0,0.0)
      call pixels (np,x1,-15.,0.0)
      call pixels (np,x1,9.0,0.0)
      call pixels (np,x1,0.0,0.0)
      GO TO 10
      END
