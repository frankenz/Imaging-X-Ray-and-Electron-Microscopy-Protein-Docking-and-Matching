C**********************************************************************
      SUBROUTINE ROTATEPOT (X,XOUT,N)
C**********************************************************************
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
      INCLUDE 'cosin.h'
C
      CHARACTER*3 CRES
      DIMENSION X(MAXPOT,3),XOUT(MAXPOT,4),A(3,3),XTT(3)
C
      A(1,1)=COPSI*COPHI-COTHE*SIPHI*SIPSI
      A(1,2)=COPSI*SIPHI+COTHE*COPHI*SIPSI
      A(1,3)=SIPSI*SITHE
      A(2,1)=-SIPSI*COPHI-COTHE*SIPHI*COPSI
      A(2,2)=-SIPSI*SIPHI+COTHE*COPHI*COPSI
      A(2,3)=COPSI*SITHE
      A(3,1)=SITHE*SIPHI
      A(3,2)=-SITHE*COPHI
      A(3,3)=COTHE
C
      DO I=1,N
        DO I1=1,3
          XTT(I1)=0.
          DO I2=1,3
            XTT(I1)=XTT(I1)+A(I1,I2)*X(I,I2)
          END DO 
	END DO         
        DO I1=1,3
          XOUT(I,I1)=XTT(I1)
        END DO
      END DO   
      RETURN
      END
