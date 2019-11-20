C**********************************************************************
      SUBROUTINE ROTXYZ (X,XOUT,N)
C**********************************************************************
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
      INCLUDE 'cosin.h'
      DIMENSION X(MAXATM,3),XOUT(MAXATM,5),A(3,3),XTT(3)
      CHARACTER*3 CRES
c
      A(1,1)=COTHE*COPHI
      A(1,2)=COTHE*SIPHI
      A(1,3)=-SITHE
      A(2,1)=SIPSI*SITHE*COPHI-COPSI*SIPHI
      A(2,2)=SIPSI*SITHE*SIPHI+COPSI*COPHI
      A(2,3)=COTHE*SIPSI
      A(3,1)=COPSI*SITHE*COPHI+SIPSI*SIPHI
      A(3,2)=COPSI*SITHE*SIPHI-SIPSI*COPHI
      A(3,3)=COTHE*COPSI
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
