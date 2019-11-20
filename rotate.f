C**********************************************************************
      SUBROUTINE ROTATE (X,XOUT,N,ITEST)
C**********************************************************************
      INCLUDE 'param.h'
      INCLUDE 'blk2_3.h'
      INCLUDE 'cosin.h'
C
      CHARACTER*3 CRES
      DIMENSION X(MAXATM,3),XOUT(MAXATM,5),A(3,3),XTT(3)
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
C CHECKING IF THERE ARE ANY SYMMETRY CONSTRAINTS AND APPLYING THEM
C
      ITEST=0
      IF (ISYMM.GT.1) THEN
        ITEST=1
        TRACE=0.
        DO IS=1,3
          TRACE=TRACE+A(IS,IS)
        END DO
        ALPHA=ACOS((TRACE-1.0)*0.5)
        IF (ALPHA.LT.0.) ALPHA=-ALPHA
c       print *,'in rotate: asymm,dsymm=',asymm,dsymm
c       print *,'alpha=',alpha
        IF (ALPHA.LE.(ASYMM+DSYMM).AND.ALPHA.GE.(ASYMM-DSYMM))ITEST=0
      END IF
c     print *,'itest=',itest
      IF (ITEST.EQ.1) RETURN 
C
c     print *,'Rotation matrix from rotate.f'
c     do i=1,3
c       print *,(a(i,j),j=1,3)
c     end do
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
c       XOUT(I,4)=X(I,4)
      END DO   
      RETURN
      END
