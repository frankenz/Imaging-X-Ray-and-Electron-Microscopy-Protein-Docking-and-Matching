C**********************************************************************
      SUBROUTINE SRTP (JP,P,NX,NY,NZ,PMIN)
C**********************************************************************
      INCLUDE 'param.h'
      DIMENSION P(MAXNPK*2),NX(MAXNPK*2),NY(MAXNPK*2),NZ(MAXNPK*2)
      DO J=1,JP
        J1=J+1
        DO K=J1,JP
          IF(P(K).GT.P(J)) THEN      
            CALL EXCHNG (P(J),P(K))
            CALL EXCHGI (NX(J),NX(K))
            CALL EXCHGI (NY(J),NY(K))
            CALL EXCHGI (NZ(J),NZ(K))
          END IF 
        END DO
      END DO  
      PMIN=P(MAXNPK+1)
      RETURN
      END
