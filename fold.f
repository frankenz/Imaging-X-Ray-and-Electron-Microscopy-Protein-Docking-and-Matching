C********************************************************
      SUBROUTINE FOLD (NP,NP2,SURF,in,jn,kn,ln)
C********************************************************
      INCLUDE 'param.h'
C
      COMPLEX*8 SURF(NP,NP,NP)
      DIMENSION CORRL(20)
C
      IL=0
      LN2=2*LN+1 
      JJZS=KN-LN
      JJZL=KN+LN
      JJYS=JN-LN
      JJYL=JN+LN
      JJXS=IN-LN
      JJXL=IN+LN
      DO K=JJZS,JJZL
        DO J=JJYS,JJYL
          do i=0,80
            k1=k+np2+1
            if(k1.le.np2) then
              kk=k1+np2
            else
              kk=k1-np2
            endif
            j1=j+np2+1
            if(j1.le.np2) then
              jj=j1+np2
            else
              jj=j1-np2
            endif
            i1=i+np2+1
            if(i1.le.np2) then
              ii=i1+np2
            else
              ii=i1-np2
            endif
            corr=real(surf(ii,jj,kk))
            IF (LN.EQ.0) THEN
              print 10,i,j,k,corr
  10          format (1x,'THE CORRELATION VALUE AT POSITION ',
     1          I4,1X,I4,1X,I4,'     IS ',F9.0)
            ELSE
              IL=IL+1
              CORRL(IL)=CORR
            END IF
          END DO
          PRINT 20,(CORRL(IL),IL=1,LN2)
  20      FORMAT (10F8.0)
          IL=0
        END DO
        print 25
  25    format ('  ')
      END DO
      RETURN
      END
