 
C**********************************************
      SUBROUTINE UNPECK (LIST,L,IN,JN,KN)
C**********************************************
      DIMENSION LIST(10000)
      NLIST=LIST(L)
      IN=INT(NLIST/1000000.)
      NLIST=NLIST-1000000*IN
      JN=INT(NLIST/1000.)
      KN=NLIST-1000*JN
      RETURN
      END
