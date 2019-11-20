C**********************************************
      SUBROUTINE REG (NP,INN,JNN,KNN,II,IEDGE,LIST)
C**********************************************
      DIMENSION LIST(10000)
      IEDGE=0
      IF (INN.EQ.1.OR.INN.EQ.NP) IEDGE=1
      IF (JNN.EQ.1.OR.JNN.EQ.NP) IEDGE=1
      IF (KNN.EQ.1.OR.KNN.EQ.NP) IEDGE=1
      LIST(II)=1000000*INN+1000*JNN+KNN
      RETURN
      END
