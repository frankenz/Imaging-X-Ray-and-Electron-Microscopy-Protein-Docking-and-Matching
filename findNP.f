C***********************************************************************
      SUBROUTINE FINDNP (NP)
C**********************************************************************
C The routine gets an approximate NP value and return the bigger 
C nearst NP value that is divided by 2,3,5

      INCLUDE 'table.h'
C
      INTEGER idx
C
      DO idx=1, 12
         IF (NP_tbl(idx) .GE. NP) THEN
            NP = NP_tbl(idx)
            RETURN
         END IF
      END DO

C If the NP calculated according to the step size is not in the range 
C of 64-128 return -1.
      NP=-1
C
      RETURN
      END
