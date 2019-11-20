C findIdx.f

C This function gets a 4 letter code, look for it in a char 
C thable and returns the index to that code.
C If the code does not apear in the table the function 
C returns -1
C--------------------------------------------------------
      SUBROUTINE findIdx (code, index)
C--------------------------------------------------------

      INCLUDE 'table.h'

      INTEGER index, i, arr_size, no_Idx
      CHARACTER*4 code
 
      no_Idx = -1
      arr_size = 33

      DO i=1, arr_size
         IF (codes(i) .EQ. code) THEN
            index = i
            RETURN
         END IF
      END DO


      IF (i .LE. arr_size) THEN
         index = i
      ELSE
         index = no_Idx
      END IF  


      RETURN
      END
