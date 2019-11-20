C-----------------------------------------------------------------       
       subroutine charToReal (cNum, rNum)
C------------------------------------------------------------------
c       
c
C This subroutine gets a charecter array and turns it into a flaot
C number. The function returns the number
c
C------------------------------------------------------------------
c
       REAL r_tmpNum, rNum, atmp
       INTEGER i, iSpace, iZero, iSign, iMinus, iDot, iflag_dot, iblank
       INTEGER iNine, iPower
       CHARACTER*1 cNum(70)
c
C initialization of variables
       r_tmpNum = 0.0
       i=1
       iSign = 1
       iflag_dot = 0
       iPower = 0
       rNum = 0.0
c
c
       iSpace = ICHAR(' ')
       iZero =  ICHAR('0')
       iNine =  ICHAR('9')
       iMinus = ICHAR('-')
       iDot =   ICHAR('.')
       iBlank = ICHAR('') 
c 
c
C Check the sign of the of the number
       IF (ICHAR(cNum(1)) .EQ. iMinus) THEN
          isign = -1
          i = i+1
       END IF
c
c
C Read a character from the character array. Substruct the ASCII value 
C of '0' from the charecter, and add this digit to the real number.
c
 10    DO WHILE (ICHAR(cNum(i)) .NE. iSpace .AND.
     1            ICHAR(cNum(i)) .NE. iBlank )    
c
C Check that the string has a number in it (isNumeric)
 20       IF ((ICHAR(cNum(i)) .NE. iDot) .AND. 
     1        ((ICHAR(cNum(i)) .LT. iZero) .OR. 
     2        (ICHAR(cNum(i)) .GT. iNine))) THEN
             print*, '*********************************'
             print*, 'ERROR, input data was not numeric'
             stop 
          END IF
c
c
C Flag myself if there is a number after the dot (matissa)
C and take care of the number according to the part it's 
C belong to
          IF (ICHAR(cNum(i)) .EQ. iDot) THEN
             iflag_dot = 1
             i = i+1
             CONTINUE
          END IF
c
          IF (iflag_dot .EQ. 0) THEN
             rNum = rNum*10
             r_tmpNum = ICHAR (cNum(i))
             r_tmpNum = r_tmpNum - iZero
          ELSE
             if (cNum(i) .EQ.  ' ')  GOTO 30
             iPower = iPower +1
             r_tmpNum = REAL(ICHAR (cNum(i)))
             r_tmpNum = r_tmpNum - iZero
             r_tmpNum =  r_tmpNum *(1.0/(10.0**iPower))
          END IF
          i = i+1
          rNum = rNum + r_tmpNum
          r_tmpNum = 0.0
 30    END DO
c
c
C Multiply the real number by it's sign (1 or -1)      
       rNum = rNum*isign      
c
c
       RETURN
       END
