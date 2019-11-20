C-----------------------------------------------------------------       
       subroutine residueID (cNum, rNum, index, cid, rid)
C------------------------------------------------------------------
c       
c
C This subroutine gets the residue parameter that the user has
c inserted. Puts the chain identifier and the residue identifier 
C in the right arrays. Makes sure that the residue num contains
C only number and send it to charToReal functions.
C The subroutine returns the residue number as real number to 
C the calling routine.
c
C------------------------------------------------------------------
c
       include 'table.h'
C
      REAL rNUM
      INTEGER index,I 
      CHARACTER*1 cNum(70), tmpNum(70), cSpace, cMinus 
      CHARACTER*1 CID(20),RID(20)
c
c
C Initialize the tmpNum string before reading the parameter
      DO i=1, 15 
         tmpNum(i)=''
      END DO
C
C Initialize variable cSpace with a blank space
      cSpace = " "
      cMinus = "-"
C
C Insert the chain identifier and into the correct global array.
C (In case there was no chain identifier insert a backslash space
C into the array so it can be compared later on in readatom - when
C omitted atoms are marked).
       CID(index) = cNUM(1)
C       
C Copy the residue number into a new char variable that will 
C contains only the number (withour the residue and chain
C identifier
       i=2
       DO WHILE (cNUM(i) .NE. cSpace)
          tmpNUM(i-1) = cNUM(i)
          i=i+1
       ENDDO

C Copy the residue identifier sign into the correct array and
C erase it from the number string.  
C (In case there was no residue identifier insert a backslash space
C into the array so it can be compared later on in readatom - when
C omitted atoms are marked).
       RID(index) = tmpNum(i-2)
       tmpNUM(i-2) = cSpace
C
C
c Call charToReal that will convet the residue number string
C into a real number.
      call charToReal(tmpNum, rNUM)
c
c
       RETURN
       END
