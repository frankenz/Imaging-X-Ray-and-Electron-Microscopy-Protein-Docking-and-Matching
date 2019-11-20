C************************************************************
      SUBROUTINE RENUM_PDB (NMOL)
C************************************************************
      CHARACTER*1 A
      DIMENSION A(80)
C     num1 is the starting residue number in the output minus 1
c      write (6,*) 'Enter number for 1st residue'
	IF(NMOL.EQ.1)THEN
	 III=110
	 IO=10
	ENDIF
	IF(NMOL.EQ.2)THEN
	 III=111
	 IO=11
	ENDIF
	num1=0
c     num1=400  
      numpre=0
   10 READ (III,12,END=100) (A(I),I=1,22),num,(a(i),i=27,72)
      if (num.ne.numpre) then
        numpre=num
        num1=num1+1
      end if
      WRITE(IO,12) (A(I),I=1,22),num1,(a(i),i=27,72)
      GO TO 10
   12 FORMAT (22A1,i4,46a1)
  100 continue
      RETURN
      END
