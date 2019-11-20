C**********************************************************************
      SUBROUTINE READDATA
C**********************************************************************
C READ INSTRUCTIONS IN FREE FORMAT
C
      INCLUDE 'blk2_3.h'
      INCLUDE 'title.h'
      INCLUDE 'main.h'
      INCLUDE 'newatm.h'
      INCLUDE 'table.h'
      INCLUDE 'param.h'
C
      CHARACTER*4 titleCode
      CHARACTER*3 CRES
      CHARACTER*1 tmpNum(15), cSpace, cNext_line, IR(76)
      CHARACTER*1 cMinus, CID(MAXOMIT),RID(MAXOMIT)
      INTEGER address_code,num_params,iSpace,iFix_params,idx
      INTEGER idx2,idx3,idx4
      INTEGER rNUM(MAXOMIT)
      REAL tmp_num_param
C
      DIMENSION AA(200)
C
   10 FORMAT(A4,76A1)
   20 FORMAT(1X,A4,76A1)
   25 FORMAT(76A1)
 26   FORMAt(1X,A4,' CONTINUE ',76A1)
C
c
c
      iSpace = ICHAR(' ')
      cSpace = ' '
      cMinus = '-'
      cNext_line  ='='
      iFix_params =1
c
C Initialize the tmpNum string before reading the parameter
      DO idx3=1, 15 
         tmpNum(idx3)=''
      END DO
c
c
c
C READ IN FREE FORMAT
C
 30   READ(5,10,END=1500) titleCode, IR
      WRITE (6,20) titleCode,IR
c
C Initialize indexes before interpreting the line
      idx1=1
      idx2=1
      idx3=1
C
C Iitialize AA real array before each iteration
      DO J=1,200
         AA(J) = 0.
      END DO
C
C Clean the arrays that hold the residues details before inserting
C details from a new data line.
      DO J=1,MAXOMIT
         CID(J) = ''
         RID(J) = ''
         rnum(J)= 0
      ENDDO

c
c
C Find the index of the 4 letter code, and use this index in order to
C find the appropriate line format and the number of parameters
C in this category. If the function findIdx returned -1, the title-
C code was not found. Them go and read the next line in the file.
 35   CALL findIdx (titleCode, idx)
c
      IF (idx .EQ. -1) THEN
         GO TO 30
      ELSE
         address_code = convert_tbl(idx,1)
         num_params = convert_tbl(idx,2)
      END IF
c
c
C Check number of parameters: 
C If it is 0, then go to the address code directly.
C If it is a number >0 then read the rest of the parameters and 
C    jump to the address code.
C If it is -1, then translate the first number in the parameter list,
C    convert it into a real number, use this number as the number of 
C    the following params. Read the other parameters and jump to the 
C    line with the appropriate address code.
c
      IF (num_params .LT. 0) THEN
         iFix_params =0
         DO WHILE (IR(idx2) .EQ. cSpace ) 
            idx2=idx2+1
         END DO
         DO WHILE (IR(idx2) .NE. cSpace)
            tmpNum(idx3)= IR(idx2)
            idx2 = idx2+1
            idx3 = idx3+1
         END DO
         CALL charToReal(tmpNum,AA(1))
         num_params = INT (AA(1))
         
C Since SITA requiers to read also the weight we had to add 1 
C to num_params.
	 IF (idx .EQ. 5 .OR. idx .EQ. 9) THEN
	     num_params = num_params + 1 
	 END IF 

         idx3=1
         GO TO 390
      ELSE IF (num_params .GT. 0) THEN
         iFix_params =1
         GOTO 390
      ELSE 
         GOTO 100
      END IF
c
C Read the parameters from the IR string.
C Read char-char from IR and translate the chars into real numbers.
 390  DO idx4=1, num_params
         DO WHILE (IR(idx2) .EQ.  cSpace) 
            idx2=idx2+1
         END DO
c
C Take care of the case the input parameters data continue in the 
C next line in the dat file
            IF (IR(idx2) .EQ. cNext_line) THEN
               READ(5,25,END=1500) IR
               WRITE (6,26),titleCode,IR
               Idx2 = 1
               DO WHILE (IR(idx2) .EQ.  cSpace) 
                  idx2=idx2+1
               END DO
            END IF
c
c
         DO WHILE (IR(idx2) .NE. cSpace)

            tmpNum(idx3)= IR(idx2)
            idx2 = idx2+1
            idx3 = idx3+1
         END DO
c
         IF (iFix_params .EQ. 0) THEN

C If OMTA or OMTB options are on than the translating of the
C residues should be done by residueID subroutine because of
C the chain and residue identifier.
C The same for SITA and SITB options
            IF ((titleCode .EQ.'OMTA') .OR. (titleCode .EQ.'OMTB').OR.
     1         (titleCode .EQ.'SITA').OR. (titleCode .EQ.'SITB')) THEN
               IF ((num_params .GT. MAXOMIT) .AND.
     1              ((titleCode .EQ.'OMTA') .OR. 
     2              (titleCode .EQ.'OMTB'))) THEN
                  print*, 'ERROR - THE MAXIMUM NUMBER OF ALLOWED',
     1                 'OMITTED RESIDUES IS 20'
                  STOP
               ENDIF
C     If it is a site option, I should take care not to read
C     the weight in residue format.
               IF ((idx4 .EQ. num_params) .AND.((titleCode .EQ.'SITA')
     1              .OR. (titleCode .EQ.'SITB'))) THEN
                  call charToReal(tmpNum, AA(idx4+1))
               ELSE
                  call residueID(tmpNum,AA(idx4+1),idx4, CID,RID)
                  rNUM(idx4) = AA(idx4+1)
               ENDIF
            ELSE
               call charToReal(tmpNum, AA(idx4+1))
            ENDIF
         ELSE
            call charToReal(tmpNum, AA(idx4))
         END IF
c
C "Clean"(initialize) the tmpNum string before reading the parameter
         DO idx3=1, 15 
            tmpNum(idx3)=''
         END DO
         idx3 = 1 
c
      END DO
c
c    
  100 GO TO(400,410,420,430,440,450,460,470,480,520,570,630,700,710,720
     1 ,730,740,750,760,800,830,840,850,860,870,880,890,900,910,920,930
     2 ,940,950,1000), idx
c
C
C OMTA
C
  400 CONTINUE

c      Print*, 'OMTA'
      IOMTA=1

C insert the residue dtails into the right arrays.
      DO J=1, num_params
         CID_omitA(J) = CID(J)
         rNUM_omitA(J) = rNUM(J)
         RID_omitA(J) = RID(J)
         
C Convert the "-" with a space " " sign so the chain id and the residue
C id will be comparable to the atom identiers that are ead from
C the ccordinates files.
         IF (CID_omitA(J) .EQ. cMinus) CID_omitA(J) = cSpace
         IF (RID_omitA(J) .EQ. cMinus) RID_omitA(J) = cSpace
      ENDDO
      GO TO 30
C
C ELEC
C
  410 CONTINUE
c      Print*, 'pot'
      IPOT=1
      PP1=AA(1)
      WF=AA(2)
      GO TO 30
C
C TRIM
C
  420 CONTINUE

c      Print*, 'TRIM'
C Initialize the number of omitted residues types and the weight
C that should be given to the omitted atoms from the AA aray
C into the correct parameters.
      KCUT=NINT(AA(1))
      XCUT=AA(2)
C Read the next line from the *.dat file in formated manner.
C Insert the residues type to be trimmed into array CRES, 
C and from where in the side chain to start cutting.
      IF (KCUT.GT.0) THEN 
         READ (5,421) (CRES(J),NCRES(J),J=1,KCUT)
         print 421,(CRES(J),NCRES(J),J=1,KCUT)
      END IF

  421 FORMAT(4X,10(A3,I3,1X))
C
      GO TO 30
C
C FMTA
C
  430 CONTINUE


c      Print*, 'FMTA'

c      DO 435 J=1,72
c  435 FRMTA(J)=IR(J)
      GO TO 30
C
C SITA
C
  440 CONTINUE
c      Print*, 'sita'
C Count the number of site options per molecule A.
C If there are more than the allowed number print a warning messege.
      NSITA = NSITA+1
      IF (NSITA .GT. MAXSITE) THEN
         PRINT*,'WARNING: ONLY 5 SITE OPTIONS PER MOLECULE ARE ALLOWED'
         GO TO 30
      ENDIF
C
C Insert the number of residues into the right variable (with the 
C right index).
      ISITA(NSITA) = NINT(AA(1))
C Insert the residues details into arrays.
      DO J=1,ISITA(NSITA)
         CID_sitA(NSITA,J)=CID(J)
         IRESITA(NSITA,J)=rNUM(J)
         RID_sitA(NSITA,J)=RID(J)
c         IRESITA(NSITA,J)=NINT(AA(J+1))
C Convert the "-" with a space " " sign so the chain id and the residue
C id will be comparable to the atom identiers that are ead from
C the ccordinates files.
         IF (CID_sitA(NSITA,J) .EQ. cMinus) CID_sitA(NSITA,J) = cSpace
         IF (RID_sitA(NSITA,J) .EQ. cMinus) RID_sitA(NSITA,J) = cSpace
      END DO
C
C Insert the weight that should be given to the residues into the
C weight array at the correct index.
      XSITA(NSITA)=AA(ISITA(NSITA)+2)
      GO TO 30
C
C STEP
C
  450 CONTINUE

c      Print*, 'step'

      ISTEP=1
      STEPSZ=AA(1)
      GO TO 30
C
C FMTB
C
  460 CONTINUE

c      Print*, 'FMTB'

      DO 465 J=1,72
  465 FRMTB(J)=IR(J)
      GO TO 30
C
C ORIG
C
  470 CONTINUE


c      Print*, 'orig'


      IORIGN=1
      DO 472 II=1,3
  472 ORIG(II)=AA(II)
      GO TO 30
C
C SITB
C
  480 CONTINUE
c      Print*, 'sitb'
      NSITB = NSITB+1
      IF (NSITB .GT. MAXSITE) THEN
         PRINT*,'WARNING: ONLY 5 SITE OPTIONS PER MOLECULE ARE ALLOWED'
         GO TO 30
      ENDIF
C
C Insert the number of residues into the right variable (with the 
C right index).
      ISITB(NSITB)=NINT(AA(1))
      DO J=1,ISITB(NSITB)
C Insert the residues details into arrays.
         CID_sitB(NSITB,J)=CID(J)
         IRESITB(NSITB,J)=rNUM(J)
         RID_sitB(NSITB,J)=RID(J)
c        IRESITB(NSITB,J)=NINT(AA(J+1))
C Convert the "-" with a space " " sign so the chain id and the residue
C id will be comparable to the atom identiers that are ead from
C the ccordinates files.
         IF (CID_sitB(NSITB,J) .EQ. cMinus) CID_sitB(NSITB,J) = cSpace
         IF (RID_sitB(NSITB,J) .EQ. cMinus) RID_sitB(NSITB,J) = cSpace
      END DO
C
C Insert the weight that should be given to the residues into the
C weight array at the correct index.
      XSITB(NSITB)=AA(ISITB(NSITB)+2)
      GO TO 30
C
C ROT
C
  520 CONTINUE

c      Print*, 'rot'


      IROT=1
      ANG1=AA(1)
      ANG2=AA(2)
      ANG3=AA(3)
      DANG=AA(4)
      GO TO 30
C
C TITL
C
  570 DO 571 J=1,76
  571 TITLE(J)=IR(J)

c      Print*, 'title'
      GO TO 30
C
C NROT
C
  630 CONTINUE

c      Print*, 'nrot'

      IIROT=NINT(AA(4))
      DO 632 II=1,3
      JROT(II)=NINT(AA(II))
  632 CONTINUE
      GO TO 30
C
C SURF
C
  700 CONTINUE

c      Print*, 'surf'

      NSURF=1
      NP=NINT(AA(1))
      if (NP .GT. 432) then
         print*, 'ERROR, THE SURFACE PARAMETER IS TOO LARGE.'
         print*, 'THE MAXIMUM SIZE ALLOWED IS 432'
         STOP
      end if
      GO TO 30
C
C PEAK
C
  710 CONTINUE

c      Print*, 'peak'

      NPEAK=NINT(AA(1))
      TPEAK=AA(2)
      GO TO 30
C
C OUTC
C
  720 CONTINUE

c      Print*, 'outc'

      IOUTC=1
      DO 725 II=1,3
      ROTT(II)=AA(II)
      ITRN(II)=NINT(AA(II+3))
  725 CONTINUE
      GO TO 30
C
C INTA
C
  730 CONTINUE

c      Print*, 'intA'

      IINTA=NINT(AA(1))
      DO II=1,IINTA
        NINTA(II)=NINT(AA(II+1))
      END DO
      GO TO 30
C
C INTB
C
  740 CONTINUE

c      Print*, 'intB'
      IINTB=NINT(AA(1))
      DO II=1,IINTB
        NINTB(II)=NINT(AA(II+1))
      END DO
      GO TO 30
C
C OMTB
C
  750 CONTINUE

c      Print*, 'OMTB'
      IOMTB=1

C insert the residue dtails into the right arrays.
      DO J=1, num_params
         CID_omitB(J) = CID(J)
         rNUM_omitB(J) = rNUM(J)
         RID_omitB(J) = RID(J)
       
C Convert the "-" with a space " " sign so the chain id and the residue
C id will be comparable to the atom identiers that are ead from
C the ccordinates files.
         IF (CID_omitB(J) .EQ. cMinus) CID_omitB(J) = cSpace
         IF (RID_omitB(J) .EQ. cMinus) RID_omitB(J) = cSpace
      ENDDO
C
      GO TO 30
C
C RHO
C
  760 CONTINUE
      RHO=AA(1)

c      print*, RHO

      GO TO 30
C
C SYMM
C
  800 CONTINUE

c      Print*, 'symm'

      ISYMM=NINT(AA(1))
      DSYMM=AA(2)
c     print *,'just read in SYMM in xdred'
c     print *,'isymm=',isymm
c     print *,'dsymm=',dsymm
      GO TO 30
C
C NTRM
C
  830 CONTINUE
C NTRIM - contains the number of the molecule that the user do
C not want to trim.
      NTRIM=NINT(AA(1))
      GO TO 30
C
C RESL
C
  840 CONTINUE
c      Print*, 'resl'

      XRESL = AA(1)
      DT = AA(2)
      GO TO 30
C
C ATOM
C
  850 CONTINUE

c      Print*, 'atom'

      NMOL=NINT(AA(1))
      INEW=NEWAT(NMOL)+1
      NEWAT(NMOL)=INEW
      LISNEW(INEW,NMOL)=NINT(AA(2))
      RVNEW(INEW,NMOL)=AA(3)
      GO TO 30
C
C HDRF
C
  860 CONTINUE
      HD=AA(1)
      IHDRF=1
      GO TO 30
C
C REFN
C
  870 CONTINUE
      IREFN = NINT(AA(1))
      GO TO 30
C
C GRID
C
  880 CONTINUE
      XXM = AA(1)
      GO TO 30
C 
C AD44
C
  890 CONTINUE
      GO TO 30
C
C REST
C
  900 IREST=1
      TREST=AA(1)
      PHREST=AA(2)
      PSREST=AA(3)

c      Print*, 'rest'


      GO TO 30
C
C ROTA
C
  910 IROTA=1
      DO IX=1,3 
        ANGSA(IX)=AA(IX)
      END DO

c      Print*, 'roTA'

      GO TO 30

C
C AD47
C
  920 CONTINUE
      GO TO 30
C
C AD48
C
  930 GO TO 30
C
C AD49
C
  940 GO TO 30
C
C AD50
C
  950 GO TO 30

C
C '    '
C
  960 GO TO 30


 1000 CONTINUE
 1100 RETURN
 1500 PRINT 1501
 1501 FORMAT (1X,'END OF INPUT CARDS')
      print*, '==============================='


      PRINT* 
      RETURN
      END
