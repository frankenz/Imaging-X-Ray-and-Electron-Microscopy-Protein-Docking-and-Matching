C-----------------------------------------------------------------------
      SUBROUTINE READATOM (INAT,NMOL,IATOM,X,IHD)
C-----------------------------------------------------------------------
C
C ROUTINE TO READ AND RECOGNIZE ATOMS BY THEIR PDB RESIDUE AND ATOM NAME
C
C
      INCLUDE 'param.h'
      INCLUDE 'title.h'
      INCLUDE 'newatm.h'
      INCLUDE 'blk2_3.h'
      INCLUDE 'table.h'
C     
c      CHARACTER*4 AA,AAA,HA,TR,IRS,IRSPRV
      CHARACTER*4 AA,AAA,HA,TR,IRSPRV
      CHARACTER*3 RESNAM,RES,CRES
      CHARACTER*1 AT1(5),AT2(18),AT3(6),ATM1,ATM2,ATM3
      CHARACTER*1 RESID,CHAINID
      INTEGER index, IRS 
      INTEGER NATHDRF
c     DIMENSION RESNAM(20),X(MAXATM,4)
      DIMENSION RESNAM(24),X(MAXATM,5),IHD(MAXATM)
      DIMENSION ICRES(20)
C     
      DATA RESNAM/'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY',
     1             'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
     2             'THR','TRP','TYR','VAL','A  ','C  ','G  ','T  '/
      DATA AT3/'1','2','3','*','P','M'/
      DATA AT2/' ','A','B','G','D','E','Z','H','X','1','2','3','4',
     1'5','6','7','8','9'/
      DATA AT1/'O','N','C','S','P'/
      DATA AA/'ATOM'/,HA/'HETA'/,TR/'TER '/
      DATA RO/1.30/,RS/2.1/
      DATA RCHn/2.00/,RCHR/1.90/
      DATA RN/1.392/
      DATA RC/1.46/
      DATA RP/2.20/
      DATA RNH/1.54/,ROH/1.50/
C
      REWIND INAT
C     !
      DO I=1,MAXATM
         X(I,5)=0.
      END DO
C
      ITR=1
      IATOM=0
      ICATOM=0
      IF (NMOL.EQ.1) NATSITA=1
      IF (NMOL.EQ.2) NATSITB=1
C
C Initialize the hydrophbicity counter.
      NATHDRF=0
C
C initialize the counter of OMITTED atoms.
      NST=0
      FOUND = 1
C
 1    IATOM=IATOM+1
      RA=1.6
C
C This format read also the chain identifier and the residue 
C identifier
      READ (INAT,15,END=3000) AAA,ILST,ATM1,ATM2,ATM3,RES,CHAINID,
     1     IRS,RESID,(X(IATOM,I),I=1,3) 
 15   FORMAT (A4,2X,I5,2X,3A1,1X,A3,1X,A1,I4,A1,3X,3F8.3)
C
C
C
C CHECK IF THIS IS AN ATOM RECORD
C
      IF (AAA.NE.AA.AND.AAA.NE.HA) THEN
         IATOM=IATOM-1
         GO TO 1
      ELSE
         IF (AAA.EQ.HA) GO TO 2970
      END IF
C     
C CHECK IF THIS RESIDUE IS ONE OF THE 'SITE' RESIDUES OF MOLECULE 1 OR 2
C     
      IF (NMOL.EQ.1 .AND. NSITA.GT.0) THEN
         DO J=1,NSITA
            DO I=1,ISITA(J)
               IF (IRS.EQ.IRESITA(J,I)) THEN
                  IF ((CID_sitA(J,I) .EQ. CHAINID) .AND.
     1                 (RID_sitA(J,I) .EQ. RESID)) THEN
                     IATSITA(NATSITA,1)=IATOM
                     IATSITA(NATSITA,2)=J
                     NATSITA=NATSITA+1
                  END IF
               END IF
            END DO
         END DO
      END IF
      IF (NMOL.EQ.2 .AND. NSITB.GT.0) THEN
         DO J=1,NSITB
            DO I=1,ISITB(J)
               IF (IRS.EQ.IRESITB(J,I)) THEN
                  IF ((CID_sitB(J,I) .EQ. CHAINID) .AND.
     1                 (RID_sitB(J,I) .EQ. RESID)) THEN
                     IATSITB(NATSITB,1)=IATOM
                     IATSITB(NATSITB,2)=J
                     NATSITB=NATSITB+1
                  END IF
               END IF
            END DO
         END DO
      END IF
C     
C IDENTIFY THE RESIDUE TYPE
C     
      IR=0
      DO I=1,24
         IF (RES.EQ.RESNAM(I)) IR=I
      END DO
      IF (IR.EQ.0) THEN
        PRINT 17,RES
   17 FORMAT (2X,'*** RESIDUE ',A3,' NOT RECOGNIZED; STANDARD VALUES
     1 APPLIED. ***')
        GO TO 2950
      END IF
C
C Check if the residue is hydrophobic (based on the residue code
C in IR). If the residue is hydrophobic flag IATOM = 1
C Currently only protein residues can be hydrophobic
       IF (IHDRF.EQ.1) THEN
          IF (IR.EQ.1.OR.IR.EQ.5.OR.IR.EQ.10.OR.IR.EQ.11.OR.IR.EQ.13
     1         .OR.IR.EQ.14.OR
     1         .IR.EQ.15.OR.IR.EQ.18.OR.IR.EQ.19.OR.IR.EQ.20) THEN
             IHD(IATOM)=1
          ELSE
             IHD(IATOM)=0
          END IF
       END IF
C     
C IDENTIFY THE ATOM TYPE
C     
       IA=0
       DO I=1,5
          IF (ATM1.EQ.AT1(I)) IA=I
       END DO
       IB=0
       DO I=1,18
          IF (ATM2.EQ.AT2(I)) THEN
             IB=I
             II=I-2
          END IF
       END DO
       IF (IA.EQ.0.OR.IB.EQ.0) THEN
         PRINT 22,ATM1,ATM2,ATM3,RES,NMOL
 22      FORMAT (2X,'ATOM ',3A1,' IN RESIDUE ',A3,' IN MOLECULE ',I1 ,
     1        ' UNKNOWN; STANDARD VALUES APPLIED.')

c         PRINT 22,ATM1,ATM2,ATM3,RES
c 22      FORMAT (2X,'ATOM ',3A1,' IN RESIDUE ',A3,' UNKNOWN; STANDARD
c     1   VALUES APPLIED.')
         GO TO 2950
      END IF
C     
C Identifying protein backbone atoms
      IF (IR.LE.20.AND.(IB.LE.2.OR.IB.EQ.9)) THEN
         IHD(IATOM)=0
         GO TO (30,40,50,60) IA
c     backbone carbonyl oxygen or CO2- terminus
 30      RA=RO
c        IF (IB.EQ.9) THEN
c           ITR=2
C     NAT1=0
c        END IF
         GO TO 2950
c     backbone NH nitrogen or N terminus (when itr.gt.0)
 40      RA=RNH
c        IF (ITR.GT.0) ITR=0
         GO TO 2950
c     backbone carbonyl carbon and C(alpha)H group
 50      IF (IB.EQ.1) THEN
c     carbonyl C
            RA=RC
         ELSE
c     C(alpha)
            RA=RCHn
         END IF
         GO TO 2950
 60      CONTINUE
      END IF
C     
C Identifying DNA backbone atoms
C
      IF (IR.GT.20.and.IB.EQ.1.OR.((IB.GE.10.and.ib.le.14)
     1   .and.(atm3.eq.at3(4).or.atm3.eq.at3(5)))) THEN
         IHD(IATOM)=0
         GO TO (31,41,51,61,71) IA
c     backbone oxygen O3',O4',O5',O1P,O2P
 31      RA=RO
         GO TO 2950
 41      continue
c     backbone carbon C1',C2',C3',C4',C5' (all are CHn)
 51      RA=RCHn
         GO TO 2950
 61      CONTINUE
c     backbone phosphate atom
 71      RA=RP
         go to 2950
      END IF
C
C     CHECK IF THIS RESIDUE IS IN THE 'TRIM' LIST
C
      IF (KCUT.GT.0) THEN
         IF (NMOL .EQ. NTRIM) THEN
            CONTINUE
         ELSE
            DO IJ=1,KCUT
               IF (RES.EQ.CRES(IJ)) THEN
                  IFLAG=0
C     CHECK IF THIS RESIDUE IS NOT IN THE LIST OF INTERIOR RESIDUES
                  IF (NMOL.EQ.1.AND.IINTA.GT.0) THEN
                     DO IJ1=1,IINTA
                        IF (IRS.EQ.NINTA(IJ1)) IFLAG=1 
                     END DO
                     if(iflag.eq.1)print*,'RESIDUE ',IRS,' NOT 
     1 TRIMMED.'
                  END IF
                  IF (NMOL.EQ.2.AND.IINTB.GT.0) THEN
                     DO IJ1=1,IINTB
                        IF (IRS.EQ.NINTB(IJ1)) IFLAG=1 
                     END DO
                     if(iflag.eq.1)print*,'RESIDUE ',IRS,' NOT', 
     1                   ' TRIMMED.'
                  END IF               
C     
                  IF ((IFLAG.EQ.0) .AND. (II.GE.NCRES(IJ))) THEN
                     X(IATOM,5)=2.0
                     ICATOM=ICATOM+1
C     IF(NMOL.EQ.2)THEN
C     WRITE(22,FRMTA) AAA,ILST,ATM1,ATM2,ATM3,RES,IRS,
c     1                         (X(IATOM,I),I=1,3)
c     ENDIF
c..   GO TO 1
                  END IF
               END IF
            END DO
         END IF
      END IF
C     
      GO TO (1000,1100,1200,1300,1400,1500,1600,1700,
     1     1800,1900,2000,2100,2200,2300,2400,2500,
     1     2600,2700,2800,2900,2920,2930,2940,2945),IR
C     
C     ALA
C     
 1000 CONTINUE
      RA=RCHn
      GO TO 2950
C     
C     ARG
C     
 1100 CONTINUE
      GO TO (1110,1110,1110,1115,1120,1115) II
C     CB, CG, CD atoms)
 1110 RA=RCHn
      GO TO 2950
 1112 RA=RCHn
      GO TO 2950
C     NE, NH1 and NH2
 1115 RA=RNH
      GO TO 2950
C     CZ
 1120 RA=RC
      GO TO 2950
C     
C     ASN
C     
 1200 GO TO (1202,1205,1210) II
 1202 RA=RCHn
      GO TO 2950
 1205 RA=RC
      GO TO 2950
 1210 IF (IA.EQ.1) THEN
c     amide carbonyl O
         RA=RO
      ELSE
c     amide NH2 group
         RA=RNH
      END IF
      GO TO 2950
C     
C     ASP
C     
 1300 GO TO (1305,1310,1320) II
 1305 RA=RCHn
      GO TO 2950
 1310 RA=RC
      GO TO 2950
 1320 RA=RO
      GO TO 2950
C     
C     CYS
C     
 1400 GO TO (1410,1420) II
 1410 RA=RCHn
      GO TO 2950
 1420 RA=RS
      GO TO 2950
C     
C     GLN
C     
 1500 CONTINUE
      GO TO (1502,1502,1505,1510) II
 1502 RA=RCHn
      GO TO 2950
 1505 RA=RC
      GO TO 2950
 1510 IF (IA.EQ.1) THEN
c     amide carbonyl oxygen
         RA=RO
      ELSE
c     amide NH2 group
         RA=RNH
      END IF
      GO TO 2950
C     
C     GLU
C     
 1600 CONTINUE
      GO TO (1602,1602,1605,1610) II
 1602 RA=RCHn
      GO TO 2950
 1605 RA=RC
      GO TO 2950
 1610 RA=RO
      GO TO 2950
C     
C     GLY
C
 1700 GO TO 2950
C     
C     HIS
C     
 1800 GO TO (1812,1815,1820,1830) II
 1812 RA=RCHn
      GO TO 2950
 1815 RA=RC
      GO TO 2950
 1820 IF (ATM1.EQ.AT1(3)) THEN
c     aromatic CH (delta)
         RA=RCHR
      ELSE
c     delta NH group
         RA=RNH
      END IF
      GO TO 2950
 1830 IF (ATM1.EQ.AT1(3)) THEN
c     aromatic CH (epsilon)
         RA=RCHR
      ELSE
c     epsilon NH group
         RA=RNH
      END IF
      GO TO 2950
C     
C     ILE
C     
 1900 GO TO (1910,1910,1910) II
 1910 RA=RCHn
      GO TO 2950
C     
C     LEU
C     
 2000 GO TO (2010,2010,2010) II
 2010 RA=RCHn
      GO TO 2950
C     
C LYS
C
 2100 GO TO (2105,2105,2105,2105,2120) II
 2105 RA=RCHn
      GO TO 2950
 2120 RA=RNH
      GO TO 2950
C
C METHIONIN - MET.
C
 2200 GO TO (2205,2205,2210,2205) II
 2205 RA=RCHn
      GO TO 2950
 2210 RA=1.74
      GO TO 2950
C
C PHE
C
 2300 GO TO (2305,2310,2315,2315,2315) II
 2305 RA=RCHn
      GO TO 2950
 2310 RA=RC
      GO TO 2950
 2315 RA=RCHR
      GO TO 2950
C     
C PRO
C
 2400 GO TO (2410,2410,2410) II
 2410 RA=RCHn
      GO TO 2950
C     
C     SER
C     
 2500 GO TO (2505,2510) II
 2505 RA=RCHn
      GO TO 2950
 2510 RA=ROH
      GO TO 2950
C
C THR
C
 2600 GO TO (2605,2610) II
 2605 RA=RCHn
      GO TO 2950
 2610 IF (ATM1.EQ.AT1(1)) THEN
         RA=ROH
      ELSE
	 RA=RCHn
      END IF
      GO TO 2950
C     
C TRP
C
 2700 GO TO (2702,2705,2710,2715,2720,2720) II
 2702 RA=RCHn
      GO TO 2950
 2705 RA=RC
      GO TO 2950
 2710 CONTINUE
c     aromatic CH group (delta)
      RA=RCHR
c     aromatic C (bridge carbon delta2)
      IF (ATM3.EQ.AT3(2)) RA=RC
      GO TO 2950
 2715 continue
c aromatic CH group epsilon 3
      RA=RCHR
c aromatic C (bridge carbon epsilon 2)
      IF (ATM1.EQ.AT1(3).AND.ATM3.EQ.AT3(2))RA=RC
c     aromatic NH group
      IF (ATM1.EQ.AT1(2))RA=RNH
      GO TO 2950
c aromatic CH groups (zeta 1,2 and h)
 2720 RA=RCHR
      GO TO 2950
C
C TYR
C
 2800 GO TO (2801,2805,2803,2803,2805,2810) II
 2801 RA=RCHn
      GO TO 2950
 2803 RA=RCHR
      GO TO 2950
 2805 RA=RC
      GO TO 2950
 2810 RA=ROH
      GO TO 2950
C
C VAL
C
 2900 GO TO (2910,2910) II
 2910 RA=RCHn
      GO TO 2950
C
C A (adenine)
C
 2920 iade=ii-7
      go to (29201,29202,29203,29204,29204,29206,29207,29208,
     1       29209) iade
29201 RA=RNH
      go to 2950
29202 ra=rchr
      go to 2950
29203 ra=rn
      go to 2950
29204 ra=rc
      go to 2950
29206 if (ia.eq.3) ra=rc
      if (ia.eq.2) ra=rnh
      go to 2950
29207 ra=rn 
      go to 2950
29208 ra=rchr
      go to 2950
29209 ra=rn
      go to 2950
C
C C (cytosine)
C
 2930 icyt=ii-7
c     print *,'icyt=',icyt
      go to (29301,29302,29303,29304,29305,29305) icyt
29301 ra=rn
      go to 2950
29302 if (ia.eq.3) ra=rc
      if (ia.eq.1) ra=ro
      go to 2950
29303 ra=rnh
      go to 2950
29304 if (ia.eq.3) ra=rc
      if (ia.eq.2) ra=rnh
      go to 2950
29305 ra=rchr
      go to 2950
C
C G (guanine)
C
 2940 igua=ii-7
      go to (29401,29402,29403,29404,29404,29406,29407,
     1       29408,29403) igua
29401 ra=rnh
      go to 2950
29402 if (ia.eq.3) ra=rc
      if (ia.eq.2) ra=rnh
      go to 2950
29403 ra=rn
      go to 2950
29404 ra=rc
      go to 2950
29406 if (ia.eq.3) ra=rc
      if (ia.eq.1) ra=ro
      go to 2950
29407 ra=rn
      go to 2950
29408 ra=rchr
      go to 2950
C
C T (thymine)
C
 2945 ithy=ii-7
      go to (29451,29452,29453,29452,29455,29456) ithy
29451 ra=rn
      go to 2950
29452 if (ia.eq.3) ra=rc
      if (ia.eq.1) ra=ro
      go to 2950
29453 ra=rnh
      go to 2950
29455 ra=rc
      if (atm3.eq.at3(6)) ra=rchn
      go to 2950
29456 ra=rchr
      go to 2950
C
 2950 CONTINUE
c     if (ir.gt.20) print *,res,atm1,atm2,atm3,ra,ir,ia,ib,ii
C For the OMTA option
C REMARK - The flag for missing atom (X(I,5)=1) must overwrite
C the flag for TRIM atom (X(I,5)=2) therefore the omit test
C must be placed after the trim test.
C Check if there is an OMTA option and if this is molecule A
      IF ((IOMTA .EQ. 1) .AND. (NMOL .EQ. 1)) THEN

C If the OMTA option is on, loop on the arrays that contain the
C details of the residue for omittiming (chain identifier, residue
C number and residue identifier) and look if the atom belongs to
C these residues. If it is mark it.
         DO index=1, MAXOMIT

C Fisrt check it the residue number of the read atom is identical
C to one of the residue numbers in the omitted residue number array.
C If it is not, than continue. If it is, check if the chain
C and residue identifiers of the atom are identicle to the omitted
C residue details. If it is than mark that atom in array X(IATOM,5).
C by giving it a value of 1.
            IF (rNUM_omitA(index) .EQ.IRS) THEN
               IF ((CID_omitA(index) .EQ. CHAINID) .AND.
     1              (RID_omitA(index) .EQ. RESID))    THEN
                  X(IATOM,5)=1.
                  NST=NST+1
                  GOTO 2970
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C     
C For the OMTB option - see for details the above OMTA option
      IF ((IOMTB .EQ. 1) .AND. (NMOL .EQ. 2)) THEN
         DO index=1, MAXOMIT
            IF (rNUM_omitB(index) .EQ.IRS) THEN
               IF ((CID_omitB(index) .EQ. CHAINID) .AND.
     1              (RID_omitB(index) .EQ. RESID))    THEN
                  X(IATOM,5)=1.
                  NST=NST+1
                  GOTO 2970
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
C
 2970 X(IATOM,4)=RA
      IF (IHD(IATOM).EQ.1) NATHDRF=NATHDRF+1
      GO TO 1
 3000 CONTINUE
      IATOM=IATOM-1
C     SITE OPTION
      IF (NMOL.EQ.1) NATSITA=NATSITA-1
      IF (NMOL.EQ.2) NATSITB=NATSITB-1
C
      PRINT 3005,IATOM,NMOL
C
C GENERAL PRINTING MESSAGES
 3005 FORMAT (I5,' ATOMS READ FOR MOLECULE ',I1)
C
C HYDROPHOBISITY PRINTING
       IF (IHDRF.EQ.1) PRINT 3004, NATHDRF,NMOL
 3004  FORMAT (I5,' HYDROPHOBIC ATOMS FOR MOLECULE ',I1,/)
C
C     SITE OPTION
      IF (NMOL.EQ.1.AND.NSITA.GT.0) THEN
         DO I=1, NSITA
            PRINT 3006,ISITA(I),NMOL,XSITA(I)
c            PRINT 3007,(IATSITA(IA,1),IA=1,NATSITA)
         END DO
         PRINT 3008, NATSITA
      END IF
      IF (NMOL.EQ.2.AND.NSITB.GT.0) THEN
         DO I=1, NSITB
            IF (NSITA.NE.0) XSITB(1)=1.0
c         IF (NSITA.EQ.0) XSITB=XSITA
            PRINT 3006,ISITB(I),NMOL,XSITB(I)
c            PRINT 3007,(IATSITB(IB,1),IB=1,NATSITB)
         END DO
         PRINT 3008, NATSITB
      END IF
 3006 FORMAT(I3,' RESIDUES IN MOL ',I1,' WERE GIVEN EXTRA WEIGHT; T='
     1     ,F5.2)
 3008 FORMAT(' THE TOTAL NUMBER OF ATOMS GIVEN EXTRA WEIGHT IS: ',I5)
c 3006 FORMAT(I3,' RESIDUES IN MOL ',I1,' WERE GIVEN EXTRA WEIGHT; T='
c     1     ,F5.2,'; THIS INCLUDES',I5,' ATOMS:')
c 3007 FORMAT(2X,10I6)
C     
      NNEW=NEWAT(NMOL)
      IF (NNEW.GT.0) THEN
         DO N=1,NNEW
            K=LISNEW(N,NMOL)
            IF (K.GT.IATOM) THEN
               PRINT 3010,K
 3010          FORMAT (1X,'***ATOM INSTRUCTION FOR ATOM ',I4,
     *              ' IGNORED.***')
            ELSE
               X(K,4)=RVNEW(N,NMOL)
               PRINT 3020,K,X(K,4)
 3020          FORMAT (1X,'RADIUS OF ATOM ',I4,
     1              ' CHANGED TO',1X,F5.2,'.')
            END IF
         END DO
      END IF
C     
C     
C     For the OMIT option
C     Print the number of omitted atoms and the details of the 
C     residues that the user inserted.
      IF ((IOMTA .EQ. 1) .AND. (NMOL .EQ. 1)) THEN
         print 3024,NST,NMOL
         DO J=1,MAXOMIT
            IF (rNUM_omitA(J) .EQ. 0) GOTO 3030
            IF (CID_omitA(J) .EQ. ' ') CID_omitA(J) = '-'
            IF (RID_omitA(J) .EQ. ' ') RID_omitA(J) = '-'
            print 3025, CID_omitA(J), rNUM_omitA(J),RID_omitA(J)
         ENDDO
      ENDIF
      IF ((IOMTB .EQ. 1).AND. (NMOL .EQ. 2)) THEN
         print 3024,NST,NMOL
         DO J=1,MAXOMIT
            IF (rNUM_omitB(J) .EQ. 0) GOTO 3030
            IF (CID_omitB(J) .EQ. ' ') CID_omitB(J) = '-'
            IF (RID_omitB(J) .EQ. ' ') RID_omitB(J) = '-'
            print 3025, CID_omitB(J), rNUM_omitB(J),RID_omitB(J)
         ENDDO
      ENDIF
       
 3024 FORMAT(I5,' ATOMS EXCLUDED FROM MOLECULE ',i1, ': CHAIN RESIDUE')
 3025 FORMAT (40X,A1,3X,I4,A1)

C
C
c 3025  FORMAT(3F8.3,1X,f3.0)

 3030 CONTINUE
C!!   
C Print a mesege if the TRIM option was on but one of the molecules
C was not trimed (because of operating option NTRM also) or else
C Print the number of trimmed atoms
      IF (KCUT.GT.0) THEN
         IF (NMOL .NE. NTRIM) THEN
            PRINT 3040, NMOL, ICATOM
         ELSE
            PRINT 3050, NMOL
         ENDIF
      ENDIF


 3040 FORMAT (1X,'THE NUMBER OF TRIMMED ATOMS IN MOLECULE ',i2,
     1     ' IS: ',I5,/)
 3050 FORMAT (' MOLECULE ',I1, ' WAS NOT TRIMMED',/)
      
      RETURN
      END
      
