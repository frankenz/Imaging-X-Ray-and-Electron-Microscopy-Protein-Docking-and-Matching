C************************************************************************
      SUBROUTINE WRTC (XB,TRNS,XYZPOT,INAT)
C************************************************************************
      INCLUDE 'param.h'
      INCLUDE 'title.h'
      INCLUDE 'blk2_3.h'
      INCLUDE 'centr.h'
C
      CHARACTER*6 AAA
      CHARACTER*3 RES,CRES
      CHARACTER*3 HT1
      CHARACTER*2 H1
      CHARACTER*1 ATM1
      CHARACTER*1 ATM2,ATM3
      REAL ONE,DENS

      DIMENSION XB(MAXATM,5),X(3),TRNS(3),XYZPOT(MAXPOT,4)
C
      REWIND INAT
c 100  FORMAT (A6,2X,I5,2X,3A1,1X,A3,2X,A4,4X,3F8.3,2X,'1.00',1X,F5.2)
 100  FORMAT (A6,I5,2X,A2,2X,A3,1X,A1,I4,4X,3F8.3,2X,F4.2,F6.1)
      DO IB=1,NB
        READ (INAT,100,END=1000) AAA,ILST,H1,HT1,CI,IRS,
     1                           (X(I),I=1,3),ONE,DENS
        DO I=1,3
          X(I)=XB(IB,I)+TRNS(I)+CENA(I)
        END DO   
        WRITE (32,100) AAA,ILST,H1,HT1,CI,IRS,
     1                           (X(I),I=1,3),ONE,DENS
      END DO
 1000 CONTINUE
      IF (IPOT.EQ.1) THEN
        INFL=62
        REWIND INFL
        DO IB=1,NPOTB
          READ (INFL,100,END=1010) AAA,ILST,ATM1,ATM2,ATM3,RES,IRS,
     1                           (X(I),I=1,3)
          DO I=1,3
            X(I)=X(I)+TRNS(I)
          END DO
          WRITE (63,301) ILST,ATM1,ATM2,ATM3,RES,ILST,
     1                           (X(I),I=1,3)
        END DO
      END IF
300   format('HETATM',i5,2x,'A1 ',1x,'POS',i6,4x,f8.3,f8.3,f8.3,
     1        2x,f10.5)
305   format('HETATM',i5,2x,'T2 ',1x,'NEG',i6,4x,f8.3,f8.3,f8.3,
     1        2x,f10.5)
301   format('HETATM',i5,2x,3A1,1x,A3,I6,4x,3f8.3,2x,f10.5)    
 1010 RETURN
      END
