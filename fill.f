C**********************************************
      SUBROUTINE FILL (NP,SURF,ISURF,RHO)
C**********************************************
      INCLUDE 'param.h'
C
      COMPLEX*8 SURF(NP,NP,NP)
      DIMENSION ISURF(NP,NP,NP),LIST(10000)
C
      DO 11 K=1,NP
      DO 11 J=1,NP
      DO 11 I=1,NP
      ISURF(I,J,K)=NINT(REAL(SURF(I,J,K)))
   11 CONTINUE
C
      DO 40 K=1,NP
      DO 20 J=1,NP
      II=0
      DO 15 I=1,NP
      IF (ISURF(I,J,K).EQ.0) II=II+1
   15 CONTINUE
      IF (II.LT.NP) GO TO 20
      DO 17 I=1,NP
   17 ISURF(I,J,K)=2
   20 CONTINUE
   40 CONTINUE
C
      NEW=0
      ILT=0
      DO 3000 K=1,NP
          DO 2000 J=1,NP
            DO 1000 I=1,NP
            KNN=K
            JNN=J
            INN=I
            II=0
            IF (ISURF(INN,JNN,KNN).GT.0) GO TO 1000
            II=II+1
            CALL REG(NP,INN,JNN,KNN,II,IEDGE,LIST)
            IF (IEDGE.EQ.1) GO TO 950
            INN=I-1
 50         IF (ISURF(INN,JNN,KNN)-1) 60,100,950
 60         II=II+1
            CALL REG(NP,INN,JNN,KNN,II,IEDGE,LIST)
            IF (IEDGE.EQ.1) GO TO 950
            INN=INN-1
            GO TO 50
 100        INN=I+1
 120        IF (ISURF(INN,JNN,KNN)-1) 130,150,950
 130        II=II+1
            CALL REG(NP,INN,JNN,KNN,II,IEDGE,LIST)
            IF (IEDGE.EQ.1) GO TO 950
            INN=INN+1
            GO TO 120
 150        INN=I
            JNN=J-1
 160        IF (ISURF(INN,JNN,KNN)-1) 170,200,950
 170        II=II+1
            CALL REG(NP,INN,JNN,KNN,II,IEDGE,LIST)
            IF (IEDGE.EQ.1) GO TO 950
            JNN=JNN-1
            GO TO 160
 200        JNN=J+1
 210        IF (ISURF(INN,JNN,KNN)-1) 220,250,950
 220        II=II+1
            CALL REG(NP,INN,JNN,KNN,II,IEDGE,LIST)
            IF (IEDGE.EQ.1) GO TO 950
            JNN=JNN+1
            GO TO 210
 250        JNN=J
            KNN=K-1
 260        IF (ISURF(INN,JNN,KNN)-1) 270,300,950
 270        II=II+1
            CALL REG(NP,INN,JNN,KNN,II,IEDGE,LIST)
            IF (IEDGE.EQ.1) GO TO 950
            KNN=KNN-1
            GO TO 260
 300        KNN=K+1
 310        IF (ISURF(INN,JNN,KNN)-1) 320,350,950
 320        II=II+1
            CALL REG(NP,INN,JNN,KNN,II,IEDGE,LIST)
            IF (IEDGE.EQ.1) GO TO 950
            KNN=KNN+1
            GO TO 310
 350        CONTINUE
C IF WE GET HERE IT COULD BE A HOLE
            IFT=2
            ILT=II
            IF (II.EQ.1) GO TO 650
 400        DO 600 IP=IFT,ILT
              CALL UNPECK (LIST,IP,IN,JN,KN)
              IF (ISURF(IN-1,JN,KN)-1) 405,430,950
 405          ILIST=1000000*(IN-1)+1000*JN+KN
                DO 410 IIP=1,ILT+NEW
                IF (ILIST.EQ.LIST(IIP)) GO TO 430
 410            CONTINUE
              NEW=NEW+1
              CALL REG(NP,IN-1,JN,KN,ILT+NEW,IEDGE,LIST)
              IF (IEDGE.EQ.1) GO TO 950
 430          IF (ISURF(IN+1,JN,KN)-1) 435,460,950
 435          ILIST=1000000*(IN+1)+1000*JN+KN
                DO 440 IIP=1,ILT+NEW
                IF (ILIST.EQ.LIST(IIP)) GO TO 460
 440            CONTINUE
              NEW=NEW+1
              CALL REG(NP,IN+1,JN,KN,ILT+NEW,IEDGE,LIST)
              IF (IEDGE.EQ.1) GO TO 950
 460          IF (ISURF(IN,JN-1,KN)-1) 465,490,950
 465          ILIST=1000000*IN+1000*(JN-1)+KN
                DO 470 IIP=1,ILT+NEW
                IF (ILIST.EQ.LIST(IIP)) GO TO 490
 470            CONTINUE
              NEW=NEW+1
              CALL REG(NP,IN,JN-1,KN,ILT+NEW,IEDGE,LIST)
              IF (IEDGE.EQ.1) GO TO 950
 490          IF (ISURF(IN,JN+1,KN)-1) 495,520,950
 495          ILIST=1000000*IN+1000*(JN+1)+KN
                DO 500 IIP=1,ILT+NEW
                IF (ILIST.EQ.LIST(IIP)) GO TO 520
 500            CONTINUE
              NEW=NEW+1
              CALL REG(NP,IN,JN+1,KN,ILT+NEW,IEDGE,LIST)
              IF (IEDGE.EQ.1) GO TO 950
 520          IF (ISURF(IN,JN,KN-1)-1) 525,550,950
 525          ILIST=1000000*IN+1000*JN+KN-1
                DO 530 IIP=1,ILT+NEW
                IF (ILIST.EQ.LIST(IIP)) GO TO 550
 530            CONTINUE
              NEW=NEW+1
              CALL REG(NP,IN,JN,KN-1,ILT+NEW,IEDGE,LIST)
              IF (IEDGE.EQ.1) GO TO 950
 550          IF (ISURF(IN,JN,KN+1)-1) 555,600,950
 555          ILIST=1000000*IN+1000*JN+KN+1
                DO 560 IIP=1,ILT+NEW
                IF (ILIST.EQ.LIST(IIP)) GO TO 600
 560            CONTINUE
              NEW=NEW+1
              CALL REG(NP,IN,JN,KN+1,ILT+NEW,IEDGE,LIST)
              IF (IEDGE.EQ.1) GO TO 950
 600        CONTINUE
C  CHECK IF NEW POINTS WERE DETECTED IN THIS HOLE
            IF (NEW.EQ.0) GO TO 650
            IFT=ILT+1
            ILT=ILT+NEW
            NEW=0
            GO TO 400
 650        CONTINUE
            DO 680 IP=1,ILT
              CALL UNPECK (LIST,IP,IN,JN,KN)
              ISURF(IN,JN,KN)=1
 680        CONTINUE
            GO TO 1000
 950        CONTINUE
            IF (ILT.EQ.0) ILT=II
            ILT=ILT+NEW
            IF (ILT.EQ.0) GO TO 1000
            DO 980 IP=1,ILT
              CALL UNPECK (LIST,IP,IN,JN,KN)
              ISURF(IN,JN,KN)=2
 980        CONTINUE
            ILT=0
            NEW=0
1000      CONTINUE
2000    CONTINUE
3000  CONTINUE
      DO 3500 I=1,NP
        DO 3500 J=1,NP
          DO 3500 K=1,NP
            TI=AIMAG(SURF(I,J,K))
3500        SURF(I,J,K)=CMPLX(0.,TI)
      DO 4000 I=1,NP
        DO 4000 J=1,NP
          DO 4000 K=1,NP
            TI=AIMAG(SURF(I,J,K))
            IF (ISURF(I,J,K).EQ.1) THEN
              SURF(I,J,K)=CMPLX(1.,TI)
            ELSE
c
              SURF(I,J,K)=CMPLX(RHO,TI)
c
            END IF
4000  CONTINUE
      RETURN
      END
