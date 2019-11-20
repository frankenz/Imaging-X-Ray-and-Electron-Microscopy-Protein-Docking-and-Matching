C------------------------------------------------------------------
       SUBROUTINE PIXELS (NP,SURF,XR,XI)
C------------------------------------------------------------------
       include 'param.h'
       COMPLEX SURF(np,np,np),TEST
C
       TEST=CMPLX(XR,XI)
       NPXL=0
       DO I1=1,NP
         DO I2=1,NP
           DO I3=1,NP
             IF (SURF(I1,I2,I3).EQ.TEST) NPXL=NPXL+1
           END DO
         END DO
       END DO
       PRINT 10,XR,XI,NPXL
  10   FORMAT(2X,'THE NUMBER OF PIXELS (',F6.2,',',F6.2,') IS ',I8)
       
C       write (26, 100) ,XR,XI,NPXL 
C 100   FORMAT(2X,'THE NUMBER OF PIXELS (',F6.2,',',F6.2,') IS ',I8)
       
       RETURN
       END

