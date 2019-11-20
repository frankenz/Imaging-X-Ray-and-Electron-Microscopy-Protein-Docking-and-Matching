C------------------------------------------------------------------
       SUBROUTINE PIXELreal(NP,SURF,XR)
C------------------------------------------------------------------
       include 'param.h'
       COMPLEX SURF(np,np,np)
C
       NPXL=0
       DO I1=1,NP
         DO I2=1,NP
           DO I3=1,NP
             IF (real(SURF(I1,I2,I3)).EQ.xr) NPXL=NPXL+1
           END DO
         END DO
       END DO
       PRINT 10,XR,NPXL
  10   FORMAT(2X,'THE NUMBER OF REAL PIXELS (',F6.2,') IS ',I8)
       RETURN
       END

