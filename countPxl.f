C Pixels subroutine changed so it will fit my tstFunction for mFFT
C------------------------------------------------------------------
       SUBROUTINE countPxl (Size, Arr, XR, XI)
C------------------------------------------------------------------



       COMPLEX Arr(Size,Size,Size),TEST
       NPXL=0
C

       print*, 'enter countPxl'
       TEST=CMPLX(XR,XI)


       print*, 'enter countPxl 2'


       DO I1=1, Size
C         print*, 'inside loops'
         DO I2=1, Size
           DO I3=1, Size

C              print*, 'inside loops', NPXL 
C              IF (REAL (Arr(I1, I2, I3)) .EQ. 1.0) NPXL=NPXL+1
C             IF (Arr(I1,I2,I3).EQ.TEST) NPXL=NPXL+1
           END DO
         END DO
       END DO
       PRINT 10,XR,XI,NPXL
  10   FORMAT(2X,'THE NUMBER OF PIXELS (',F6.2,',',F6.2,') IS ',I8)
       RETURN
       END


